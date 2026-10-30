// Guards the half-band optimisations against a naive scalar reference.
//
// HalfBand2x no longer computes the textbook sum. The delay line is mirrored
// (2L slots) so the FIR window is contiguous, and the taps are evaluated in
// symmetric-folded form with four accumulators, which changes both the
// summation grouping and the order. None of that is supposed to move the
// output, so this file keeps an obviously-correct scalar implementation of the
// same filter and compares against it.
//
// The comparison is exact at float precision, because that is the precision the
// oversampler actually delivers: the double-level difference is a couple of
// ulps (~5e-16), which is roughly 5e-9 of a float ulp, so it rounds away.
// "folded FIR stays within a few double ulps of the scalar sum" pins that
// relationship down directly, so a future change that widens it fails here
// rather than silently drifting.

#include "Dsp/Oversampling.h"
#include "Tests/TestHarness.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Deterministic input. A fixed seed keeps the corpus identical across runs and
// platforms, so a failure is reproducible rather than a floating-point lottery.
// ---------------------------------------------------------------------------
struct Rng
{
    uint64_t state = 0x9E3779B97F4A7C15ull;

    uint64_t next()
    {
        uint64_t x = state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        state = x;
        return x * 0x2545F4914F6CDD1Dull;
    }

    double uniform() { return static_cast<double>(next() >> 11) / 9007199254740992.0; }
    double bipolar() { return uniform() * 2.0 - 1.0; }
};

// ---------------------------------------------------------------------------
// Naive reference: one delay line of exactly L slots, modulo indexing, plain
// running sum, no folding and no mirroring. This is the form the header had
// before any of the optimisations, and the form a reader can check by eye.
// ---------------------------------------------------------------------------
template <typename Taps>
struct RefHalfBand2x
{
    static constexpr int L = sizeof(Taps::hE) / sizeof(double);
    static constexpr int eCap = 64;

    std::vector<double> dl;
    double eBuf[eCap] = {};
    int w = 0;
    int ew = 0;

    void init()
    {
        dl.assign(static_cast<size_t>(L), 0.0);
        reset();
    }

    void reset()
    {
        std::fill(dl.begin(), dl.end(), 0.0);
        std::fill(std::begin(eBuf), std::end(eBuf), 0.0);
        w = ew = 0;
    }

    void push(double x)
    {
        dl[static_cast<size_t>(w)] = x;
        w = (w + 1) % L;
    }

    double at(int age) const
    {
        return dl[static_cast<size_t>((w - 1 - age + 2 * L) % L)];
    }

    double eDelayed(double e)
    {
        eBuf[ew] = e;
        const double v = 0.5 * eBuf[(ew - Taps::p + eCap) % eCap];
        ew = (ew + 1) % eCap;
        return v;
    }

    void up(double x, double& e, double& o)
    {
        push(x);
        e = at(Taps::p + 1);
        double acc = 0.0;
        for (int i = 0; i < L; ++i)
            acc += Taps::hE[i] * at(i);
        o = 2.0 * acc;
    }

    double down(double e, double o)
    {
        push(o);
        const double ev = eDelayed(e);
        double acc = 0.0;
        for (int j = 0; j < L; ++j)
            acc += Taps::hE[j] * at(j);
        return ev + acc;
    }
};

// Same cascade as eon::Oversampler, built from the naive stages.
struct RefOversampler
{
    RefHalfBand2x<eon::HbTaps71> s0, down0;
    RefHalfBand2x<eon::HbTaps47> s1, down1;
    RefHalfBand2x<eon::HbTaps35> s2, down2;
    RefHalfBand2x<eon::HbTaps51> s0e, down0e;

    int stages = 2;
    int capacity = 0;
    bool ecoStage0 = false;
    std::vector<float> tmp, tmp2;

    void setStages(int s) { stages = std::clamp(s, 1, 3); }
    int factor() const { return 1 << stages; }

    void prepare(int maxBlock)
    {
        s0.init(); down0.init();
        s0e.init(); down0e.init();
        s1.init(); down1.init();
        s2.init(); down2.init();
        capacity = std::max(0, maxBlock);
        tmp.assign(static_cast<size_t>(capacity) * 8u + 8u, 0.0f);
        tmp2.assign(static_cast<size_t>(capacity) * 8u + 8u, 0.0f);
    }

    void reset()
    {
        s0.reset(); down0.reset();
        s0e.reset(); down0e.reset();
        s1.reset(); down1.reset();
        s2.reset(); down2.reset();
    }

    void upStage(const float* in, int n, float* out, int stage)
    {
        for (int i = 0; i < n; ++i)
        {
            double e, o;
            if (stage == 0)      { if (ecoStage0) s0e.up(in[i], e, o); else s0.up(in[i], e, o); }
            else if (stage == 1) s1.up(in[i], e, o);
            else                 s2.up(in[i], e, o);
            out[2 * i] = static_cast<float>(e);
            out[2 * i + 1] = static_cast<float>(o);
        }
    }

    void up(const float* in, int n, float* out)
    {
        n = std::clamp(n, 0, capacity);
        if (stages == 1)
        {
            upStage(in, n, out, 0);
            return;
        }
        upStage(in, n, tmp.data(), 0);
        const float* src = tmp.data();
        for (int s = 1; s < stages; ++s)
        {
            const bool last = (s == stages - 1);
            upStage(src, n << s, last ? out : tmp2.data(), s);
            src = tmp2.data();
        }
    }

    void downStage(const float* in, float* out, int n, int stage)
    {
        for (int i = 0; i < n; ++i)
        {
            const double y = stage == 0
                ? (ecoStage0 ? down0e.down(in[2 * i], in[2 * i + 1])
                             : down0.down(in[2 * i], in[2 * i + 1]))
                : stage == 1 ? down1.down(in[2 * i], in[2 * i + 1])
                             : down2.down(in[2 * i], in[2 * i + 1]);
            out[i] = static_cast<float>(y);
        }
    }

    void down(const float* in, float* out, int n)
    {
        n = std::clamp(n, 0, capacity);
        if (stages == 1)
        {
            for (int i = 0; i < n; ++i)
                out[i] = static_cast<float>(ecoStage0 ? down0e.down(in[2 * i], in[2 * i + 1])
                                                      : down0.down(in[2 * i], in[2 * i + 1]));
            return;
        }
        downStage(in, tmp.data(), n << (stages - 1), stages - 1);
        for (int s = stages - 2; s > 0; --s)
            downStage(tmp.data(), tmp.data(), n << s, s);
        downStage(tmp.data(), out, n, 0);
    }
};

std::string describe(int stages, int eco, const char* what, int index,
                     float actual, float expected)
{
    char buffer[256];
    std::snprintf(buffer, sizeof buffer,
                  "stages=%d eco=%d %s[%d]: optimised=%.9g reference=%.9g",
                  stages, eco, what, index,
                  static_cast<double>(actual), static_cast<double>(expected));
    return std::string(buffer);
}

} // namespace

// ---------------------------------------------------------------------------
// End to end: the shipping oversampler must reproduce the naive cascade exactly
// at float precision, across every stage count, both stage-0 tap sets, random
// block sizes and a mid-run reset.
// ---------------------------------------------------------------------------
EON_TEST_CASE("optimised oversampler matches the naive cascade at float precision")
{
    constexpr int max_block = 256;
    constexpr int total_samples = 20000;

    for (int stages = 1; stages <= 3; ++stages)
    {
        for (int eco = 0; eco <= 1; ++eco)
        {
            eon::Oversampler optimised;
            optimised.setStages(stages);
            optimised.ecoStage0 = (eco != 0);
            optimised.prepare(max_block);

            RefOversampler reference;
            reference.setStages(stages);
            reference.ecoStage0 = (eco != 0);
            reference.prepare(max_block);

            std::vector<float> in(static_cast<size_t>(max_block));
            std::vector<float> up_opt(static_cast<size_t>(max_block) * 8u);
            std::vector<float> up_ref(static_cast<size_t>(max_block) * 8u);
            std::vector<float> down_opt(static_cast<size_t>(max_block));
            std::vector<float> down_ref(static_cast<size_t>(max_block));

            Rng rng;
            int done = 0;
            bool reset_done = false;

            while (done < total_samples)
            {
                int n = 1 + static_cast<int>(rng.next() % 64u);
                if (done + n > total_samples)
                    n = total_samples - done;

                for (int i = 0; i < n; ++i)
                    in[static_cast<size_t>(i)] = static_cast<float>(rng.bipolar());

                optimised.up(in.data(), n, up_opt.data());
                reference.up(in.data(), n, up_ref.data());

                const int up_count = n * optimised.factor();
                for (int i = 0; i < up_count; ++i)
                {
                    if (up_opt[static_cast<size_t>(i)] != up_ref[static_cast<size_t>(i)])
                        eon::test::fail("up output matches the naive cascade", __FILE__, __LINE__,
                                        describe(stages, eco, "up", i,
                                                 up_opt[static_cast<size_t>(i)],
                                                 up_ref[static_cast<size_t>(i)]));
                }

                optimised.down(up_opt.data(), down_opt.data(), n);
                reference.down(up_ref.data(), down_ref.data(), n);

                for (int i = 0; i < n; ++i)
                {
                    if (down_opt[static_cast<size_t>(i)] != down_ref[static_cast<size_t>(i)])
                        eon::test::fail("down output matches the naive cascade", __FILE__, __LINE__,
                                        describe(stages, eco, "down", i,
                                                 down_opt[static_cast<size_t>(i)],
                                                 down_ref[static_cast<size_t>(i)]));
                }

                done += n;

                if (!reset_done && done >= total_samples / 2)
                {
                    optimised.reset();
                    reference.reset();
                    reset_done = true;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The structural half of the same claim: folding and the four accumulators only
// regroup the sum, so the double-precision result may move by a couple of ulps
// but not more. Measured worst case on this corpus is ~5.6e-16 for inputs in
// [-1, 1]; 1e-14 leaves room for other compilers without hiding a real bug.
// ---------------------------------------------------------------------------
EON_TEST_CASE("folded FIR stays within a few double ulps of the scalar sum")
{
    constexpr double tolerance = 1e-14;
    constexpr int iterations = 20000;

    auto probe = [&](auto tag)
    {
        using Taps = decltype(tag);
        eon::HalfBand2x<Taps> stage;
        stage.init();

        Rng rng;
        double worst = 0.0;

        for (int it = 0; it < iterations; ++it)
        {
            stage.push(rng.bipolar());

            double scalar = 0.0;
            for (int i = 0; i < eon::HalfBand2x<Taps>::L; ++i)
                scalar += Taps::hE[i] * stage.at(i);

            worst = std::max(worst, std::abs(stage.folded() - scalar));
        }

        EON_CHECK_NEAR(worst, 0.0, tolerance);
    };

    probe(eon::HbTaps71{});
    probe(eon::HbTaps47{});
    probe(eon::HbTaps35{});
    probe(eon::HbTaps51{});
}
