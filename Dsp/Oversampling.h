// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <vector>
#include <cstddef>
#include <cassert>
#include <algorithm>

namespace eon {

// ---------------------------------------------------------------------------
// Pure-DSP half-band polyphase oversampler — no framework dependency, so the
// same code serves JUCE plugins, iPlug2 projects and headless measurement
// tools. Kaiser-windowed half-band FIRs, coefficients embedded:
//
//   stage 0 (fs <-> 2fs):  N=71  -> -113 dB stopband, flat to 0.43*pi
//   stage 1 (2fs <-> 4fs): N=47  -> -102 dB
//   stage 2 (4fs <-> 8fs): N=35  ->  -90 dB
//
// The strongest filter sits at the lowest rate where aliasing would fold
// straight back into the audio band. Polyphase decomposition (verified to
// machine precision against direct convolution):
//
//   up  : out[2n]   = x[n - p - 1]
//         out[2n+1] = 2 * sum_i hE[i] * x[n - i]
//   down: y[m]      = 0.5 * e[m - p] + sum_j hE[j] * o[m - j]
//
// where hE = even-indexed taps of the half-band kernel (the odd taps are all
// zero except the centre 0.5), p = (N-3)/4.
// ---------------------------------------------------------------------------

struct HbTaps71
{
    static constexpr int p = 17;
    static constexpr double hE[36] = {
        -1.2477987897e-06,  1.1227480081e-05, -4.2213195345e-05,  1.1705190499e-04,
        -2.7257107958e-04,  5.6356400147e-04, -1.0666100923e-03,  1.8837546556e-03,
        -3.1467476030e-03,  5.0238835479e-03, -7.7341787388e-03,  1.1579539945e-02,
        -1.7020470787e-02,  2.4864705852e-02, -3.6794020275e-02,  5.7171141772e-02,
        -1.0208649494e-01,  3.1694998446e-01,  3.1694998446e-01, -1.0208649494e-01,
         5.7171141772e-02, -3.6794020275e-02,  2.4864705852e-02, -1.7020470787e-02,
         1.1579539945e-02, -7.7341787388e-03,  5.0238835479e-03, -3.1467476030e-03,
         1.8837546556e-03, -1.0666100923e-03,  5.6356400147e-04, -2.7257107958e-04,
         1.1705190499e-04, -4.2213195345e-05,  1.1227480081e-05, -1.2477987897e-06 };
};

struct HbTaps47
{
    static constexpr int p = 11;
    static constexpr double hE[24] = {
        -4.9151111759e-06,  6.5122378319e-05, -2.8715247150e-04,  8.7813478063e-04,
        -2.1786964615e-03,  4.6960948225e-03, -9.1540676449e-03,  1.6629580781e-02,
        -2.9012452245e-02,  5.0748369748e-02, -9.7845787349e-02,  3.1546742160e-01,
         3.1546742160e-01, -9.7845787349e-02,  5.0748369748e-02, -2.9012452245e-02,
         1.6629580781e-02, -9.1540676449e-03,  4.6960948225e-03, -2.1786964615e-03,
         8.7813478063e-04, -2.8715247150e-04,  6.5122378319e-05, -4.9151111759e-06 };
};

struct HbTaps35
{
    static constexpr int p = 8;
    static constexpr double hE[18] = {
         1.7121717578e-05, -2.6899266591e-04,  1.2540227767e-03, -3.9265122295e-03,
         9.8371917699e-03, -2.1471419812e-02,  4.3766261984e-02, -9.2881796171e-02,
         3.1366794053e-01,  3.1366794053e-01, -9.2881796171e-02,  4.3766261984e-02,
        -2.1471419812e-02,  9.8371917699e-03, -3.9265122295e-03,  1.2540227767e-03,
        -2.6899266591e-04,  1.7121717578e-05 };
};

// One 2x stage. Each instance owns delay history for one processing direction.
template <typename Taps>
struct HalfBand2x
{
    static constexpr int L = sizeof (Taps::hE) / sizeof (double);
    std::vector<double> dl;      // shared delay line
    int w = 0;                   // write cursor (circular)

    void reset()
    {
        std::fill (dl.begin(), dl.end(), 0.0);
        std::fill (std::begin (eBuf), std::end (eBuf), 0.0);
        w = ew = 0;
    }
    void init() { dl.assign (L, 0.0); reset(); }

    inline void push (double x)
    {
        dl[(size_t) w] = x;
        w = (w + 1) % L;
    }
    inline double at (int age) const        // age 0 = newest
    {
        return dl[(size_t) ((w - 1 - age + 2 * L) % L)];
    }

    // consume one base-rate sample, emit a pair
    inline void up (double x, double& e, double& o)
    {
        push (x);
        e = at (Taps::p + 1);
        double acc = 0.0;
        for (int i = 0; i < L; ++i) acc += Taps::hE[i] * at (i);
        o = 2.0 * acc;
    }

    // consume an (even, odd) pair at the doubled rate, emit one base sample
    inline double down (double e, double o)
    {
        push (o);                            // odd samples run the FIR branch
        const double ev = eDelayed (e);
        double acc = 0.0;
        for (int j = 0; j < L; ++j) acc += Taps::hE[j] * at (j);
        return ev + acc;
    }

    // even branch needs its own p-deep delay (small fixed line; p <= 17)
    static constexpr int eCap = 64;
    double eBuf[eCap] = {};
    int ew = 0;
    inline double eDelayed (double e)
    {
        eBuf[ew] = e;
        const double v = 0.5 * eBuf[(ew - Taps::p + eCap) % eCap];
        ew = (ew + 1) % eCap;
        return v;
    }
};

// Cascaded 2x stages -> 2x / 4x / 8x. Stage order matters: index 0 is the
// strongest filter (lowest rate). Upsample applies 0->N, downsample N->0.
struct Oversampler
{
    HalfBand2x<HbTaps71> s0;
    HalfBand2x<HbTaps47> s1;
    HalfBand2x<HbTaps35> s2;
    int stages = 2;                        // default 4x
    int capacity = 0;                      // largest host-rate block accepted
    std::vector<float> tmp, tmp2;          // ping-pong scratch for mid rates
    HalfBand2x<HbTaps71> down0;
    HalfBand2x<HbTaps47> down1;
    HalfBand2x<HbTaps35> down2;

    void setStages (int s) { stages = std::clamp (s, 1, 3); }
    int  factor() const { return 1 << stages; }

    void prepare (int maxBlock)
    {
        s0.init(); down0.init();
        s1.init(); down1.init();
        s2.init(); down2.init();
        capacity = std::max (0, maxBlock);
        tmp.assign  ((size_t) maxBlock * 8u + 8u, 0.0f);
        tmp2.assign ((size_t) maxBlock * 8u + 8u, 0.0f);
    }
    void reset()
    {
        s0.reset(); down0.reset();
        s1.reset(); down1.reset();
        s2.reset(); down2.reset();
    }

    // in[n] -> out[n * factor()]; `out` must hold n * factor() floats and
    // must NOT alias `in`.
    void up (const float* in, int n, float* out)
    {
        assert (n <= capacity && "Oversampler::prepare(maxBlock) was not called for this block size");
        if (stages == 1)
        {
            upStage (in, n, out, 0);
            return;
        }
        upStage (in, n, tmp.data(), 0);
        const float* src = tmp.data();      // only used when stages == 3
        for (int s = 1; s < stages; ++s)
        {
            const bool last = (s == stages - 1);
            upStage (src, n << s, last ? out : tmp2.data(), s);
            src = tmp2.data();
        }
    }

    void upStage (const float* in, int n, float* out, int stage)
    {
        for (int i = 0; i < n; ++i)
        {
            double e, o;
            if (stage == 0)      s0.up (in[i], e, o);
            else if (stage == 1) s1.up (in[i], e, o);
            else                 s2.up (in[i], e, o);
            out[2 * i] = (float) e; out[2 * i + 1] = (float) o;
        }
    }

    // in[n * factor()] -> out[n]
    void down (const float* in, float* out, int n)
    {
        assert (n <= capacity && "Oversampler::prepare(maxBlock) was not called for this block size");
        if (stages == 1)
        {
            for (int i = 0; i < n; ++i)
                out[i] = (float) down0.down (in[2 * i], in[2 * i + 1]);
            return;
        }
        // downsample highest stage first: stage (stages-1) ... stage 0
        downStage (in, tmp.data(), n << (stages - 1), stages - 1);
        for (int s = stages - 2; s > 0; --s)
            downStage (tmp.data(), tmp.data(), n << s, s);
        downStage (tmp.data(), out, n, 0);
    }

    void downStage (const float* in, float* out, int n, int stage)
    {
        for (int i = 0; i < n; ++i)
        {
            const double y = stage == 0 ? down0.down (in[2 * i], in[2 * i + 1])
                           : stage == 1 ? down1.down (in[2 * i], in[2 * i + 1])
                                        : down2.down (in[2 * i], in[2 * i + 1]);
            out[i] = (float) y;
        }
    }

    // Roundtrip latency (up + down) in base-rate samples.
    double latencySamples() const
    {
        double l = 0.0;
        for (int s = 0; s < stages; ++s)
        {
            const int N = s == 0 ? 71 : s == 1 ? 47 : 35;
            l += (double) (N - 1) / (double) (1 << (s + 1));
        }
        return l;
    }
};

} // namespace eon
