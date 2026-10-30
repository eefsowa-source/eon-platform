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
// tools. Least-squares half-band FIRs (tools/design_hb_ls.py, firls with
// stopband weight 256, even-offset taps snapped to the half-band form).
// Passband flat to 0.43*pi kernel-rate (~20.6 kHz at 48 kHz), stopband above
// 0.57*pi:
//
//   stage 0 (fs <-> 2fs):  N=71  -> -65 dB stopband
//   stage 1 (2fs <-> 4fs): N=47  -> -42 dB
//   stage 2 (4fs <-> 8fs): N=35  -> -31 dB
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
// zero except the centre 0.5), p = (N-3)/4. Both FIR sums are evaluated in
// symmetric-folded form (see HalfBand2x::folded): the same sum up to float
// rounding, at half the multiplies.
// ---------------------------------------------------------------------------

struct HbTaps71
{
    static constexpr int p = 17;
    static constexpr double hE[36] = {
         3.3289602328e-05,  3.4356458805e-07, -7.8164908403e-05,  2.4337014826e-04,
        -5.4402521075e-04,  1.0401463111e-03, -1.8052237108e-03,  2.9281700804e-03,
        -4.5172283335e-03,  6.7082677504e-03, -9.6822483189e-03,  1.3702185294e-02,
        -1.9194495770e-02,  2.6943212850e-02, -3.8621329087e-02,  5.8600319255e-02,
        -1.0299698567e-01,  3.1726266868e-01,  3.1726266868e-01, -1.0299698567e-01,
         5.8600319255e-02, -3.8621329087e-02,  2.6943212850e-02, -1.9194495770e-02,
         1.3702185294e-02, -9.6822483189e-03,  6.7082677504e-03, -4.5172283335e-03,
         2.9281700804e-03, -1.8052237108e-03,  1.0401463111e-03, -5.4402521075e-04,
         2.4337014826e-04, -7.8164908403e-05,  3.4356458805e-07,  3.3289602328e-05 };
};

struct HbTaps47
{
    static constexpr int p = 11;
    static constexpr double hE[24] = {
         4.2156740182e-04,  2.9450201696e-04, -1.2900127880e-03,  2.9445653780e-03,
        -5.5192454860e-03,  9.3367347080e-03, -1.4866523379e-02,  2.2915362402e-02,
        -3.5155388620e-02,  5.5932624902e-02, -1.0131548269e-01,  3.1668826443e-01,
         3.1668826443e-01, -1.0131548269e-01,  5.5932624902e-02, -3.5155388620e-02,
         2.2915362402e-02, -1.4866523379e-02,  9.3367347080e-03, -5.5192454860e-03,
         2.9445653780e-03, -1.2900127880e-03,  2.9450201696e-04,  4.2156740182e-04 };
};

struct HbTaps35
{
    static constexpr int p = 8;
    static constexpr double hE[18] = {
        -1.5356513996e-03, -1.7148253779e-03,  5.3182484742e-03, -1.0759138535e-02,
         1.8981958009e-02, -3.1692856342e-02,  5.3222964522e-02, -9.9589022715e-02,
         3.1609537124e-01,  3.1609537124e-01, -9.9589022715e-02,  5.3222964522e-02,
        -3.1692856342e-02,  1.8981958009e-02, -1.0759138535e-02,  5.3182484742e-03,
        -1.7148253779e-03, -1.5356513996e-03 };
};

// Eco stage-0 set: same LS spec (0.43*pi edge) at 51 taps. -53 dB ripple for
// 13 MACs per output sample vs HbTaps71's -65 dB for 18 -- both counts are the
// folded polyphase loop below (one multiply per tap pair) -- for hosts that
// trade a little image rejection for a cheaper 2x path.
struct HbTaps51
{
    static constexpr int p = 12;
    static constexpr double hE[26] = {
         3.5315496922e-04, -8.7769728054e-04,  1.7554353133e-03, -3.1152835401e-03,
         5.1154070483e-03, -7.9584281393e-03,  1.1926000692e-02, -1.7457920319e-02,
         2.5344539270e-02, -3.7256989657e-02,  5.7556431367e-02, -1.0234153397e-01,
         3.1703918947e-01,  3.1703918947e-01, -1.0234153397e-01,  5.7556431367e-02,
        -3.7256989657e-02,  2.5344539270e-02, -1.7457920319e-02,  1.1926000692e-02,
        -7.9584281393e-03,  5.1154070483e-03, -3.1152835401e-03,  1.7554353133e-03,
        -8.7769728054e-04,  3.5315496922e-04 };
};

// One 2x stage. Each instance owns delay history for one processing direction.
//
// The delay line is indexed by hand instead of with `%`, and the half-band taps
// are folded: hE[i] == hE[L-1-i], so
//
//   sum_{i<L} hE[i] * at(i) == sum_{i<L/2} hE[i] * (at(i) + at(L-1-i))
//
// which halves the multiplies and lets the FIR loop walk two cursors that wrap
// on a compare. The two sums differ only by float rounding.
namespace detail {

// The folded FIR equals the plain convolution only when hE is symmetric, so
// every tap set is checked at compile time.
template <typename Taps>
constexpr bool hbTapsSymmetric()
{
    constexpr int L = sizeof (Taps::hE) / sizeof (double);
    for (int i = 0; i < L / 2; ++i)
        if (Taps::hE[i] != Taps::hE[L - 1 - i])
            return false;
    return true;
}

} // namespace detail

template <typename Taps>
struct HalfBand2x
{
    static constexpr int L = sizeof (Taps::hE) / sizeof (double);
    static_assert (L % 2 == 0, "half-band hE must hold an even number of taps");
    static_assert (Taps::p + 2 < L, "even-branch delay must fit the odd-branch line");
    static_assert (detail::hbTapsSymmetric<Taps>(),
                   "half-band hE must be symmetric for folded evaluation");

    std::vector<double> dl;      // shared delay line
    int w = 0;                   // write cursor (circular)

    void reset()
    {
        std::fill (dl.begin(), dl.end(), 0.0);
        std::fill (std::begin (eBuf), std::end (eBuf), 0.0);
        w = ew = 0;
    }
    // Mirrored delay line: every push writes the low and high copies, so the
    // FIR window is a contiguous read with no modulo per tap.
    void init() { dl.assign (2 * L, 0.0); reset(); }

    inline void push (double x)
    {
        dl[(size_t) w] = x;
        dl[(size_t) (w + L)] = x;
        w = (w + 1 < L) ? w + 1 : 0;
    }

    // Mirrored line: the low and high copies always hold the same sample, so
    // every age below L + w is a valid index and no modulo is needed. Ages at
    // or above that would wrap past the array, which is why the callers stay
    // inside L.
    inline double at (int age) const        // age 0 = newest
    {
        return dl[(size_t) (w - 1 - age + L)];
    }

    // Symmetric-folded FIR over the mirrored window. Ages (i, L-1-i) land at
    // dl[w+i] and dl[w-1+L-i], so the two cursors become contiguous reads
    // walking in opposite directions and the loop carries no wrap branches.
    inline double folded() const
    {
        const double* fwd = dl.data() + (size_t) w;              // age L-1 -> 0
        const double* bwd = dl.data() + (size_t) (w - 1 + L);    // age 0 -> L-1
        // Four accumulators break the sequential add chain; the pairwise
        // reduction at the end keeps the result identical to the scalar sum.
        double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
        int i = 0;
        for (; i + 3 < L / 2; i += 4)
        {
            a0 += Taps::hE[i]     * (fwd[i]     + bwd[-i]);
            a1 += Taps::hE[i + 1] * (fwd[i + 1] + bwd[-(i + 1)]);
            a2 += Taps::hE[i + 2] * (fwd[i + 2] + bwd[-(i + 2)]);
            a3 += Taps::hE[i + 3] * (fwd[i + 3] + bwd[-(i + 3)]);
        }
        double acc = (a0 + a1) + (a2 + a3);
        for (; i < L / 2; ++i)
            acc += Taps::hE[i] * (fwd[i] + bwd[-i]);
        return acc;
    }

    // consume one base-rate sample, emit a pair
    inline void up (double x, double& e, double& o)
    {
        push (x);
        e = at (Taps::p + 1);
        o = 2.0 * folded();
    }

    // consume an (even, odd) pair at the doubled rate, emit one base sample
    inline double down (double e, double o)
    {
        push (o);                            // odd samples run the FIR branch
        return eDelayed (e) + folded();
    }

    // even branch needs its own p-deep delay (small fixed line; p <= 17)
    static constexpr int eCap = 64;
    static_assert ((eCap & (eCap - 1)) == 0, "eCap must be a power of two");
    double eBuf[eCap] = {};
    int ew = 0;
    inline double eDelayed (double e)
    {
        eBuf[ew] = e;
        const double v = 0.5 * eBuf[(ew - Taps::p + eCap) & (eCap - 1)];
        ew = (ew + 1) & (eCap - 1);
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
    HalfBand2x<HbTaps51> s0e;              // used when ecoStage0 is set
    int stages = 2;                        // default 4x
    int capacity = 0;                      // largest host-rate block accepted
    std::vector<float> tmp, tmp2;          // ping-pong scratch for mid rates
    HalfBand2x<HbTaps71> down0;
    HalfBand2x<HbTaps47> down1;
    HalfBand2x<HbTaps35> down2;
    HalfBand2x<HbTaps51> down0e;

    // Swap stage 0 (the strongest, most expensive filter) for the 51-tap eco
    // set in both directions. Leave false for the full -65 dB stage.
    //
    // The eco and full stage-0 filters keep SEPARATE delay histories and have
    // different latency. Change this flag only while the stream is stopped,
    // followed by reset() (and a latency re-report): flipping it mid-stream
    // switches to an empty history and produces a click.
    bool ecoStage0 = false;

    void setStages (int s) { stages = std::clamp (s, 1, 3); }
    int  factor() const { return 1 << stages; }

    void prepare (int maxBlock)
    {
        s0.init(); down0.init();
        s0e.init(); down0e.init();
        s1.init(); down1.init();
        s2.init(); down2.init();
        capacity = std::max (0, maxBlock);
        tmp.assign  ((size_t) capacity * 8u + 8u, 0.0f);
        tmp2.assign ((size_t) capacity * 8u + 8u, 0.0f);
    }
    void reset()
    {
        s0.reset(); down0.reset();
        s0e.reset(); down0e.reset();
        s1.reset(); down1.reset();
        s2.reset(); down2.reset();
    }

    // in[n] -> out[n * factor()]; `out` must hold n * factor() floats and
    // must NOT alias `in`.
    void up (const float* in, int n, float* out)
    {
        assert (n <= capacity && "Oversampler::prepare(maxBlock) was not called for this block size");
        n = std::clamp (n, 0, capacity);    // release builds: clamp instead of overrunning tmp/tmp2
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
            if (stage == 0)      { if (ecoStage0) s0e.up (in[i], e, o); else s0.up (in[i], e, o); }
            else if (stage == 1) s1.up (in[i], e, o);
            else                 s2.up (in[i], e, o);
            out[2 * i] = (float) e; out[2 * i + 1] = (float) o;
        }
    }

    // in[n * factor()] -> out[n]
    void down (const float* in, float* out, int n)
    {
        assert (n <= capacity && "Oversampler::prepare(maxBlock) was not called for this block size");
        n = std::clamp (n, 0, capacity);
        if (stages == 1)
        {
            for (int i = 0; i < n; ++i)
                out[i] = (float) (ecoStage0 ? down0e.down (in[2 * i], in[2 * i + 1])
                                           : down0.down  (in[2 * i], in[2 * i + 1]));
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
            const double y = stage == 0 ? (ecoStage0 ? down0e.down (in[2 * i], in[2 * i + 1])
                                                     : down0.down  (in[2 * i], in[2 * i + 1]))
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
            const int N = s == 0 ? (ecoStage0 ? 51 : 71) : s == 1 ? 47 : 35;
            l += (double) (N - 1) / (double) (1 << (s + 1));
        }
        return l;
    }
};

} // namespace eon
