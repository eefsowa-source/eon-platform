// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <numbers>
#include <vector>
#include <algorithm>

namespace eon {
namespace measure {

// ---------------------------------------------------------------------------
// Shared measurement primitives — the same tools every plugin's test harness
// needs, so evidence stays comparable across the portfolio.
// ---------------------------------------------------------------------------

// Single-frequency DFT magnitude via direct sum (Goertzel-equivalent accuracy
// for non-bin-aligned frequencies). Phase-locked reference — robust for THD.
inline double binMag (const float* buf, size_t N, double freq, double sr)
{
    const double w = 2.0 * std::numbers::pi * freq / sr;
    double re = 0.0, im = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        re += (double) buf[i] * std::cos (w * (double) i);
        im -= (double) buf[i] * std::sin (w * (double) i);
    }
    return std::sqrt (re * re + im * im) * (2.0 / (double) N);
}

// Harmonic amplitudes H1..maxH of a recorded sine response.
inline std::vector<double> harmonics (const float* buf, size_t N,
                                      double f0, double sr, int maxH = 6)
{
    std::vector<double> h ((size_t) maxH);
    for (int k = 1; k <= maxH; ++k)
        h[(size_t) k - 1] = binMag (buf, N, f0 * (double) k, sr);
    return h;
}

// THD% = 100 * sqrt(H2^2+..+HmaxH^2) / H1
inline double thdPercent (const float* buf, size_t N, double f0, double sr,
                          int maxH = 6)
{
    const auto h = harmonics (buf, N, f0, sr, maxH);
    double p2 = 0.0;
    for (size_t k = 1; k < h.size(); ++k) p2 += h[k] * h[k];
    return 100.0 * std::sqrt (p2) / h[0];
}

// Static transfer curve: feed constant x, settle, record steady output.
// `proc` is any callable float -> float that carries state; reset between
// points is the caller's choice (continuous sweep keeps hysteresis visible).
template <typename F>
inline std::vector<std::pair<double, double>>
transferCurve (F&& proc, double xLo, double xHi, int points = 101, int settle = 64)
{
    std::vector<std::pair<double, double>> curve;
    curve.reserve ((size_t) points);
    for (int i = 0; i < points; ++i)
    {
        const double x = xLo + (xHi - xLo) * (double) i / (points - 1);
        double y = 0.0;
        for (int s = 0; s < settle; ++s) y = (double) proc ((float) x);
        curve.push_back ({ x, y });
    }
    return curve;
}

// Sine-fill helper: fill buf with `amp` sine at `freq`, `sr`, after `settle`
// warm-up samples processed identically (so states settle before capture).
template <typename F>
inline void renderSine (F&& proc, float* buf, size_t N, double freq,
                        double amp, double sr, int settle)
{
    const double w = 2.0 * std::numbers::pi * freq / sr;
    for (int i = 0; i < settle; ++i)
        proc ((float) (amp * std::sin (w * (double) i)));
    for (size_t i = 0; i < N; ++i)
        buf[i] = proc ((float) (amp * std::sin (w * (double) (i + settle))));
}

// Magnitude response at `freq` for an LTI-ish path: render sine, measure.
template <typename F>
inline double magAt (F&& proc, double freq, double amp, double sr,
                     size_t N = 8192, int settle = 4096)
{
    std::vector<float> buf (N);
    renderSine (proc, buf.data(), N, freq, amp, sr, settle);
    return binMag (buf.data(), N, freq, sr) / amp;
}

// Alias estimate: drive the processor with a sine at `freq` (meant to live
// above the audio band once nonlinearities spread it), return the magnitude
// of the component that folded back to `aliasFreq`.
template <typename F>
inline double foldedMagAt (F&& proc, double freq, double aliasFreq,
                           double amp, double sr,
                           size_t N = 8192, int settle = 4096)
{
    std::vector<float> buf (N);
    renderSine (proc, buf.data(), N, freq, amp, sr, settle);
    return binMag (buf.data(), N, aliasFreq, sr);
}

} // namespace measure
} // namespace eon
