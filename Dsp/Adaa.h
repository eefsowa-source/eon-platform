#pragma once
#include <cmath>

namespace eon {

// ---------------------------------------------------------------------------
// ADAA — Antiderivative Anti-Aliasing (Parker/Bilbao)
//
// ADAA1: y[n] = (F1(x0) - F1(x1)) / (x0 - x1),      F1' = f
// ADAA2: y[n] = 2 * (D01 - D12) / (x0 - x2),        Dij = (F2(xi) - F2(xj))/(xi - xj), F2' = F1
//
// State and all antiderivative arithmetic are DOUBLE. This is not optional:
// the second divided difference subtracts F2 values that differ by ~1e-6
// over ~1e-5 intervals — in float32 the quantization noise is the same order
// as the result and the output oscillates wildly near extrema (verified:
// +-0.7 error on a sine peak). With doubles the noise floor is irrelevant.
//
// Requires *exact* antiderivatives. When consecutive samples are too close,
// fall back to midpoint evaluations (F2' = F1) instead of dividing by ~0.
// ---------------------------------------------------------------------------

struct ADAA1
{
    double x1 = 0.0, F1x1 = 0.0;

    void reset() { x1 = 0.0; F1x1 = 0.0; }

    template <typename F, typename F1>
    inline float process (float x0, F f, F1 f1)
    {
        const double F1x0 = f1 (x0);
        const double dx = (double) x0 - x1;
        const double y = std::abs (dx) < 1e-9 ? f (0.5 * ((double) x0 + x1))
                                              : (F1x0 - F1x1) / dx;
        x1 = x0; F1x1 = F1x0;
        return (float) y;
    }
};

struct ADAA2
{
    double x1 = 0.0, x2 = 0.0;
    double F2x1 = 0.0, F2x2 = 0.0;

    void reset() { x1 = x2 = 0.0; F2x1 = F2x2 = 0.0; }

    template <typename F, typename F1, typename F2>
    inline float process (float x0, F f, F1 f1, F2 f2)
    {
        constexpr double eps = 1e-9;
        const double F2x0 = f2 (x0);
        const double dx01 = (double) x0 - x1;
        const double dx12 = x1 - x2;
        const double dx02 = (double) x0 - x2;

        double y;
        if (std::abs (dx02) < eps)
        {
            y = f (x0); // x0 ~= x1 ~= x2: converged
        }
        else
        {
            const double d01 = std::abs (dx01) < eps ? f1 (0.5 * ((double) x0 + x1))
                                                     : (F2x0 - F2x1) / dx01;
            const double d12 = std::abs (dx12) < eps ? f1 (0.5 * (x1 + x2))
                                                     : (F2x1 - F2x2) / dx12;
            y = 2.0 * (d01 - d12) / dx02;
        }

        x2 = x1; x1 = x0;
        F2x2 = F2x1; F2x1 = F2x0;
        return (float) y;
    }
};

// ---------------------------------------------------------------------------
// tanh: exact first antiderivative exists -> use with ADAA1.
//   F1(v) = ln(cosh v) = |v| + ln(1 + e^-2|v|) - ln 2   (numerically stable)
// (F2 for tanh has no elementary form — do NOT fake it with polynomials that
//  are discontinuous; the original project's piecewise F2 had a jump at |v|=2.)
// ---------------------------------------------------------------------------

inline double tanhF1 (double v)
{
    const double a = std::abs (v);
    return a + std::log1p (std::exp (-2.0 * a)) - 0.6931471805599453;
}

struct TanhSat
{
    ADAA1 adaa;
    void reset() { adaa.reset(); }
    inline float process (float x)
    {
        return adaa.process (x,
            [] (double v) { return std::tanh (v); },
            [] (double v) { return tanhF1 (v); });
    }
};

// ---------------------------------------------------------------------------
// Cubic soft clip: f(x) = x - x^3/3 for |x| <= 1, else sign(x)*2/3.
// Exact, C1-continuous F1 and F2 -> true ADAA2 with no discontinuities.
//   F1(x) = x^2/2 - x^4/12            (|x|<=1), even extension beyond
//   F2(x) = x^3/6 - x^5/60            (|x|<=1), odd extension beyond
// ---------------------------------------------------------------------------

struct SoftClip
{
    static inline double f (double x)
    {
        const double a = std::abs (x);
        return a <= 1.0 ? x - x * x * x / 3.0
                        : std::copysign (2.0 / 3.0, x);
    }

    static inline double F1 (double x)
    {
        const double a = std::abs (x);
        if (a <= 1.0)
            return x * x / 2.0 - x * x * x * x / 12.0;
        // value at 1: 1/2 - 1/12 = 5/12; slope beyond = f(1) = 2/3
        return 5.0 / 12.0 + (2.0 / 3.0) * (a - 1.0);
    }

    static inline double F2 (double x)
    {
        const double a = std::abs (x);
        if (a <= 1.0)
            return x * x * x / 6.0 - x * x * x * x * x / 60.0;
        // value at 1: 1/6 - 1/60 = 3/20; slope = F1(1) = 5/12
        const double t = a - 1.0;
        return std::copysign (t * t / 3.0 + (5.0 / 12.0) * t + 3.0 / 20.0, x);
    }
};

struct SoftClipSat
{
    ADAA2 adaa;
    void reset() { adaa.reset(); }
    inline float process (float x)
    {
        return adaa.process (x, SoftClip::f, SoftClip::F1, SoftClip::F2);
    }
};

} // namespace eon
