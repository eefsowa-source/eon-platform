// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <algorithm>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace eon {

// ---------------------------------------------------------------------------
// Solver toolkit for circuit models.
//
// lambertW0      — principal branch of W(x) for x > 0, solved in log space so
//                  diode-scale arguments (R*Is/Vt * exp((R*Is+a)/Vt), which
//                  overflow double at a > ~0.7V) never materialise.
// newtonScalar   — damped Newton with a maintained sign-change bracket;
//                  falls back to bisection whenever a step leaves the bracket,
//                  so it converges on any continuous monotone-ish residual.
// solveDense     — Gaussian elimination with partial pivoting for the small
//                  dense systems that appear in MNA / DK-method solves (N <= 8).
// ---------------------------------------------------------------------------

// W0(x) for x > 0 via w + ln(w) = ln(x). ~5 iterations to 1e-14.
inline double lambertW0log (double logX)
{
    if (std::isnan (logX)) return logX;
    if (logX == -std::numeric_limits<double>::infinity()) return 0.0;
    if (logX == std::numeric_limits<double>::infinity()) return logX;

    // log(denorm_min) for IEEE-754 binary64, kept constant to avoid a
    // function-local runtime initializer on the audio thread.
    constexpr double smallestLog = -744.440071921381218089;
    if (logX < smallestLog) return 0.0;

    // Initial guess: asymptotic L - ln L for large logs, exp for small.
    double w = logX > 1.0 ? logX - std::log (logX) : std::exp (logX);
    if (! (w > 0.0) || ! std::isfinite (w)) return 0.0;

    for (int i = 0; i < 40; ++i)
    {
        const double f = w + std::log (w) - logX;
        if (f == 0.0) break;

        // The reciprocal derivative w/(w+1) avoids 1/w overflowing for
        // subnormal initial guesses. Backtracking keeps Newton in its domain.
        const double step = f * (w / (w + 1.0));
        double candidate = w - step;
        if (! (candidate > 0.0) || ! std::isfinite (candidate)
            || std::abs (candidate + std::log (candidate) - logX) >= std::abs (f))
        {
            double scale = 0.5;
            candidate = w;
            for (int backtrack = 0; backtrack < 64; ++backtrack)
            {
                const double damped = w - scale * step;
                if (damped > 0.0 && std::isfinite (damped)
                    && std::abs (damped + std::log (damped) - logX) < std::abs (f))
                {
                    candidate = damped;
                    break;
                }
                scale *= 0.5;
            }
        }

        if (candidate == w) break;
        const double delta = std::abs (candidate - w);
        w = candidate;
        if (delta < 1e-14 * std::max (1.0, std::abs (w))) break;
    }
    return w;
}

namespace detail {

// Sign test that cannot underflow. NaN operands report "no sign change" so a
// non-finite sample can never manufacture a bracket.
inline bool signDiffers (double a, double b)
{
    // Strict opposite signs: an exact zero endpoint is not a sign change.
    return std::isfinite (a) && std::isfinite (b)
        && ((a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0));
}

} // namespace detail

inline double lambertW0 (double x)
{
    if (x <= 0.0) return 0.0;
    return lambertW0log (std::log (x));
}

// Scalar Newton on f(x) = 0 inside [lo, hi]. f must be evaluable; df may be
// nullptr -> secant estimate. Keeps a bracket once a sign change is seen and
// bisects whenever Newton leaves it — guaranteed progress on sane inputs.
//
// Scale invariant: the bracket is detected by sign bits, not by f(lo)*f(hi)
// (that product underflows to zero for residuals below ~1e-154), and any
// non-zero finite derivative is used however small it is (an absolute
// |df| < 1e-20 cutoff silently disabled Newton for scaled residuals).
template <typename F, typename DF>
inline double newtonScalar (F&& f, DF&& df, double x0, double lo, double hi,
                            int maxIter = 20, double tol = 1e-12)
{
    double x = std::clamp (x0, lo, hi);
    double flo = f (lo);
    const double fhi = f (hi);
    const bool bracketed = detail::signDiffers (flo, fhi);

    // Convergence tolerance is relative to the residual scale of the bracket, so
    // f scaled by 1e-200 converges on the same root as f scaled by 1. An
    // absolute |fx| < tol cutoff would accept x0 immediately in that case.
    const double scaleA = std::isfinite (flo) ? std::abs (flo) : 0.0;
    const double scaleB = std::isfinite (fhi) ? std::abs (fhi) : 0.0;
    const double ftol = tol * std::max (scaleA, scaleB);

    for (int i = 0; i < maxIter; ++i)
    {
        const double fx = f (x);
        if (fx == 0.0 || std::abs (fx) < ftol) return x;

        if (bracketed && std::isfinite (fx))
        {
            if (! detail::signDiffers (fx, flo)) { lo = x; flo = fx; }
            else                                 { hi = x; }
        }

        double d;
        if constexpr (! std::is_same_v<std::decay_t<DF>, std::nullptr_t>)
        {
            d = df (x);
        }
        else
        {
            const double h = 1e-6 * std::max (1.0, std::abs (x));
            d = (f (x + h) - f (x - h)) / (2.0 * h);
        }

        // Bisection is the safe default; a Newton step is taken only when it
        // lands inside the current bracket.
        double xn = 0.5 * (lo + hi);
        if (std::isfinite (d) && d != 0.0)
        {
            const double step = fx / d;
            const double candidate = x - step;
            if (std::isfinite (candidate) && candidate > lo && candidate < hi)
                xn = candidate;
        }
        if (std::abs (xn - x) < tol * std::max (1.0, std::abs (xn))) return xn;
        x = xn;
    }
    return x;
}

// Small dense solve A x = b, N <= 8. Partial pivoting; returns false on
// singularity (caller decides the fallback — usually "use last iterate").
template <int N>
inline bool solveDense (double (&A)[N][N], double (&b)[N], double (&x)[N])
{
    double M[N][N + 1];
    for (int r = 0; r < N; ++r)
    {
        for (int c = 0; c < N; ++c) M[r][c] = A[r][c];
        M[r][N] = b[r];
    }
    for (int col = 0; col < N; ++col)
    {
        int piv = col;
        for (int r = col + 1; r < N; ++r)
            if (std::abs (M[r][col]) > std::abs (M[piv][col])) piv = r;
        if (std::abs (M[piv][col]) < 1e-18) return false;
        if (piv != col) for (int c = col; c <= N; ++c) std::swap (M[piv][c], M[col][c]);
        const double d = M[col][col];
        for (int r = col + 1; r < N; ++r)
        {
            const double m = M[r][col] / d;
            for (int c = col; c <= N; ++c) M[r][c] -= m * M[col][c];
        }
    }
    for (int r = N - 1; r >= 0; --r)
    {
        double s = M[r][N];
        for (int c = r + 1; c < N; ++c) s -= M[r][c] * x[c];
        x[r] = s / M[r][r];
    }
    return true;
}

} // namespace eon
