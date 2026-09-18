#pragma once
#include <cmath>
#include <algorithm>
#include <cstddef>
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
    // initial guess: asymptotic L - ln L for large logs, exp for small
    double w = logX > 1.0 ? logX - std::log (std::max (logX, 1e-10)) : std::exp (logX);
    if (! (w > 0.0)) w = 1e-10;
    for (int i = 0; i < 40; ++i)
    {
        const double f  = w + std::log (w) - logX;
        const double dw = f / (1.0 + 1.0 / w);      // Newton on g(w)=w+ln w
        w -= dw;
        if (std::abs (dw) < 1e-14 * std::max (1.0, std::abs (w))) break;
    }
    return w;
}

inline double lambertW0 (double x)
{
    if (x <= 0.0) return 0.0;
    return lambertW0log (std::log (x));
}

// Scalar Newton on f(x) = 0 inside [lo, hi]. f must be evaluable; df may be
// nullptr -> secant estimate. Keeps a bracket once a sign change is seen and
// bisects whenever Newton leaves it — guaranteed progress on sane inputs.
template <typename F, typename DF>
inline double newtonScalar (F&& f, DF&& df, double x0, double lo, double hi,
                            int maxIter = 20, double tol = 1e-12)
{
    double x = std::clamp (x0, lo, hi);
    double flo = f (lo), fhi = f (hi);
    bool bracketed = (flo * fhi < 0.0);

    for (int i = 0; i < maxIter; ++i)
    {
        const double fx = f (x);
        if (std::abs (fx) < tol) return x;

        if (bracketed)
        {
            if (fx * flo > 0.0) { lo = x; flo = fx; }
            else                { hi = x; fhi = fx; }
        }

        double d;
        if constexpr (! std::is_same_v<std::decay_t<DF>, std::nullptr_t>)
        {
            d = df (x);
            if (! std::isfinite (d) || std::abs (d) < 1e-20) d = 0.0;
        }
        else
        {
            const double h = 1e-6 * std::max (1.0, std::abs (x));
            d = (f (x + h) - f (x - h)) / (2.0 * h);
        }

        double xn = (d != 0.0) ? x - fx / d : x;
        if (! std::isfinite (xn) || xn <= lo || xn >= hi || (bracketed && (xn <= lo || xn >= hi)))
            xn = 0.5 * (lo + hi);                       // bisection fallback
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
