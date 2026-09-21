#pragma once
#include <cmath>
#include <algorithm>
#include <limits>
#include "Solvers.h"

namespace eon {

// ---------------------------------------------------------------------------
// Wave Digital Filter core (Fettweis / Yeh formulation).
//
// Wave variables per port:  a = v + R*i (incident),  b = v - R*i (reflected).
// Elements are leaves; series/parallel adaptors compose them into a tree with
// ONE root — the element that terminates the delay-free loop.
//
// Rules that keep the math correct (both verified against references):
//  * Adaptors pick their parent port resistance so they are reflection-free:
//      series   R0 = R1 + R2 + ...
//      parallel G0 = G1 + G2 + ...        (R0 = parallel combination)
//  * An IDEAL voltage source works only as the ROOT (b = 2Vs - a preserves
//    the feedback path — a resistive source at the root collapses the pole
//    to backward-Euler, measured pole 0.955 vs bilinear 0.910).
//  * As a LEAF the source must not depend on its incident wave, so use the
//    resistive form (emitted = Vs, internal Rs folds into its port R).
//  * Nonlinear elements (diode, diode pair) sit at the root and solve
//    v/i in the wave domain — no iteration inside the tree.
//
// Voltage read-out at any element: v = 0.5*(a + b) with the element's own
// sign convention — for a cap across the output, read 0.5*(incident+emitted).
// ---------------------------------------------------------------------------

struct WdfPort
{
    double R = 1.0;                               // port resistance [ohm]
    virtual ~WdfPort() = default;
    virtual double emitted() const = 0;           // b: wave toward parent
    virtual void   incident (double a) = 0;       // a: wave received from parent
    virtual void   reset() {}
};

// ------------------------------- linear leaves ------------------------------

struct WdfResistor : WdfPort
{
    explicit WdfResistor (double r) { R = r; }
    double emitted() const override { return 0.0; }
    void   incident (double) override {}
};

struct WdfCapacitor : WdfPort
{
    double C, state = 0.0, voltageState = 0.0;
    explicit WdfCapacitor (double c) : C (c) {}
    void setSampleRate (double fs) { R = 1.0 / (2.0 * C * fs); }
    void reset() override { state = voltageState = 0.0; }
    double emitted() const override { return state; }      // b[n] = a[n-1]
    void   incident (double a) override
    {
        voltageState = 0.5 * (a + state);
        state = a;
    }
    double voltage() const { return voltageState; }
};

struct WdfInductor : WdfPort
{
    double L, state = 0.0;
    explicit WdfInductor (double l) : L (l) {}
    void setSampleRate (double fs) { R = 2.0 * L * fs; }
    void reset() override { state = 0.0; }
    double emitted() const override { return -state; }     // b[n] = -a[n-1]
    void   incident (double a) override { state = a; }
};

// Resistive voltage source — LEAF ONLY. emitted wave ignores the incident,
// so it can never close a loop. Rs sets its port resistance.
struct WdfVSourceRes : WdfPort
{
    double Vs = 0.0;
    explicit WdfVSourceRes (double rs) { R = rs; }
    double emitted() const override { return Vs; }
    void   incident (double) override {}
};

// Ideal voltage source — ROOT ONLY (b = 2Vs - a). Using it as a leaf creates
// an algebraic loop; using a resistive source at the root breaks the pole.
struct WdfVSourceIdeal : WdfPort
{
    double Vs = 0.0;
    double aIn = 0.0;
    double emitted() const override { return 2.0 * Vs - aIn; }
    void   incident (double a) override { aIn = a; }
};

// Resistive current source — leaf. i = Is fixed -> b = a - 2R*Is.
struct WdfISourceRes : WdfPort
{
    double Is = 0.0;
    double aIn = 0.0;
    explicit WdfISourceRes (double rs) { R = rs; }
    double emitted() const override { return aIn - 2.0 * R * Is; }
    void   incident (double a) override { aIn = a; }
    void   reset() override { aIn = 0.0; }
};

// ------------------------------- adaptors -----------------------------------

// Series junction with two children. Parent port R = c1.R + c2.R (adapted).
struct WdfSeries2 : WdfPort
{
    WdfPort& c1; WdfPort& c2;
    WdfSeries2 (WdfPort& a, WdfPort& b) : c1 (a), c2 (b) { R = c1.R + c2.R; }
    double emitted() const override { return -(c1.emitted() + c2.emitted()); }
    void incident (double a0) override
    {
        const double a1 = c1.emitted(), a2 = c2.emitted();
        const double s  = a0 + a1 + a2;
        c1.incident (a1 - (c1.R / R) * s);
        c2.incident (a2 - (c2.R / R) * s);
    }
    void reset() override { c1.reset(); c2.reset(); }
};

struct WdfSeries3 : WdfPort
{
    WdfPort& c1; WdfPort& c2; WdfPort& c3;
    WdfSeries3 (WdfPort& a, WdfPort& b, WdfPort& c) : c1 (a), c2 (b), c3 (c)
    {
        R = c1.R + c2.R + c3.R;
    }
    double emitted() const override
    {
        return -(c1.emitted() + c2.emitted() + c3.emitted());
    }
    void incident (double a0) override
    {
        const double a1 = c1.emitted(), a2 = c2.emitted(), a3 = c3.emitted();
        const double s  = a0 + a1 + a2 + a3;
        c1.incident (a1 - (c1.R / R) * s);
        c2.incident (a2 - (c2.R / R) * s);
        c3.incident (a3 - (c3.R / R) * s);
    }
    void reset() override { c1.reset(); c2.reset(); c3.reset(); }
};

// Parallel junction with two children. Parent port R = 1/(1/c1.R + 1/c2.R).
struct WdfParallel2 : WdfPort
{
    WdfPort& c1; WdfPort& c2;
    double G0;
    WdfParallel2 (WdfPort& a, WdfPort& b) : c1 (a), c2 (b)
    {
        const double G = 1.0 / c1.R + 1.0 / c2.R;
        R = 1.0 / G;  G0 = G;
    }
    double emitted() const override                       // b0 = sum(Gc/G0)*ac
    {
        return (c1.emitted() / c1.R + c2.emitted() / c2.R) / G0;
    }
    void incident (double a0) override
    {
        const double a1 = c1.emitted(), a2 = c2.emitted();
        const double v  = (a0 / R + a1 / c1.R + a2 / c2.R) / (1.0 / R + G0);
        c1.incident (2.0 * v - a1);
        c2.incident (2.0 * v - a2);
    }
    void reset() override { c1.reset(); c2.reset(); }
};

// --------------------------- nonlinear root devices -------------------------

// Single diode (Shockley) as tree ROOT. Closed form via Lambert W:
//   b = a + 2R*Is - 2Vt*W( (R*Is/Vt) * exp((R*Is + a)/Vt) )
// evaluated in log space — the exp() argument overflows double for a > ~0.7V.
struct WdfDiode : WdfPort
{
    double Is = 2.52e-9;      // 1N914-ish
    double Vt = 25.85e-3;
    double aIn = 0.0;
    void   incident (double a) override { aIn = a; }
    double emitted() const override
    {
        const double logZ = std::log (R * Is / Vt) + (R * Is + aIn) / Vt;
        return aIn + 2.0 * R * Is - 2.0 * Vt * lambertW0log (logZ);
    }
};

// Antiparallel diode pair (clipper staple) as root.
//   i(v) = 2*Is*sinh(v/Vt) = (a - v)/R   -> safeguarded solve, then b = 2v - a.
struct WdfDiodePair : WdfPort
{
    double Is = 2.52e-9;
    double Vt = 25.85e-3;
    int    iterations = 8;
    double aIn = 0.0;
    mutable double vPrev = 0.0;              // warm start for the Newton loop
    mutable bool solveSucceeded = false;
    void   incident (double a) override { aIn = a; solveSucceeded = false; }
    void   reset() override { aIn = vPrev = 0.0; solveSucceeded = false; }
    double emitted() const override
    {
        solveSucceeded = false;
        if (! std::isfinite (aIn))
        {
            vPrev = 0.0;
            return 0.0;
        }
        constexpr double logMaxDouble = 709.782712893383973096;
        constexpr double logMinSubnormal = -744.440071921381218089;
        const double maxDouble = std::numeric_limits<double>::max();
        if (! std::isfinite (R) || ! (R > 0.0)
            || ! std::isfinite (Is) || ! (Is > 0.0)
            || ! std::isfinite (Vt) || ! (Vt > 0.0))
        {
            vPrev = 0.0;
            return -aIn;
        }

        if (iterations <= 0)
        {
            if (! std::isfinite (vPrev)) vPrev = 0.0;
            return vPrev - (aIn - vPrev);
        }

        if (aIn == 0.0)
        {
            vPrev = 0.0;
            solveSucceeded = true;
            return 0.0;
        }

        const double absA = std::abs (aIn);
        const double logA = std::log (absA);
        const double logR = std::log (R);
        const double logIs = std::log (Is);
        // Keep both Is*exp(q) and R*Is*exp(q) representable inside the bracket.
        const double qMaxCurrent = logMaxDouble - logIs;
        const double qMaxScaledCurrent = logMaxDouble - logR - logIs;
        const double qMax = std::min (qMaxCurrent, qMaxScaledCurrent);
        if (! (qMax > 0.0) || ! std::isfinite (qMax))
        {
            vPrev = 0.0;
            return -aIn;
        }

        // Four log-current e-folds above the load current give the upper end
        // a genuine positive residual when representable.
        const double qNeeded = logA - logR - logIs + 4.0;
        const double qBound = std::min (qMax, std::max (1.0, qNeeded));
        const double rawVoltageBound = Vt * qBound;
        const double voltageBound = std::isfinite (rawVoltageBound)
                                  ? rawVoltageBound : maxDouble;
        const double bound = std::min (absA, voltageBound);
        if (! (bound > 0.0) || ! std::isfinite (bound))
        {
            vPrev = 0.0;
            return -aIn;
        }

        const auto scaledDiodeCurrent = [logR, logIs, maxDouble](double q)
        {
            const double absQ = std::abs (q);
            double logMagnitude;
            double sign;
            if (absQ < 20.0)
            {
                const double factor = 2.0 * std::sinh (q);
                if (factor == 0.0) return 0.0;
                logMagnitude = logR + logIs + std::log (std::abs (factor));
                sign = factor;
            }
            else
            {
                // For large |q|, 2*sinh(q) is sign(q)*exp(|q|) to double precision.
                logMagnitude = logR + logIs + absQ;
                sign = q;
            }

            if (logMagnitude < logMinSubnormal) return std::copysign (0.0, sign);
            if (logMagnitude >= logMaxDouble) return std::copysign (maxDouble, sign);
            return std::copysign (std::exp (logMagnitude), sign);
        };
        const auto scaledDerivative = [this, logR, logIs](double q)
        {
            const double absQ = std::abs (q);
            double logMagnitude;
            if (absQ < 20.0)
            {
                const double factor = 2.0 * std::cosh (q);
                logMagnitude = logR + logIs + std::log (factor) - std::log (Vt);
            }
            else
            {
                logMagnitude = logR + logIs + absQ - std::log (Vt);
            }

            if (logMagnitude < logMinSubnormal) return 0.0;
            if (logMagnitude >= logMaxDouble) return std::numeric_limits<double>::infinity();
            return std::exp (logMagnitude);
        };
        const auto residual = [this, &scaledDiodeCurrent](double v)
        {
            return (v - aIn) + scaledDiodeCurrent (v / Vt);
        };
        const auto derivative = [this, &scaledDerivative](double v)
        {
            return 1.0 + scaledDerivative (v / Vt);
        };

        double lo = -bound, hi = bound;
        double fLo = residual (lo);
        const double fHi = residual (hi);
        if (! (fLo <= 0.0 && fHi >= 0.0))
        {
            vPrev = 0.0;
            return -aIn;
        }

        double initialVoltage;
        if (absA <= Vt)
        {
            const double logSlope = std::log (2.0) + logR + logIs - std::log (Vt);
            const double slope = logSlope >= logMaxDouble
                               ? std::numeric_limits<double>::infinity()
                               : (logSlope < logMinSubnormal ? 0.0 : std::exp (logSlope));
            initialVoltage = aIn / (1.0 + slope);
        }
        else
        {
            const double logRatio = logA - std::log (2.0) - logR - logIs;
            const double qGuess = logRatio > 20.0 ? logRatio + std::log (2.0)
                                : (logRatio < -20.0 ? std::exp (logRatio)
                                                   : std::asinh (std::exp (logRatio)));
            const double rawGuess = Vt * qGuess;
            initialVoltage = std::copysign (std::isfinite (rawGuess)
                                           ? std::min (bound, rawGuess) : bound, aIn);
        }
        initialVoltage = std::clamp (initialVoltage, lo, hi);
        double v = initialVoltage;
        if (std::isfinite (vPrev) && vPrev != 0.0
            && std::signbit (vPrev) == std::signbit (aIn)
            && vPrev >= lo && vPrev <= hi)
        {
            const double candidateResidual = residual (initialVoltage);
            const double warmResidual = residual (vPrev);
            if (std::isfinite (warmResidual)
                && (! std::isfinite (candidateResidual)
                    || std::abs (warmResidual) < std::abs (candidateResidual)))
                v = vPrev;
        }

        const double residualScale = std::max (absA, std::abs (v));
        const double residualTolerance = std::max (
            64.0 * std::numeric_limits<double>::epsilon() * residualScale,
            std::numeric_limits<double>::denorm_min());

        for (int i = 0; i < iterations; ++i)
        {
            const double f = residual (v);
            if (std::isfinite (f) && std::abs (f) <= residualTolerance)
            {
                solveSucceeded = true;
                break;
            }
            if (! std::isfinite (f))
            {
                if (f < 0.0) { lo = v; fLo = f; }
                else          { hi = v; }
                v = 0.5 * lo + 0.5 * hi;
                continue;
            }

            if ((f < 0.0) == (fLo < 0.0)) { lo = v; fLo = f; }
            else                           { hi = v; }

            const double df = derivative (v);
            double next = std::isfinite (df) && df > 0.0 ? v - f / df : v;
            if (! std::isfinite (next) || next <= lo || next >= hi)
                next = 0.5 * lo + 0.5 * hi;
            const double voltageScale = std::max ({std::abs (v), std::abs (next),
                                                   std::abs (initialVoltage)});
            const double voltageTolerance = std::max (
                8.0 * std::numeric_limits<double>::epsilon() * voltageScale,
                std::numeric_limits<double>::denorm_min());
            if (std::abs (next - v) <= voltageTolerance)
            {
                v = next;
                solveSucceeded = true;
                break;
            }
            v = next;
        }

        v = std::clamp (v, -bound, bound);
        if (! std::isfinite (v)) v = 0.0;
        vPrev = v;
        return v - (aIn - v);
    }
};

} // namespace eon
