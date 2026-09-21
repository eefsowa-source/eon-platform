#pragma once
#include <cmath>
#include <algorithm>
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
//   i(v) = 2*Is*sinh(v/Vt) = (a - v)/R   -> Newton on v, then b = 2v - a.
struct WdfDiodePair : WdfPort
{
    double Is = 2.52e-9;
    double Vt = 25.85e-3;
    int    iterations = 8;
    double aIn = 0.0;
    mutable double vPrev = 0.0;              // warm start for the Newton loop
    void   incident (double a) override { aIn = a; }
    void   reset() override { aIn = vPrev = 0.0; }
    double emitted() const override
    {
        if (aIn == 0.0)
        {
            vPrev = 0.0;
            return 0.0;
        }

        const double bound = std::abs (aIn);
        double lo = -bound, hi = bound;
        constexpr double maxSinhArgument = 700.0;
        const auto residual = [this, maxSinhArgument](double v)
        {
            const double x = std::clamp (v / Vt, -maxSinhArgument, maxSinhArgument);
            return 2.0 * Is * std::sinh (x) - (aIn - v) / R;
        };
        const auto derivative = [this, maxSinhArgument](double v)
        {
            const double x = std::clamp (v / Vt, -maxSinhArgument, maxSinhArgument);
            return (2.0 * Is / Vt) * std::cosh (x) + 1.0 / R;
        };

        double fLo = residual (lo);
        double v = std::isfinite (vPrev) ? std::clamp (vPrev, lo, hi) : 0.0;
        if (v <= lo || v >= hi) v = 0.0;
        const double tolerance = 1e-12 * std::max (std::abs (aIn / R), 1e-12);

        for (int i = 0; i < std::max (iterations, 1); ++i)
        {
            const double f = residual (v);
            if (std::abs (f) <= tolerance) break;

            if ((f < 0.0) == (fLo < 0.0)) { lo = v; fLo = f; }
            else                           hi = v;

            const double df = derivative (v);
            double next = std::isfinite (df) && df > 0.0 ? v - f / df : v;
            const double scaledVoltage = v / Vt;
            if (std::abs (scaledVoltage) >= maxSinhArgument
                || ! std::isfinite (next) || next <= lo || next >= hi)
                next = 0.5 * lo + 0.5 * hi;
            if (next == v) break;
            v = next;
        }

        // Newton is fast near the root, while bisection guarantees progress
        // when a large drive starts outside the unclamped sinh range.
        for (int i = 0; i < 64; ++i)
        {
            const double f = residual (v);
            if (std::abs (f) <= tolerance) break;
            if ((f < 0.0) == (fLo < 0.0)) { lo = v; fLo = f; }
            else                           hi = v;
            const double next = 0.5 * lo + 0.5 * hi;
            if (next == v) break;
            v = next;
        }

        v = std::clamp (v, -bound, bound);
        if (! std::isfinite (v)) v = 0.0;
        vPrev = v;
        return v - (aIn - v);
    }
};

} // namespace eon
