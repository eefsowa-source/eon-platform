#pragma once
#include <cmath>
#include <algorithm>

namespace eon {

// ---------------------------------------------------------------------------
// Zero-delay-feedback (TPT) filters — Zavalishin's formulation.
// These discretise analog filters without the direct-form delay-free loop
// problem: the integrator carries the state, the instantaneous feedback is
// solved exactly. Correct cutoff tracking under per-sample modulation and
// stable saturation inside the loop — the standard for synth/VA filters.
//
// All coefficients use the trapezoidal prewarp: g = tan(pi * fc / fs).
// ---------------------------------------------------------------------------

// Single TPT one-pole lowpass. v = G*(x - s); y = v + s; s += 2v.
struct OnePoleTPT
{
    double g = 0.0, s = 0.0;
    void setCutoff (double fc, double fs)
    {
        g = std::tan (M_PI * std::min (fc, 0.49 * fs) / fs);
    }
    void reset() { s = 0.0; }
    inline double process (double x)
    {
        const double G = g / (1.0 + g);
        const double v = (x - s) * G;
        const double y = v + s;
        s = y + v;
        return y;
    }
};

// TPT state-variable filter: simultaneous lp / bp / hp / notch.
// R = 1/(2Q) damping. Q = 0.5 is butterworth-ish.
struct SvfTPT
{
    double g = 0.0, R = 1.0, s1 = 0.0, s2 = 0.0;
    void setParams (double fc, double q, double fs)
    {
        g = std::tan (M_PI * std::min (fc, 0.49 * fs) / fs);
        R = 1.0 / (2.0 * std::max (0.05, q));
    }
    void reset() { s1 = s2 = 0.0; }
    struct Out { double lp, bp, hp; };
    inline Out process (double x)
    {
        const double den = 1.0 + 2.0 * R * g + g * g;
        const double hp = (x - (2.0 * R + g) * s1 - s2) / den;
        const double bp = g * hp + s1;
        const double lp = g * bp + s2;
        s1 = 2.0 * bp - s1;
        s2 = 2.0 * lp - s2;
        return { lp, bp, hp };
    }
};

// 4-pole resonant ladder (Moog-style) with ZDF feedback.
//   y4 = [G^4 * u0 + (1-G) * S] / (1 + k * G^4),  S = G^3 s0 + G^2 s1 + G s2 + s3
// solved exactly for the instantaneous feedback, so resonance stays accurate
// at any cutoff — no delay-compensation hacks.
struct Ladder4
{
    double g = 0.0, k = 0.0;      // k = resonance 0..4
    double s[4] = {0, 0, 0, 0};
    double drive = 1.0;           // pre-gain into the saturating input
    bool   saturate = true;       // tanh on the feedback input

    void setCutoff (double fc, double fs)
    {
        g = std::tan (M_PI * std::min (fc, 0.49 * fs) / fs);
    }
    void setResonance (double r) { k = std::clamp (r, 0.0, 4.0); }
    void reset() { s[0] = s[1] = s[2] = s[3] = 0.0; }

    inline double process (double x)
    {
        const double G  = g / (1.0 + g);
        const double G4 = G * G * G * G;
        const double S  = G * G * G * s[0] + G * G * s[1] + G * s[2] + s[3];

        double u0 = (x * drive - k * (1.0 - G) * S) / (1.0 + k * G4);
        if (saturate) u0 = std::tanh (u0);

        double u = u0;
        for (int i = 0; i < 4; ++i)
        {
            const double v = (u - s[i]) * G;
            const double y = v + s[i];
            s[i] = y + v;
            u = y;
        }
        return u;
    }
};

} // namespace eon
