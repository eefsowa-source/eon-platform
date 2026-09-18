#pragma once
#include <cmath>
#include <algorithm>
#include "Rng.h"

namespace eon {

// ---------------------------------------------------------------------------
// Jiles-Atherton hysteresis — textbook form.
//
//   He   = H + alpha*M
//   Man  = Ms * (coth(He/a) - a/He)                     (Langevin)
//   dMirr/dH = (Man - Mirr) / (k*delta - alpha*(Man - Mirr)),  delta = sign(dH/dt)
//   M    = Mirr + c*(Man - Mirr)                        (reversible component)
//
// The only slow memory is the hysteresis state itself — no free-running
// "iron memory" accumulator injecting DC wander.
// ---------------------------------------------------------------------------

struct JilesAtherton
{
    float Ms    = 1.10f;  // saturation magnetisation
    float a     = 0.22f;  // domain wall density / knee sharpness
    float alpha = 0.001f; // interdomain coupling
    float k     = 0.30f;  // coercivity (hysteresis width)
    float c     = 0.18f;  // reversibility
    float drive = 1.9f;   // input -> H scaling
    float tolerance = 1.f; // unit-to-unit component variation

    double Mirr = 0.0;
    double M    = 0.0;
    double H    = 0.0;

    Rng rng;

    void setTolerance (float t) { tolerance = t; }
    void reset() { Mirr = M = H = 0.0; }

    inline double langevin (double He) const
    {
        const double q = He / a;
        if (std::abs (q) < 1e-3)
            return Ms * q / 3.0;
        return Ms * (1.0 / std::tanh (q) - 1.0 / q);
    }

    inline float process (float x)
    {
        const double Hn    = (double) x * (drive * tolerance);
        const double dH    = Hn - H;
        const double delta = dH >= 0.0 ? 1.0 : -1.0;

        const double He  = Hn + alpha * M;
        const double Man = langevin (He);
        const double num = Man - Mirr;
        const double den = k * delta - alpha * num;

        double dMirr_dH = 0.0;
        if (std::abs (den) > 1e-9)
            dMirr_dH = num / den;

        Mirr += std::clamp (dMirr_dH * dH, -1.0, 1.0);
        Mirr  = std::clamp (Mirr, -1.3 * (double) Ms, 1.3 * (double) Ms);
        M     = Mirr + c * num;
        H     = Hn;

        const double B = std::tanh (M * 0.92 + Hn * 0.08);
        const double activity = std::min (0.001, std::abs (dMirr_dH) * 1e-5);
        const double bark = (rng.next() - 0.5) * activity;
        return (float) (B + bark);
    }
};

} // namespace eon
