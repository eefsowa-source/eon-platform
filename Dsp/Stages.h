#pragma once
#include <cmath>
#include <algorithm>
#include "Rng.h"

namespace eon {

// ---------------------------------------------------------------------------
// DC blocker: y[n] = x[n] - x[n-1] + R*y[n-1]
// ---------------------------------------------------------------------------
struct DCBlocker
{
    double x1 = 0.0, y1 = 0.0;
    double R = 0.9997;
    void reset() { x1 = y1 = 0.0; }
    inline float process (float x)
    {
        const double y = (double) x - x1 + R * y1;
        x1 = x; y1 = y;
        return (float) y;
    }
};

// ---------------------------------------------------------------------------
// Inductor resonator (honest version of the old "WDFInductorGod"):
// tanh saturation whose depth follows frequency + a resonant lowpass tail.
// Not a wave-digital model — just a saturable inductor with energy storage.
// ---------------------------------------------------------------------------
struct InductorResonator
{
    double lp1 = 0.0, lp2 = 0.0;
    float tolerance = 1.f;
    void setTolerance (float t) { tolerance = t; }
    void reset() { lp1 = lp2 = 0.0; }

    inline float process (float x, float freqHz, float q, float drive)
    {
        const double wd   = std::min (1.0, freqHz / 8000.0);
        const double in   = (double) x * drive * tolerance;
        const double satD = 1.0 + 0.35 / (0.18 + wd);
        const double out  = std::tanh (in * (0.55 + q * 0.32) * satD);
        // two leaky integrators give the inductor's stored-energy tail
        lp1 = out * (0.18 + q * 0.08) + lp1 * (0.80 - wd * 0.06);
        lp2 = lp1 * 0.14 + lp2 * 0.86;
        return (float) (out * (0.82 + q * 0.22) + lp1 * 0.14 + lp2 * 0.05);
    }
};

// ---------------------------------------------------------------------------
// Analog "air" noise floor — pink-ish noise at ~-90 dBFS.
// Per-instance RNG (the original shared getSystemRandom per-sample, per-stage,
// which is unnecessary contention).
// ---------------------------------------------------------------------------
struct AnalogAir
{
    Rng rng;
    double pink = 0.0, lp = 0.0;
    void reset() { pink = lp = 0.0; }
    inline float process()
    {
        const double w = rng.next() * 2.0 - 1.0;
        lp   = lp * 0.90 + w * 0.10;
        pink = pink * 0.985 + lp * 0.015;
        return (float) (pink * 5e-5);
    }
};

// ---------------------------------------------------------------------------
// Power supply sag + thermal drift — BOUNDED.
// The original's `drift += rand*2e-7` was an unbounded random walk that would
// wander the output DC level forever. Here drift is an Ornstein-Uhlenbeck
// process (mean-reverting) and sag/thermal are leaky envelope followers.
// ---------------------------------------------------------------------------
struct SagThermal
{
    double sag = 0.0, thermal = 0.0, drift = 0.0;
    Rng rng;
    void reset() { sag = thermal = drift = 0.0; }

    inline void tick (double env)
    {
        sag     = sag * 0.9970  + env * 0.0030;
        thermal = thermal * 0.99990 + env * 0.00010;
        // OU: mean-reverting -> bounded, still gives slow random "breathing"
        drift = drift * 0.99998 + (rng.next() - 0.5) * 4e-6;
        drift = std::clamp (drift, -0.001, 0.001);
    }
};

// ---------------------------------------------------------------------------
// Avalon-style Class-A path: asymmetric even-harmonic generator +
// dielectric-absorption memory that is AC-COUPLED (original summed raw
// integrator outputs into the audio -> guaranteed DC buildup).
// Here each DA cap is servoed to its own mean so only the *variation* is added.
// ---------------------------------------------------------------------------
struct ClassAStage
{
    double cap1 = 0.0, cap2 = 0.0, cap3 = 0.0;
    double cap1m = 0.0, cap2m = 0.0, cap3m = 0.0; // running means (DC servo)
    float tolerance = 1.f;
    SagThermal sagTh;
    void setTolerance (float t) { tolerance = t; }
    void reset() { cap1 = cap2 = cap3 = cap1m = cap2m = cap3m = 0.0; sagTh.reset(); }

    inline float process (float x)
    {
        const double in = (double) x * (0.40 * tolerance);

        cap1 = cap1 * 0.9993  + in * 0.0007;
        cap2 = cap2 * 0.99996 + in * 0.00004;
        cap3 = cap3 * 0.999998 + in * 0.000002;
        cap1m = cap1m * 0.99999 + cap1 * 0.00001;
        cap2m = cap2m * 0.99999 + cap2 * 0.00001;
        cap3m = cap3m * 0.999999 + cap3 * 0.000001;

        // only the AC deviation of each cap couples back in
        const double da = in + (cap1 - cap1m) * 2.5
                             + (cap2 - cap2m) * 0.8
                             + (cap3 - cap3m) * 0.2;

        const double a2  = da * da;
        double out = da + 0.115 * a2 * (da >= 0 ? 1.0 : -0.60) - 0.0035 * a2 * da;

        sagTh.tick (std::abs (out));
        out *= (1.0 - sagTh.sag * 0.045 - sagTh.thermal * 0.01);

        if (std::abs (out) > 1.5)
            out = 1.5 * std::tanh (out / 1.5);

        return (float) (out + sagTh.drift);
    }
};

} // namespace eon
