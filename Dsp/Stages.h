#pragma once
#include <cmath>
#include <algorithm>
#include <numbers>
#include "Rng.h"
#include "Rate.h"

namespace eon {

// ---------------------------------------------------------------------------
// DC blocker: y[n] = x[n] - x[n-1] + R*y[n-1]
//
// R is derived from a cutoff in hertz and the processing rate, so the blocker
// keeps its physical corner under any oversampling factor. The unprepared
// default (R = 0.9997) is retained for callers that never prepare.
// ---------------------------------------------------------------------------
struct DCBlocker
{
    double x1 = 0.0, y1 = 0.0;
    double R = 0.9997;
    void prepare (double sampleRate, double cutoffHz = 18.0) noexcept
    {
        R = std::exp (-2.0 * std::numbers::pi * cutoffHz / sampleRate);
    }
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

    // Prepared per-sample coefficients. The defaults are the authoring-rate
    // constants, so a stage that never calls prepare() behaves as before.
    double pole1 = 0.80, pole2 = 0.86;
    double feed1 = 0.18, feed2 = 0.14, feed3 = 0.05;
    double satScale = 0.55, outScale = 0.82;
    double rateHz = 0.0;
    float cachedFreq = -1.f, cachedQ = -1.f;

    void setTolerance (float t) { tolerance = t; }
    void reset() { lp1 = lp2 = 0.0; }

    // rateHz = processing rate — pass the oversampled rate when this stage runs
    // inside an oversampling loop. Call before reset().
    void prepare (double newRateHz)
    {
        rateHz = newRateHz;
        cachedFreq = -1.f;              // force a refresh on the next sample
        refresh();
    }

    // freqHz/q are block-constant in every caller, so the pole rescale (a pow)
    // is cached rather than evaluated per sample.
    inline void refresh()
    {
        const double wd = std::min (1.0, (double) cachedFreq / 8000.0);
        const double q  = (double) cachedQ;

        const double pole1Orig = 0.80 - wd * 0.06;
        const double pole2Orig = 0.86;
        pole1 = detail::rescalePole (pole1Orig, rateHz);
        pole2 = detail::rescalePole (pole2Orig, rateHz);

        // A one-pole pair is only rate independent when the feed tracks the
        // pole: keeping the time constant but not the DC gain would make the
        // inductor's low-frequency gain scale with the oversampling factor.
        const double dc1 = (0.18 + q * 0.08) / (1.0 - pole1Orig);
        const double dc2 = 0.14 / (1.0 - pole2Orig);
        feed1 = dc1 * (1.0 - pole1);
        feed2 = dc2 * (1.0 - pole2);

        satScale = 0.55 + q * 0.32;
        outScale = 0.82 + q * 0.22;
    }

    inline float process (float x, float freqHz, float q, float drive)
    {
        if (freqHz != cachedFreq || q != cachedQ)
        {
            cachedFreq = freqHz; cachedQ = q;
            refresh();
        }

        const double wd   = std::min (1.0, (double) freqHz / 8000.0);
        const double in   = (double) x * drive * tolerance;
        const double satD = 1.0 + 0.35 / (0.18 + wd);
        const double out  = std::tanh (in * satScale * satD);
        // two leaky integrators give the inductor's stored-energy tail
        lp1 = out * feed1 + lp1 * pole1;
        lp2 = lp1 * feed2 + lp2 * pole2;
        return (float) (out * outScale + lp1 * feed2 + lp2 * feed3);
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

    // Prepared constants (defaults = authoring-rate values).
    double poleLp = 0.90, polePink = 0.985, amplitude = 5e-5;

    void reset() { pink = lp = 0.0; }

    // rateHz = processing rate. Both the shaping poles and the per-sample
    // amplitude are rate-compensated so the noise floor does not move when the
    // oversampling mode changes.
    void prepare (double rateHz)
    {
        poleLp    = detail::rescalePole  (0.90,  rateHz);
        polePink  = detail::rescalePole  (0.985, rateHz);
        amplitude = detail::rescaleNoise (5e-5,  rateHz);
    }

    inline float process()
    {
        const double w = rng.next() * 2.0 - 1.0;
        lp   = lp * poleLp + w * (1.0 - poleLp);
        pink = pink * polePink + lp * (1.0 - polePink);
        return (float) (pink * amplitude);
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

    // Prepared constants (defaults = authoring-rate values).
    double poleSag = 0.9970, poleThermal = 0.99990, poleDrift = 0.99998;
    double driftStep = 4e-6;

    void reset() { sag = thermal = drift = 0.0; }

    // rateHz = processing rate. PSU sag and thermal drift are physical time
    // constants; they must not shorten when the oversampling factor rises.
    void prepare (double rateHz)
    {
        poleSag     = detail::rescalePole (0.9970,  rateHz);
        poleThermal = detail::rescalePole (0.99990, rateHz);
        poleDrift   = detail::rescalePole (0.99998, rateHz);
        // The random walk's per-step variance must also be rate-compensated,
        // otherwise the drift excursion grows as sqrt(rate) between updates.
        driftStep = detail::rescaleNoise (4e-6, rateHz);
    }

    inline void tick (double env)
    {
        sag     = sag * poleSag + env * (1.0 - poleSag);
        thermal = thermal * poleThermal + env * (1.0 - poleThermal);
        // OU: mean-reverting -> bounded, still gives slow random "breathing"
        drift = drift * poleDrift + (rng.next() - 0.5) * driftStep;
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

    // Prepared constants (defaults = authoring-rate values).
    double pole1 = 0.9993, pole2 = 0.99996, pole3 = 0.999998;
    double mean1 = 0.99999, mean2 = 0.99999, mean3 = 0.999999;

    void setTolerance (float t) { tolerance = t; }
    void reset() { cap1 = cap2 = cap3 = cap1m = cap2m = cap3m = 0.0; sagTh.reset(); }

    // rateHz = processing rate. The dielectric-absorption time constants are
    // physical (milliseconds to seconds), so they are rescaled to the actual
    // processing rate instead of being replayed per sample at 8x-32x speed.
    void prepare (double rateHz)
    {
        pole1 = detail::rescalePole (0.9993,   rateHz);
        pole2 = detail::rescalePole (0.99996,  rateHz);
        pole3 = detail::rescalePole (0.999998, rateHz);
        mean1 = detail::rescalePole (0.99999,  rateHz);
        mean2 = detail::rescalePole (0.99999,  rateHz);
        mean3 = detail::rescalePole (0.999999, rateHz);
        sagTh.prepare (rateHz);
    }

    inline float process (float x)
    {
        const double in = (double) x * (0.40 * tolerance);

        cap1 = cap1 * pole1 + in * (1.0 - pole1);
        cap2 = cap2 * pole2 + in * (1.0 - pole2);
        cap3 = cap3 * pole3 + in * (1.0 - pole3);
        cap1m = cap1m * mean1 + cap1 * (1.0 - mean1);
        cap2m = cap2m * mean2 + cap2 * (1.0 - mean2);
        cap3m = cap3m * mean3 + cap3 * (1.0 - mean3);

        // only the AC deviation of each cap couples back in
        const double da = in + (cap1 - cap1m) * 2.5
                             + (cap2 - cap2m) * 0.8
                             + (cap3 - cap3m) * 0.2;

        const double a2  = da * da;
        double out = da + 0.115 * a2 * (da >= 0 ? 1.0 : -0.60) - 0.0035 * a2 * da;

        sagTh.tick (std::abs (out));
        out *= (1.0 - sagTh.sag * 0.045 - sagTh.thermal * 0.01);

        out = 1.5 * std::tanh (out / 1.5);

        return (float) (out + sagTh.drift);
    }
};

} // namespace eon
