// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>

namespace eon {
namespace detail {

// ---------------------------------------------------------------------------
// Rate independence for per-sample one-pole coefficients.
//
// The poles in the emulation stages were authored for a 48 kHz processing rate
// — the rate the voicing was tuned at. Running the same per-sample constants at
// an oversampled rate (384 kHz … 1.536 MHz) moves every pole up by exactly that
// factor, so the analog "memory" (inductor tail, dielectric absorption, PSU sag,
// noise shaping) audibly changes when the oversampling mode changes. A hardware
// emulation's time constants must not depend on the host's internal rate.
//
// rescalePole() keeps the physical time constant fixed:
//
//     a(rate) = a(48k) ^ (48k / rate)
//
// Pass the *processing* rate — the oversampled rate when the stage runs inside
// an oversampling loop. A stage that never calls prepare() keeps the raw
// authoring-rate constants, so existing callers and tests are unaffected.
// ---------------------------------------------------------------------------
inline constexpr double kAuthoringRate = 48000.0;

inline double rescalePole (double poleAtAuthoringRate, double rate) noexcept
{
    if (! (poleAtAuthoringRate > 0.0) || poleAtAuthoringRate >= 1.0 || ! (rate > 0.0))
        return poleAtAuthoringRate;
    return std::pow (poleAtAuthoringRate, kAuthoringRate / rate);
}

// A per-sample white-noise source has PSD ~ sigma^2 / rate, so holding the
// in-band noise power (and the random-walk variance) fixed across rates
// requires sigma ~ sqrt(rate).
inline double rescaleNoise (double sigmaAtAuthoringRate, double rate) noexcept
{
    if (! (rate > 0.0))
        return sigmaAtAuthoringRate;
    return sigmaAtAuthoringRate * std::sqrt (rate / kAuthoringRate);
}

} // namespace detail
} // namespace eon
