// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <cstdint>
#if defined (__SSE__)
  #include <immintrin.h>
#endif

namespace eon {

// ---------------------------------------------------------------------------
// Real-time hygiene helpers.
//
// ScopedDenormalsOff — flush denormals to zero for the scope. Denormals cost
// ~10-100x normal FP cost and appear exactly when the audio decays (reverb
// tails, IIR states), i.e. the worst possible time. Plugin processors should
// put one of these at the top of processBlock.
// ---------------------------------------------------------------------------

struct ScopedDenormalsOff
{
#if defined (__aarch64__) || defined (__arm64__) || defined (__arm64e__)
    uint64_t oldFpcr;
    ScopedDenormalsOff()
    {
        __asm__ volatile ("mrs %0, fpcr" : "=r" (oldFpcr));
        __asm__ volatile ("msr fpcr, %0" :: "r" (oldFpcr | (1ull << 24))); // FZ
    }
    ~ScopedDenormalsOff()
    {
        __asm__ volatile ("msr fpcr, %0" :: "r" (oldFpcr));
    }
#elif defined (__SSE__)
    unsigned int oldMxcsr;
    ScopedDenormalsOff()
    {
        oldMxcsr = _mm_getcsr();
        _mm_setcsr (oldMxcsr | 0x8040);     // FTZ (bit15) + DAZ (bit6)
    }
    ~ScopedDenormalsOff() { _mm_setcsr (oldMxcsr); }
#else
    ScopedDenormalsOff() {}
#endif
    ScopedDenormalsOff (const ScopedDenormalsOff&) = delete;
    ScopedDenormalsOff& operator= (const ScopedDenormalsOff&) = delete;
};

// Cheap flush-to-zero for a single value inside a feedback path.
inline float ftz (float x)
{
    constexpr float tiny = 1e-20f;          // ~-400 dB, above FLT_MIN normal
    return std::abs (x) < tiny ? 0.0f : x;
}

} // namespace eon
