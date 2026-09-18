#pragma once
#include <cstdint>

namespace eon {

// Tiny deterministic PRNG (splitmix64 -> xorshift128+). Keeps the DSP headers
// free of any framework dependency and makes noise sources reproducible
// per instance — no global RNG contention, no thread surprises.
struct Rng
{
    uint64_t s0 = 0x9E3779B97F4A7C15ull, s1 = 0xBF58476D1CE4E5B9ull;

    explicit Rng (uint64_t seed = 0x243F6A8885A308D3ull)
    {
        // splitmix64 to fill the state
        for (int i = 0; i < 2; ++i)
        {
            seed += 0x9E3779B97F4A7C15ull;
            uint64_t z = seed;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
            z ^= z >> 31;
            (i == 0 ? s0 : s1) = z;
        }
    }

    inline uint64_t nextU64()
    {
        uint64_t x = s0, y = s1;
        s0 = y;
        x ^= x << 23;
        s1 = x ^ y ^ (x >> 17) ^ (y >> 26);
        return s1 + y;
    }

    // uniform in [0, 1)
    inline double next()
    {
        return (nextU64() >> 11) * (1.0 / 9007199254740992.0);
    }
};

} // namespace eon
