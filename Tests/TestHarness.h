#pragma once

#include <cmath>
#include <cstdio>

namespace eon::test {

inline int failures = 0;

inline void check(bool condition, const char* message)
{
    std::printf("%s  %s\n", condition ? "ok  " : "FAIL", message);
    if (!condition)
        ++failures;
}

inline bool near(double actual, double expected, double tolerance)
{
    return std::isfinite(actual) && std::abs(actual - expected) <= tolerance;
}

} // namespace eon::test
