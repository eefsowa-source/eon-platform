#include "Dsp/Adaa.h"
#include "Tests/TestHarness.h"

#include <cmath>

EON_TEST_CASE("nonlinear suite smoke test")
{
    eon::SoftClipSat saturator;
    constexpr int sample_count = 4;
    const float input[sample_count] = {-0.5f, 0.25f, 1.2f, -0.1f};
    float reference[sample_count] = {};

    for (int i = 0; i < sample_count; ++i)
    {
        reference[i] = saturator.process(input[i]);
        EON_CHECK(std::isfinite(reference[i]));
    }

    saturator.reset();
    for (int i = 0; i < sample_count; ++i)
        EON_CHECK_NEAR(saturator.process(input[i]), reference[i], 0.0);
}
