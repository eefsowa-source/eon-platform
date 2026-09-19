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

EON_TEST_CASE("ADAA2 handles repeated outer nodes")
{
    eon::SoftClipSat saturator;
    saturator.process(0.5f);
    saturator.process(-0.5f);
    const float actual = saturator.process(0.5f);

    constexpr double a = 0.5;
    constexpr double b = -0.5;
    const double difference = a - b;
    const double expected = 2.0 * (eon::SoftClip::F1(a) * difference
                                 - eon::SoftClip::F2(a)
                                 + eon::SoftClip::F2(b))
                          / (difference * difference);
    EON_CHECK_NEAR(actual, expected, 1e-6);
}

EON_TEST_CASE("ADAA2 remains finite for near-equal outer nodes")
{
    eon::SoftClipSat saturator;
    saturator.process(0.5f);
    saturator.process(-0.5f);
    const float near_outer = std::nextafter(0.5f, 1.0f);

    EON_CHECK(std::isfinite(saturator.process(near_outer)));
}

EON_TEST_CASE("ADAA2 reaches the constant-input limit")
{
    eon::SoftClipSat saturator;
    saturator.process(0.5f);
    saturator.process(0.5f);
    const float actual = saturator.process(0.5f);

    EON_CHECK_NEAR(actual, eon::SoftClip::f(0.5), 1e-6);
}
