#include "Dsp/Adaa.h"
#include "Tests/TestHarness.h"

#include <cmath>
#include <limits>

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

EON_TEST_CASE("ADAA2 is stable when the middle node nearly repeats the outer nodes")
{
    eon::SoftClipSat saturator;
    saturator.process(1.0f);
    saturator.process(std::nextafter(1.0f, 0.0f));
    const float actual = saturator.process(1.0f);

    EON_CHECK_NEAR(actual, 2.0 / 3.0, 1e-6);
}

EON_TEST_CASE("ADAA2 remains accurate across repeated-middle ULP distances")
{
    struct SweepCase
    {
        float outer;
        int ulps;
    };
    constexpr SweepCase cases[] = {
        {1.0f, 10}, {0.5f, 10}, {0.005f, 10},
        {1.0f, 1000}, {0.5f, 1000}, {0.005f, 1000},
        {1.0f, 5000}, {1.0f, 5791}, {1.0f, 5793}, {1.0f, 10000}
    };

    for (const auto& sample : cases)
    {
        float middle = sample.outer;
        for (int i = 0; i < sample.ulps; ++i)
            middle = std::nextafter(middle, 0.0f);

        eon::SoftClipSat saturator;
        saturator.process(sample.outer);
        saturator.process(middle);
        const float actual = saturator.process(sample.outer);
        const double center = (2.0 * sample.outer + middle) / 3.0;
        EON_CHECK_NEAR(actual, eon::SoftClip::f(center), 1e-6);
    }
}

EON_TEST_CASE("ADAA2 remains accurate through the middle-node switch")
{
    constexpr float outers[] = {0.005f, 0.5f, 1.0f};
    constexpr double offsets[] = {0.9, 1.1};
    const double switch_scale = std::sqrt(std::numeric_limits<float>::epsilon());

    for (const float outer : outers)
    {
        const double scale = std::fmax(1.0, std::abs(static_cast<double>(outer)));
        const double switch_distance = switch_scale * scale;
        for (const double offset : offsets)
        {
            const float middle = static_cast<float>(outer - offset * switch_distance);
            eon::SoftClipSat saturator;
            saturator.process(outer);
            saturator.process(middle);
            const float actual = saturator.process(outer);
            const double center = (2.0 * outer + middle) / 3.0;
            EON_CHECK_NEAR(actual, eon::SoftClip::f(center), 1e-6);
        }
    }
}

EON_TEST_CASE("ADAA2 matches the repeated-outer limit across the node threshold")
{
    constexpr float a = 0.005f;
    constexpr float b = -0.5f;
    const float inside_threshold = std::nextafter(a, 1.0f);
    const float outside_threshold = std::nextafter(
        std::nextafter(inside_threshold, 1.0f), 1.0f);
    EON_CHECK(std::abs(static_cast<double>(inside_threshold) - a) < 1e-9);
    EON_CHECK(std::abs(static_cast<double>(outside_threshold) - a) > 1e-9);

    const double difference = static_cast<double>(a) - b;
    const double expected = 2.0 * (eon::SoftClip::F1(a) * difference
                                 - eon::SoftClip::F2(a)
                                 + eon::SoftClip::F2(b))
                          / (difference * difference);
    const auto process = [](float final_outer)
    {
        eon::SoftClipSat saturator;
        saturator.process(a);
        saturator.process(b);
        return saturator.process(final_outer);
    };
    const float actual_inside = process(inside_threshold);
    const float actual_outside = process(outside_threshold);

    EON_CHECK_NEAR(actual_inside, expected, 1e-6);
    EON_CHECK_NEAR(actual_outside, expected, 1e-6);
    EON_CHECK_NEAR(actual_inside, actual_outside, 1e-6);
}

EON_TEST_CASE("ADAA2 keeps small repeated nodes accurate")
{
    eon::SoftClipSat saturator;
    constexpr float a = 0.005f;
    const float b = std::nextafter(a, 0.0f);
    saturator.process(a);
    saturator.process(b);
    const float actual = saturator.process(a);

    const double center = static_cast<double>(a)
                        + (static_cast<double>(b) - a) / 3.0;
    EON_CHECK_NEAR(actual, eon::SoftClip::f(center), 1e-8);
}

EON_TEST_CASE("ADAA2 uses the three-node centroid when all nodes are close")
{
    eon::SoftClipSat saturator;
    constexpr float a = 4.9e-10f;
    saturator.process(a);
    saturator.process(-a);
    const float actual = saturator.process(a);
    const double center = (static_cast<double>(a) - a + a) / 3.0;

    EON_CHECK_NEAR(actual, eon::SoftClip::f(center), 2e-17);
}

EON_TEST_CASE("ADAA2 reaches the constant-input limit")
{
    eon::SoftClipSat saturator;
    saturator.process(0.5f);
    saturator.process(0.5f);
    const float actual = saturator.process(0.5f);

    EON_CHECK_NEAR(actual, eon::SoftClip::f(0.5), 1e-6);
}
