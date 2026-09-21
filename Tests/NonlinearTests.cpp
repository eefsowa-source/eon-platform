#include "Dsp/Adaa.h"
#include "Dsp/Stages.h"
#include "Dsp/Triode.h"
#include "Tests/TestHarness.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

static void checkTriodeSelfBiasResetStability(float cathodeRk)
{
    eon::TriodeStage stage;
    stage.iterations = 2;
    stage.setCathode(cathodeRk, 25e-6f, 48000.0 * 32.0);
    stage.reset();

    EON_CHECK(stage.iterations == 2);
    EON_CHECK_NEAR(stage.VpPrev, stage.Vq, 0.0f);
    const float expected_out_scale = 1.0f / std::max(40.0f, stage.Bplus - stage.Vq);
    EON_CHECK_NEAR(stage.outScale, expected_out_scale, 0.0f);
    const double reset_cathode_vk = stage.cathodeVk;
    const float reset_last_ip = stage.lastIp;
    const float reset_vq = stage.Vq;
    const float reset_plate_ip = stage.korenIp(
        stage.biasVg - static_cast<float>(reset_cathode_vk),
        reset_vq - static_cast<float>(reset_cathode_vk));
    const double cathode_residual = reset_last_ip * cathodeRk - reset_cathode_vk;
    const double plate_residual = (stage.Bplus - reset_vq) / stage.Rload
                                - reset_plate_ip;
    const float immediate_output = stage.process(0.0f);
    float maximum_output = std::abs(immediate_output);
    for (int sample = 0; sample < 500000; ++sample)
        maximum_output = std::max(maximum_output, std::abs(stage.process(0.0f)));

    EON_CHECK_NEAR(cathode_residual, 0.0, 1e-3);
    EON_CHECK_NEAR(reset_last_ip, reset_plate_ip, 1e-8f);
    EON_CHECK_NEAR(plate_residual, 0.0, 1e-8);
    EON_CHECK(std::abs(immediate_output) < 1e-4f);
    EON_CHECK(maximum_output < 1e-4f);
}

EON_TEST_CASE("Triode 2700-ohm self-bias reset converges and remains stable")
{
    checkTriodeSelfBiasResetStability(2700.0f);
}

EON_TEST_CASE("Triode 3510-ohm self-bias reset converges and remains stable")
{
    checkTriodeSelfBiasResetStability(3510.0f);
}

EON_TEST_CASE("Triode grounded-cathode reset starts at its zero-input operating point")
{
    eon::TriodeStage stage;
    stage.iterations = 2;
    stage.reset();

    EON_CHECK(stage.iterations == 2);
    EON_CHECK_NEAR(stage.VpPrev, stage.Vq, 0.0f);
    const float expected_out_scale = 1.0f / std::max(40.0f, stage.Bplus - stage.Vq);
    EON_CHECK_NEAR(stage.outScale, expected_out_scale, 0.0f);
    const float reset_last_ip = stage.lastIp;
    const float reset_plate_ip = stage.korenIp(stage.biasVg, stage.Vq);
    const float output = stage.process(0.0f);
    const float process_plate_ip = stage.korenIp(stage.biasVg, stage.VpPrev);
    // A 1e-5 normalized bound leaves ample room for float Newton rounding.
    EON_CHECK(std::abs(output) < 1e-5f);
    EON_CHECK_NEAR(stage.lastIp, process_plate_ip, 1e-8f);

    const double plate_residual = (stage.Bplus - stage.Vq) / stage.Rload
                                - reset_plate_ip;
    EON_CHECK_NEAR(reset_last_ip, reset_plate_ip, 1e-8f);
    EON_CHECK_NEAR(plate_residual, 0.0, 1e-8);
}

EON_TEST_CASE("Triode self-bias reset stays at its zero-input operating point")
{
    eon::TriodeStage stage;
    stage.iterations = 2;
    stage.setCathode(1500.0f, 25e-6f, 48000.0 * 32.0);
    stage.reset();

    EON_CHECK(stage.iterations == 2);
    EON_CHECK_NEAR(stage.VpPrev, stage.Vq, 0.0f);
    const float expected_out_scale = 1.0f / std::max(40.0f, stage.Bplus - stage.Vq);
    EON_CHECK_NEAR(stage.outScale, expected_out_scale, 0.0f);
    const double reset_cathode_vk = stage.cathodeVk;
    const float reset_last_ip = stage.lastIp;
    const float reset_plate_ip = stage.korenIp(
        stage.biasVg - static_cast<float>(reset_cathode_vk),
        stage.Vq - static_cast<float>(reset_cathode_vk));
    const float process_vk = static_cast<float>(reset_cathode_vk);
    const float immediate_output = stage.process(0.0f);
    const float immediate_last_ip = stage.lastIp;
    const float immediate_plate_ip = stage.korenIp(
        stage.biasVg - process_vk, stage.VpPrev - process_vk);
    float maximum_output = std::abs(immediate_output);
    for (int sample = 0; sample < 500000; ++sample)
        maximum_output = std::max(maximum_output, std::abs(stage.process(0.0f)));

    // Repeated float state updates stay within this 1e-4 normalized bound.
    EON_CHECK(maximum_output < 1e-4f);
    EON_CHECK(std::abs(immediate_output) < 1e-5f);
    EON_CHECK_NEAR(immediate_last_ip, immediate_plate_ip, 1e-8f);
    const double plate_residual = (stage.Bplus - stage.Vq) / stage.Rload
                                - reset_plate_ip;
    EON_CHECK_NEAR(reset_last_ip, reset_plate_ip, 1e-8f);
    EON_CHECK_NEAR(plate_residual, 0.0, 1e-8);
    EON_CHECK_NEAR(reset_cathode_vk, reset_last_ip * stage.cathodeRk, 1e-3);
}

EON_TEST_CASE("DC blocker prepare keeps its cutoff in hertz across rates")
{
    const eon::DCBlocker unprepared;
    EON_CHECK_NEAR(unprepared.R, 0.9997, 0.0);

    constexpr double rates[] = {48000.0 * 8.0, 48000.0 * 16.0, 48000.0 * 32.0};
    for (const double sample_rate : rates)
    {
        eon::DCBlocker blocker;
        blocker.prepare(sample_rate, 18.0);
        const double recovered_cutoff = -std::log(blocker.R) * sample_rate
                                      / (2.0 * std::numbers::pi);
        EON_CHECK_NEAR(recovered_cutoff, 18.0, 1e-10);
    }
}

EON_TEST_CASE("DC blocker reset clears prepared impulse history")
{
    eon::DCBlocker blocker;
    blocker.prepare(48000.0, 18.0);
    EON_CHECK_NEAR(blocker.process(1.0f), 1.0, 0.0);
    EON_CHECK(blocker.process(0.0f) < 0.0f);

    blocker.reset();
    EON_CHECK_NEAR(blocker.process(0.0f), 0.0, 0.0);
}

EON_TEST_CASE("ClassA positive limiting is continuous around its saturation boundary")
{
    constexpr float below_boundary = 3.27077f;
    constexpr float above_boundary = 3.27078f;
    const auto first_output = [](float input)
    {
        eon::ClassAStage stage;
        return stage.process(input);
    };

    const float positive_below = first_output(below_boundary);
    const float positive_above = first_output(above_boundary);

    EON_CHECK_NEAR(std::abs(positive_above - positive_below), 0.0, 1e-3);
}

EON_TEST_CASE("ClassA negative limiting is continuous around its saturation boundary")
{
    constexpr float below_boundary = -3.43981f;
    constexpr float above_boundary = -3.43982f;
    const auto first_output = [](float input)
    {
        eon::ClassAStage stage;
        return stage.process(input);
    };

    const float negative_below = first_output(below_boundary);
    const float negative_above = first_output(above_boundary);

    EON_CHECK_NEAR(std::abs(negative_above - negative_below), 0.0, 1e-3);
}

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
