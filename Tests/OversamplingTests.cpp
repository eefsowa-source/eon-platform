#include "Dsp/Oversampling.h"
#include "Tests/TestHarness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

constexpr int sine_sample_count = 8192;
constexpr int warmup_sample_count = 256;
constexpr int largest_block_size = 256;

void render_round_trip(eon::Oversampler& oversampler,
                       const std::vector<float>& input,
                       int block_size,
                       std::vector<float>& upsampled,
                       std::vector<float>& output)
{
    oversampler.reset();
    for (int offset = 0; offset < static_cast<int>(input.size()); offset += block_size)
    {
        const int count = std::min(block_size, static_cast<int>(input.size()) - offset);
        oversampler.up(input.data() + offset, count, upsampled.data());
        oversampler.down(upsampled.data(), output.data() + offset, count);
    }
}

double peak_difference_after_warmup(const std::vector<float>& first,
                                    const std::vector<float>& second)
{
    double peak = 0.0;
    for (int i = warmup_sample_count; i < static_cast<int>(first.size()); ++i)
        peak = std::max(peak, std::abs(static_cast<double>(first[static_cast<size_t>(i)]) -
                                      static_cast<double>(second[static_cast<size_t>(i)])));
    return peak;
}

} // namespace

EON_TEST_CASE("oversampling suite smoke test")
{
    eon::Oversampler oversampler;
    oversampler.setStages(1);
    EON_CHECK(oversampler.factor() == 2);

    constexpr int sample_count = 4;
    const float input[sample_count] = {0.25f, -0.5f, 0.125f, 0.0f};
    float upsampled[sample_count * 2] = {};
    float output[sample_count] = {};
    float reference[sample_count] = {};
    oversampler.prepare(sample_count);

    const auto render = [&]()
    {
        oversampler.up(input, sample_count, upsampled);
        for (const float sample : upsampled)
            EON_CHECK(std::isfinite(sample));

        oversampler.down(upsampled, output, sample_count);
        for (const float sample : output)
            EON_CHECK(std::isfinite(sample));
    };

    render();
    for (int i = 0; i < sample_count; ++i)
        reference[i] = output[i];

    oversampler.reset();
    render();
    for (int i = 0; i < sample_count; ++i)
        EON_CHECK_NEAR(output[i], reference[i], 0.0);
}

EON_TEST_CASE("oversampling round trip is independent of block partition")
{
    std::vector<float> input(sine_sample_count);
    std::vector<float> upsampled(static_cast<size_t>(largest_block_size) * 8u);
    std::vector<float> output_small(sine_sample_count);
    std::vector<float> output_large(sine_sample_count);
    for (int i = 0; i < sine_sample_count; ++i)
        input[static_cast<size_t>(i)] = static_cast<float>(0.5 * std::sin(
            2.0 * std::acos(-1.0) * 1000.0 * static_cast<double>(i) / 48000.0));

    double peak_differences[3] = {};
    int factors[3] = {};
    for (int stages = 1; stages <= 3; ++stages)
    {
        eon::Oversampler oversampler;
        oversampler.setStages(stages);
        oversampler.prepare(largest_block_size);
        factors[stages - 1] = oversampler.factor();

        render_round_trip(oversampler, input, 16, upsampled, output_small);
        render_round_trip(oversampler, input, largest_block_size, upsampled, output_large);
        peak_differences[stages - 1] = peak_difference_after_warmup(output_small, output_large);
    }

    for (int i = 0; i < 3; ++i)
        std::cout << "[INFO] oversampling factor " << factors[i]
                  << " block partition peak difference after warmup: "
                  << std::setprecision(12) << peak_differences[i] << '\n';

    for (double peak_difference : peak_differences)
        EON_CHECK(peak_difference < 1.0e-6);
}

EON_TEST_CASE("oversampling reset clears interpolation and decimation state at every factor")
{
    constexpr int sample_count = 64;
    std::vector<float> input(sample_count);
    std::vector<float> downsample_input(static_cast<size_t>(sample_count) * 8u);
    std::vector<float> actual_up(static_cast<size_t>(sample_count) * 8u);
    std::vector<float> expected_up(static_cast<size_t>(sample_count) * 8u);
    std::vector<float> actual_down(sample_count);
    std::vector<float> expected_down(sample_count);
    for (int i = 0; i < sample_count; ++i)
        input[static_cast<size_t>(i)] = static_cast<float>(std::sin(0.17 * static_cast<double>(i)));
    for (size_t i = 0; i < downsample_input.size(); ++i)
        downsample_input[i] = static_cast<float>(std::cos(0.11 * static_cast<double>(i)));

    for (int stages = 1; stages <= 3; ++stages)
    {
        eon::Oversampler oversampler;
        eon::Oversampler fresh_reference;
        oversampler.setStages(stages);
        fresh_reference.setStages(stages);
        oversampler.prepare(sample_count);
        fresh_reference.prepare(sample_count);

        oversampler.up(input.data(), sample_count, actual_up.data());
        oversampler.down(downsample_input.data(), actual_down.data(), sample_count);
        oversampler.reset();

        oversampler.up(input.data(), sample_count, actual_up.data());
        fresh_reference.up(input.data(), sample_count, expected_up.data());
        for (int i = 0; i < sample_count * oversampler.factor(); ++i)
            EON_CHECK_NEAR(actual_up[static_cast<size_t>(i)], expected_up[static_cast<size_t>(i)], 0.0);

        oversampler.down(downsample_input.data(), actual_down.data(), sample_count);
        fresh_reference.down(downsample_input.data(), expected_down.data(), sample_count);
        for (int i = 0; i < sample_count; ++i)
            EON_CHECK_NEAR(actual_down[static_cast<size_t>(i)], expected_down[static_cast<size_t>(i)], 0.0);
    }
}
