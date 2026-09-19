#include "Dsp/Oversampling.h"
#include "Tests/TestHarness.h"

#include <cmath>

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
