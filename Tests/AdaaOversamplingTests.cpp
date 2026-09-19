#include "TestHarness.h"

#include "Dsp/Adaa.h"
#include "Dsp/Oversampling.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

namespace {

std::vector<float> renderRoundTrip(int blockSize, int stages)
{
    constexpr int total = 8192;
    eon::Oversampler oversampler;
    oversampler.setStages(stages);
    oversampler.prepare(blockSize);

    std::vector<float> input(static_cast<std::size_t>(blockSize));
    std::vector<float> up(static_cast<std::size_t>(blockSize) *
                          static_cast<std::size_t>(oversampler.factor()));
    std::vector<float> output(static_cast<std::size_t>(total));
    for (int offset = 0; offset < total; offset += blockSize)
    {
        for (int i = 0; i < blockSize; ++i)
        {
            input[static_cast<std::size_t>(i)] = static_cast<float>(
                0.5 * std::sin(2.0 * std::numbers::pi * 1000.0 * (offset + i) / 48000.0));
        }
        oversampler.up(input.data(), blockSize, up.data());
        oversampler.down(up.data(), output.data() + offset, blockSize);
    }
    return output;
}

} // namespace

void runAdaaOversamplingTests()
{
    {
        eon::SoftClipSat sat;
        sat.reset();
        sat.process(0.5f);
        sat.process(-0.5f);
        const double actual = sat.process(0.5f);
        const double a = 0.5;
        const double b = -0.5;
        const double expected =
            2.0 * (eon::SoftClip::F1(a) * (a - b) - eon::SoftClip::F2(a) + eon::SoftClip::F2(b)) /
            ((a - b) * (a - b));
        eon::test::check(eon::test::near(actual, expected, 1.0e-6),
                         "ADAA2 uses the repeated-node limit when x0 equals x2");
    }

    {
        eon::SoftClipSat sat;
        for (float x : {0.5f, -0.5f, std::nextafter(0.5f, 1.0f)})
        {
            eon::test::check(std::isfinite(sat.process(x)),
                             "ADAA2 stays finite around the repeated-node branch");
        }
    }

    for (const int stages : {1, 2, 3})
    {
        const auto small = renderRoundTrip(16, stages);
        const auto large = renderRoundTrip(256, stages);
        double peak = 0.0;
        for (std::size_t i = 256; i < small.size(); ++i)
            peak = std::max(peak, std::abs(static_cast<double>(small[i]) - large[i]));
        eon::test::check(peak < 1.0e-6,
                         "oversampling output is independent of host block partition");
    }
}
