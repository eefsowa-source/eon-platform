#include "TestHarness.h"

#include "Dsp/Stages.h"
#include "Dsp/Triode.h"

#include <cmath>

void runStageTriodeTests()
{
    for (double rate : {48000.0 * 8.0, 48000.0 * 16.0, 48000.0 * 32.0})
    {
        eon::DCBlocker blocker;
        blocker.prepare(rate, 18.0);
        const double measured = -std::log(blocker.R) * rate / (2.0 * M_PI);
        eon::test::check(eon::test::near(measured, 18.0, 1.0e-9),
                         "DC blocker cutoff remains 18 Hz at every processing rate");
    }

    {
        eon::ClassAStage below;
        eon::ClassAStage above;
        const double y0 = below.process(3.27077f);
        const double y1 = above.process(3.27078f);
        eon::test::check(std::abs(y1 - y0) < 1.0e-3,
                         "Class-A limiter is continuous across its former threshold");
    }

    {
        eon::TriodeStage triode;
        triode.iterations = 2;
        triode.reset();
        eon::test::check(std::abs(triode.process(0.0f)) < 1.0e-5,
                         "grounded-cathode triode starts at zero output");
    }

    {
        eon::TriodeStage triode;
        triode.iterations = 2;
        triode.setCathode(1500.0f, 25.0e-6f, 48000.0 * 32.0);
        triode.reset();
        float output = 0.0f;
        for (int i = 0; i < 500000; ++i)
            output = triode.process(0.0f);
        eon::test::check(std::abs(output) < 1.0e-4,
                         "self-biased triode remains at its reset operating point");
    }
}
