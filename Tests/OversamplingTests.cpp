#include "Dsp/Oversampling.h"
#include "Tests/TestHarness.h"

EON_TEST_CASE("oversampling suite smoke test")
{
    eon::Oversampler oversampler;
    oversampler.setStages(1);
    EON_CHECK(oversampler.factor() == 2);
}
