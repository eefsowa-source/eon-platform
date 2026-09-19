#include "Dsp/Adaa.h"
#include "Tests/TestHarness.h"

EON_TEST_CASE("nonlinear suite smoke test")
{
    EON_CHECK(eon::SoftClip::f(0.5) > 0.0);
}
