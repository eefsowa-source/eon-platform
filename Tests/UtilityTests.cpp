#include "Dsp/Rng.h"
#include "Tests/TestHarness.h"

EON_TEST_CASE("utility suite smoke test")
{
    eon::Rng first(1234);
    eon::Rng second(1234);
    EON_CHECK(first.nextU64() == second.nextU64());
    EON_CHECK_NEAR(first.next(), second.next(), 0.0);
    EON_CHECK_NEAR(0.1 + 0.2, 0.3, 1e-15);
}
