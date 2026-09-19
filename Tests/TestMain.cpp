#include "TestHarness.h"

void runAdaaOversamplingTests();
void runStageTriodeTests();
void runSolverWdfZdfTests();

int main()
{
    runAdaaOversamplingTests();
    runStageTriodeTests();
    runSolverWdfZdfTests();
    std::printf("\n%d failure(s)\n", eon::test::failures);
    return eon::test::failures == 0 ? 0 : 1;
}
