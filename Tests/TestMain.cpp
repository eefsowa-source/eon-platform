#include "Tests/TestHarness.h"

#include <cstddef>
#include <exception>
#include <iostream>

int main()
{
    std::size_t failures = 0;
    for (const auto& test : eon::test::registry())
    {
        try
        {
            test.run();
            std::cout << "[PASS] " << test.name << '\n';
        }
        catch (const std::exception& error)
        {
            ++failures;
            std::cerr << "[FAIL] " << test.name << ": " << error.what() << '\n';
        }
        catch (...)
        {
            ++failures;
            std::cerr << "[FAIL] " << test.name << ": unknown exception\n";
        }
    }

    const std::size_t total = eon::test::registry().size();
    std::cout << total - failures << '/' << total << " tests passed";
    if (failures != 0)
        std::cout << ", " << failures << " failed";
    std::cout << '\n';
    return failures == 0 ? 0 : 1;
}
