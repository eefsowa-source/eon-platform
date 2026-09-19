#pragma once

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace eon::test {

struct TestCase
{
    std::string_view name;
    void (*run)();
};

inline std::vector<TestCase>& registry()
{
    static std::vector<TestCase> tests;
    return tests;
}

class Registrar
{
public:
    Registrar(std::string_view name, void (*run)())
    {
        registry().push_back({name, run});
    }
};

class AssertionFailure : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

[[noreturn]] inline void fail(std::string_view expression,
                              std::string_view file,
                              int line,
                              std::string_view detail = {})
{
    std::ostringstream message;
    message << file << ':' << line << ": assertion '" << expression << "' failed";
    if (!detail.empty())
        message << " (" << detail << ')';
    throw AssertionFailure(message.str());
}

inline void check(bool condition,
                  std::string_view expression,
                  std::string_view file,
                  int line)
{
    if (!condition)
        fail(expression, file, line);
}

inline void check_near(double actual,
                       double expected,
                       double tolerance,
                       std::string_view expression,
                       std::string_view file,
                       int line)
{
    if (!(std::abs(actual - expected) <= tolerance))
    {
        std::ostringstream detail;
        detail.precision(17);
        detail << "actual=" << actual << ", expected=" << expected
               << ", tolerance=" << tolerance;
        fail(expression, file, line, detail.str());
    }
}

} // namespace eon::test

#define EON_TEST_DETAIL_JOIN_IMPL(a, b) a##b
#define EON_TEST_DETAIL_JOIN(a, b) EON_TEST_DETAIL_JOIN_IMPL(a, b)
#define EON_TEST_DETAIL_CASE(name, unique_id) \
    static void EON_TEST_DETAIL_JOIN(eon_test_case_, unique_id)(); \
    static const ::eon::test::Registrar EON_TEST_DETAIL_JOIN(eon_test_registration_, unique_id)( \
        name, &EON_TEST_DETAIL_JOIN(eon_test_case_, unique_id)); \
    static void EON_TEST_DETAIL_JOIN(eon_test_case_, unique_id)()
// Put each EON_TEST_CASE invocation on a distinct source line; __LINE__ is its unique ID.
#define EON_TEST_CASE(name) EON_TEST_DETAIL_CASE(name, __LINE__)

#define EON_CHECK(condition) \
    ::eon::test::check(static_cast<bool>(condition), #condition, __FILE__, __LINE__)
#define EON_CHECK_NEAR(actual, expected, tolerance) \
    ::eon::test::check_near(static_cast<double>(actual), \
                            static_cast<double>(expected), \
                            static_cast<double>(tolerance), \
                            #actual " ~= " #expected, __FILE__, __LINE__)
