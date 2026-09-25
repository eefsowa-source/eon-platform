cmake_minimum_required(VERSION 3.16)

project(eon_dsp LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Build quality for a numerical DSP core.
#
# The core is header-only and framework-free, so it must not impose its own
# optimisation flags on consumers. What it *does* own is the correctness
# policy: the tests are compiled with a strict floating-point model and with
# warnings that catch the failure modes this library actually has (implicit
# narrowing, sign changes, uninitialised state, shadowed members).
#
# -ffast-math is deliberately NOT used anywhere. It assumes no NaN/Inf, allows
# reassociation, and makes std::isfinite() unreliable — this library's solvers
# (WDF diode, Lambert W, triode Newton) rely on those checks to stay finite.
# ---------------------------------------------------------------------------
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(EON_DSP_SANITIZERS "Build tests with AddressSanitizer + UndefinedBehaviorSanitizer" OFF)
option(EON_DSP_FAST_MATH "Allow -ffast-math in tests (parity checks only, not recommended)" OFF)

add_library(eon_dsp INTERFACE)
target_include_directories(eon_dsp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(eon_dsp INTERFACE cxx_std_20)

# Warning set shared by every test target.
#
# Deliberately excludes -Wconversion/-Wsign-conversion: the DSP code moves
# between float and double at every stage boundary by design, so those two
# would emit hundreds of lines of noise and bury the warnings that matter.
add_library(eon_dsp_warnings INTERFACE)
if(MSVC)
    target_compile_options(eon_dsp_warnings INTERFACE /W4 /permissive-)
else()
    target_compile_options(eon_dsp_warnings INTERFACE
        -Wall -Wextra -Wpedantic -Wshadow -Wcast-qual -Wnon-virtual-dtor)
endif()

# Deterministic floating point: no contraction, no reassociation. The library
# promises block-partition invariance and reproducible solver limits, which
# requires the compiler not to fuse or reorder the arithmetic underneath it.
add_library(eon_dsp_strict_fp INTERFACE)
if(MSVC)
    target_compile_options(eon_dsp_strict_fp INTERFACE /fp:strict)
else()
    target_compile_options(eon_dsp_strict_fp INTERFACE -ffp-contract=off -fno-fast-math)
endif()

add_library(eon_dsp_sanitizers INTERFACE)
if(EON_DSP_SANITIZERS)
    if(MSVC)
        target_compile_options(eon_dsp_sanitizers INTERFACE /fsanitize=address)
        target_link_options(eon_dsp_sanitizers INTERFACE /fsanitize=address)
    else()
        target_compile_options(eon_dsp_sanitizers INTERFACE
            -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(eon_dsp_sanitizers INTERFACE -fsanitize=address,undefined)
    endif()
endif()

if(EON_DSP_FAST_MATH AND NOT MSVC)
    target_compile_options(eon_dsp INTERFACE -ffast-math)
endif()

include(CTest)

if(BUILD_TESTING)
    function(eon_add_test_target target)
        target_link_libraries(${target} PRIVATE
            eon_dsp eon_dsp_warnings eon_dsp_strict_fp eon_dsp_sanitizers)
    endfunction()

    add_executable(eon_dsp_tests
        Tests/TestMain.cpp
        Tests/NonlinearTests.cpp
        Tests/OversamplingTests.cpp
        Tests/RateIndependenceTests.cpp
        Tests/UtilityTests.cpp
    )
    eon_add_test_target(eon_dsp_tests)

    add_executable(eon_dsp_empty_tests Tests/TestMain.cpp)
    eon_add_test_target(eon_dsp_empty_tests)

    add_test(NAME eon_dsp_tests COMMAND eon_dsp_tests)
    add_test(NAME eon_dsp_empty_registry COMMAND eon_dsp_empty_tests)
    set_tests_properties(eon_dsp_empty_registry PROPERTIES WILL_FAIL TRUE)

    # Every public header must compile standalone (no hidden include order).
    add_library(eon_dsp_header_checks OBJECT
        Tests/HeaderCompile/Adaa.cpp
        Tests/HeaderCompile/Measure.cpp
        Tests/HeaderCompile/Oversampling.cpp
        Tests/HeaderCompile/Rate.cpp
        Tests/HeaderCompile/Rng.cpp
        Tests/HeaderCompile/RtGuard.cpp
        Tests/HeaderCompile/Solvers.cpp
        Tests/HeaderCompile/Stages.cpp
        Tests/HeaderCompile/Transformer.cpp
        Tests/HeaderCompile/Triode.cpp
        Tests/HeaderCompile/Wdf.cpp
        Tests/HeaderCompile/Zdf.cpp
    )
    eon_add_test_target(eon_dsp_header_checks)
endif()
