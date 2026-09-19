cmake_minimum_required(VERSION 3.16)

project(eon_dsp LANGUAGES CXX)

add_library(eon_dsp INTERFACE)
target_include_directories(eon_dsp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(eon_dsp INTERFACE cxx_std_20)

include(CTest)

if(BUILD_TESTING)
    function(eon_enable_warnings target)
        if(MSVC)
            target_compile_options(${target} PRIVATE /W4 /permissive-)
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
            target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        endif()
    endfunction()

    add_executable(eon_dsp_tests
        Tests/TestMain.cpp
        Tests/NonlinearTests.cpp
        Tests/OversamplingTests.cpp
        Tests/UtilityTests.cpp
    )
    target_link_libraries(eon_dsp_tests PRIVATE eon_dsp)
    eon_enable_warnings(eon_dsp_tests)

    add_executable(eon_dsp_empty_tests Tests/TestMain.cpp)
    target_link_libraries(eon_dsp_empty_tests PRIVATE eon_dsp)
    eon_enable_warnings(eon_dsp_empty_tests)

    add_test(NAME eon_dsp_tests COMMAND eon_dsp_tests)
    add_test(NAME eon_dsp_empty_registry COMMAND eon_dsp_empty_tests)
    set_tests_properties(eon_dsp_empty_registry PROPERTIES WILL_FAIL TRUE)

    add_library(eon_dsp_header_checks OBJECT
        Tests/HeaderCompile/Adaa.cpp
        Tests/HeaderCompile/Measure.cpp
        Tests/HeaderCompile/Oversampling.cpp
        Tests/HeaderCompile/Rng.cpp
        Tests/HeaderCompile/RtGuard.cpp
        Tests/HeaderCompile/Solvers.cpp
        Tests/HeaderCompile/Stages.cpp
        Tests/HeaderCompile/Transformer.cpp
        Tests/HeaderCompile/Triode.cpp
        Tests/HeaderCompile/Wdf.cpp
        Tests/HeaderCompile/Zdf.cpp
    )
    target_link_libraries(eon_dsp_header_checks PRIVATE eon_dsp)
    eon_enable_warnings(eon_dsp_header_checks)
endif()
