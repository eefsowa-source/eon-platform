cmake_minimum_required(VERSION 3.22)
project(eon_dsp LANGUAGES CXX)

add_library(eon_dsp INTERFACE)
add_library(eon::dsp ALIAS eon_dsp)
target_include_directories(eon_dsp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(eon_dsp INTERFACE cxx_std_20)

include(CTest)
if(BUILD_TESTING)
    add_executable(eon_dsp_tests
        Tests/TestMain.cpp
        Tests/AdaaOversamplingTests.cpp
        Tests/StageTriodeTests.cpp
        Tests/SolverWdfZdfTests.cpp
    )
    target_link_libraries(eon_dsp_tests PRIVATE eon::dsp)
    if(MSVC)
        target_compile_options(eon_dsp_tests PRIVATE /W4)
    else()
        target_compile_options(eon_dsp_tests PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    add_test(NAME eon_dsp_core COMMAND eon_dsp_tests)
endif()
