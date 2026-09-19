cmake_minimum_required(VERSION 3.16)

project(eon_dsp LANGUAGES CXX)

add_library(eon_dsp INTERFACE)
target_include_directories(eon_dsp INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(eon_dsp INTERFACE cxx_std_20)

include(CTest)

if(BUILD_TESTING)
    add_executable(eon_dsp_tests
        Tests/TestMain.cpp
        Tests/NonlinearTests.cpp
        Tests/OversamplingTests.cpp
        Tests/UtilityTests.cpp
    )
    target_link_libraries(eon_dsp_tests PRIVATE eon_dsp)

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(eon_dsp_tests PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_test(NAME eon_dsp_tests COMMAND eon_dsp_tests)
endif()
