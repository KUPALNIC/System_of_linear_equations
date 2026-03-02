cmake_minimum_required(VERSION 3.20)
project(SystemOfLinearEquations)

option(WITH_TESTS "enable tests" OFF)

add_library(${PROJECT_NAME})

add_subdirectory(src)

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()