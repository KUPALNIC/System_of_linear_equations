cmake_minimum_required(VERSION 3.20)
project(SystemOfLinearEquations)

option(WITH_TESTS "enable tests" OFF)

add_library(${PROJECT_NAME})

add_subdirectory(src)

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

add_executable(benchmark_run theory/comp_CSR_matrix.cpp src/CSR/csr.cpp src/matrix/matrix.cpp)
target_include_directories(benchmark_run PRIVATE src/)
