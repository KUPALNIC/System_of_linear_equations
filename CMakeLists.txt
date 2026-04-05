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
add_executable(benchmark_4
    theory/Solvers/solvers_bench.cpp 
    src/CSR/csr.cpp 
    src/matrix/matrix.cpp
    src/Solvers/GaussSeidelSolver.cpp
    src/Solvers/JacobiSolver.cpp
    src/Solvers/SimpleIterationSolver.cpp
   
)
add_executable(test_1 theory/test_1.cpp 
    src/matrix/matrix.cpp
    src/CSR/csr.cpp
    src/Solvers/JacobiSolver.cpp)

add_executable(bench_simm theory/simmetric/bench_sim.cpp 
    src/CSR/csr.cpp 
    src/matrix/matrix.cpp
    src/Solvers/GaussSeidelSolver.cpp
    src/Solvers/JacobiSolver.cpp
    src/Solvers/SimpleIterationSolver.cpp
    src/Solvers/SymmetricGaussSeidel.cpp
   
)
target_include_directories(benchmark_run PRIVATE src/)
target_include_directories(benchmark_4 PRIVATE src/)
