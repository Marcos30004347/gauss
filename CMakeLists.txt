cmake_minimum_required(VERSION 3.10)

project (math) 


add_library(
    math
    src/Core/AST/AST.cpp
    src/Core/Debug/Assert.cpp
    src/Core/Algebra/Algebra.cpp
    src/Core/Algebra/List.cpp
    src/Core/Algebra/Set.cpp
    src/Core/Polynomial/Polynomial.cpp
    src/Core/Trigonometry/Trigonometry.cpp
    src/Core/Calculus/Calculus.cpp
    src/Core/Rational/Rational.cpp
    src/Core/Primes/Primes.cpp
    src/Core/Expand/Expand.cpp
    src/Core/Expand/Multiplication.cpp
    src/Core/Expand/Polynomial.cpp
    src/Core/Expand/Division.cpp
    src/Core/Simplification/Rationals.cpp
    src/Core/Simplification/Power.cpp
    src/Core/Simplification/Multiplication.cpp
    src/Core/Simplification/Addition.cpp
    src/Core/Simplification/Subtraction.cpp
    src/Core/Simplification/Division.cpp
    src/Core/Simplification/Simplification.cpp
    src/Core/Simplification/Factorial.cpp
    src/Core/Simplification/Trigonometry.cpp
    src/Core/Simplification/Exponential.cpp
    src/Core/Polynomial/Zp.cpp
    src/Core/Polynomial/Factorization.cpp
)

target_include_directories(
    math
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

enable_testing()
add_subdirectory(tests)


