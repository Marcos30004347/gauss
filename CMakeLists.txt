cmake_minimum_required(VERSION 3.10)

project (math)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(
    math
    src/Core/AST/AST.cpp
    src/Core/AST/Integer.cpp
    src/Core/AST/Int.cpp
    src/Core/Debug/Assert.cpp
    src/Core/Algebra/Algebra.cpp
    src/Core/Algebra/List.cpp
    src/Core/Algebra/Set.cpp
    src/Core/Trigonometry/Trigonometry.cpp
    src/Core/Exponential/Exponential.cpp
    src/Core/Calculus/Derivative.cpp
    src/Core/Calculus/Integral.cpp
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
    #src/Core/Polynomial/Zp.cpp
    #src/Core/GaloisField/GaloisField.cpp
    #src/Core/Factorization/SquareFree.cpp
    #src/Core/Factorization/Hensel.cpp
    #src/Core/Factorization/Utils.cpp
    #src/Core/Factorization/Wang.cpp
    #src/Core/Factorization/Zassenhaus.cpp
    #src/Core/Factorization/Berlekamp.cpp
    #src/Core/Polynomial/Factorization.cpp
    src/Core/Polynomial/Polynomial.cpp
    src/Core/Polynomial/Algebraic.cpp
    #src/Core/Polynomial/Hensel.cpp
    src/Core/Polynomial/Resultant.cpp
)

target_include_directories(
    math
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

enable_testing()
add_subdirectory(tests)
