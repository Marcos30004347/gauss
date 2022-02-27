cmake_minimum_required(VERSION 3.10)

project (math)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3")

set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++0x")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


add_library(
  math
  src/MathSystem/Algebra/Expression.cpp
  src/MathSystem/Algebra/Integer.cpp
  src/MathSystem/Algebra/Int.cpp
    # src/MathSystem/Trigonometry/Trigonometry.cpp
    # src/MathSystem/Exponential/Exponential.cpp
    src/MathSystem/Calculus/Derivative.cpp
    # src/MathSystem/Calculus/Integral.cpp
    # src/MathSystem/Rational/Rational.cpp
    src/MathSystem/Primes/Primes.cpp
    # src/MathSystem/Expand/Expand.cpp
    # src/MathSystem/Expand/Multiplication.cpp
    # src/MathSystem/Expand/Polynomial.cpp
    # src/MathSystem/Expand/Division.cpp
    # src/MathSystem/Simplification/Rationals.cpp
    # src/MathSystem/Simplification/Power.cpp
    # src/MathSystem/Simplification/Multiplication.cpp
    # src/MathSystem/Simplification/Addition.cpp
    # src/MathSystem/Simplification/Subtraction.cpp
    # src/MathSystem/Simplification/Division.cpp
    # src/MathSystem/Simplification/Simplification.cpp
    # src/MathSystem/Simplification/Factorial.cpp
    # src/MathSystem/Simplification/Trigonometry.cpp
    # src/MathSystem/Simplification/Exponential.cpp

    src/MathSystem/Factorization/Utils.cpp
    src/MathSystem/Factorization/Hensel.cpp
		src/MathSystem/Factorization/Berlekamp.cpp
		src/MathSystem/Factorization/SquareFree.cpp
    src/MathSystem/Factorization/Zassenhaus.cpp
    src/MathSystem/Factorization/Wang.cpp

		# #src/MathSystem/Polynomial/Factorization.cpp

    src/MathSystem/GaloisField/GaloisField.cpp
		src/MathSystem/Polynomial/Resultant.cpp
		src/MathSystem/Polynomial/Polynomial.cpp

		# src/MathSystem/Polynomial/Algebraic.cpp

		# #src/MathSystem/Polynomial/Hensel.cpp
		src/TypeSystem/Types.cpp
		src/TypeSystem/SymbolTable.cpp
)

target_include_directories(
    math
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

enable_testing()
add_subdirectory(tests)
