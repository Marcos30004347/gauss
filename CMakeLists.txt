cmake_minimum_required(VERSION 3.10)

project (math) 

add_library(
    math
    src/Core/AST/AST.cpp
    src/Core/Algebra/Algebra.cpp
    src/Core/Polynomial/Polynomial.cpp
    src/Core/Expand/Expand.cpp
    src/Core/Expand/Multiplication.cpp
    src/Core/Expand/Polynomial.cpp
    src/Core/Expand/Division.cpp
    src/Core/Reduce/Rationals.cpp
    src/Core/Reduce/Power.cpp
    src/Core/Reduce/Multiplication.cpp
    src/Core/Reduce/Addition.cpp
    src/Core/Reduce/Subtraction.cpp
    src/Core/Reduce/Division.cpp
    src/Core/Reduce/Reduce.cpp
    src/Core/Reduce/Factorial.cpp

    # src/Core/.Reduce/Rationals.cpp
    # src/Core/.Reduce/Power.cpp
    # src/Core/.Reduce/Multiplication.cpp
    # src/Core/.Reduce/Addition.cpp
    # src/Core/.Reduce/Subtraction.cpp
    # src/Core/.Reduce/Division.cpp
    # src/Core/.Reduce/Reduce.cpp
    # src/Core/.Reduce/Factorial.cpp
)

target_include_directories(
    math
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)

enable_testing()
add_subdirectory(tests)


