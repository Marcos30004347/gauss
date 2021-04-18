cmake_minimum_required(VERSION 3.10)

project (math) 

add_library(
    math
    src/core/expression.cpp
    src/core/power.cpp
    src/core/rational.cpp
    src/ordering/ordering.cpp
    src/evaluation/evaluate.cpp
    src/simplification/rationals.cpp
    src/simplification/summations.cpp
    src/simplification/quotient.cpp
    src/simplification/powers.cpp
    src/simplification/difference.cpp
    src/simplification/products.cpp
)


target_include_directories(math PUBLIC src)

enable_testing()
add_subdirectory(tests)


