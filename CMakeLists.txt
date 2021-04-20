cmake_minimum_required(VERSION 3.10)

project (math) 

add_library(
    math
    src/core/memory.cpp
    src/algebra/expression.cpp
    src/algebra/power.cpp
    src/algebra/rational.cpp
    src/algebra/product.cpp
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


