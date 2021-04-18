cmake_minimum_required(VERSION 3.10)

include(CTest)

# project(core_test)
# add_executable(core_test core/core.cpp)
# target_link_libraries(core_test math)
# add_test(NAME core COMMAND core_test)

# project(evaluation_test)
# add_executable(evaluation_test evaluation/evaluate.cpp)
# target_link_libraries(evaluation_test math)
# add_test(NAME evaluation COMMAND evaluation_test)

# project(simplification_rational_test)
# add_executable(simplification_rational_test simplification/rational.cpp)
# target_link_libraries(simplification_rational_test math)
# add_test(NAME simplification_rational COMMAND simplification_rational_test)

# project(ordering_test)
# add_executable(ordering_test ordering/ordering.cpp)
# target_link_libraries(ordering_test math)
# add_test(NAME ordering COMMAND ordering_test)

project(simplification_summation_test)
add_executable(simplification_summation_test simplification/summation.cpp)
target_link_libraries(simplification_summation_test math)
add_test(NAME simplification_summation COMMAND simplification_summation_test)


project(simplification_product_test)
add_executable(simplification_product_test simplification/product.cpp)
target_link_libraries(simplification_product_test math)
add_test(NAME simplification_product COMMAND simplification_product_test)
