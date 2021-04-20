cmake_minimum_required(VERSION 3.10)

include(CTest)

project(containers_test)
add_executable(containers_test core/containers.cpp)
target_link_libraries(containers_test math)
add_test(NAME containers COMMAND containers_test)


# project(algebra_test)
# add_executable(algebra_test algebra/algebra.cpp)
# target_link_libraries(algebra_test math)
# add_test(NAME algebra COMMAND algebra_test)

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

# project(simplification_summation_test)
# add_executable(simplification_summation_test simplification/summation.cpp)
# target_link_libraries(simplification_summation_test math)
# add_test(NAME simplification_summation COMMAND simplification_summation_test)


# project(simplification_product_test)
# add_executable(simplification_product_test simplification/product.cpp)
# target_link_libraries(simplification_product_test math)
# add_test(NAME simplification_product COMMAND simplification_product_test)
