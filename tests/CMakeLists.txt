cmake_minimum_required(VERSION 3.10)

include(CTest)

project(ASTTests)
add_executable(ASTTests Core/AST/AST.cpp)
target_link_libraries(ASTTests math)
add_test(NAME ASTTests COMMAND ASTTests)

project(AlgebraTests)
add_executable(AlgebraTests Core/Algebra/Algebra.cpp)
target_link_libraries(AlgebraTests math)
add_test(NAME AlgebraTests COMMAND AlgebraTests)

project(ListTests)
add_executable(ListTests Core/Algebra/List.cpp)
target_link_libraries(ListTests math)
add_test(NAME ListTests COMMAND ListTests)

project(SetTests)
add_executable(SetTests Core/Algebra/Set.cpp)
target_link_libraries(SetTests math)
add_test(NAME SetTests COMMAND SetTests)

project(PolynomialTests)
add_executable(PolynomialTests Core/Polynomial/Polynomial.cpp)
target_link_libraries(PolynomialTests math)
add_test(NAME PolynomialTests COMMAND PolynomialTests)

# project(PolynomialZnTests)
# add_executable(PolynomialZnTests Core/Polynomial/Zn.cpp)
# target_link_libraries(PolynomialZnTests math)
# add_test(NAME PolynomialZnTests COMMAND PolynomialZnTests)

project(PolynomialFactorizationTests)
add_executable(PolynomialFactorizationTests Core/Polynomial/Factorization.cpp)
target_link_libraries(PolynomialFactorizationTests math)
add_test(NAME PolynomialFactorizationTests COMMAND PolynomialFactorizationTests)

project(RationalTests)
add_executable(RationalTests Core/Rational/Rational.cpp)
target_link_libraries(RationalTests math)
add_test(NAME RationalTests COMMAND RationalTests)

project(SimplificationRationals)
add_executable(SimplificationRationals Core/Simplification/Rationals.cpp)
target_link_libraries(SimplificationRationals math)
add_test(NAME SimplificationRationals COMMAND SimplificationRationals)

project(SimplificationPowers)
add_executable(SimplificationPowers Core/Simplification/Power.cpp)
target_link_libraries(SimplificationPowers math)
add_test(NAME SimplificationPowers COMMAND SimplificationPowers)

project(SimplificationAdditions)
add_executable(SimplificationAdditions Core/Simplification/Addition.cpp)
target_link_libraries(SimplificationAdditions math)
add_test(NAME SimplificationAdditions COMMAND SimplificationAdditions)

project(SimplificationMultiplications)
add_executable(SimplificationMultiplications Core/Simplification/Multiplication.cpp)
target_link_libraries(SimplificationMultiplications math)
add_test(NAME SimplificationMultiplications COMMAND SimplificationMultiplications)

project(SimplificationDivisions)
add_executable(SimplificationDivisions Core/Simplification/Division.cpp)
target_link_libraries(SimplificationDivisions math)
add_test(NAME SimplificationDivisions COMMAND SimplificationDivisions)

project(SimplificationSubtractions)
add_executable(SimplificationSubtractions Core/Simplification/Subtraction.cpp)
target_link_libraries(SimplificationSubtractions math)
add_test(NAME SimplificationSubtractions COMMAND SimplificationSubtractions)

project(ExpandTests)
add_executable(ExpandTests Core/Expand/Expand.cpp)
target_link_libraries(ExpandTests math)
add_test(NAME ExpandTests COMMAND ExpandTests)

project(TrigonometryTests)
add_executable(TrigonometryTests Core/Trigonometry/Trigonometry.cpp)
target_link_libraries(TrigonometryTests math)
add_test(NAME TrigonometryTests COMMAND TrigonometryTests)

project(ExponentialTests)
add_executable(ExponentialTests Core/Exponential/Exponential.cpp)
target_link_libraries(ExponentialTests math)
add_test(NAME ExponentialTests COMMAND ExponentialTests)
