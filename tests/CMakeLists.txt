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

project(PolynomialTests)
add_executable(PolynomialTests Core/Polynomial/Polynomial.cpp)
target_link_libraries(PolynomialTests math)
add_test(NAME PolynomialTests COMMAND PolynomialTests)


project(ExpandTests)
add_executable(ExpandTests Core/Expand/Expand.cpp)
target_link_libraries(ExpandTests math)
add_test(NAME ExpandTests COMMAND ExpandTests)

project(SimplificationRationals)
add_executable(SimplificationRationals Core/Reduce/Rationals.cpp)
target_link_libraries(SimplificationRationals math)
add_test(NAME SimplificationRationals COMMAND SimplificationRationals)

project(SimplificationPowers)
add_executable(SimplificationPowers Core/Reduce/Power.cpp)
target_link_libraries(SimplificationPowers math)
add_test(NAME SimplificationPowers COMMAND SimplificationPowers)

project(SimplificationAdditions)
add_executable(SimplificationAdditions Core/Reduce/Addition.cpp)
target_link_libraries(SimplificationAdditions math)
add_test(NAME SimplificationAdditions COMMAND SimplificationAdditions)

project(SimplificationMultiplications)
add_executable(SimplificationMultiplications Core/Reduce/Multiplication.cpp)
target_link_libraries(SimplificationMultiplications math)
add_test(NAME SimplificationMultiplications COMMAND SimplificationMultiplications)

project(SimplificationDivisions)
add_executable(SimplificationDivisions Core/Reduce/Division.cpp)
target_link_libraries(SimplificationDivisions math)
add_test(NAME SimplificationDivisions COMMAND SimplificationDivisions)

project(SimplificationSubtractions)
add_executable(SimplificationSubtractions Core/Reduce/Subtraction.cpp)
target_link_libraries(SimplificationSubtractions math)
add_test(NAME SimplificationSubtractions COMMAND SimplificationSubtractions)

