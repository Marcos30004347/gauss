cmake_minimum_required(VERSION 3.10)

include(CTest)


project(ASTTests)
add_executable(ASTTests Core/AbstractSyntaxTree/AST.cpp)
target_link_libraries(ASTTests math)
add_test(NAME ASTTests COMMAND ASTTests)

project(AlgebraTests)
add_executable(AlgebraTests Core/Algebra/Algebra.cpp)
target_link_libraries(AlgebraTests math)
add_test(NAME AlgebraTests COMMAND AlgebraTests)

project(SimplificationRationals)
add_executable(SimplificationRationals Core/Reduce/Rationals.cpp)
target_link_libraries(SimplificationRationals math)
add_test(NAME SimplificationRationals COMMAND SimplificationRationals)

project(SimplificationPowers)
add_executable(SimplificationPowers Core/Reduce/Power.cpp)
target_link_libraries(SimplificationPowers math)
add_test(NAME SimplificationPowers COMMAND SimplificationPowers)

