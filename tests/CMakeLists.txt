cmake_minimum_required(VERSION 3.10)

include(CTest)

project(IntegerTests)
add_executable(IntegerTests MathSystem/Algebra/Integer.cpp)
target_link_libraries(IntegerTests math)
target_include_directories(IntegerTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME IntegerTests COMMAND IntegerTests)

project(IntTests)
add_executable(IntTests MathSystem/Algebra/Int.cpp)
target_link_libraries(IntTests math)
target_include_directories(IntTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME IntTests COMMAND IntTests)

project(ExpressionTests)
add_executable(ExpressionTests MathSystem/Algebra/Expression.cpp)
target_link_libraries(ExpressionTests math)
target_include_directories(ExpressionTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME ExpressionTests COMMAND ExpressionTests)

project(MatrixTests)
add_executable(MatrixTests MathSystem/Algebra/Matrix.cpp)
target_link_libraries(MatrixTests math)
target_include_directories(MatrixTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME MatrixTests COMMAND MatrixTests)

project(LexerTests)
add_executable(LexerTests MathSystem/Algebra/Lexer.cpp)
target_link_libraries(LexerTests math)
target_include_directories(LexerTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME LexerTests COMMAND LexerTests)

project(ParserTests)
add_executable(ParserTests MathSystem/Algebra/Parser.cpp)
target_link_libraries(ParserTests math)
target_include_directories(ParserTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME ParserTests COMMAND ParserTests)


# project(AlgebraTests)
# add_executable(AlgebraTests MathSystem/Algebra/Algebra.cpp)
# target_link_libraries(AlgebraTests math)
# target_include_directories(AlgebraTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME AlgebraTests COMMAND AlgebraTests)

# project(ListTests)
# add_executable(ListTests MathSystem/Algebra/List.cpp)
# target_link_libraries(ListTests math)
# target_include_directories(ListTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME ListTests COMMAND ListTests)

# project(SetTests)
# add_executable(SetTests MathSystem/Algebra/Set.cpp)
# target_link_libraries(SetTests math)
# target_include_directories(SetTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME SetTests COMMAND SetTests)

# project(RationalTests)
# add_executable(RationalTests MathSystem/Rational/Rational.cpp)
# target_link_libraries(RationalTests math)
# target_include_directories(RationalTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME RationalTests COMMAND RationalTests)

# project(SimplificationRationals)
# add_executable(SimplificationRationals MathSystem/Simplification/Rationals.cpp)
# target_link_libraries(SimplificationRationals math)
# target_include_directories(SimplificationRationals PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME SimplificationRationals COMMAND SimplificationRationals)

# project(SimplificationPowers)
# add_executable(SimplificationPowers MathSystem/Simplification/Power.cpp)
# target_link_libraries(SimplificationPowers math)
# target_include_directories(SimplificationPowers PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME SimplificationPowers COMMAND SimplificationPowers)

# project(SimplificationAdditions)
# add_executable(SimplificationAdditions MathSystem/Simplification/Addition.cpp)
# target_link_libraries(SimplificationAdditions math)
# target_include_directories(SimplificationAdditions PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME SimplificationAdditions COMMAND SimplificationAdditions)

# project(SimplificationMultiplications)
# add_executable(SimplificationMultiplications MathSystem/Simplification/Multiplication.cpp)
# target_link_libraries(SimplificationMultiplications math)
# target_include_directories(SimplificationMultiplications PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME SimplificationMultiplications COMMAND SimplificationMultiplications)

# project(SimplificationDivisions)
# add_executable(SimplificationDivisions MathSystem/Simplification/Division.cpp)
# target_link_libraries(SimplificationDivisions math)
# target_include_directories(SimplificationDivisions PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME SimplificationDivisions COMMAND SimplificationDivisions)

# project(SimplificationSubtractions)
# add_executable(SimplificationSubtractions MathSystem/Simplification/Subtraction.cpp)
# target_link_libraries(SimplificationSubtractions math)
# target_include_directories(SimplificationSubtractions PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# add_test(NAME SimplificationSubtractions COMMAND SimplificationSubtractions)

project(PolynomialTests)
add_executable(PolynomialTests MathSystem/Polynomial/Polynomial.cpp)
target_link_libraries(PolynomialTests math)
target_include_directories(PolynomialTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME PolynomialTests COMMAND PolynomialTests)

project(ResultantTests)
add_executable(ResultantTests MathSystem/Polynomial/Resultant.cpp)
target_link_libraries(ResultantTests math)
target_include_directories(ResultantTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME ResultantTests COMMAND ResultantTests)

project(GaloisFieldTests)
add_executable(GaloisFieldTests MathSystem/GaloisField/GaloisField.cpp)
target_link_libraries(GaloisFieldTests math)
target_include_directories(GaloisFieldTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME GaloisFieldTests COMMAND GaloisFieldTests)

# # project(TrigonometryTests)
# # add_executable(TrigonometryTests MathSystem/Trigonometry/Trigonometry.cpp)
# # target_link_libraries(TrigonometryTests math)
# # add_test(NAME TrigonometryTests COMMAND TrigonometryTests)

# # project(ExponentialTests)
# # add_executable(ExponentialTests MathSystem/Exponential/Exponential.cpp)
# # target_link_libraries(ExponentialTests math)
# # add_test(NAME ExponentialTests COMMAND ExponentialTests)

project(CalculusTests)
add_executable(CalculusTests MathSystem/Calculus/Calculus.cpp)
target_link_libraries(CalculusTests math)
target_include_directories(CalculusTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME CalculusTests COMMAND CalculusTests)

project(FactorizationSquareFreeTests)
add_executable(FactorizationSquareFreeTests MathSystem/Factorization/SquareFree.cpp)
target_link_libraries(FactorizationSquareFreeTests math)
target_include_directories(FactorizationSquareFreeTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME FactorizationSquareFreeTests COMMAND FactorizationSquareFreeTests)

project(FactorizationHenselTests)
add_executable(FactorizationHenselTests MathSystem/Factorization/Hensel.cpp)
target_link_libraries(FactorizationHenselTests math)
target_include_directories(FactorizationHenselTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME FactorizationHenselTests COMMAND FactorizationHenselTests)

###############################################################################################
# project(FactorizationBerlekampTests)																											  #
# add_executable(FactorizationBerlekampTests MathSystem/Factorization/Berlekamp.cpp)				  #
# target_link_libraries(FactorizationBerlekampTests math)																		  #
# target_include_directories(FactorizationBerlekampTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}) #
# add_test(NAME FactorizationBerlekampTests COMMAND FactorizationBerlekampTests)						  #
###############################################################################################

project(FactorizationZassenhausTests)
add_executable(FactorizationZassenhausTests MathSystem/Factorization/Zassenhaus.cpp)
target_link_libraries(FactorizationZassenhausTests math)
target_include_directories(FactorizationZassenhausTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME FactorizationZassenhausTests COMMAND FactorizationZassenhausTests)

project(FactorizationWangTests)
add_executable(FactorizationWangTests MathSystem/Factorization/Wang.cpp)
target_link_libraries(FactorizationWangTests math)
target_include_directories(FactorizationWangTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME FactorizationWangTests COMMAND FactorizationWangTests)


project(SymbolTableTests)
add_executable(SymbolTableTests TypeSystem/SymbolTable.cpp)
target_link_libraries(SymbolTableTests math)
target_include_directories(SymbolTableTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME SymbolTableTests COMMAND SymbolTableTests)

project(ASTTests)
		add_executable(ASTTests TypeSystem/AST.cpp)
target_link_libraries(ASTTests math)
target_include_directories(ASTTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME ASTTests COMMAND ASTTests)
