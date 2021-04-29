#ifndef MATH_ALGEBRA_POLYNOMIALS_H
#define MATH_ALGEBRA_POLYNOMIALS_H

#include "Core/Algebra/Algebra.hpp"

#include <vector>
#include <utility> 

namespace polynomial {
	
/**
 * Return the variable parts, so if thos variables are removed
 * from u, only rational coefficients are left.
 */
std::vector<ast::AST*> variables(ast::AST* u);

/**
 * Return if u is a polynomial in the variables defined
 * in vars.
 */
bool isPolynomialGPE(ast::AST* u, std::vector<ast::AST*> vars);

/**
 * Returns the biggest degree of x in u, by default.
 * The degree of 0 monomial is -infinity, so if no 
 * x is found in u, -infinity will be returned
 */
ast::AST* degreeGPE(ast::AST* u, ast::AST* x);

/**
 * Return the sum of the coefficients of x^j in u
 * 
 * EXAMPLE: coefficientGPE(ax^2 + bx^2, x^2) = a + b
 */
ast::AST* coefficientGPE(ast::AST* u, ast::AST* x);

/**
 * Returns the coefficient of x with the biggest degree in u.
 */
ast::AST* leadingCoefficientGPE(ast::AST* u, ast::AST* x);

/**
 * Divide the polynomial u by the polynomial v using the x variable.
 * The result is a pair where the first members is the quotient ant 
 * the second member is the remainder
 */
std::pair<ast::AST*, ast::AST*> divideGPE(ast::AST* u, ast::AST* v, ast::AST* x);

/**
 * Return the quotient of the division of u by v using the 
 * x variable.
 */
ast::AST* quotientGPE(ast::AST* u, ast::AST* v, ast::AST* x);

/**
 * Return the remainder of the division of u by v using the 
 * x variable.
 */
ast::AST* remainderGPE(ast::AST* u, ast::AST* v, ast::AST* x);

ast::AST* expandGPE(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* t);
}

#endif
