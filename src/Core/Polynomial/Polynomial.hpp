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

/**
 * Expand u in terms of v, and substitute the ocurrences of v by t.
 */
ast::AST* expandGPE(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* t);

/**
 * Calculate the greatest commom divisor betwwen gpes u and v 
 * using x;
 */
ast::AST* gcdGPE(ast::AST* u, ast::AST* v, ast::AST* x);

/**
 * Extended Euclidean Algorithm for gpes.
 */
std::vector<ast::AST*> extendedEuclideanAlgGPE(ast::AST* u, ast::AST* v, ast::AST* x);

// 
// Those functions use expression p as side relation 
//
ast::AST* algCoeffSimp(ast::AST* u, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algMulInverseAST(ast::AST* v, ast::AST* p, ast::AST* a);
ast::AST* algDivideAST(ast::AST* u, ast::AST* v, ast::AST* p, ast::AST* a);
std::vector<ast::AST*> algPolynomialDivisionAST(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algPolynomialRemainderAST(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algPolynomialQuotientAST(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algPolynomialGCDAST(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algMonicAST(ast::AST* u,ast::AST* x, ast::AST* p, ast::AST* a);

// TEST THIS
//
// Multivariable
//

// TODO TEST
ast::AST* polynomialContent(ast::AST* u, ast::AST* x, ast::AST* R, ast::AST* K);

// TODO TEST
// u = 5x⁴y³ + 3xy + 2,    v = 2x³y + 2x + 3
// p1 = 5y³x, 			s1 = −10x²y³ − 15xy³ + 6xy² + 4y
// termitates with σ = 1
// r = -20x²y⁴ - 30xy⁴ + 12xy³ + 8y²
ast::AST* pseudoDivision(ast::AST* u, ast::AST* v, ast::AST* x);
ast::AST* pseudoQuotient(ast::AST* u, ast::AST* v, ast::AST* x);
ast::AST* pseudoRemainder(ast::AST* u, ast::AST* v, ast::AST* x);

// TODO TEST
// normalizePoly((2 y + 1) x + (6 y + 3), [x,y], Q) = (y + 1/2)x + (3 y + 3/2).
ast::AST* normalizePoly(ast::AST* u, ast::AST* L, ast::AST* K);

// TODO TEST
ast::AST* mvPolyGCD(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);

ast::AST* recPolyDiv(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);
ast::AST* recQuotient(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);
ast::AST* recRemainder(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);

ast::AST* leadingMonomial(ast::AST* u, ast::AST* L);

ast::AST* monomialPolyDiv(ast::AST* u, ast::AST* v, ast::AST* L);
ast::AST* monomialPolyRem(ast::AST* u, ast::AST* v, ast::AST* L);
ast::AST* monomialPolyQuo(ast::AST* u, ast::AST* v, ast::AST* L);

// TODO TEST
// monomialPolyExpansion(a^2*b + 2*a*b^2 + b^3 + 2*a + 2*b + 3, a+b, [a, b], t) -> b*t^2 + 2*t + 3
ast::AST* monomialPolyExpansion(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* t);

}

#endif
