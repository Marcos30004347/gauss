#ifndef MATH_ALGEBRA_POLYNOMIALS_H
#define MATH_ALGEBRA_POLYNOMIALS_H

#include "Core/Algebra/Algebra.hpp"

#include <vector>

namespace polynomial {
	

/**
 * Return true if u is a generalized monomial expression in the
 * generalized variable v or in the set of variables v
 */
bool isGeneralMonomial(ast::AST* u, ast::AST* v);

/**
 * Return true if u is a generalized polynomial expression in the
 * generalized variable v or in the set of variables v
 */
bool isGerenalPolynomial(ast::AST* u, ast::AST* v);

/**
 * Return the variable parts, so if thos variables are removed
 * from u, only rational coefficients are left.
 */
ast::AST* variables(ast::AST* u);

/**
 * Return a list with the coefficient and variable parts,
 * the variables symbols are given in the set S.
 */
ast::AST* coeffVarMonomial(ast::AST* u, ast::AST* S);

/**
 * Return u with the terms of the variables S collected.
 * EXAMPLE: collectTerms(ax + bc + c + d, {x}) -> (a + b)x + c + d;
 */
ast::AST* collectTerms(ast::AST* u, ast::AST* S);

/**
 * Expands the expression u and all its
 * childs recursivelly
 */
ast::AST* algebraicExpand(ast::AST* u);

/**
 * Expands only the expression u and lets they
 * operands as they are
 */
ast::AST* algebraicExpandRoot(ast::AST* u);

/**
 * Returns the biggest degree of x in u, by default.
 * The degree of 0 monomial is -infinity, so if no 
 * x is found in u, -infinity will be returned
 */
ast::AST* degreeGPE(ast::AST* u, ast::AST* x);

/**
 * Return the sum of the coefficients of x^j in u
 * 
 * EXAMPLE: coefficientGPE(ax^2 + bx^2, x, 2) = a + b
 */
ast::AST* coefficientGPE(ast::AST* u, ast::AST* x, ast::AST* j);

/**
 * Returns the coefficient of x with the biggest degree in u.
 */
ast::AST* leadingCoefficientGPE(ast::AST* u, ast::AST* x);

/**
 * Divide the polynomial u by the polynomial v using the x variable.
 * The result is a pair where the first members is the quotient ant 
 * the second member is the remainder
 */
ast::AST* divideGPE(ast::AST* u, ast::AST* v, ast::AST* x);

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
ast::AST* extendedEuclideanAlgGPE(ast::AST* u, ast::AST* v, ast::AST* x);

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

//
// Multivariable
//

// TODO TEST
ast::AST* polynomialContent(ast::AST* u, ast::AST* x, ast::AST* R, ast::AST* K);

ast::AST* pseudoDivision(ast::AST* u, ast::AST* v, ast::AST* x);
ast::AST* pseudoQuotient(ast::AST* u, ast::AST* v, ast::AST* x);
ast::AST* pseudoRemainder(ast::AST* u, ast::AST* v, ast::AST* x);

ast::AST* normalizePoly(ast::AST* u, ast::AST* L, ast::AST* K);

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
