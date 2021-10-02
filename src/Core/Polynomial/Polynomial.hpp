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
ast::AST* degree(ast::AST* u, ast::AST* x);

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
ast::AST* algPolynomialDivisionAST(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algPolynomialRemainderAST(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algPolynomialQuotientAST(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algPolynomialGCDAST(ast::AST* u, ast::AST* v, ast::AST* x, ast::AST* p, ast::AST* a);
ast::AST* algMonicAST(ast::AST* u,ast::AST* x, ast::AST* p, ast::AST* a);

//
// Multivariable
//
/**
 * Calculate the content of the multivariable polynomial u using
 * x as main variable, R is a list containing all variables of the
 * u polynomial and K is either Z, if u in a polynomial in Z[R...],
 * or Q if u is a polynomial in Q[R...].
 */
ast::AST* polynomialContent(ast::AST* u, ast::AST* x, ast::AST* R, ast::AST* K);

/**
 * Calculate the sub resultant content of the multivariable polynomial
 * u using x as main variable, R is a list containing all variables of the
 * u polynomial and K is either Z, if u in a polynomial in Z[R...],
 * or Q if u is a polynomial in Q[R...]. Note that the list L should
 * be empty if u is univariate in x.
 */
ast::AST* polynomialContentSubResultant(ast::AST* u, ast::AST* x, ast::AST* R, ast::AST* K);

/**
 * Pseudo division of the multivariable polynomial u by the 
 * multivariable polynomial v using x.
 */
ast::AST* pseudoDivision(ast::AST* u, ast::AST* v, ast::AST* x);

/**
 * Pseudo quotient of the multivariable polynomial u 
 * divided by the multivariable polynomial v using x.
 */
ast::AST* pseudoQuotient(ast::AST* u, ast::AST* v, ast::AST* x);

/**
 * Pseudo remainder of the multivariable polynomial u 
 * divided by the multivariable polynomial v using x .
 */
ast::AST* pseudoRemainder(ast::AST* u, ast::AST* v, ast::AST* x);

/**
 * Normalize the multivariable polynomial u in the field K
 * with variables defined inside the list L.
 */
ast::AST* normalizePoly(ast::AST* u, ast::AST* L, ast::AST* K);


/**
 * Return the GCD between the multivariable polynomials u and v
 * with variables defined in the list L in the field K
 */
ast::AST* mvPolyGCD(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);


/**
 * Return the GCD between the multivariable polynomials u and v
 * with variables defined in the list L in the field K using sub
 * resultant.
 */
ast::AST* mvSubResultantGCD(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);

/**
 * Return the GCD between the multivariable polynomials u and v
 * with variables defined in the list L in the field K using
 * the sub resultant content
 */
ast::AST* subResultantGCDRec(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);

/**
 * Recursive polynomial divisition between the multivariable 
 * polynomials u and v with variables defined in the list L in the field K
 */
ast::AST* recPolyDiv(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);

/**
 * Recursive polynomial quotient between the division of the multivariable 
 * polynomials u and v with variables defined in the list L in the field K
 */
ast::AST* recQuotient(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);

/**
 * Recursive polynomial remainder between the division of the multivariable 
 * polynomials u and v with variables defined in the list L in the field K
 */
ast::AST* recRemainder(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);

/**
 * Return the leading monomial of the multivariable polynomial
 * u with variables defined on the list L
 */
ast::AST* leadingMonomial(ast::AST* u, ast::AST* L);

/**
 * Return the quotient and remainder of the divisision between
 * the multivariable monomials u and v with variables defined
 * in the list L.
 */
ast::AST* monomialPolyDiv(ast::AST* u, ast::AST* v, ast::AST* L);

/**
 * Return the remainder of the divisision between
 * the multivariable monomials u and v with variables defined
 * in the list L.
 */
ast::AST* monomialPolyRem(ast::AST* u, ast::AST* v, ast::AST* L);

/**
 * Return the quotient of the divisision between
 * the multivariable monomials u and v with variables defined
 * in the list L.
 */
ast::AST* monomialPolyQuo(ast::AST* u, ast::AST* v, ast::AST* L);

/**
 * Return u in terms of v and replace the v part with t,
 * both u and v are multivariable polynomials with variables
 * defined in L.
 * u will be represented as:
 * 	 u = d[k]*v^k + d[k-1]*v^(k-1) + ... + d[0]
 * where d are also polynomials in Q[L...]
 * 
 * 
 * monomialBasedPolyExpansion can be used to rewrite the default operations
 * like degree and coefficientGPE for multivariable polynomials
 * reducing them to single variable polynomials
 * and them using degree and coefficientGPE to make the query
 */
ast::AST* monomialBasedPolyExpansion(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* t);

// TODO: refactor this to use polynomialContent
ast::AST* cont(ast::AST* u, ast::AST* x);

ast::AST* pdiv(ast::AST* u, ast::AST* v, ast::AST* x);


ast::AST* mulPoly(ast::AST* p1, ast::AST* p2);
ast::AST* addPoly(ast::AST* p1, ast::AST* p2);
ast::AST* subPoly(ast::AST* p1, ast::AST* p2);
}

#endif
