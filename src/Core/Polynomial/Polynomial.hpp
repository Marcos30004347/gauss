#ifndef MATH_ALGEBRA_POLYNOMIALS_H
#define MATH_ALGEBRA_POLYNOMIALS_H

#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"

#include <vector>

namespace polynomial {

/**
 * Return true if u is a generalized monomial expression in the
 * generalized variable v or in the set of variables v
 */
bool isGeneralMonomial(ast::Expr u, ast::Expr v);

/**
 * Return true if u is a generalized polynomial expression in the
 * generalized variable v or in the set of variables v
 */
bool isGerenalPolynomial(ast::Expr u, ast::Expr v);

/**
 * Return all the variables in u ordered by degree,
 * function calls and symbols are considered symbols.
 */
ast::Expr variables(ast::Expr u);

/**
 * Return a list with the coeff and variable parts,
 * the variables symbols are given in the set S.
 */
ast::Expr coeffVarMonomial(ast::Expr u, ast::Expr S);

/**
 * Return u with the terms of the variables S collected.
 * EXAMPLE: collectTerms(ax + bc + c + d, {x}) -> (a + b)x + c + d;
 */
ast::Expr collectTerms(ast::Expr u, ast::Expr S);

/**
 * Expands the expression u and all its
 * childs recursivelly
 */
ast::Expr algebraicExpand(ast::Expr u);

/**
 * Expands only the expression u and lets they
 * operands as they are
 */
ast::Expr algebraicExpandRoot(ast::Expr u);

/**
 * Returns the biggest degree of x in u, by default.
 * The degree of 0 is -infinity, so if
 * x is not found in u, -infinity will be returned
 */
ast::Expr degree(ast::Expr u, ast::Expr x);

/**
 * Return the sum of the coefficients of x^j in u
 *
 * EXAMPLE: coeff(ax^2 + bx^2, x, 2) = a + b
 */
ast::Expr coeff(ast::Expr u, ast::Expr x, ast::Expr j);

/**
 * Returns the coeff of x with the biggest degree in u.
 */
ast::Expr leadCoeff(ast::Expr u, ast::Expr x);

/**
 * Divide the polynomial u by the polynomial v using the x variable.
 * The result is a pair where the first members is the quotient ant
 * the second member is the remainder
 */
ast::Expr divideGPE(ast::Expr u, ast::Expr v, ast::Expr x);

/**
 * Return the quotient of the division of u by v using the
 * x variable.
 */
ast::Expr quotientGPE(ast::Expr u, ast::Expr v, ast::Expr x);

/**
 * Return the remainder of the division of u by v using the
 * x variable.
 */
ast::Expr remainderGPE(ast::Expr u, ast::Expr v, ast::Expr x);

/**
 * Expand u in terms of v, and substitute the ocurrences of v by t.
 */
ast::Expr expandGPE(ast::Expr u, ast::Expr v, ast::Expr x, ast::Expr t);

/**
 * Calculate the greatest commom divisor betwwen gpes u and v
 * using x;
 */
ast::Expr gcdGPE(ast::Expr u, ast::Expr v, ast::Expr x);

/**
 * Extended Euclidean Algorithm for gpes.
 */
ast::Expr extendedEuclideanAlgGPE(ast::Expr u, ast::Expr v, ast::Expr x);

/**
 * @brief Computes the content of a multivariate polynomial u(L...) in K[L...]
 *
 * @param u A polynomial in K[L...]
 *
 * @param L The list of variables of u
 *
 * @param K The field that us is defined, either Z or Q
 *
 * @return The content of the polynomial u
 */
ast::Expr cont(ast::Expr u, ast::Expr L, ast::Expr K);

/**
 * @brief Computes the primitive part of a multivariate polynomial u(L...) in
 * K[L...]
 *
 * @param u A polynomial in K[L...]
 *
 * @param L The list of variables of u
 *
 * @param K The field that us is defined, either Z or Q
 *
 * @return The content of the polynomial u
 */
ast::Expr pp(ast::Expr u, ast::Expr L, ast::Expr K);

/**
 * @brief Computes the primitive part of a multivariate polynomial u(L...) in
 * K[L...]
 *
 * @param u A polynomial in K[L...]
 *
 * @param c The content of u
 *
 * @param L The list of variables of u
 *
 * @param K The field that us is defined, either Z or Q
 *
 * @return The content of the polynomial u
 */
ast::Expr pp(ast::Expr u, ast::Expr c, ast::Expr L, ast::Expr K);

/**
 * Pseudo division of the multivariable polynomial u by the
 * multivariable polynomial v using x.
 */
ast::Expr pseudoDivision(ast::Expr u, ast::Expr v, ast::Expr x);

/**
 * Pseudo quotient of the multivariable polynomial u
 * divided by the multivariable polynomial v using x.
 */
ast::Expr pseudoQuotient(ast::Expr u, ast::Expr v, ast::Expr x);

/**
 * Pseudo remainder of the multivariable polynomial u
 * divided by the multivariable polynomial v using x .
 */
ast::Expr pseudoRemainder(ast::Expr u, ast::Expr v, ast::Expr x);

/**
 * Normalize the multivariable polynomial u in the field K
 * with variables defined inside the list L.
 */
ast::Expr normalizePoly(ast::Expr u, ast::Expr L, ast::Expr K);

/**
 * Return the GCD between the multivariable polynomials u and v
 * with variables defined in the list L in the field K
 */
ast::Expr mvPolyGCD(ast::Expr u, ast::Expr v, ast::Expr L, ast::Expr K);

/**
 * Return the GCD between the multivariable polynomials u and v
 * with variables defined in the list L in the field K using sub
 * resultant.
 */
ast::Expr mvSubResultantGCD(ast::Expr u, ast::Expr v, ast::Expr L, ast::Expr K);

/**
 * Return the GCD between the multivariable polynomials u and v
 * with variables defined in the list L in the field K using
 * the sub resultant content
 */
ast::Expr subResultantGCDRec(ast::Expr u, ast::Expr v, ast::Expr L,
                             ast::Expr K);

/**
 * Recursive polynomial divisition between the multivariable
 * polynomials u and v with variables defined in the list L in the field K
 */
ast::Expr recPolyDiv(ast::Expr u, ast::Expr v, ast::Expr L, ast::Expr K);

/**
 * Recursive polynomial quotient between the division of the multivariable
 * polynomials u and v with variables defined in the list L in the field K
 */
ast::Expr recQuotient(ast::Expr u, ast::Expr v, ast::Expr L, ast::Expr K);

/**
 * Recursive polynomial remainder between the division of the multivariable
 * polynomials u and v with variables defined in the list L in the field K
 */
ast::Expr recRemainder(ast::Expr u, ast::Expr v, ast::Expr L, ast::Expr K);

/**
 * Return the leading monomial of the multivariable polynomial
 * u with variables defined on the list L
 */
ast::Expr leadMonomial(ast::Expr u, ast::Expr L);

/**
 * Return the quotient and remainder of the divisision between
 * the multivariable monomials u and v with variables defined
 * in the list L.
 */
ast::Expr monomialPolyDiv(ast::Expr u, ast::Expr v, ast::Expr L);

/**
 * Return the remainder of the divisision between
 * the multivariable monomials u and v with variables defined
 * in the list L.
 */
ast::Expr monomialPolyRem(ast::Expr u, ast::Expr v, ast::Expr L);

/**
 * Return the quotient of the divisision between
 * the multivariable monomials u and v with variables defined
 * in the list L.
 */
ast::Expr monomialPolyQuo(ast::Expr u, ast::Expr v, ast::Expr L);

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
 * like degree and coeff for multivariable polynomials
 * reducing them to single variable polynomials
 * and them using degree and coeff to make the query
 */
ast::Expr monomialBasedPolyExpansion(ast::Expr u, ast::Expr v, ast::Expr L,
                                     ast::Expr t);

// TODO: refactor this to use polynomialContent
ast::Expr cont(ast::Expr u, ast::Expr x);

ast::Expr pdiv(ast::Expr u, ast::Expr v, ast::Expr x);

ast::Expr collect(ast::Expr &&u, ast::Expr &&L);
ast::Expr collect(ast::Expr &u, ast::Expr &L);
ast::Expr collect(ast::Expr &&u, ast::Expr &L);
ast::Expr collect(ast::Expr &u, ast::Expr &&L);

ast::Expr mulPoly(ast::Expr p1, ast::Expr p2);
ast::Expr addPoly(ast::Expr p1, ast::Expr p2);
ast::Expr subPoly(ast::Expr p1, ast::Expr p2);
ast::Expr raisePoly(ast::Expr f, long n);

bool isZeroColPoly(ast::Expr& u);

ast::Expr mulColPoly(ast::Expr& p1, ast::Expr& p2);
ast::Expr addColPoly(ast::Expr& p1, ast::Expr& p2);
ast::Expr subColPoly(ast::Expr& p1, ast::Expr& p2);
ast::Expr powColPoly(ast::Expr& p1, Int n);

ast::Expr mulColPoly(ast::Expr&& p1, ast::Expr&& p2);
ast::Expr addColPoly(ast::Expr&& p1, ast::Expr&& p2);
ast::Expr subColPoly(ast::Expr&& p1, ast::Expr&& p2);
ast::Expr powColPoly(ast::Expr&& p1, Int n);
ast::Expr recColPolyDiv(ast::Expr u, ast::Expr v, ast::Expr L, ast::Expr K);
} // namespace polynomial

#endif
