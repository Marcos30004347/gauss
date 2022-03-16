#ifndef MATH_ALGEBRA_POLYNOMIALS_H
#define MATH_ALGEBRA_POLYNOMIALS_H

#include "gauss/Algebra/Expression.hpp"
//#include "gauss/Algebra/Algebra.hpp"

#include <vector>

namespace poly {
// /**
//  * Return true if u is a generalized monomial expression in the
//  * generalized variable v or in the set of variables v
//  */
// bool isGeneralMonomial(alg::expr u, alg::expr v);

// /**
//  * Return true if u is a generalized polynomial expression in the
//  * generalized variable v or in the set of variables v
//  */
// bool isGerenalPolynomial(alg::expr u, alg::expr v);

// /**
//  * Return all the variables in u ordered by degree,
//  * function calls and symbols are considered symbols.
//  */
// alg::set variables(alg::expr u);

// /**
//  * Return a list with the coeff and variable parts,
//  * the variables symbols are given in the set S.
//  */
// alg::list coeffVarMonomial(alg::expr u, alg::set S);

// /**
//  * Return u with the terms of the variables S collected.
//  * EXAMPLE: collectTerms(ax + bc + c + d, {x}) -> (a + b)x + c + d;
//  */
// alg::expr collectTerms(alg::expr u, alg::set S);

// /**
//  * Expands only the expression u and lets they
//  * operands as they are
//  */
// alg::expr algebraicExpandRoot(alg::expr u);

// /**
//  * Returns the biggest degree of x in u, by default.
//  * The degree of 0 is -infinity, so if
//  * x is not found in u, -infinity will be returned
//  */
// alg::expr degree(alg::expr u, alg::expr x);

// /**
//  * Return the sum of the coefficients of x^j in u
//  *
//  * EXAMPLE: coeff(ax^2 + bx^2, x, 2) = a + b
//  */
// alg::expr coeff(alg::expr u, alg::expr x, alg::expr j);

// /**
//  * Returns the coeff of x with the biggest degree in u.
//  */
// alg::expr leadCoeff(alg::expr u, alg::expr x);

// /**
//  * Divide the polynomial u by the polynomial v using the x variable.
//  * The result is a pair where the first members is the quotient ant
//  * the second member is the remainder
//  */
// alg::expr divideGPE(alg::expr u, alg::expr v, alg::expr x);

// /**
//  * Return the quotient of the division of u by v using the
//  * x variable.
//  */
// alg::expr quotientGPE(alg::expr u, alg::expr v, alg::expr x);

// /**
//  * Return the remainder of the division of u by v using the
//  * x variable.
//  */
// alg::expr remainderGPE(alg::expr u, alg::expr v, alg::expr x);

// /**
//  * Expand u in terms of v, and substitute the ocurrences of v by t.
//  */
// alg::expr expandGPE(alg::expr u, alg::expr v, alg::expr x, alg::expr t);

// /**
//  * Calculate the greatest commom divisor betwwen gpes u and v
//  * using x;
//  */
// alg::expr gcdGPE(alg::expr u, alg::expr v, alg::expr x);

// /**
//  * Extended Euclidean Algorithm for gpes.
//  */
// alg::expr extendedEuclideanAlgGPE(alg::expr u, alg::expr v, alg::expr x);

// /**
//  * @brief Computes the content of a multivariate polynomial u(L...) in K[L...]
//  *
//  * @param u A polynomial in K[L...]
//  *
//  * @param L The list of variables of u
//  *
//  * @param K The field that us is defined, either Z or Q
//  *
//  * @return The content of the polynomial u
//  */
// alg::expr cont(alg::expr u, alg::expr L, alg::expr K);

// /**
//  * @brief Computes the primitive part of a multivariate polynomial u(L...) in
//  * K[L...]
//  *
//  * @param u A polynomial in K[L...]
//  *
//  * @param L The list of variables of u
//  *
//  * @param K The field that us is defined, either Z or Q
//  *
//  * @return The content of the polynomial u
//  */
// alg::expr pp(alg::expr u, alg::expr L, alg::expr K);

// /**
//  * @brief Computes the primitive part of a multivariate polynomial u(L...) in
//  * K[L...]
//  *
//  * @param u A polynomial in K[L...]
//  *
//  * @param c The content of u
//  *
//  * @param L The list of variables of u
//  *
//  * @param K The field that us is defined, either Z or Q
//  *
//  * @return The content of the polynomial u
//  */
// alg::expr pp(alg::expr u, alg::expr c, alg::expr L, alg::expr K);

// /**
//  * Pseudo division of the multivariable polynomial u by the
//  * multivariable polynomial v using x.
//  */
// alg::expr pseudoDivision(alg::expr u, alg::expr v, alg::expr x);

// /**
//  * Pseudo quotient of the multivariable polynomial u
//  * divided by the multivariable polynomial v using x.
//  */
// alg::expr pseudoQuotient(alg::expr u, alg::expr v, alg::expr x);

// /**
//  * Pseudo remainder of the multivariable polynomial u
//  * divided by the multivariable polynomial v using x .
//  */
// alg::expr pseudoRemainder(alg::expr u, alg::expr v, alg::expr x);

// /**
//  * Normalize the multivariable polynomial u in the field K
//  * with variables defined inside the list L.
//  */
// alg::expr normalizePoly(alg::expr u, alg::expr L, alg::expr K);

// /**
//  * Return the GCD between the multivariable polynomials u and v
//  * with variables defined in the list L in the field K
//  */
// alg::expr mvPolyGCD(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

// /**
//  * Return the GCD between the multivariable polynomials u and v
//  * with variables defined in the list L in the field K using sub
//  * resultant.
//  */
// alg::expr mvSubResultantGCD(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

// /**
//  * Return the GCD between the multivariable polynomials u and v
//  * with variables defined in the list L in the field K using
//  * the sub resultant content
//  */
// alg::expr subResultantGCDRec(alg::expr u, alg::expr v, alg::expr L,
//                              alg::expr K);

// /**
//  * Recursive polynomial divisition between the multivariable
//  * polynomials u and v with variables defined in the list L in the field K
//  */
// alg::expr recPolyDiv(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

// /**
//  * Recursive polynomial quotient between the division of the multivariable
//  * polynomials u and v with variables defined in the list L in the field K
//  */
// alg::expr recQuotient(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

// /**
//  * Recursive polynomial remainder between the division of the multivariable
//  * polynomials u and v with variables defined in the list L in the field K
//  */
// alg::expr recRemainder(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

// /**
//  * Return the leading monomial of the multivariable polynomial
//  * u with variables defined on the list L
//  */
// alg::expr leadMonomial(alg::expr u, alg::expr L);

// /**
//  * Return the quotient and remainder of the divisision between
//  * the multivariable monomials u and v with variables defined
//  * in the list L.
//  */
// alg::expr monomialPolyDiv(alg::expr u, alg::expr v, alg::expr L);

// /**
//  * Return the remainder of the divisision between
//  * the multivariable monomials u and v with variables defined
//  * in the list L.
//  */
// alg::expr monomialPolyRem(alg::expr u, alg::expr v, alg::expr L);

// /**
//  * Return the quotient of the divisision between
//  * the multivariable monomials u and v with variables defined
//  * in the list L.
//  */
// alg::expr monomialPolyQuo(alg::expr u, alg::expr v, alg::expr L);

// /**
//  * Return u in terms of v and replace the v part with t,
//  * both u and v are multivariable polynomials with variables
//  * defined in L.
//  * u will be represented as:
//  * 	 u = d[k]*v^k + d[k-1]*v^(k-1) + ... + d[0]
//  * where d are also polynomials in Q[L...]
//  *
//  *
//  * monomialBasedPolyExpansion can be used to rewrite the default operations
//  * like degree and coeff for multivariable polynomials
//  * reducing them to single variable polynomials
//  * and them using degree and coeff to make the query
//  */
// alg::expr monomialBasedPolyExpansion(alg::expr u, alg::expr v, alg::expr L,
//                                      alg::expr t);

// // TODO: refactor this to use polynomialContent
// alg::expr cont(alg::expr u, alg::expr x);

// alg::expr pdiv(alg::expr u, alg::expr v, alg::expr x);


/**
 * @brief Get the coefficient of x^d on u.
 * @param[in] u A expression.
 * @param[in] x A symbol.
 * @param[in] d A integer degree.
 * @return The coefficient of x^d on u.
 */
alg::expr coeff(alg::expr u, alg::expr x, alg::expr d);

/**
 * @brief Get the greatest degree of u on v.
 * @param[in] u A expression.
 * @param[in] v A symbol.
 * @return The degree of v on u.
 */
alg::expr degree(alg::expr u, alg::expr v);

// alg::expr mulPoly(alg::expr p1, alg::expr p2);
// alg::expr addPoly(alg::expr p1, alg::expr p2);
// alg::expr subPoly(alg::expr p1, alg::expr p2);
// alg::expr raisePoly(alg::expr f, long n);


/**
 * @brief Computes the gcd between two multivariate polynomials.
 * @param u polynomial in K[L]
 * @param v polynomial in K[L]
 * @param L list of variables of u and v
 * @param K the field of u and v
 * @return the gcd between u and v
 */
// alg::expr gcdPoly(alg::expr u, alg::expr v, alg::expr L, alg::expr K);


/**
 * @brief Computes the gcd between two multivariate polynomials using heuristic methods.
 * @param u polynomial in K[L]
 * @param v polynomial in K[L]
 * @param L list of variables of u and v
 * @param K the field of u and v
 * @return the gcd between u and v or fail
 */
// alg::expr heuristicGcdPoly(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

/**
 * @brief put an expanded and reduced expression into a polynomial collected form.
 *
 * @param u expanded and reduced expression
 * @param L list of symbols in u
 * @return alg::expr
 */
alg::expr polyExpr(alg::expr &&u, alg::expr &&L);
alg::expr polyExpr(alg::expr &&u, alg::expr &L);
alg::expr polyExpr(alg::expr &u, alg::expr &L);
alg::expr polyExpr(alg::expr &u, alg::expr &&L);

/**
 * @brief Computes the gcd between two multivariate poly expressions using heuristic methods.
 * @param u polynomial expression in Z[L]
 * @param v polynomial expression in Z[L]
 * @param L list of variables of u and v
 * @param K the field of u and v, only Z is allowed
 * @return the gcd between u and v or fail
 */
alg::expr heuristicGcdPolyExpr(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

/**
 * @brief Remove the denominators from all the coefficients of the polynomial and return the lcm and the polynomial with the denominators removed.
 * @param u A polynomial in K[L]
 * @param L the list of symbols of u
 * @param K the field of u, Z or Q
 * @return a list with lcm and the polynomial 'v' such that u = v/lcm
 */
// alg::expr removeDenominatorsPoly(alg::expr u, alg::expr L, alg::expr K);


/**
 * @brief Remove the denominators from all the coefficients of a polynomial expression and return the lcm and the polynomial expression with the denominators removed.
 * @param u A polynomial expression in K[L]
 * @param L the list of symbols of u
 * @param K the field of u, Z or Q
 * @return a list with lcm and the polynomial expression 'v' such that u = v/lcm
 */
alg::expr removeDenominatorsPolyExpr(alg::expr u, alg::expr L, alg::expr K);


/**
 * @brief Return if a given expression is a zero
 *
 * @param u An polynomial expression
 * @return true if u is zero
 * @return false if u is not zero
 */
bool isZeroPolyExpr(alg::expr& u);


/**
 * @brief Return true if u is a constant equal to v
 *
 * @param u A polynomial expression
 * @param v A const Integer or Fraction
 * @return true if u is equal to v
 * @return false if u is not equal to v
 */
bool isConstantPolyExpr(alg::expr& u, alg::expr v);

/**
 * @brief Return true if u is a Integer of a fraction
 *
 * @param u A polynomial expression
 * @return true if u is a constant
 * @return false if u is not a constant
 */
bool isConstantPolyExpr(alg::expr& u);

/**
 * @brief Add two polynomial expressions, both have to
 * be transformed with the same set of symbols.
 *
 * @param p1 A polynomial expression
 * @param p2 A polynomial expression
 * @return alg::expr sum of p1 and p2
 */
alg::expr addPolyExpr(alg::expr&& p1, alg::expr&& p2);
alg::expr addPolyExpr(alg::expr& p1, alg::expr& p2);

/**
 * @brief Multiply two polynomial expressions, both have to
 * be transformed with the same set of symbols.
 *
 * @param p1 A polynomial expression
 * @param p2 A polynomial expression
 * @return alg::expr product of p1 and p2
 */
alg::expr mulPolyExpr(alg::expr& p1, alg::expr& p2);
alg::expr mulPolyExpr(alg::expr&& p1, alg::expr&& p2);

/**
 * @brief Subtract two polynomial expressions, both have to
 * be transformed with the same set of symbols.
 *
 * @param p1 A polynomial expression
 * @param p2 A polynomial expression
 * @return alg::expr difference of p1 and p2
 */
alg::expr subPolyExpr(alg::expr& p1, alg::expr& p2);
alg::expr subPolyExpr(alg::expr&& p1, alg::expr&& p2);

/**
 * @brief elevate a polynomial expression to the n'th power
 *
 * @param p1 A polynomial expression
 * @param p2 A polynomial expression
 * @return alg::expr p1^n
 */
alg::expr powPolyExpr(alg::expr& p1, Int n);
alg::expr powPolyExpr(alg::expr&& p1, Int n);

/**
 * @brief Divide two polynomial expressions, both have to
 * be transformed with the same set of symbols.
 *
 * @param u A polynomial expression
 * @param v A polynomial expression
 * @param L list of symbols in u and v
 * @param K Field where the division should take place, Z or Q are the valid fields
 * @return alg::expr list with the quotient and remainder of the division of p1 by p2
 */
alg::expr divPolyExpr(alg::expr&& u, alg::expr&& v, alg::expr& L, alg::expr& K);
alg::expr divPolyExpr(alg::expr& u, alg::expr& v, alg::expr& L, alg::expr& K);

/**
 * @brief Divide two polynomial expressions and return the quotient,
 * both have to be transformed with the same set of symbols.
 *
 * @param u A polynomial expression
 * @param v A polynomial expression
 * @param L list of symbols in u and v
 * @param K Field where the division should take place, Z or Q are the valid fields
 * @return alg::expr quotient of the division of p1 by p2
 */
alg::expr quoPolyExpr(alg::expr& u, alg::expr& v, alg::expr& L, alg::expr& K);

/**
 * @brief Divide two polynomial expressions and return the remainder,
 * both have to be transformed with the same set of symbols.
 *
 * @param u A polynomial expression
 * @param v A polynomial expression
 * @param L list of symbols in u and v
 * @param K Field where the division should take place, Z or Q are the valid fields
 * @return alg::expr remainder of the division of p1 by p2
 */
alg::expr remPolyExpr(alg::expr& u, alg::expr& v, alg::expr& L, alg::expr& K);

/**
 * @brief Pseudo Divide two polynomial expressions, both have to
 * be transformed with the same set of symbols.
 *
 * @param u A polynomial expression
 * @param v A polynomial expression
 * @param L list of symbols in u and v
 * @return alg::expr quotient and remainder of the division of p1 by p2
 */
alg::expr pseudoDivPolyExpr(alg::expr& u, alg::expr& v, alg::expr& L);
alg::expr pseudoDivPolyExpr(alg::expr&& u, alg::expr&& v, alg::expr& L);

/**
 * @brief Pseudo Divide two polynomial expressions and return the quotient,
 * both have to be transformed with the same set of symbols.
 *
 * @param u A polynomial expression
 * @param v A polynomial expression
 * @param L list of symbols in u and v
 * @return alg::expr quotient of the division of p1 by p2
 */
alg::expr pseudoQuoPolyExpr(alg::expr& u, alg::expr& v, alg::expr& L);

/**
 * @brief Pseudo Divide two polynomial expressions and return the remainder,
 * both have to be transformed with the same set of symbols.
 *
 * @param u A polynomial expression
 * @param v A polynomial expression
 * @param L list of symbols in u and v
 * @return alg::expr remainder of the division of p1 by p2
 */
alg::expr pseudoRemPolyExpr(alg::expr& u, alg::expr& v, alg::expr& L);

/**
 * @brief Get the leading coefficient of u
 *
 * @param u A polynomial expression
 * @return alg::expr leading coefficient of u
 */
alg::expr leadCoeffPolyExpr(alg::expr& u);

/**
 * @brief Get the highest degree of u
 *
 * @param u A polynomial expression
 * @return alg::expr the highest degree of u
 */
alg::expr degreePolyExpr(alg::expr& u);


/**
 * @brief Get the highest degree of u on x
 *
 * @param u A polynomial expression
 * @param x A symbol
 * @return alg::expr the highest degree of u
 */
alg::expr degreePolyExpr(alg::expr u, alg::expr x);

/**
 * @brief return a polynimal expression of x^e on the set of symbols L
 *
 * @param x A Symbol or FunctionCall
 * @param e A Integer
 * @param L a list of symbols
 * @return alg::expr a Polomial Expression
 */
alg::expr raiseToPolyExpr(alg::expr& x, Int e, alg::expr& L);
alg::expr raiseToPolyExpr(alg::expr&& x, Int e,alg::expr&& L);

/**
 * @brief return a polynimal expression of u*x^e
 *
 * @param u A Polynomial Expression
 * @param e A Integer
 * @param x a Symbol or FunctionCall
 * @return alg::expr a Polomial Expression
 */
alg::expr raisePolyExpr(alg::expr& u, Int exp, alg::expr& x);
alg::expr raisePolyExpr(alg::expr&& u, Int exp, alg::expr& x);

/**
 * @brief Get the gcd of two polynomial Expressions
 *
 * @param u A polynomial expression
 * @param v A polynomial expression
 * @param L A list of variables both in u and v
 * @param K A field, Z or Q are the only options
 * @return alg::expr the GCD of u and v
 */
alg::expr gcdPolyExpr(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

/**
 * @brief Return the normalized form of a polynomial expression
 *
 * @param u A polynomial Expression
 * @param L A list of variables
 * @param K A field, either Z or Q
 * @return alg::expr The normalized polynomial form of u
 */
alg::expr normalizePolyExpr(alg::expr& u, alg::expr& L, alg::expr& K);

/**
 * @brief Return the content of a polynomial Expression
 *
 * @param u A polynomial Expression
 * @param L The list of variables in u
 * @param K A field, either Z or Q
 * @return alg::expr the content of u
 */
alg::expr contPolyExpr(alg::expr& u, alg::expr& L, alg::expr& K);
alg::expr contPolyExpr(alg::expr&& u, alg::expr& L, alg::expr& K);

/**
 * @brief Return the primitive part of a polynomial Expression
 *
 * @param u A polynomial Expression
 * @param L The list of variables in u
 * @param K A field, either Z or Q
 * @return alg::expr the primitive part of u
 */
alg::expr ppPolyExpr(alg::expr& u, alg::expr& L, alg::expr& K);
alg::expr ppPolyExpr(alg::expr&& u, alg::expr& L, alg::expr& K);

/**
 * @brief Return the content and primitive part of a polynomial Expression
 *
 * @param u A polynomial Expression
 * @param L The list of variables in u
 * @param K A field, either Z or Q
 * @return alg::expr list with content and primitive part of u
 */
alg::expr contAndPpPolyExpr(alg::expr& u, alg::expr& L, alg::expr& K);
alg::expr contAndPpPolyExpr(alg::expr&& u, alg::expr& L, alg::expr& K);


/**
 * @brief Transform a poly expr into a normal expression
 *
 * @param u A polynomial expression
 * @return A expression that is the expansion of u
 */
alg::expr expandPolyExpr(alg::expr &u);
alg::expr expandPolyExpr(alg::expr &&u);


/**
 * @brief Differentiate a poly expr on the x variable
 *
 * @param u A polynomial expression
 * @param x the variable where u is to be differentiated
 * @return alg::expr the polynomial expression that is the derivative of u on x
 */
alg::expr diffPolyExpr(alg::expr& u, alg::expr& x);

alg::expr groundLeadCoeffPolyExpr(alg::expr u);

alg::expr groundContPolyExpr(alg::expr f);

alg::expr groundPPPolyExpr(alg::expr F);

alg::expr groundInvertPolyExpr(alg::expr p);

alg::expr evalPolyExpr(alg::expr u, alg::expr x, Int c);
alg::expr evalPolyExpr(alg::expr u, alg::expr x, alg::expr c);
alg::expr evalTailPolyExpr(alg::expr u, alg::expr L, alg::expr A, Int from = 0);

alg::expr insertSymbolPolyExpr(alg::expr u, alg::expr x, Int d, Int level = 0);
alg::expr insertSymbolsPolyExpr(alg::expr u, alg::expr L, Int d, Int level = 0);
alg::expr groundDivPolyExpr(alg::expr u, alg::expr v);
alg::expr groundMulPolyExpr(alg::expr u, Int v);

alg::expr raiseToExpression(alg::expr c, alg::expr u);


/**
 * @brief Return the list of variables in the expression 'a'.
 * @details The returned list is sorted by order of degree in 'a'.
 *
 * @param[in] a A expression.
 *
 * @return List of free variables of 'a'.
 */
alg::list getVariableListForPolyExpr(alg::expr a);

bool isPolynomial(alg::expr& a);


/**
 * @brief Compute the factorized form of P
 * @param[in] P A polynomial expression
 * @param[in] L A list of symbols of P
 * @param[in] K The field identifier, either Z or Q
 * @return The factors of P.
 */
alg::expr factorPolyExprAndExpand(alg::expr P, alg::expr L, alg::expr K);

/**
 * @brief Given two expressions, return two polynomial
 * expressions in the same set of variables.
 * @param[in] u A expression.
 * @param[in] v A expression.
 * @return list where the first element is the set of variables,
 * and the following two elements are the expressions u and v on
 * polynomial form;
 */
alg::expr normalizeToPolyExprs(alg::expr u, alg::expr v);

/**
 * @brief Compute the least commom multiple of two polynomials.
 * @param[in] u A polynomial expression.
 * @param[in] v A polynomial expression.
 * @param[in] L A list of symbols of P
 * @param[in] K The field identifier, either Z or Q
 * @return A polynomial expression equal to the lcm of u and v.
 */
alg::expr lcmPolyExpr(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

} // namespace polynomial

#endif
