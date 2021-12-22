#ifndef FACTORIZATION_SQUARE_FREE_H
#define FACTORIZATION_SQUARE_FREE_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief Given a primitive polynomial a(x) in R[x],
 * 				computes the square-free factorization of a(x).
 *
 * @param ax 	A plynomial in GF(q)[x].
 * @param x 	A symbol.
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 		that f[i] is square free.
 */
ast::Expr squareFreeFactorization(ast::Expr ax, ast::Expr x);

/**
 * @brief Given a primitive polynomial a(x) in R[x],
 * 				computes the square-free factorization of a(x).
 *
 * @param ax 	A plynomial in GF(q)[x].
 * @param x 	A symbol.
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 		that f[i] is square free.
 */
ast::Expr squareFreeFactorization2(ast::Expr ax, ast::Expr x);


/**
 * @brief Given a primitive polynomial a(x) in GF(q)[x],
 * 				with GF(q) a Galois field of order q = p^m,
 * 				computes the square-free factorization of a(x)
 *
 * @param ax 	A plynomial in Zq[x]
 * @param x 	A symbol
 * @param q 	a power p^m
 * @param sym true if the result should be in symmetric representation over Zq[x]
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 		that f[i] is square free
 */
ast::Expr squareFreeFactorizationFiniteField(ast::Expr ax, ast::Expr x, Int q, bool sym = true);

/**
 * @brief Computes the square free part of a polynomial in K[L...].
 *
 * @param f A multivariate polynomial in K[L...]
 * @param L The list of variable symbols of f
 * @param K Either Z or Q
 * @return The square free part of f
 */
ast::Expr squareFreePart(ast::Expr f, ast::Expr L, ast::Expr K);

/**
 * @brief Computes if f is square free in Zp[x]
 *
 * @param f A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A integer
 * @param sym true if f is in symmetric form over Zp[x], false otherwise.
 * @return true if f is quare free, false otherwise
 */
bool isSquareFreeInZp(ast::Expr f, ast::Expr x, long p, bool sym = true);

/**
 * @brief Computes if f is square free in K[x]
 *
 * @param f A polynomial in K[x]
 * @param x The symbol x
 * @return true if f is quare free, false otherwise
 */
bool isSquareFree(ast::Expr f, ast::Expr x, ast::Expr K);

}

#endif
