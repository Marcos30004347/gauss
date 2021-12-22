#ifndef FACTORIZATION_SQUARE_FREE_H
#define FACTORIZATION_SQUARE_FREE_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief Given a primitive polynomial a(x) in R[x],
 * 				computes the square-free factorization of a(x).
 *
 * @param ax 	A plynomial.
 * @param x 	A symbol.
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 		that f[i] is square free.
 */
ast::Expr squareFreeFactorization(ast::Expr ax, ast::Expr x);

/**
 * @brief Given a primitive polynomial expression a(x) in Z[x],
 * 				computes the square-free factorization of a(x).
 *
 * @param ax A plynomial expression.
 * @param L The list of symbols of ax, this list can have size at most 1.
 * @param K The field of az, Only Z is allowed
 * @return The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 that f[i] is square free.
 */
ast::Expr squareFreeFactorizationPolyExpr(ast::Expr ax, ast::Expr L, ast::Expr K);

/**
 * @brief Given a primitive polynomial a(x) in Z[x],
 * 				computes the square-free factorization of a(x).
 *
 * @param ax 	A plynomial in Z[x].
 * @param x 	A symbol.
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 		that f[i] is square free.
 */
ast::Expr squareFreeFactorization2(ast::Expr ax, ast::Expr x);

/**
 * @brief Given a primitive polynomial expression a(x) in Z[x],
 * 				computes the square-free factorization of a(x).
 *
 * @param ax A plynomial in Z[x].
 * @param L The list of symbols of ax, this list can have size at most 1.
 * @param Z The field of a(x), only Z is allowed
 * @return The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 that f[i] is square free.
 */
ast::Expr squareFreeFactorizationPolyExpr2(ast::Expr ax, ast::Expr L, ast::Expr Z);

/**
 * @brief Given a primitive polynomial a(x) in GF(q)[x],
 * 				with GF(q) a Galois field of order q = p^m,
 * 				computes the square-free factorization of a(x)
 *
 * @param ax 	A plynomial in Zq[x]
 * @param x 	A symbol
 * @param q 	a prime integer
 * @param sym true if the result should be in symmetric representation over Zq[x]
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 		that f[i] is square free
 */
ast::Expr squareFreeFactorizationFiniteField(ast::Expr ax, ast::Expr x, Int q, bool sym = true);

/**
 * @brief Given a primitive polynomial expression a(x) in Zq[x],
 * computes its square free factorization.
 *
 * @param ax 	A plynomial expression in Zq[x]
 * @param L 	The list of symbols in a(x), this list can have at most 1
element
 * @param q 	a prime integer
 * @param K the field of a(x), only Z is allowed
 * @param sym true if the result should be in symmetric representation over Zq[x]
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 		that f[i] is square free
 */
ast::Expr squareFreeFactorizationFiniteFieldPolyExpr(ast::Expr ax, ast::Expr L, ast::Expr K, Int q, bool sym = true);


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
 * @brief Computes the square free part of a polynomial expression in K[L...].
 *
 * @param f A multivariate polynomial expression in K[L...]
 * @param L The list of variable symbols of f
 * @param K Either Z or Q
 * @return The square free part of f
 */
ast::Expr squareFreePartPolyExpr(ast::Expr f, ast::Expr L, ast::Expr K);

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
