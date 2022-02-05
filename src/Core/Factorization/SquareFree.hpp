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
alg::expr squareFreeFactorization(alg::expr ax, alg::expr x);

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
alg::expr squareFreeFactorizationPolyExpr(alg::expr ax, alg::expr L, alg::expr K);

/**
 * @brief Given a primitive polynomial a(x) in Z[x],
 * 				computes the square-free factorization of a(x).
 *
 * @param ax 	A plynomial in Z[x].
 * @param x 	A symbol.
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such
 * 				 		that f[i] is square free.
 */
alg::expr squareFreeFactorization2(alg::expr ax, alg::expr x);

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
alg::expr squareFreeFactorizationPolyExpr2(alg::expr ax, alg::expr L, alg::expr Z);

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
alg::expr squareFreeFactorizationFiniteField(alg::expr ax, alg::expr x, Int q, bool sym = true);

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
alg::expr squareFreeFactorizationFiniteFieldPolyExpr(alg::expr ax, alg::expr L, alg::expr K, Int q, bool sym = true);


/**
 * @brief Computes the square free part of a polynomial in K[L...].
 *
 * @param f A multivariate polynomial in K[L...]
 * @param L The list of variable symbols of f
 * @param K Either Z or Q
 * @return The square free part of f
 */
alg::expr squareFreePart(alg::expr f, alg::expr L, alg::expr K);


/**
 * @brief Computes the square free part of a polynomial expression in K[L...].
 *
 * @param f A multivariate polynomial expression in K[L...]
 * @param L The list of variable symbols of f
 * @param K Either Z or Q
 * @return The square free part of f
 */
alg::expr squareFreePartPolyExpr(alg::expr f, alg::expr L, alg::expr K);

/**
 * @brief Computes if f is square free in Zp[x]
 *
 * @param f A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A integer
 * @param sym true if f is in symmetric form over Zp[x], false otherwise.
 * @return true if f is quare free, false otherwise
 */
bool isSquareFreeInZp(alg::expr f, alg::expr x, long p, bool sym = true);

/**
 * @brief Computes if f is square free in K[x]
 *
 * @param f A polynomial in K[x]
 * @param x The symbol x
 * @return true if f is quare free, false otherwise
 */
bool isSquareFree(alg::expr f, alg::expr x, alg::expr K);


/**
 * @brief Computes if f is square free in K[x]
 *
 * @param f A polynomial expression in K[x]
 * @param x The symbol x
 * @return true if f is quare free, false otherwise
 */
bool isSquareFreePolyExpr(alg::expr f, alg::expr x, alg::expr K);

}

#endif
