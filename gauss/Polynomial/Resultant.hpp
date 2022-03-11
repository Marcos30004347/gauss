#ifndef MATH_ALGEBRA_RESULTANT_H
#define MATH_ALGEBRA_RESULTANT_H

#include "gauss/Algebra/Expression.hpp"
#include "Polynomial.hpp"

#include <vector>
#include <utility>

namespace poly {

/**
 * @brief Euclidean algorithm that obtains the resultant in F[x].
 *
 * @param u A non-zero polynomial u(x) in F[x] where all field in F are
 * obtained with automatic simplification.
 * @param v A non-zero polynomial v(x) in F[x] where all field in F are
 * obtained with automatic simplification.
 * @param x The main variable of u(x) and v(x)
 * @return alg::expr
 */
// alg::expr univariateResultant(alg::expr u, alg::expr v, alg::expr x);

/**
 * @brief Euclidean algorithm that obtains the resultant for multivariate
 * polynomials.
 *
 * @param u A non-zero multivariate polynomial in symbols in L
 * with coefficients in Z or Q.
 * @param v A non-zero multivariate polynomial in symbols in L
 * with coefficients in Z or Q.
 * @param L A list of symbols with main variable being the first element of the list.
 * @param K The symbol Z or Q.
 * @return alg::expr
 */
// alg::expr multivariateResultant(alg::expr u, alg::expr v, alg::expr L, alg::expr K);

/**
 * @brief Recursive procedure that obtains the resultant for multivariate
 * polynomials using the sub-resultant remainder sequence. This methods is
 * more efficient than the "multivariateResultant".
 *
 * @param u A non-zero multivariate polynomial in symbols in L
 * with coefficients in Z or Q.
 * @param v A non-zero multivariate polynomial in symbols in L
 * with coefficients in Z or Q.
 * @param L A list of symbols with main variable being the first element of the list.
 * @param K The symbol Z or Q.
 * @return alg::expr
 */
// alg::expr polynomialResultant(alg::expr	u, alg::expr v, alg::expr L, alg::expr K);

/**
 * @brief Computes the polynomial remainder sequence between
 * polynomials F1(L...) and F2(L...) and return the gcd and resultant
 * between the two.
 *
 * @param F1 A nonzero polynomial defined in K[L...].
 * @param F2 A nonzero polynomial defined in K[L...].
 * @param L The list of symbols that define the variables of F1 and F2,
 * the remainder sequence will use the first symbol in L as main variable.
 * @param K The field of the polynomials, should be either Z or Q.
 * @return A list [gcd, res]
 */
// alg::expr polyRemSeq(alg::expr F1, alg::expr F2, alg::expr L, alg::expr K);

alg::expr remSeqPolyExpr(alg::expr F1, alg::expr F2, alg::expr L, alg::expr K);
alg::expr resultantPolyExpr(alg::expr u, alg::expr v, alg::expr L, alg::expr K);


}

#endif
