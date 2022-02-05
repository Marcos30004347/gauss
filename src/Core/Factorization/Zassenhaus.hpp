#ifndef FACTORIZATION_ZASSENHAUS_H
#define FACTORIZATION_ZASSENHAUS_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief Computes the irreductible factors f[i] in Z[x] of f
 *
 * @param f A polynomial in Z[x]
 * @param x The symbol x
 * @param The field of f(x), only Z is accepted
 * @return A list with all the factors of f in Z[x]
 */
alg::expr zassenhaus(alg::expr f, alg::expr x, alg::expr K);

/**
 * @brief Computes the irreductible factors f[i] in Z[x] of f
 *
 * @param f A polynomial expression
 * @param The list of symbols of f, this list can have at most one element
 * @param The field of f(x), only Z is accepted
 * @return A list with all the factors of f in Z[x]
 */
alg::expr zassenhausPolyExpr(alg::expr f, alg::expr L, alg::expr K);


/**
 * @brief Given a square-free polynomial a(x) in Zp[x],
 * 				computes its partial distinct degree factorization
 *
 * @param a A square free polynomial in Zp[x]
 * @param x The symbol x
 * @param q A prime integer p
 * @return A list of tuples, the first element of the tuple is a factor, and the second is its degree
 */
alg::expr cantorZassenhausDDF(alg::expr a, alg::expr x, Int p);


/**
 * @brief Given a square-free polynomial expression a(x) in Zp[x],
 * 				computes its partial distinct degree factorization
 *
 * @param a A square free polynomial in Zp[x]
 * @param L The list nof symbols in a, this list should have at most 1 element
 * @param q A prime integer p
 * @return A list of tuples, the first element of the tuple is a factor, and the second is its degree
 */
alg::expr cantorZassenhausDDFPolyExpr(alg::expr a, alg::expr L, Int p);

/**
 * @brief Given a square-free polynomial a(x) and a integer n in Zp[x],
 * 				computes the equal degree factorization of a(x)^n
 *
 * @param a A square free polynomial in Zp[x]
 * @param x The symbol x
 * @param n An integer
 * @param p A prime integer
 * @return A list of equal degree factors
 */
alg::expr cantorZassenhausEDF(alg::expr a, alg::expr x, Int n, Int p);

/**
 * @brief Given a square-free polynomial a(x) and a integer n in Zp[x],
 * 				computes the equal degree factorization of a(x)^n
 *
 * @param a A square free polynomial in Zp[x]
 * @param L The list of symbols of a, this list can have at most 1 symbol
 * @param n An integer
 * @param p A prime integer
 * @return A list of equal degree factors
 */
alg::expr cantorZassenhausEDFPolyExpr(alg::expr a, alg::expr L, Int n, Int p);



/**
 * @brief Given an univariate polynomial in Zp[x], Compute its unique factorization in Zp[x]
 *
 * @param u A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @return A product expression with all the irreductible factors of u(x) in Zp[x] as operands
 */
alg::expr cantorZassenhaus(alg::expr u, alg::expr x, Int p);

/**
 * @brief Given an univariate polynomial expression in Zp[x], Compute its unique factorization in Zp[x]
 *
 * @param u A polynomial in Zp[x]
 * @param L The list of symbols in u(x), this list can have at most one element
 * @param p A prime integer
 * @return A product expression with all the irreductible factors of u(x) in Zp[x] as operands
 */
alg::expr cantorZassenhausPolyExpr(alg::expr u, alg::expr L, Int p);



}

#endif
