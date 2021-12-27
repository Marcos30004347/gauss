#ifndef FACTORIZATION_ZASSENHAUS_H
#define FACTORIZATION_ZASSENHAUS_H

#include "Core/AST/AST.hpp"
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
	ast::Expr zassenhaus(ast::Expr f, ast::Expr x, ast::Expr K);

/**
 * @brief Given a square-free polynomial a(x) in Zp[x],
 * 				computes its partial distinct degree factorization
 *
 * @param a A square free polynomial in Zp[x]
 * @param x The symbol x
 * @param q A prime integer p
 * @return A list of tuples, the first element of the tuple is a factor, and the second is its degree
 */
ast::Expr cantorZassenhausDDF(ast::Expr a, ast::Expr x, Int p);

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
ast::Expr cantorZassenhausEDF(ast::Expr a, ast::Expr x, Int n, Int p);

/**
 * @brief Given an univariate polynomial in Zp[x], Compute its unique factorization in Zp[x]
 *
 * @param u A polynomial in Zp[x]
 * @param x The symbol x
 * @param p A prime integer
 * @return A product expression with all the irreductible factors of u(x) in Zp[x] as operands
 */
ast::Expr cantorZassenhaus(ast::Expr u, ast::Expr x, Int p);

}

#endif
