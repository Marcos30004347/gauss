#ifndef FACTORIZATION_UTILS_H
#define FACTORIZATION_UTILS_H

#include "Core/Algebra/Algebra.hpp"
#include <limits>

namespace factorization {

/**
 * @brief Computes the Landau Mignotte Bound for u(x) in Z[x].
 * 				such that if g(x) in Z[x] is a factor of u(x), then
 * 				every coefficieng of g(x), lets say g[i], obey
 * 				|g[i]| <= |'Landau-Mignotte-Bound'|.
 *
 * @param u A polynomial u(x) in Z[x].
 *
 * @param x A symbol.
 *
 * @return The Landau-Mignotte Bound B, such that if
 * 				 g(x) is a factor if u(x), every coeff
 * 				 g[i] of g(x) is smaller in magnitude or
 * 				 equal than |B|.
 */
Int landauMignotteBound(ast::Expr u, ast::Expr x);


/**
 * @brief Computes the infinite norm of the polynomial.
 * 				It is equal to the bigger absolute coefficient
 * 				of u in K[L...]
 *
 * @param u A polynomial in K[L...]
 * @param L The list of symbols
 * @param K The field, either K or Q
 * @param i The index of the first variable, defaults to zero
 * @return The magnitude of the largest coefficient in u
 */
Int norm(ast::Expr u, ast::Expr L, ast::Expr K, long long i = 0);

/**
 * @brief Computes the infinite norm of the polynomial.
 * 				It is equal to the bigger absolute coefficient
 * 				of u in Z[x]
 *
 * @param u A polynomial in Z[x]
 * @param x The symbol x
 * @return The magnitude of the largest coefficient in u
 */
Int norm(ast::Expr u, ast::Expr x);

/**
 * @brief Computes the infinite norm of a polynomial expression.
 * 				It is equal to the bigger absolute coefficient
 * 				of u in Z[x]
 *
 * @param u A univariate polynomial expression
 * @return The magnitude of the largest coefficient in u
 */
Int normPolyExpr(ast::Expr u);


/**
 * @brief Computes the L1 norm of the polynomial.
 * 				It is equal to sum of the absolute value
 * 				of the coefficients of u in K[L...]
 *
 * @param u A polynomial in K[L...]
 * @param L The list of symbols
 * @param K The field, either K or Q
 * @param i The index of the first variable, defaults to zero
 * @return The magnitude of the largest coefficient in u
 */
Int l1norm(ast::Expr u, ast::Expr L, ast::Expr K, long long i = 0);

/**
 * @brief Computes the L1 norm of the polynomial.
 * 				It is equal to sum of the absolute value
 * 				of the coefficients of u in Z[x]
 *
 * @param u A polynomial in Z[x]
 * @param x The symbol x
 * @return The magnitude of the largest coefficient in u
 */
Int l1norm(ast::Expr u, ast::Expr x);

/**
 * @brief Computes the L1 norm of a polynomial expression.
 * 				It is equal to sum of the absolute value
 * 				of the coefficients of u in Z[x]
 *
 * @param u A univariate polynomial expression
 * @return The magnitude of the largest coefficient in u
 */
Int l1normPolyExpr(ast::Expr u);


/**
 * @brief Pick a random number between min and max
 *
 * @param min A integer
 * @param max A integer
 * @return A random number between [min, max]
 */
Int random(long long min = std::numeric_limits<long long>::min(), long long max = std::numeric_limits<long long>::max());



/**
 * @brief Sort Operands on F
 *
 * @param F a expression
 * @return A expression of the same type of F with the operands sorted
 */
ast::Expr sortTerms(ast::Expr& F);

}


#endif
