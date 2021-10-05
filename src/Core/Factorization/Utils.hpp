#ifndef FACTORIZATION_UTILS_H
#define FACTORIZATION_UTILS_H

#include "Core/Algebra/Algebra.hpp"

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
long landauMignotteBound(ast::AST* u, ast::AST* x);


ast::AST* repeatSquaring(ast::AST* a, ast::AST* x, long n, long m);

}

#endif
