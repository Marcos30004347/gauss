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
ast::AST* squareFreeFactorization(ast::AST* ax, ast::AST* x);

/**
 * @brief Given a primitive polynomial a(x) in R[x], 
 * 				computes the square-free factorization of a(x).
 * 
 * @param ax 	A plynomial in GF(q)[x].
 * @param x 	A symbol.
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such 
 * 				 		that f[i] is square free.
 */
ast::AST* squareFreeFactorization2(ast::AST* ax, ast::AST* x);


/**
 * @brief Given a primitive polynomial a(x) in GF(q)[x], 
 * 				with GF(q) a Galois field of order q = p^m, 
 * 				computes the square-free factorization of a(x).
 * 
 * @param ax 	A plynomial in GF(q)[x].
 * @param x 	A symbol.
 * @param q 	a power p^m.
 * @return 		The factorization of a(x) = f[1]*f[2]*...*f[n] such 
 * 				 		that f[i] is square free.
 */
ast::AST* squareFreeFactorizationFiniteField(ast::AST* ax, ast::AST* x, ast::AST* q);

}

#endif
