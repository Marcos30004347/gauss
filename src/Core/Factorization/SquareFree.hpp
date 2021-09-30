#ifndef FACTORIZATION_SQUARE_FREE_H
#define FACTORIZATION_SQUARE_FREE_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * Algorithm 8.1 from "Algorithms for computer algebra"
 * 
 * Given a primitive polynomial a(x) in R[x], with
 * characteristic zero, computes the square-free
 * factorization of a(x).
 */
ast::AST* squareFreeFactorization(ast::AST* ax, ast::AST* x);

/**
 * Algorithm 8.2 from "Algorithms for computer algebra"
 * 
 * Given a primitive polynomial a(x) in R[x], with
 * characteristic zero, computes the square-free
 * factorization of a(x).
 */
ast::AST* squareFreeFactorization2(ast::AST* ax, ast::AST* x);

/**
 * Algorithm 8.3 from "Algorithms for computer algebra"
 * 
 * Given a primitive polynomial a(x) in GF(q)[x], with
 * GF(q) a Galois field of order q = p^m, 
 * computes the square-free factorization of a(x).
 */
ast::AST* squareFreeFactorizationFiniteField(ast::AST* ax, ast::AST* x, ast::AST* q);

}

#endif
