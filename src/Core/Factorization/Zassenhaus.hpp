#ifndef FACTORIZATION_ZASSENHAUS_H
#define FACTORIZATION_ZASSENHAUS_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief Computes the irreductible factors f[i] in Z[L[0]] of f
 * 
 * @param f A polynomial in Z[L[0]]
 * @param L A list that should have just one element,
 * 					that is the only variable of f
 * @param K The field, should be Z
 * @return A list with all the factors of f in Z[L[0]]
 */
ast::AST* zassenhaus(ast::AST* f, ast::AST* L, ast::AST* K);

ast::AST* distinctDegreeFactorization(ast::AST* v, ast::AST* L, ast::AST* K, ast::AST* q);
ast::AST* equalDegreeFactorization(ast::AST* a, ast::AST* x, long n, long p);

}

#endif
