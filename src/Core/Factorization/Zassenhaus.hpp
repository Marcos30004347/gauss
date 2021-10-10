#ifndef FACTORIZATION_ZASSENHAUS_H
#define FACTORIZATION_ZASSENHAUS_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief Computes the irreductible factors f[i] in Z[x] of f
 * 
 * @param f A polynomial in Z[x]
 * @param x The symbol x
 * @return A list with all the factors of f in Z[x]
 */
ast::AST* zassenhaus(ast::AST* f, ast::AST* x);

/**
 * @brief Given a square-free polynomial a(x) in Zp[x], 
 * 				computes its partial distinct degree factorization
 * 
 * @param a A square free polynomial in Zp[x]
 * @param x The symbol x
 * @param q A prime integer p
 * @return A list of tuples, the first element of the tuple is a factor, and the second is its degree
 */
ast::AST* cantorZassenhausDDF(ast::AST* a, ast::AST* x, long p);

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
ast::AST* cantorZassenhausEDF(ast::AST* a, ast::AST* x, long n, long p);

}

#endif
