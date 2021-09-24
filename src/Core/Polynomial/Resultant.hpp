#ifndef MATH_ALGEBRA_HENSEL_H
#define MATH_ALGEBRA_HENSEL_H

#include "Polynomial.hpp"
#include "Core/Algebra/List.hpp"

#include <vector>
#include <utility> 

namespace polynomial {

/**
 * @brief Euclidean algorithm that obtains the resultant in F[x].
 * 
 * @param u A non-zero polynomial u(x) in F[x] where all field in F are 
 * obtained with automatic simplification.
 * @param v A non-zero polynomial v(x) in F[x] where all field in F are 
 * obtained with automatic simplification.
 * @param x The main variable of u(x) and v(x)
 * @return ast::AST* 
 */
ast::AST* univariateResultant(ast::AST* u, ast::AST* v, ast::AST* x);

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
 * @return ast::AST* 
 */
ast::AST* multivariateResultant(ast::AST* u, ast::AST* v, ast::AST* L, ast::AST* K);

/**
 * @brief Recursive procedure that obtains the resultant for multivariate 
 * polynomials using the sub-resultant remainder sequence.
 * 
 * @param u A non-zero multivariate polynomial in symbols in L
 * with coefficients in Z or Q.
 * @param v A non-zero multivariate polynomial in symbols in L
 * with coefficients in Z or Q.
 * @param L A list of symbols with main variable being the first element of the list.
 * @param K The symbol Z or Q.
 * @return ast::AST* 
 */
ast::AST* srPolynomialResultant(ast::AST*	u, ast::AST* v, ast::AST* L, ast::AST* K);

}

#endif
