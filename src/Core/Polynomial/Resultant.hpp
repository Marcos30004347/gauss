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
 * polynomials using the sub-resultant remainder sequence. This methods is
 * more efficient than the "multivariateResultant".
 * 
 * This method is currently not running for the polynomials:
 * 
 * u = x^3y^2 + gx^4y + 9x^5 + 4x^2y^2 + 24x^y + 36x^4 + 5xy^2 + 30yx^2
 * + 45x^3 + 2y^2 + 12yx + 18x^2
 * 
 * v = x^5y^2 8x^4y + 16x^3 + 12x^4y^2 + 96x^3y + 192x^2 + 45x^3y^2
 * + 360yx^2 + 720x + 50x^2y^2 + 400yx + 800
 * 
 * in Z[x,y]
 * 
 * The resultant between those two should be 0(zero)
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
