#ifndef MATH_ALGEBRA_HENSEL_H
#define MATH_ALGEBRA_HENSEL_H

#include "Polynomial.hpp"

#include <vector>
#include <utility> 

namespace polynomial {

/**
 * @brief 
 * 
 * @param ax - A primitive Polynomial a(x) defined in Z[x]
 * @param p - A prime integer p which does not divide leadCoeff(a(x))
 * @param ux - A polynomial u(x) relative prime with w(x) such that a(x) = u(x)*w(x) (mod p)
 * @param wx - A polynomial w(x) relative prime with u(x) such that a(x) = u(x)*w(x) (mod p)
 * @param B - A integer which bounds the magnitudes of all integers coefficients appearing in a(x)
 * and in any of its possible factors with degrees not exceeding max{deg(u(x)), deg(w(x))}, this
 * value limits the ammount of iterations taken by the algorithm.
 * @param gamma - Optionally, an integer gamma that exists in Z which is known to be a multiple of
 * leadCoeff(u(x)), where u(x) is one of the factors of a(x) in Z[x] to be computed
 * @return ast::AST* 
 */
ast::AST* univariateHensel(ast::AST* ax, ast::AST* x, ast::AST* p, ast::AST* ux, ast::AST* wx, ast::AST* B, ast::AST* gamma);

ast::AST* leadCoeffReplace(ast::AST* ux, ast::AST* x, ast::AST* c);

}

#endif
