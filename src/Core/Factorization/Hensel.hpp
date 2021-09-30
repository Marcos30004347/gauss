#ifndef FACTORIZATION_HENSEL_H
#define FACTORIZATION_HENSEL_H

#include "Core/Algebra/List.hpp"
#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief Lift the factors u(x) and w(x) in Zp[x] to Z[x] using Hensel Lift process.
 * 
 * @param ax		A primitive Polynomial a(x) defined in Z[x].
 * 
 * @param p			A prime integer p which does not divide leadCoeff(a(x)).
 * 
 * @param ux		A polynomial u(x) relative prime with w(x) such that 
 * 							a(x) = u(x)*w(x) (mod p).
 * 
 * @param wx		A polynomial w(x) relative prime with u(x) such that 
 * 							a(x) = u(x)*w(x) (mod p).
 * 
 * @param B			A integer which bounds the magnitudes of all integers coefficients 
 * 							appearing in a(x) and in any of its possible factors with degrees
 * 							not exceeding max{deg(u(x)), deg(w(x))}, this value limits the ammount 
 * 							of iterations taken by the algorithm.
 * 
 * @param gamma	Optionally, an integer gamma that exists in Z which is known to be 
 * 							a multiple of leadCoeff(u(x)), where u(x) is one of the factors of 
 * 							a(x) in Z[x] to be computed.
 * 
 * @return 			The list with the two lifted factors.
 */
ast::AST* univariateHensel(ast::AST* ax, ast::AST* x, ast::AST* p, ast::AST* ux, ast::AST* wx, ast::AST* B, ast::AST* gamma);

}

#endif
