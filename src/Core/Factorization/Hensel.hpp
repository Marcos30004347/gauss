#ifndef FACTORIZATION_HENSEL_H
#define FACTORIZATION_HENSEL_H

#include "Core/Algebra/List.hpp"
#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

// /**
//  * @brief Lift the factors u(x) and w(x) in Zp[x] to Z[x] using Hensel Lift process.
//  *
//  * @param ax		A primitive Polynomial a(x) defined in Z[x].
//  * @param p			A prime integer p which does not divide leadCoeff(a(x)).
//  * @param ux		A polynomial u(x) relative prime with w(x) such that
//  * 							a(x) = u(x)*w(x) (mod p).
//  * @param wx		A polynomial w(x) relative prime with u(x) such that
//  * 							a(x) = u(x)*w(x) (mod p).
//  * @param B			A integer which bounds the magnitudes of all integers coefficients
//  * 							appearing in a(x) and in any of its possible factors with degrees
//  * 							not exceeding max{deg(u(x)), deg(w(x))}, this value limits the ammount
//  * 							of iterations taken by the algorithm.
//  * @param gamma	Optionally, an integer gamma that exists in Z which is known to be
//  * 							a multiple of leadCoeff(u(x)), where u(x) is one of the factors of
//  * 							a(x) in Z[x] to be computed.
//  * @param sym True if Zp[x] is in symmetric representation, false otherwise
//  * @return 			The list with the two lifted factors.
//  */
// ast::Expr univariateHensel(ast::Expr ax, ast::Expr x, ast::Expr p, ast::Expr ux, ast::Expr wx, ast::Expr B, ast::Expr gamma, bool sym = true);

/**
 * @brief Computes polynomials G, H, S, T in K[x] such that
 * 				f = G*H mod (m * m) and S*G + T*H = 1 mod (m * m)
 *
 * @param f A polynomial in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param g A polynomial in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param h A polynomial in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param s A polynomial in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param t A polynomial in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param x The symbol x
 * @param m A integer such that f = g*h mod m and s*g t*h = 1 mod m
 * @param sym True if Zp[x] is in symmetric representation, false otherwise
 * @return G, H, S, T such that
 * 				 f = G*H mod (m * m) and S*G + T*H = 1 mod (m * m)
 */
ast::Expr henselSep(ast::Expr f, ast::Expr g, ast::Expr h, ast::Expr s, ast::Expr t, ast::Expr x, Int m, bool sym = true);


/**
 * @brief Compute monic polynomials F[0],...,F[r - 1] in Zp[x]
 * 				with f = lc(f) * F[0] * ... * F[r-1] mod (p^l) and
 * 				F[i] = f[i] mod p for all 0 <= i < r
 *
 * @param f A polynomial in Zp[x]
 * @param H A list of factors to be lifted
 * @param x The symbol x
 * @param p A integer
 * @param l A integer
 * @param sym True if Zp[x] is in symmetric representation, false otherwise
 * @return The list of lifted factors
 */
ast::Expr multifactorHenselLifting(ast::Expr f, ast::Expr H, ast::Expr x, Int p, Int l, bool sym = true);

}

#endif
