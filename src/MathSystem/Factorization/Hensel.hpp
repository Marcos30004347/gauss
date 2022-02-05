#ifndef FACTORIZATION_HENSEL_H
#define FACTORIZATION_HENSEL_H

#include "MathSystem/Polynomial/Polynomial.hpp"

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
// alg::expr univariateHensel(alg::expr ax, alg::expr x, alg::expr p, alg::expr ux, alg::expr wx, alg::expr B, alg::expr gamma, bool sym = true);

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
alg::expr henselSep(alg::expr f, alg::expr g, alg::expr h, alg::expr s, alg::expr t, alg::expr x, Int m, bool sym = true);


/**
 * @brief Computes polynomial expressions G, H, S, T in K[x] such that
 * 				f = G*H mod (m * m) and S*G + T*H = 1 mod (m * m)
 *
 * @param f A polynomial expressions in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param g A polynomial expressions in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param h A polynomial expressions in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param s A polynomial expressions in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param t A polynomial expressions in K[x] such that f = g*h mod m and s*g t*h = 1 mod m
 * @param L The list symbols L, that list can have at most one element
 * @param m A integer such that f = g*h mod m and s*g t*h = 1 mod m
 * @param sym True if Zp[x] is in symmetric representation, false otherwise
 * @return G, H, S, T such that
 * 				 f = G*H mod (m * m) and S*G + T*H = 1 mod (m * m)
 */
alg::expr henselSepPolyExpr(alg::expr f, alg::expr g, alg::expr h, alg::expr s, alg::expr t, alg::expr L, Int m, bool sym = true);

/**
 * @brief Compute monic polynomials F[0],...,F[r - 1] in Zp[x]
 * 				with f = lc(f) * F[0] * ... * F[r-1] mod (p^l)
 * and F[i] = f[i] mod p for all 0 <= i < r
 *
 * @param f A polynomial in Zp[x]
 * @param H A list of factors to be lifted
 * @param x The symbol x
 * @param p A integer
 * @param l A integer
 * @param sym True if Zp[x] is in symmetric representation, false otherwise
 * @return The list of lifted factors
 */
alg::expr multifactorHenselLifting(alg::expr f, alg::expr H, alg::expr x, Int p,
                                   Int l, bool sym = true);

/**
 * @brief Compute monic polynomials F[0],...,F[r - 1] in Zp[x]
 * 				with f = lc(f) * F[0] * ... * F[r-1] mod (p^l)
 * and F[i] = f[i] mod p for all 0 <= i < r
 *
 * @param f A polynomial in Zp[x]
 * @param H A list of factors to be lifted
 * @param L The list of symbols of f, this list should have size one
 * @param p A integer
 * @param l A integer
 * @param sym True if Zp[x] is in symmetric representation, false otherwise
 * @return The list of lifted factors
 */
alg::expr multifactorHenselLiftingPolyExpr(alg::expr f, alg::expr H,
                                           alg::expr L, Int p, Int l,
                                           bool sym = true);

} // namespace factorization

#endif
