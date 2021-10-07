#ifndef FACTORIZATION_HENSEL_H
#define FACTORIZATION_HENSEL_H

#include "Core/Algebra/List.hpp"
#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief Lift the factors u(x) and w(x) in Zp[x] to Z[x] using Hensel Lift process.
 * 
 * @param ax		A primitive Polynomial a(x) defined in Z[x].
 * @param p			A prime integer p which does not divide leadCoeff(a(x)).
 * @param ux		A polynomial u(x) relative prime with w(x) such that 
 * 							a(x) = u(x)*w(x) (mod p).
 * @param wx		A polynomial w(x) relative prime with u(x) such that 
 * 							a(x) = u(x)*w(x) (mod p).
 * @param B			A integer which bounds the magnitudes of all integers coefficients 
 * 							appearing in a(x) and in any of its possible factors with degrees
 * 							not exceeding max{deg(u(x)), deg(w(x))}, this value limits the ammount 
 * 							of iterations taken by the algorithm.
 * @param gamma	Optionally, an integer gamma that exists in Z which is known to be 
 * 							a multiple of leadCoeff(u(x)), where u(x) is one of the factors of 
 * 							a(x) in Z[x] to be computed.
 * @return 			The list with the two lifted factors.
 */
ast::AST* univariateHensel(ast::AST* ax, ast::AST* x, ast::AST* p, ast::AST* ux, ast::AST* wx, ast::AST* B, ast::AST* gamma);

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
 * @return G, H, S, T such that
 * 				 f = G*H mod (m * m) and S*G + T*H = 1 mod (m * m) 
 */
ast::AST* henselSep(ast::AST* f, ast::AST* g, ast::AST* h, ast::AST* s, ast::AST* t, ast::AST* x, long m);


/**
 * @brief Compute monic polynomials F[0],...,F[r - 1] in K[x] 
 * 				with f = lc(f) * F[0] * ... * F[r-1] mod (p^l) and
 * 				F[i] = f[i] mod p for all 0 <= i < r
 * 				
 * @param f A polynomial in K[x]
 * @param H A list of factors to be lifted
 * @param L The list of variables, this list should have just the symbol x on it
 * @param K The field symbol, either Z or Q
 * @param p A integer
 * @param l A integer
 * @return The list of lifted factors
 */
ast::AST* multifactorHenselLifting(ast::AST* f, ast::AST* H, ast::AST* L, ast::AST* K, long p, long l);

}

#endif
