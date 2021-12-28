#ifndef FACTORIZATION_WANG_H
#define FACTORIZATION_WANG_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

/**
 * @brief This is a helper function that is used to test
 * 				evaluation point suitability.
 *
 * @param G An integer that multiplies F[0]*...*F[i]
 *
 * @param F The list of factors F[i] that will be tested.
 *
 * @param c The content of U[0](x)
 *
 * @param L A list containing the variablex x[i]
 *
 * @param K The field K[L...] of F[i](L...)
 *
 *  @return The set of integers d[i] defined as above,
 * 					this array must be pre-allocated before
 * 					the call.
 *
 */
ast::Expr nondivisors(Int G, ast::Expr F, Int c, ast::Expr L, ast::Expr K);

ast::Expr factors(ast::Expr f, ast::Expr L, ast::Expr K);

ast::Expr getEvaluationPoints(ast::Expr f, ast::Expr G, ast::Expr F, ast::Expr L, ast::Expr K, Int p, ast::Expr S);

ast::Expr factorsWang(ast::Expr f, ast::Expr L, ast::Expr K);
ast::Expr trialDivision(ast::Expr f, ast::Expr F, ast::Expr L, ast::Expr K);

ast::Expr groundLeadCoeff(ast::Expr f, ast::Expr L);
ast::Expr groundCont(ast::Expr f, ast::Expr L, ast::Expr K);
ast::Expr groundPP(ast::Expr f, ast::Expr L, ast::Expr K);
ast::Expr groundPP(ast::Expr f, ast::Expr c, ast::Expr L, ast::Expr K);

/**
 * @brief Computes s, t such that s*a + t*b = 1 mod p^k
 * 		    with deg(s) < deg(b) and deg(t) < deg(a).
 * 				Assumption: GCD(a mod p, b mod p) = 1 in Zp[x]
 *
 * @param a A list a of r > 1 polynomials in the domain Zp^k[L]
 * @param b A polynomial in x
 * @param x The symbol x
 * @param p A prime integer
 * @param k A integer
 * @return Expr
 */
ast::Expr EEAlift(ast::Expr a, ast::Expr b, ast::Expr x, Int p, Int k);

/**
 * @brief Computes s1,...,sr such that
 * 				s1*b1 + ... + sr*br = 1 mod p^k
 * 				with deg(sj) < deg(aj) where, in terms
 * 				of the given list of polynomials a1, ..., ar
 * 				the polynomial bi are defined by:
 *
 * 				b[i] = a[1] * ... * a[i-1] * a[i+1] * ... * a[r], i = 1, ..., r.
 *
 * 				Contitions: p must not divide lcoeff(a[i]), i = 1, ..., r
 * 				a[i] mod p, i = 1,...,r, must be pairwise relatively prime
 * 				in Zp[x]
 * @param a A list a of r > 1 polynomials in the domain Zp^k[L]
 * @param L The list variables of a
 * @param p A prime integer p
 * @param k A positive integer k specifying that the coefficient arithmetic is to be performed modulo p^k
 * @return ast::Expr
 */
ast::Expr multiTermEEAlift(ast::Expr a, ast::Expr L, Int p, Int k);

/**
 * @brief Solve in the domain Zp^k[L] the multivariate polynomial diophantine
 * 				equation sig[1]*b[1] + ... + sig[r]*b[r] = x^m mod <I^(d+1), p^k>.
 *
 * @param a A list a of r > 1 polynomials in the domain Zp^k[L]
 * @param c A polynomial in Zp^k[L]
 * @param L The list variables of a and c
 * @param I list of numbers
 * @param d A nonnegative integer specifying the maximum total degree with respect to x2, ..., xv
 * @param p A prime integer p
 * @param k A positive integer k specifying that the coefficient arithmetic is to be performed modulo p^k
 * @return ast::Expr
 */
ast::Expr multivariateDiophant(ast::Expr a, ast::Expr c, ast::Expr L, ast::Expr I, Int d, Int p, Int k);



/**
 * @brief Solve in the domain Zp^k[x] the univariate polynomial diophantine
 * 				equation sig[1]*b[1] + ... + sig[r]*b[r] = x^m mod p^k.
 *
 * @param a A list a of r > 1 polynomials in the domain Zp^k[x]
 * @param L The symbol x
 * @param m A integer
 * @param p A prime integer p
 * @param k A positive integer k specifying that the coefficient arithmetic is to be performed modulo p^k
 * @return ast::Expr
 */
ast::Expr univariateDiophant(ast::Expr a, ast::Expr x, Int m, Int p, Int k);

Int mignotteBound(ast::Expr f, ast::Expr L, ast::Expr K);
Int mignoteExpoent(ast::Expr f, ast::Expr L, ast::Expr K, Int p);

ast::Expr groundPP(ast::Expr f, ast::Expr c, ast::Expr L, ast::Expr K);
ast::Expr groundCont(ast::Expr f, ast::Expr L, ast::Expr K);

ast::Expr sqfFactors(ast::Expr f, ast::Expr x, ast::Expr K);
ast::Expr wangLeadingCoeff(ast::Expr f, ast::Expr delta, ast::Expr u, ast::Expr F, ast::Expr sF, ast::Expr a, ast::Expr L, ast::Expr K);
ast::Expr wangEEZ(ast::Expr f, ast::Expr u, ast::Expr lc, ast::Expr a, Int p, ast::Expr L, ast::Expr K);
}

#endif
