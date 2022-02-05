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
alg::expr nondivisors(Int G, alg::expr F, Int c, alg::expr L, alg::expr K);
alg::expr nondivisorsPolyExpr(Int G, alg::expr F, Int c, alg::expr L, alg::expr K);

alg::expr getEvaluationPoints(alg::expr& f, alg::expr& G, alg::expr& F, alg::expr& L, alg::expr K, Int p, alg::expr S);
alg::expr getEvaluationPointsPolyExpr(alg::expr& f, alg::expr& G, alg::expr& F, alg::expr& L, alg::expr K, Int p, alg::expr S);

alg::expr factorsWang(alg::expr& f, alg::expr& L, alg::expr K);
alg::expr factorsWangPolyExpr(alg::expr& f, alg::expr& L, alg::expr K);

alg::expr trialDivision(alg::expr& f, alg::expr& F, alg::expr& L, alg::expr K);
alg::expr trialDivisionPolyExpr(alg::expr& f, alg::expr& F, alg::expr& L, alg::expr K);

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
alg::expr EEAlift(alg::expr& a, alg::expr& b, alg::expr& x, Int p, Int k);
alg::expr EEAliftPolyExpr(alg::expr& a, alg::expr& b, alg::expr& x, Int p, Int k);

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
 * @return alg::expr
 */
alg::expr multiTermEEAlift(alg::expr& a, alg::expr& L, Int p, Int k);
alg::expr multiTermEEAliftPolyExpr(alg::expr& a, alg::expr& L, Int p, Int k);

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
 * @return alg::expr
 */
alg::expr multivariateDiophant(alg::expr& a, alg::expr& c, alg::expr& L, alg::expr& I, Int d, Int p, Int k);
alg::expr multivariateDiophantPolyExpr(alg::expr& a, alg::expr& c, alg::expr& L, alg::expr& I, Int d, Int p, Int k, alg::expr K);



/**
 * @brief Solve in the domain Zp^k[x] the univariate polynomial diophantine
 * 				equation sig[1]*b[1] + ... + sig[r]*b[r] = x^m mod p^k.
 *
 * @param a A list a of r > 1 polynomials in the domain Zp^k[x]
 * @param L The symbol x
 * @param m A integer
 * @param p A prime integer p
 * @param k A positive integer k specifying that the coefficient arithmetic is to be performed modulo p^k
 * @return alg::expr
 */
alg::expr univariateDiophant(alg::expr a, alg::expr x, Int m, Int p, Int k);
alg::expr univariateDiophantPolyExpr(alg::expr a, alg::expr L, Int m, Int p, Int k, alg::expr K);

Int mignotteBound(alg::expr f, alg::expr L, alg::expr K);
Int mignotteBoundPolyExpr(alg::expr f, alg::expr L, alg::expr K);
Int mignoteExpoentPolyExpr(alg::expr f, alg::expr L, alg::expr K, Int p);

// alg::expr groundPP(alg::expr f, alg::expr c, alg::expr L, alg::expr K);
// alg::expr groundCont(alg::expr f, alg::expr L, alg::expr K);

alg::expr sqfFactors(alg::expr f, alg::expr x, alg::expr K);
alg::expr sqfFactorsPolyExpr(alg::expr f, alg::expr x, alg::expr K);


alg::expr wangLeadingCoeff(alg::expr f, alg::expr delta, alg::expr u, alg::expr F, alg::expr sF, alg::expr a, alg::expr L, alg::expr K);
alg::expr wangLeadingCoeffPolyExpr(alg::expr f, alg::expr delta, alg::expr u, alg::expr F, alg::expr sF, alg::expr a, alg::expr L, alg::expr K);


alg::expr wangEEZ(alg::expr f, alg::expr u, alg::expr lc, alg::expr a, Int p, alg::expr L, alg::expr K);
alg::expr wangEEZPolyExpr(alg::expr f, alg::expr u, alg::expr lc, alg::expr a, Int p, alg::expr L, alg::expr K);

alg::expr factors(alg::expr f, alg::expr L, alg::expr K);
alg::expr factorsPolyExpr(alg::expr f, alg::expr L, alg::expr K);

}

#endif
