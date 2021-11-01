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
ast::AST* nondivisors(Int G, ast::AST* F, Int c, ast::AST* L, ast::AST* K);

// bool getEvaluationPoint(ast::AST* F, ast::AST* L, ast::AST* K, long mod = 3);
ast::AST* factors(ast::AST* f, ast::AST* L, ast::AST* K);
ast::AST* factorsWang(ast::AST* f, ast::AST* L, ast::AST* K);
ast::AST* trialDivision(ast::AST* f, ast::AST* F, ast::AST* L, ast::AST* K);
ast::AST* groundLeadCoeff(ast::AST* f, ast::AST* L);



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
 * @return AST* 
 */
ast::AST* EEAlift(ast::AST* a, ast::AST* b, ast::AST* x, Int p, Int k);

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
 * @return ast::AST* 
 */
ast::AST* multiTermEEAlift(ast::AST* a, ast::AST* L, Int p, Int k);

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
 * @return ast::AST* 
 */
ast::AST* multivariateDiophant(ast::AST* a, ast::AST* c, ast::AST* L, ast::AST* I, Int d, Int p, Int k);



/**
 * @brief Solve in the domain Zp^k[x] the univariate polynomial diophantine
 * 				equation sig[1]*b[1] + ... + sig[r]*b[r] = x^m mod p^k.
 *
 * @param a A list a of r > 1 polynomials in the domain Zp^k[x]
 * @param L The symbol x
 * @param m A integer
 * @param p A prime integer p
 * @param k A positive integer k specifying that the coefficient arithmetic is to be performed modulo p^k
 * @return ast::AST* 
 */
ast::AST* univariateDiophant(ast::AST* a, ast::AST* x, Int m, Int p, Int k);

Int mignotteBound(ast::AST* f, ast::AST* L, ast::AST* K);
Int mignoteExpoent(ast::AST* f, ast::AST* L, ast::AST* K, Int p);

ast::AST* sqfFactors(ast::AST* f, ast::AST* x, ast::AST* K);
ast::AST* wangLeadingCoeff(ast::AST* f, ast::AST* delta, ast::AST* u, ast::AST* F, ast::AST* sF, ast::AST* a, ast::AST* L, ast::AST* K);

}

#endif
