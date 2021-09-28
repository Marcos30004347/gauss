#ifndef POLYNOMIAL_FACTORIZATION_H
#define POLYNOMIAL_FACTORIZATION_H

#include "Polynomial.hpp"

namespace polynomial {


void RMatrix(ast::AST* u, ast::AST* x, ast::AST* n_, int p);
void destroyRMatrix(int n);
int getRMatrixValue(int i, int j);

ast::AST* berlekampFactor(ast::AST* u, ast::AST* x, int p);

ast::AST* genExtendSigmaP(ast::AST* V, ast::AST* x, unsigned p);
ast::AST* genExtendRP(ast::AST* V, ast::AST* S, ast::AST* F, ast::AST* x, unsigned p);

ast::AST* polynomialHeight_Z(ast::AST* u, ast::AST* x);

// TODO: TEST
ast::AST* trueFactors(ast::AST* u, ast::AST* l, ast::AST* x, ast::AST* p, ast::AST* k);
ast::AST* henselLift(ast::AST* u, ast::AST* S, ast::AST* x, int p, int k);
ast::AST* irreducibleFactor(ast::AST* u, ast::AST* x, ast::AST* y);
ast::AST* squareFreeFactor(ast::AST* u, ast::AST* x);

/**
 * Algorithm 8.1 from "Algorithms for computer algebra"
 * 
 * Given a primitive polynomial a(x) in R[x], with
 * characteristic zero, computes the square-free
 * factorization of a(x).
 */
ast::AST* squareFreeFactorization(ast::AST* ax, ast::AST* x);

/**
 * Algorithm 8.2 from "Algorithms for computer algebra"
 * 
 * Given a primitive polynomial a(x) in R[x], with
 * characteristic zero, computes the square-free
 * factorization of a(x).
 */
ast::AST* squareFreeFactorization2(ast::AST* ax, ast::AST* x);

/**
 * Algorithm 8.3 from "Algorithms for computer algebra"
 * 
 * Given a primitive polynomial a(x) in GF(q)[x], with
 * GF(q) a Galois field of order q = p^m, 
 * computes the square-free factorization of a(x).
 */
ast::AST* squareFreeFactorizationFiniteField(ast::AST* ax, ast::AST* x, ast::AST* q);

ast::AST* formMatrixQ(ast::AST* ax, ast::AST* x, ast::AST* q);
ast::AST* formMatrixQBinary(ast::AST* ax, ast::AST* x, ast::AST* q);
ast::AST* berlekamp(ast::AST* ax, ast::AST* x, ast::AST* q);

std::pair<ast::AST*, ast::AST*> getPolynomialInZ(ast::AST* ax, ast::AST* x);
ast::AST* irreductibleFactors(ast::AST* ux, ast::AST* x, ast::AST* y);

ast::AST* res(ast::AST* p1, ast::AST* p2, ast::AST* x);

ast::AST* algFactorization(ast::AST* az, ast::AST* z, ast::AST* mx, ast::AST* x, ast::AST* a, ast::AST* y);

}

#endif
