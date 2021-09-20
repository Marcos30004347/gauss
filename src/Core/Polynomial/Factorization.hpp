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
ast::AST* trueFactors(ast::AST* u, ast::AST* l, ast::AST* x, ast::AST* p, ast::AST* k);

ast::AST* squareFreeFactor(ast::AST* u, ast::AST* x);

}

#endif
