#ifndef FACTORIZATION_BERLEKAMP_H
#define FACTORIZATION_BERLEKAMP_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

ast::AST* nullSpace_Zp(ast::AST* M, long q);

ast::AST* buildBerkelampBasisMatrix(ast::AST* ax, ast::AST* x, ast::AST* p);

ast::AST* berlekamp(ast::AST* ax, ast::AST* x, ast::AST* q);

ast::AST* berlekampFactors(ast::AST* u, ast::AST* x, int p);
}

#endif
