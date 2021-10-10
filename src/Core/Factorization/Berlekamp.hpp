#ifndef FACTORIZATION_BERLEKAMP_H
#define FACTORIZATION_BERLEKAMP_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

ast::AST* buildBerlekampBasis(ast::AST* M, ast::AST* q, bool sym);

ast::AST* buildBerkelampMatrix(ast::AST* ax, ast::AST* x, ast::AST* p, bool sym);

ast::AST* berlekampFactors(ast::AST* sfx, ast::AST* x, ast::AST* p, bool sym = false);

}

#endif
