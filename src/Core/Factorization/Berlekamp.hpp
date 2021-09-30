#ifndef FACTORIZATION_BERLEKAMP_H
#define FACTORIZATION_BERLEKAMP_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

// ast::AST* berlekamp(ast::AST* ax, ast::AST* x, ast::AST* q);

ast::AST* berlekampFactors(ast::AST* u, ast::AST* x, int p);
}

#endif
