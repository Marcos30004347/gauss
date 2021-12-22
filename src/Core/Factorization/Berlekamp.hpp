#ifndef FACTORIZATION_BERLEKAMP_H
#define FACTORIZATION_BERLEKAMP_H

#include "Core/Polynomial/Polynomial.hpp"

namespace factorization {

ast::Expr buildBerlekampBasis(ast::Expr M, ast::Expr q, bool sym);

ast::Expr buildBerkelampMatrix(ast::Expr ax, ast::Expr x, ast::Expr p, bool sym);

ast::Expr berlekampFactors(ast::Expr sfx, ast::Expr x, ast::Expr p, bool sym = false);

}

#endif
