#ifndef POLYNOMIAL_FACTORIZATION_H
#define POLYNOMIAL_FACTORIZATION_H

#include "Polynomial.hpp"

namespace polynomial {

std::vector<ast::AST*> trueFactors(ast::AST* u, ast::AST* l, ast::AST* x, ast::AST* p, ast::AST* k);

}

#endif
