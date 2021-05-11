#ifndef MATH_ALGEBRA_RATIONAL_EXPRESSION_H
#define MATH_ALGEBRA_RATIONAL_EXPRESSION_H

#include "Core/Algebra/Algebra.hpp"

#include <vector>

namespace rational {

bool isRationalExpression(ast::AST* u, ast::AST* S);

ast::AST* rationalVariables(ast::AST* u);

ast::AST* numerator(ast::AST* u);

ast::AST* denominator(ast::AST* u);

ast::AST* rationalize(ast::AST* u);

ast::AST* expandRational(ast::AST* u);

}

#endif
