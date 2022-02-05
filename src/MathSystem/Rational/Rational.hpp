#ifndef MATH_ALGEBRA_RATIONAL_EXPRESSION_H
#define MATH_ALGEBRA_RATIONAL_EXPRESSION_H

#include "MathSystem/Algebra/Algebra.hpp"

#include <vector>

namespace rational {

bool isRationalExpression(ast::Expr u, ast::Expr S);

ast::Expr rationalVariables(ast::Expr u);

ast::Expr numerator(ast::Expr u);

ast::Expr denominator(ast::Expr u);

ast::Expr rationalize(ast::Expr u);

ast::Expr expandRational(ast::Expr u);

}

#endif
