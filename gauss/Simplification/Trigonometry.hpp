#ifndef SIMPLIFICATION_TRIGONOMETRIC_H
#define SIMPLIFICATION_TRIGONOMETRIC_H

#include "gauss/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceTrigonometricAST(ast::Expr& u);
ast::Expr reduceTrigonometricAST(ast::Expr&& u);

}

#endif