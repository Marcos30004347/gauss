#ifndef SIMPLIFICATION_TRIGONOMETRIC_H
#define SIMPLIFICATION_TRIGONOMETRIC_H

#include "MathSystem/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceTrigonometricAST(ast::Expr& u);
ast::Expr reduceTrigonometricAST(ast::Expr&& u);

}

#endif
