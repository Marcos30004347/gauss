#ifndef SIMPLIFICATION_ADDITION_H
#define SIMPLIFICATION_ADDITION_H

#include "MathSystem/AST/AST.hpp"
#include "MathSystem/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceAdditionAST(ast::Expr u);
ast::Expr reduceAdditionExpr(ast::Expr&& u);
ast::Expr reduceAdditionExpr(ast::Expr& u);
}

#endif
