#ifndef SIMPLIFICATION_ADDITION_H
#define SIMPLIFICATION_ADDITION_H

#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceAdditionAST(ast::Expr u);
ast::Expr reduceAdditionExpr(ast::Expr u);
}

#endif
