#ifndef SIMPLIFICATION_MULTIPLICATION_H
#define SIMPLIFICATION_MULTIPLICATION_H

#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include <vector>

namespace simplification {

ast::Expr reduceMultiplicationAST(ast::Expr u);

ast::Expr reduceMultiplicationExpr(ast::Expr&& u);
ast::Expr reduceMultiplicationExpr(ast::Expr& u);
}

#endif
