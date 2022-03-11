#ifndef SIMPLIFICATION_MULTIPLICATION_H
#define SIMPLIFICATION_MULTIPLICATION_H

#include "gauss/AST/AST.hpp"
#include "gauss/Algebra/Algebra.hpp"
#include <vector>

namespace simplification {

ast::Expr reduceMultiplicationAST(ast::Expr u);

ast::Expr reduceMultiplicationExpr(ast::Expr&& u);
ast::Expr reduceMultiplicationExpr(ast::Expr& u);
}

#endif
