#ifndef SIMPLIFICATION_H
#define SIMPLIFICATION_H

#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceAST(ast::Expr&& u);
ast::Expr reduceAST(ast::Expr& u);

}

#endif
