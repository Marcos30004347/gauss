#ifndef SIMPLIFICATION_H
#define SIMPLIFICATION_H

#include "MathSystem/AST/AST.hpp"
#include "MathSystem/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceAST(ast::Expr&& u);
ast::Expr reduceAST(ast::Expr& u);

}

#endif
