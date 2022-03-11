#ifndef SIMPLIFICATION_SUBTRACTION_H
#define SIMPLIFICATION_SUBTRACTION_H

#include "gauss/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceSubtractionAST(ast::Expr u);

ast::Expr reduceSubtractionExpr(ast::Expr& u);
ast::Expr reduceSubtractionExpr(ast::Expr&& u);

}

#endif
