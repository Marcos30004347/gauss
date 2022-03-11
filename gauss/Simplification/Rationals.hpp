#ifndef SIMPLIFICATION_RATIONAL_H
#define SIMPLIFICATION_RATIONAL_H

#include "gauss/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceRNEAST(ast::Expr& u);
ast::Expr reduceRNEAST(ast::Expr&& u);

}

#endif