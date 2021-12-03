#ifndef SIMPLIFICATION_DIVISION_H
#define SIMPLIFICATION_DIVISION_H

#include "Core/Algebra/Algebra.hpp"

namespace simplification {

ast::Expr reduceDivisionAST(ast::Expr& u);
ast::Expr reduceDivisionAST(ast::Expr&& u);

}

#endif
