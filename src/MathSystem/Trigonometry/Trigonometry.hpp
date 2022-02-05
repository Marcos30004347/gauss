#ifndef TRIGONOMETRY_H
#define TRIGONOMETRY_H

#include "MathSystem/Algebra/Algebra.hpp"

namespace trigonometry {

ast::Expr substituteTrig(ast::Expr u);
ast::Expr expandTrig(ast::Expr u);
ast::Expr contractTrig(ast::Expr u);
}

#endif
