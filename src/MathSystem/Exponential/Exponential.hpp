#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

#include "MathSystem/Algebra/Algebra.hpp"

namespace exponential {

ast::Expr expandExponential(ast::Expr u);
ast::Expr contractExponential(ast::Expr u);
}

#endif
