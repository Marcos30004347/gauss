#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

#include "Core/Algebra/Algebra.hpp"

namespace exponential {

ast::Expr expandExponential(ast::Expr u);
ast::Expr contractExponential(ast::Expr u);
}

#endif
