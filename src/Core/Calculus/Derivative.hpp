#ifndef CALCULUS_H
#define CALCULUS_H

#include "Core/Algebra/Algebra.hpp"

namespace calculus {

ast::Expr derivative(ast::Expr u, ast::Expr x);
ast::Expr derivate(ast::Expr u, ast::Expr x);

}

#endif
