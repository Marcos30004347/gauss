#ifndef CALCULUS_H
#define CALCULUS_H

#include "Core/Algebra/Algebra.hpp"

namespace calculus {

ast::AST* derivative(ast::AST* u, ast::AST* x);
ast::AST* derivate(ast::AST* u, ast::AST* x);

}

#endif
