#ifndef EXPONENTIAL_H
#define EXPONENTIAL_H

#include "Core/Algebra/Algebra.hpp"

namespace exponential {

ast::AST* expandExponential(ast::AST* u);
ast::AST* contractExponential(ast::AST* u);
}

#endif
