#ifndef TRIGONOMETRY_H
#define TRIGONOMETRY_H

#include "Core/Algebra/Algebra.hpp"

namespace trigonometry {

ast::AST* substituteTrig(ast::AST* u);
ast::AST* expandTrig(ast::AST* u);
ast::AST* contractTrig(ast::AST* u);
}

#endif
