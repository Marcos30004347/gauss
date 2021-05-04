#ifndef MATH_ALGEBRA_LIST_H
#define MATH_ALGEBRA_LIST_H

#include "Core/AST/AST.hpp"

namespace algebra {

// Returns a list with the elements in the vector e
ast::AST* list(std::vector<ast::AST*> e);

// return a new list with the elements of L that are not in M
ast::AST* listDifference(ast::AST* L, ast::AST* M);

// Returns a new list with members of L followed by the members of M
ast::AST* listJoin(ast::AST* L, ast::AST* M);

// Returns a new list with x adjoined to the beginning of L
ast::AST* listAdjoin(ast::AST* x, ast::AST* L);

// Return a copy of the First element of L
ast::AST* listFirst(ast::AST* L);

// Return a new list with the copy from the i'th element of L onwards
ast::AST* listRest(ast::AST* L, int i = 1);

}

#endif
