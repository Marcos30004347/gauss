#ifndef MATH_ALGEBRA_LIST_H
#define MATH_ALGEBRA_LIST_H

#include "Core/AST/AST.hpp"

namespace algebra {

// Returns a list with the elements in the vector e
ast::AST* list(std::vector<ast::AST*> e);

// return a new list with the elements of L that are not in M
ast::AST* remove(ast::AST* L, ast::AST* M);

// Returns a new list with members of L followed by the members of M
ast::AST* join(ast::AST* L, ast::AST* M);

// get the result from f([x, L[0]]) and return a new list with the result appended
// to the beginning of  L, if the result of f is a list, the elements of it will be
// appended to a new List followed be the elements of L. If f is not provided a new
// list is returned with x appended to the beginning of L.
ast::AST* adjoin(ast::AST* x, ast::AST* L, ast::AST* (*f)(ast::AST* const) = nullptr);

// Return a copy of the First element of L
ast::AST* first(ast::AST* L);

// Return a new list with the copy from the i'th element of L onwards
ast::AST* rest(ast::AST* L, int i = 1);

}

#endif
