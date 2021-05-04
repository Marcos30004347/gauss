#ifndef MATH_ALGEBRA_SET_H
#define MATH_ALGEBRA_SET_H

#include "Core/AST/AST.hpp"

namespace algebra {

// Returns a set with the elements in the vector e
ast::AST* set(std::vector<ast::AST*> e);

// return a new set with the elements of L that are not in M
ast::AST* setDifference(ast::AST* L, ast::AST* M);

// Returns a new set with members of L followed by the members of M
ast::AST* setUnion(ast::AST* L, ast::AST* M);

// Returns the intersection of the elements in L and M
ast::AST* setIntersection(ast::AST* L, ast::AST* M);

// Return true if the element e exists on the set L
bool setExists(ast::AST* L, ast::AST* e);

// Return a set with all combination with size k of the elements of n
ast::AST* combination(ast::AST* n, int k);

}

#endif
