#ifndef MATH_ALGEBRA_SET_H
#define MATH_ALGEBRA_SET_H

#include "MathSystem/AST/AST.hpp"

namespace algebra {

// Returns a set with the elements in the vector e
ast::Expr set(std::vector<ast::Expr> e);

// return a new set with the elements of L that are not in M
ast::Expr difference(ast::Expr L, ast::Expr M);

// Returns a new set with members of L followed by the members of M
ast::Expr unification(ast::Expr L, ast::Expr M);

// Returns the intersection of the elements in L and M
ast::Expr intersection(ast::Expr L, ast::Expr M);

// Return true if the element e exists on the set L
bool exists(ast::Expr L, ast::Expr e);

// Return a set with all combination with size k of the elements of n
ast::Expr combination(ast::Expr n, ast::Expr k);

// Cleanup receive a set of sets C and a set t, and return a new set
// with all sets of C that dont have any element in t
ast::Expr cleanUp(ast::Expr C, ast::Expr t);
}

#endif
