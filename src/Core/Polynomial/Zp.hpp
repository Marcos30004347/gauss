#ifndef MATH_ALGEBRA_POLYNOMIALS_ZP_H
#define MATH_ALGEBRA_POLYNOMIALS_ZP_H

#include "Polynomial.hpp"

#include <vector>
#include <utility> 

namespace polynomial {

int modInverse_p(int a, int p);
int division_Zp(int s, int t, int p);
int division_Sp(int s, int t, int p) ;
int mul_Zp(int s, int t, int p);
int mul_Sp(int s, int t, int p);
int S(int b, int m);

// Non negative representation of Zm
ast::AST* Tnn(ast::AST* u, ast::AST* x, int m);
// Symetric representation of Zm
ast::AST* Ts(ast::AST* u, ast::AST* x, int m);

ast::AST* divideGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* remainderGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* quotientGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* gcdGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* extendedEuclideanAlgGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);

ast::AST* divideGPE_Sp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* remainderGPE_Sp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* quotientGPE_Sp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* gcdGPE_Sp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* extendedEuclideanAlgGPE_Sp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* nullSpace_Sp(ast::AST* M, signed long q);


}
#endif
