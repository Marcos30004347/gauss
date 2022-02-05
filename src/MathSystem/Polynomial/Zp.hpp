#ifndef MATH_ALGEBRA_POLYNOMIALS_ZP_H
#define MATH_ALGEBRA_POLYNOMIALS_ZP_H

#include "Polynomial.hpp"

#include <vector>
#include <utility> 

namespace polynomial {

long modInverse(long a, long p);
long division_Zp(long s, long t, long p);
long division_sZp(long s, long t, long p) ;
long mul_Zp(long s, long t, long p);
long mul_sZp(long s, long t, long p);
long sZp(long b, long m);
long Zp(long b, long m);
long modInverse_sZp(long a, long p);

// Non negative representation of Zm
ast::AST* Zp(ast::AST* u, ast::AST* x, int m);
// Symetric representation of Zm
ast::AST* sZp(ast::AST* u, ast::AST* x, int m);

ast::AST* divideGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* remainderGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* quotientGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* gcdGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* extendedEuclideanAlgGPE_Zp(ast::AST* u, ast::AST* v, ast::AST* x, int p);

ast::AST* divideGPE_sZp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* remainderGPE_sZp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* quotientGPE_sZp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* gcdGPE_sZp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* extendedEuclideanAlgGPE_sZp(ast::AST* u, ast::AST* v, ast::AST* x, int p);
ast::AST* nullSpace_sZp(ast::AST* M, signed long q);

ast::AST* monic_Zp(ast::AST* f, ast::AST* x, long p);
ast::AST* monic_sZp(ast::AST* f, ast::AST* x, long p);
ast::AST* extendedGCDGf(ast::AST* f, ast::AST* g, ast::AST* x, long p);

}
#endif
