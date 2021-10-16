#ifndef MATH_ALGEBRA_H
#define MATH_ALGEBRA_H

#include "Core/AST/AST.hpp"

namespace algebra {

bool isConstant(ast::AST* u);

bool isTerminal(ast::AST* u);

bool isRNE(ast::AST* u);

bool orderRelation(ast::AST* a, ast::AST* b);

bool isLessEqZero(ast::AST* u);
bool isLessZero(ast::AST* u);
bool isEqZero(ast::AST* u);
bool isGreaterEqZero(ast::AST* u);
bool isGreaterZero(ast::AST* u);
bool isDivisionByZero(ast::AST* k);

bool isGreatherThan(ast::AST* a, ast::AST* b);
bool isLessThan(ast::AST* a, ast::AST* b);
bool isGreatherOrEqualThan(ast::AST* a, ast::AST* b);
bool isLessOrEqualThan(ast::AST* a, ast::AST* b);

ast::AST* completeSubExpressions(ast::AST* u);

/**
 * Return a and b of u such that u = a*x + b
 */
std::pair<ast::AST*, ast::AST*> linearForm(ast::AST* u, ast::AST* x);

ast::AST* integer(signed long val);

ast::AST* symbol(const char* identifier);

ast::AST* add(std::vector<ast::AST*>);

ast::AST* sub(std::vector<ast::AST*>);

ast::AST* mul(std::vector<ast::AST*>);

ast::AST* div(ast::AST* numerator, ast::AST* denominator);

ast::AST* power(ast::AST* bas, ast::AST* expoent);

ast::AST* fraction(signed long n, signed long d);

ast::AST* fraction(ast::AST* n, ast::AST* d);

ast::AST* factorial(ast::AST* u);

ast::AST* base(ast::AST* u);

ast::AST* expoent(ast::AST* u);

ast::AST* binomial(signed long n, std::vector<signed long> ks);

ast::AST* funCall(const char* id, std::vector<ast::AST*> args);

ast::AST* integerGCD(ast::AST* a, ast::AST* b);

ast::AST* min(ast::AST* a, ast::AST* b);

ast::AST* max(ast::AST* a, ast::AST* b);

ast::AST* undefined();
ast::AST* fail();

ast::AST* leastCommomMultiple(ast::AST* op);

ast::AST* leastCommomMultiple(ast::AST* a, ast::AST* b);

ast::AST* sinh(ast::AST* x);
ast::AST* cosh(ast::AST* x);
ast::AST* tanh(ast::AST* x);
ast::AST* exp(ast::AST* x);
ast::AST* cos(ast::AST* x);
ast::AST* sin(ast::AST* x);
ast::AST* tan(ast::AST* x);
ast::AST* csc(ast::AST* x);
ast::AST* cot(ast::AST* x);
ast::AST* log(ast::AST* x);
ast::AST* ln(ast::AST* x);
ast::AST* sec(ast::AST* x);
ast::AST* coth(ast::AST* x);
ast::AST* sech(ast::AST* x);
ast::AST* csch(ast::AST* x);
ast::AST* abs(ast::AST* x);
ast::AST* arccos(ast::AST* x);
ast::AST* arcsin(ast::AST* x);
ast::AST* arctan(ast::AST* x);
ast::AST* arccot(ast::AST* x);
ast::AST* arcsec(ast::AST* x);
ast::AST* arccsc(ast::AST* x);
ast::AST* arccosh(ast::AST* x);
ast::AST* arctanh(ast::AST* x);


ast::AST* matrix(ast::AST* rows, ast::AST* cols);
ast::AST* matrix(std::vector<ast::AST*> M);

ast::AST* getSymbols(ast::AST* u);

int mod(int a, int b);
long gcd(long a, long b);
long fat(long a);
} // algebra

#endif
