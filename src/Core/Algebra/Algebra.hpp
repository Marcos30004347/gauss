#ifndef MATH_ALGEBRA_H
#define MATH_ALGEBRA_H

#include "Core/AST/AST.hpp"

namespace algebra {

bool isConstant(ast::Expr u);
bool isTerminal(ast::Expr u);
bool isRNE(ast::Expr u);
bool orderRelation(ast::Expr a, ast::Expr b);
bool isDivisionByZero(ast::Expr k);
ast::Expr completeSubExpressions(ast::Expr u);
std::pair<ast::Expr, ast::Expr> linearForm(ast::Expr u, ast::Expr x);
ast::Expr integer(Int val);
ast::Expr symbol(const char* identifier);
ast::Expr add(std::vector<ast::Expr>);
ast::Expr sub(std::vector<ast::Expr>);
ast::Expr mul(std::vector<ast::Expr>);
ast::Expr div(ast::Expr numerator, ast::Expr denominator);
ast::Expr power(ast::Expr bas, ast::Expr expoent);
ast::Expr fraction(ast::Expr n, ast::Expr d);
ast::Expr factorial(ast::Expr u);
ast::Expr base(ast::Expr u);
ast::Expr expoent(ast::Expr u);
ast::Expr binomial(Int n, std::vector<Int> ks);
ast::Expr funCall(const char* id, std::vector<ast::Expr> args);
ast::Expr integerGCD(ast::Expr a, ast::Expr b);
ast::Expr min(ast::Expr a, ast::Expr b);
ast::Expr max(ast::Expr a, ast::Expr b);
ast::Expr leastCommomMultiple(ast::Expr op);
ast::Expr leastCommomMultiple(ast::Expr a, ast::Expr b);
ast::Expr sinh(ast::Expr x);
ast::Expr cosh(ast::Expr x);
ast::Expr tanh(ast::Expr x);
ast::Expr exp(ast::Expr x);
ast::Expr cos(ast::Expr x);
ast::Expr sin(ast::Expr x);
ast::Expr tan(ast::Expr x);
ast::Expr csc(ast::Expr x);
ast::Expr cot(ast::Expr x);
ast::Expr log(ast::Expr x);
ast::Expr ln(ast::Expr x);
ast::Expr sec(ast::Expr x);
ast::Expr coth(ast::Expr x);
ast::Expr sech(ast::Expr x);
ast::Expr csch(ast::Expr x);
ast::Expr abs(ast::Expr x);
ast::Expr arccos(ast::Expr x);
ast::Expr arcsin(ast::Expr x);
ast::Expr arctan(ast::Expr x);
ast::Expr arccot(ast::Expr x);
ast::Expr arcsec(ast::Expr x);
ast::Expr arccsc(ast::Expr x);
ast::Expr arccosh(ast::Expr x);
ast::Expr arctanh(ast::Expr x);
ast::Expr matrix(ast::Expr rows, ast::Expr cols);
ast::Expr matrix(std::vector<ast::Expr> M);
ast::Expr getSymbols(ast::Expr u);

int mod(int a, int b);
long gcd(long a, long b);
long fat(long a);
} // algebra

#endif
