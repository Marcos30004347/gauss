#include <assert.h>
#include <cstdio>

#include "Core/Simplification/Power.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_powers() {
	Expr exp0 = power(integer(2), integer(2));
	Expr exp1 = power(fraction(1, 2), integer(2));
	Expr exp2 = power(integer(2), integer(0));
	Expr exp3 = power(integer(2), integer(1));
	Expr exp4 = power(symbol("x"), integer(0));
	Expr exp5 = power(integer(0), integer(0));
	Expr exp6 = power(integer(0), integer(3));

	Expr res_exp0 = reducePowerAST(exp0);
	Expr res_exp1 = reducePowerAST(exp1);
	Expr res_exp2 = reducePowerAST(exp2);
	Expr res_exp3 = reducePowerAST(exp3);
	Expr res_exp4 = reducePowerAST(exp4);
	Expr res_exp5 = reducePowerAST(exp5);
	Expr res_exp6 = reducePowerAST(exp6);

	assert(res_exp0 == 4);
	assert(res_exp1 == fraction(1, 4));
	assert(res_exp2 == 1);
	assert(res_exp3 == 2);
	assert(res_exp4 == 1);
	assert(res_exp5 == undefined());
}


int main() {
	should_simplify_powers();
	return 0;
}
