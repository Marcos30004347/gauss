#include <assert.h>
#include <cstdio>

#include "Core/Simplification/Power.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_powers() {
	AST* exp0 = power(integer(2), integer(2));
	AST* exp1 = power(fraction(1, 2), integer(2));
	AST* exp2 = power(integer(2), integer(0));
	AST* exp3 = power(integer(2), integer(1));
	AST* exp4 = power(symbol("x"), integer(0));
	AST* exp5 = power(integer(0), integer(0));
	AST* exp6 = power(integer(0), integer(3));

	AST* res_exp0 = reducePowerAST(exp0);
	AST* res_exp1 = reducePowerAST(exp1);
	AST* res_exp2 = reducePowerAST(exp2);
	AST* res_exp3 = reducePowerAST(exp3);
	AST* res_exp4 = reducePowerAST(exp4);
	AST* res_exp5 = reducePowerAST(exp5);
	AST* res_exp6 = reducePowerAST(exp6);

	assert(res_exp0->kind() == Kind::Integer);
	assert(res_exp0->value() == 4);
	assert(res_exp1->kind() == Kind::Fraction);
	assert(res_exp1->operand(0)->kind() == Kind::Integer);
	assert(res_exp1->operand(1)->kind() == Kind::Integer);
	assert(res_exp1->operand(0)->value() == 1);
	assert(res_exp1->operand(1)->value() == 4);
	assert(res_exp2->kind() == Kind::Integer);
	assert(res_exp2->value() == 1);
	assert(res_exp3->kind() == Kind::Integer);
	assert(res_exp3->value() == 2);
	assert(res_exp4->kind() == Kind::Integer);
	assert(res_exp4->value() == 1);
	assert(res_exp5->kind() == Kind::Undefined);
	// assert(res_exp6->kind() == Kind::Integer);
	// assert(res_exp6->value() == 0);

	delete exp0;
	delete exp1;
	delete exp2;
	delete exp3;
	delete exp4;
	delete exp5;
	delete exp6;
	delete res_exp0;
	delete res_exp1;
	delete res_exp2;
	delete res_exp3;
	delete res_exp4;
	delete res_exp5;
	delete res_exp6;
}


int main() {
	should_simplify_powers();
	return 0;
}
