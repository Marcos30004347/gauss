#include <assert.h>
#include <cstdio>

#include "Core/Simplification/Multiplication.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_products() {
	AST* exp0 = mul({ inte(2), inte(2) });
	AST* exp1 = mul({ mul({ inte(2), inte(3) }), mul({ inte(4), inte(5) }) });
	AST* exp2 = mul({ symb("x"), symb("x")});
	AST* exp3 = mul({ mul({inte(2), symb("x")}), mul({inte(2), symb("x")})});

	AST* res_exp0 = reduceMultiplicationAST(exp0);
	AST* res_exp1 = reduceMultiplicationAST(exp1);
	AST* res_exp2 = reduceMultiplicationAST(exp2);
	AST* res_exp3 = reduceMultiplicationAST(exp3);

	assert(res_exp0->kind() == Kind::Integer);
	assert(res_exp0->value() == 4);

	assert(res_exp1->kind() == Kind::Integer);
	assert(res_exp1->value() == 120);

	assert(res_exp2->kind() == Kind::Power);
	assert(res_exp2->operand(0)->kind() == Kind::Symbol);
	assert(res_exp2->operand(0)->identifier() == "x");
	assert(res_exp2->operand(1)->kind() == Kind::Integer);
	assert(res_exp2->operand(1)->value() == 2);

	assert(res_exp3->kind() == Kind::Multiplication);
	assert(res_exp3->operand(0)->kind() == Kind::Integer);
	assert(res_exp3->operand(0)->value() == 4);
	assert(res_exp3->operand(1)->kind() == Kind::Power);
	assert(res_exp3->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(res_exp3->operand(1)->operand(0)->identifier() == "x");
	assert(res_exp3->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(res_exp3->operand(1)->operand(1)->value() == 2);

	delete exp0;
	delete exp1;
	delete exp2;
	delete exp3;

	delete res_exp0;
	delete res_exp1;
	delete res_exp2;
	delete res_exp3;
}

int main() {
	should_simplify_products();
	return 0;
}
