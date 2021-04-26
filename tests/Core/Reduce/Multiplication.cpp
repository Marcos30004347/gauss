#include <assert.h>
#include <cstdio>

#include "Core/Reduce/Multiplication.hpp"

using namespace ast;
using namespace reduce;
using namespace algebra;

void should_simplify_products() {
	AST* exp0 = mul({ inte(2), inte(2) });
	AST* exp1 = mul({ add({ inte(2), inte(4) }), add({ inte(5), inte(6) }) });
	AST* exp2 = mul({ inte(3), add({ inte(5), inte(6) }) });
	AST* exp3 = mul({ add({ inte(5), inte(6) }), inte(3) });
	AST* exp4 = mul({ mul({ inte(2), inte(3) }), mul({ inte(4), inte(5) }) });
	AST* exp5 = mul({ symb("x"), symb("x")});
	AST* exp6 = mul({ mul({inte(2), symb("x")}), mul({inte(2), symb("x")})});

	AST* res_exp0 = reduceMultiplicationAST(exp0);
	AST* res_exp1 = reduceMultiplicationAST(exp1);
	AST* res_exp2 = reduceMultiplicationAST(exp2);
	AST* res_exp3 = reduceMultiplicationAST(exp3);
	AST* res_exp4 = reduceMultiplicationAST(exp4);
	AST* res_exp5 = reduceMultiplicationAST(exp5);
	AST* res_exp6 = reduceMultiplicationAST(exp6);

	assert(res_exp0->kind() == Kind::Integer);
	assert(res_exp0->value() == 4);

	assert(res_exp1->kind() == Kind::Integer);
	assert(res_exp1->value() == 66);

	assert(res_exp2->kind() == Kind::Integer);
	assert(res_exp2->value() == 33);

	assert(res_exp3->kind() == Kind::Integer);
	assert(res_exp3->value() == 33);

	assert(res_exp4->kind() == Kind::Integer);
	assert(res_exp4->value() == 120);

	assert(res_exp5->kind() == Kind::Power);
	assert(res_exp5->operand(0)->kind() == Kind::Symbol);
	assert(res_exp5->operand(0)->identifier() == "x");
	assert(res_exp5->operand(1)->kind() == Kind::Integer);
	assert(res_exp5->operand(1)->value() == 2);

	assert(res_exp6->kind() == Kind::Multiplication);
	assert(res_exp6->operand(0)->kind() == Kind::Integer);
	assert(res_exp6->operand(0)->value() == 4);
	assert(res_exp6->operand(1)->kind() == Kind::Power);
	assert(res_exp6->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(res_exp6->operand(1)->operand(0)->identifier() == "x");
	assert(res_exp6->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(res_exp6->operand(1)->operand(1)->value() == 2);

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
	should_simplify_products();
	return 0;
}
