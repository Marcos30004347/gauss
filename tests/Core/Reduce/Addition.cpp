#include <assert.h>
#include <cstdio>

#include "Core/Reduce/Addition.hpp"

using namespace ast;
using namespace reduce;
using namespace algebra;

void should_simplify_additions() {
	AST* exp0 = add({ inte(2), inte(2) });
	AST* exp1 = add({ add({ inte(2), inte(4) }), add({ inte(5), inte(6) }) });
	AST* exp2 = add({ inte(3), add({ inte(5), inte(6) }) });
	AST* exp3 = add({ add({ inte(5), inte(6) }), inte(3) });
	AST* exp4 = add({ add({ inte(2), inte(3) }), add({ inte(4), inte(5) }) });
	AST* exp5 = add({ symb("x"), symb("x")});
	AST* exp6 = add({ add({inte(2), symb("x")}), add({inte(2), symb("x")})});

	AST* res_exp0 = reduceAdditionAST(exp0);
	AST* res_exp1 = reduceAdditionAST(exp1);
	AST* res_exp2 = reduceAdditionAST(exp2);
	AST* res_exp3 = reduceAdditionAST(exp3);
	AST* res_exp4 = reduceAdditionAST(exp4);
	AST* res_exp5 = reduceAdditionAST(exp5);
	AST* res_exp6 = reduceAdditionAST(exp6);

	assert(res_exp0->kind() == Kind::Integer);
	assert(res_exp0->value() == 4);

	assert(res_exp1->kind() == Kind::Integer);
	assert(res_exp1->value() == 17);

	assert(res_exp2->kind() == Kind::Integer);
	assert(res_exp2->value() == 14);

	assert(res_exp3->kind() == Kind::Integer);
	assert(res_exp3->value() == 14);

	assert(res_exp4->kind() == Kind::Integer);
	assert(res_exp4->value() == 14);

	assert(res_exp5->kind() == Kind::Multiplication);
	assert(res_exp5->operand(0)->kind() == Kind::Integer);
	assert(res_exp5->operand(0)->value() == 2);
	assert(res_exp5->operand(1)->kind() == Kind::Symbol);
	assert(res_exp5->operand(1)->identifier() == "x");

	assert(res_exp6->kind() == Kind::Addition);
	assert(res_exp6->operand(0)->kind() == Kind::Integer);
	assert(res_exp6->operand(0)->value() == 4);
	assert(res_exp6->operand(1)->kind() == Kind::Multiplication);
	assert(res_exp6->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(res_exp6->operand(1)->operand(0)->value() == 2);
	assert(res_exp6->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp6->operand(1)->operand(1)->identifier() == "x");

	delete exp0;
	delete res_exp0;
	delete exp1;
	delete res_exp1;
	delete exp2;
	delete res_exp2;
	delete exp3;
	delete res_exp3;
	delete exp4;
	delete res_exp4;
	delete exp5;
	delete res_exp5;
	delete exp6;
	delete res_exp6;
}

int main() {
	should_simplify_additions();
	return 0;
}
