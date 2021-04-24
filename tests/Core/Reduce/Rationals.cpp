#include <assert.h>
#include <cstdio>

#include "Core/Reduce/Rationals.hpp"

using namespace ast;
using namespace reduce;
using namespace algebra;

void should_simplify_rational_number_expressions() {
	AST* exp0 = add({ inte(1), inte(2) });
	AST* exp1 = sub({ inte(1), inte(2) });
	AST* exp2 = sub({ inte(2), inte(1) });
	AST* exp3 = mul({ inte(2), inte(2) });
	AST* exp4 = div(inte(4), inte(2));
	AST* exp5 = pow(inte(4), inte(2));
	AST* exp6 = add({frac(1,2), frac(1,2)});
	AST* exp7 = sub({frac(1,2), frac(1,2)});
	AST* exp8 = mul({frac(1,2), frac(1,2)});
	AST* exp9 = div(frac(1,2), frac(1,2));
	AST* exp10 = pow(frac(1,2), inte(2));
	AST* exp11 = add({inte(3), frac(1,2)});
	AST* exp12 = add({frac(1,2), inte(3)});
	AST* exp13 = sub({inte(3), frac(1,2)});
	AST* exp14 = sub({frac(1,2), inte(3)});

	AST* res_exp0 = reduceRNEAST(exp0);
	AST* res_exp1 = reduceRNEAST(exp1);
	AST* res_exp2 = reduceRNEAST(exp2);
	AST* res_exp3 = reduceRNEAST(exp3);
	AST* res_exp4 = reduceRNEAST(exp4);
	AST* res_exp5 = reduceRNEAST(exp5);
	AST* res_exp6 = reduceRNEAST(exp6);
	AST* res_exp7 = reduceRNEAST(exp7);
	AST* res_exp8 = reduceRNEAST(exp8);
	AST* res_exp9 = reduceRNEAST(exp9);
	AST* res_exp10 = reduceRNEAST(exp10);
	AST* res_exp11 = reduceRNEAST(exp11);
	AST* res_exp12 = reduceRNEAST(exp12);
	AST* res_exp13 = reduceRNEAST(exp13);
	AST* res_exp14 = reduceRNEAST(exp14);

	assert(res_exp0->kind() == Kind::Integer);
	assert(res_exp0->value() == 3);
	assert(res_exp1->kind() == Kind::Integer);
	assert(res_exp1->value() == -1);
	assert(res_exp2->kind() == Kind::Integer);
	assert(res_exp2->value() == 1);
	assert(res_exp3->kind() == Kind::Integer);
	assert(res_exp3->value() == 4);
	assert(res_exp4->kind() == Kind::Integer);
	assert(res_exp4->value() == 2);
	assert(res_exp5->kind() == Kind::Integer);
	assert(res_exp5->value() == 16);
	assert(res_exp6->kind() == Kind::Integer);
	assert(res_exp6->value() == 1);
	assert(res_exp7->kind() == Kind::Integer);
	assert(res_exp7->value() == 0);
	assert(res_exp8->kind() == Kind::Fraction);
	assert(res_exp8->operand(0)->kind() == Kind::Integer);
	assert(res_exp8->operand(1)->kind() == Kind::Integer);
	assert(res_exp8->operand(0)->value() == 1);
	assert(res_exp8->operand(1)->value() == 4);
	assert(res_exp9->kind() == Kind::Integer);
	assert(res_exp9->value() == 1);
	assert(res_exp10->kind() == Kind::Fraction);
	assert(res_exp10->operand(0)->kind() == Kind::Integer);
	assert(res_exp10->operand(1)->kind() == Kind::Integer);
	assert(res_exp10->operand(0)->value() == 1);
	assert(res_exp10->operand(1)->value() == 4);
	assert(res_exp11->kind() == Kind::Fraction);
	assert(res_exp11->operand(0)->kind() == Kind::Integer);
	assert(res_exp11->operand(1)->kind() == Kind::Integer);
	assert(res_exp11->operand(0)->value() == 7);
	assert(res_exp11->operand(1)->value() == 2);
	assert(res_exp12->kind() == Kind::Fraction);
	assert(res_exp12->operand(0)->kind() == Kind::Integer);
	assert(res_exp12->operand(1)->kind() == Kind::Integer);
	assert(res_exp12->operand(0)->value() == 7);
	assert(res_exp12->operand(1)->value() == 2);
	assert(res_exp13->kind() == Kind::Fraction);
	assert(res_exp13->operand(0)->kind() == Kind::Integer);
	assert(res_exp13->operand(1)->kind() == Kind::Integer);
	assert(res_exp13->operand(0)->value() == 5);
	assert(res_exp13->operand(1)->value() == 2);
	assert(res_exp14->kind() == Kind::Fraction);
	assert(res_exp14->operand(0)->kind() == Kind::Integer);
	assert(res_exp14->operand(1)->kind() == Kind::Integer);
	assert(res_exp14->operand(0)->value() == -5);
	assert(res_exp14->operand(1)->value() == 2);

	delete exp0;
	delete exp1;
	delete exp2;
	delete exp3;
	delete exp4;
	delete exp5;
	delete exp6;
	delete exp7;
	delete exp8;
	delete exp9;
	delete exp10;
	delete exp11;
	delete exp12;
	delete exp13;
	delete exp14;
	delete res_exp0;
	delete res_exp1;
	delete res_exp2;
	delete res_exp3;
	delete res_exp4;
	delete res_exp5;
	delete res_exp6;
	delete res_exp7;
	delete res_exp8;
	delete res_exp9;
	delete res_exp10;
	delete res_exp11;
	delete res_exp12;
	delete res_exp13;
	delete res_exp14;
}

int main() {
	should_simplify_rational_number_expressions();
	return 0;
}
