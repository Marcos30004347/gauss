#include "Core/Algebra/Algebra.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;

void should_create_algebraic_expressions() {
	AST* e0 = inte(3);
	AST* e1 = add({ inte(1), inte(2), inte(3) });
	AST* e2 = sub({ inte(1), inte(2), inte(3) });
	AST* e3 = symb("x");
	AST* e4 = frac(1, 2);
	AST* e5 = pow(symb("x"), inte(4));
	AST* e6 = div(inte(1), inte(4));
	AST* e7 = mul({inte(2), inte(5)});

	assert(e0->kind() == Kind::Integer);
	assert(e0->value() == 3);

	assert(e1->kind() == Kind::Addition);
	assert(e1->operand(0)->kind() == Kind::Integer);
	assert(e1->operand(1)->kind() == Kind::Integer);
	assert(e1->operand(2)->kind() == Kind::Integer);
	assert(e1->operand(0)->value() == 1);
	assert(e1->operand(1)->value() == 2);
	assert(e1->operand(2)->value() == 3);

	assert(e2->kind() == Kind::Subtraction);
	assert(e2->operand(0)->kind() == Kind::Integer);
	assert(e2->operand(1)->kind() == Kind::Integer);
	assert(e2->operand(2)->kind() == Kind::Integer);
	assert(e2->operand(0)->value() == 1);
	assert(e2->operand(1)->value() == 2);
	assert(e2->operand(2)->value() == 3);

	assert(e3->kind() == Kind::Symbol);
	assert(e3->identifier() == "x");

	assert(e4->kind() == Kind::Fraction);
	assert(e4->operand(0)->kind() == Kind::Integer);
	assert(e4->operand(1)->kind() == Kind::Integer);
	assert(e4->operand(0)->value() == 1);
	assert(e4->operand(1)->value() == 2);

	assert(e5->kind() == Kind::Power);
	assert(e5->operand(0)->kind() == Kind::Symbol);
	assert(e5->operand(1)->kind() == Kind::Integer);
	assert(e5->operand(0)->identifier() == "x");
	assert(e5->operand(1)->value() == 4);

	assert(e6->kind() == Kind::Division);
	assert(e6->operand(0)->kind() == Kind::Integer);
	assert(e6->operand(1)->kind() == Kind::Integer);
	assert(e6->operand(0)->value() == 1);
	assert(e6->operand(1)->value() == 4);

	assert(e7->kind() == Kind::Multiplication);
	assert(e7->operand(0)->kind() == Kind::Integer);
	assert(e7->operand(1)->kind() == Kind::Integer);
	assert(e7->operand(0)->value() == 2);
	assert(e7->operand(1)->value() == 5);

	delete e0;
	delete e1;
	delete e2;
	delete e3;
	delete e4;
	delete e5;
	delete e6;
	delete e7;
}

void should_get_info_of_algebraic_expressions() {
	AST* exp0 = inte(3);
	AST* exp1 = frac(4, 2);
	AST* exp2 = pow(inte(5), inte(2));
	
	AST* exp3 = base(exp0);
	AST* exp4 = exp(exp0);
	AST* exp5 = base(exp2);
	AST* exp6 = exp(exp2);
	AST* exp7 = num(exp0);
	AST* exp8 = den(exp0);
	AST* exp9 = num(exp1);
	AST* exp10 = den(exp1);

	assert(exp3->kind() == Kind::Integer);
	assert(exp3->value() == 3);
	assert(exp4->kind() == Kind::Integer);
	assert(exp4->value() == 1);
	assert(exp5->kind() == Kind::Integer);
	assert(exp5->value() == 5);
	assert(exp6->kind() == Kind::Integer);
	assert(exp6->value() == 2);
	assert(exp7->kind() == Kind::Integer);
	assert(exp7->value() == 3);
	assert(exp8->kind() == Kind::Integer);
	assert(exp8->value() == 1);
	assert(exp9->kind() == Kind::Integer);
	assert(exp9->value() == 4);
	assert(exp10->kind() == Kind::Integer);
	assert(exp10->value() == 2);

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

}

void should_order_realate_expressions() {
	AST* exp0 = inte(1);
	AST* exp1 = inte(2);
	AST* exp2 = symb("x");
	AST* exp3 = pow(symb("x"), inte(2));
	AST* exp4 = add({symb("x"), pow(symb("x"), inte(2))});
	AST* exp5 = add({symb("x"), pow(symb("x"), inte(3))});
	AST* exp6 = mul({symb("x"), pow(symb("x"), inte(2))});
	AST* exp7 = mul({symb("x"), pow(symb("x"), inte(3))});

	assert(orderRelation(exp0, exp1));
	assert(orderRelation(exp1, exp2));
	assert(orderRelation(exp2, exp3));
	assert(orderRelation(exp4, exp5));
	assert(orderRelation(exp6, exp7));

	delete exp0;
	delete exp1;
	delete exp2;
	delete exp3;
	delete exp4;
	delete exp5;
	delete exp6;
	delete exp7;
}

int main() {
	should_create_algebraic_expressions();
	should_get_info_of_algebraic_expressions();
	should_order_realate_expressions();
	return 0;
}
