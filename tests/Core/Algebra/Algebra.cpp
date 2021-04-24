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

void should_get_correct_info_of_algebraic_expressions() {

}

int main() {
	should_create_algebraic_expressions();
	should_get_correct_info_of_algebraic_expressions();
	return 0;
}
