#include "Core/Algebra/Algebra.hpp"
#include "Core/Expand/Expand.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;
using namespace expand;

void should_expand_polynomials() {
	AST* exp0 = pow(add({symb("x"), inte(2)}), inte(3));
	AST* exp1 = pow(sub({symb("x"), inte(2)}), inte(3));
	AST* res_exp0 = expandAST(exp0);
	AST* res_exp1 = expandAST(exp1);
	
	printf("%s\n", res_exp0->toString().c_str());
	printf("%s\n", res_exp1->toString().c_str());
	
	// assert(res_exp0->kind() == Kind::Addition);
	// assert(res_exp0->operand(0)->kind() == Kind::Integer);
	// assert(res_exp0->operand(0)->value() == 8);

	// assert(res_exp0->operand(1)->kind() == Kind::Multiplication);
	// assert(res_exp0->operand(1)->operand(0)->kind() == Kind::Integer);
	// assert(res_exp0->operand(1)->operand(0)->value() == 12);
	// assert(res_exp0->operand(1)->operand(1)->kind() == Kind::Symbol);
	// assert(res_exp0->operand(1)->operand(1)->identifier() == "x");

	// assert(res_exp0->operand(2)->kind() == Kind::Multiplication);
	// assert(res_exp0->operand(2)->operand(0)->kind() == Kind::Integer);
	// assert(res_exp0->operand(2)->operand(0)->value() == 6);
	// assert(res_exp0->operand(2)->operand(1)->kind() == Kind::Power);
	// assert(res_exp0->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	// assert(res_exp0->operand(2)->operand(1)->operand(0)->identifier() =="x");
	// assert(res_exp0->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	// assert(res_exp0->operand(2)->operand(1)->operand(1)->value() == 2);
	
	// assert(res_exp0->operand(3)->kind() == Kind::Power);
	// assert(res_exp0->operand(3)->operand(0)->kind() == Kind::Symbol);
	// assert(res_exp0->operand(3)->operand(0)->identifier() =="x");
	// assert(res_exp0->operand(3)->operand(1)->kind() == Kind::Integer);
	// assert(res_exp0->operand(3)->operand(1)->value() == 3);


	// assert(res_exp1->kind() == Kind::Addition);
	// assert(res_exp1->operand(0)->kind() == Kind::Integer);
	// assert(res_exp1->operand(0)->value() == -8);

	// assert(res_exp1->operand(1)->kind() == Kind::Multiplication);
	// assert(res_exp1->operand(1)->operand(0)->kind() == Kind::Integer);
	// assert(res_exp1->operand(1)->operand(0)->value() == 12);
	// assert(res_exp1->operand(1)->operand(1)->kind() == Kind::Symbol);
	// assert(res_exp1->operand(1)->operand(1)->identifier() == "x");

	// assert(res_exp1->operand(2)->kind() == Kind::Multiplication);
	// assert(res_exp1->operand(2)->operand(0)->kind() == Kind::Integer);
	// assert(res_exp1->operand(2)->operand(0)->value() == -6);
	// assert(res_exp1->operand(2)->operand(1)->kind() == Kind::Power);
	// assert(res_exp1->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	// assert(res_exp1->operand(2)->operand(1)->operand(0)->identifier() =="x");
	// assert(res_exp1->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	// assert(res_exp1->operand(2)->operand(1)->operand(1)->value() == 2);
	
	// assert(res_exp1->operand(3)->kind() == Kind::Power);
	// assert(res_exp1->operand(3)->operand(0)->kind() == Kind::Symbol);
	// assert(res_exp1->operand(3)->operand(0)->identifier() =="x");
	// assert(res_exp1->operand(3)->operand(1)->kind() == Kind::Integer);
	// assert(res_exp1->operand(3)->operand(1)->value() == 3);
	
	delete exp0;
	delete exp1;
	delete res_exp0;
	delete res_exp1;
}

void should_expand_products() {
	AST* exp0 = mul({add({symb("x"), inte(4)}), add({symb("x"), inte(3)})});
	AST* exp1 = mul({add({symb("x"), inte(4)}), add({symb("x"), inte(3)}), add({symb("x"), inte(5)})});
	
	AST* res_exp0 = expandAST(exp0);
	AST* res_exp1 = expandAST(exp1);

	assert(res_exp0->kind() == Kind::Addition);
	assert(res_exp0->operand(0)->kind() == Kind::Integer);
	assert(res_exp0->operand(0)->value() == 12);
	assert(res_exp0->operand(1)->kind() == Kind::Multiplication);
	assert(res_exp0->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(res_exp0->operand(1)->operand(0)->value() == 7);
	assert(res_exp0->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp0->operand(1)->operand(1)->identifier() == "x");
	assert(res_exp0->operand(2)->kind() == Kind::Power);
	assert(res_exp0->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(res_exp0->operand(2)->operand(0)->identifier() == "x");
	assert(res_exp0->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(res_exp0->operand(2)->operand(1)->value() == 2);
	
	assert(res_exp1->kind() == Kind::Addition);
	assert(res_exp1->operand(0)->kind() == Kind::Integer);
	assert(res_exp1->operand(0)->value() == 60);
	assert(res_exp1->operand(1)->kind() == Kind::Multiplication);
	assert(res_exp1->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(res_exp1->operand(1)->operand(0)->value() == 47);
	assert(res_exp1->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp1->operand(1)->operand(1)->identifier() == "x");
	assert(res_exp1->operand(2)->kind() == Kind::Multiplication);
	assert(res_exp1->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(res_exp1->operand(2)->operand(0)->value() == 12);
	assert(res_exp1->operand(2)->operand(1)->kind() == Kind::Power);
	assert(res_exp1->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(res_exp1->operand(2)->operand(1)->operand(0)->identifier() == "x");
	assert(res_exp1->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(res_exp1->operand(2)->operand(1)->operand(1)->value() == 2);
	assert(res_exp1->operand(3)->kind() == Kind::Power);
	assert(res_exp1->operand(3)->operand(0)->kind() == Kind::Symbol);
	assert(res_exp1->operand(3)->operand(0)->identifier() == "x");
	assert(res_exp1->operand(3)->operand(1)->kind() == Kind::Integer);
	assert(res_exp1->operand(3)->operand(1)->value() == 3);

	delete exp1;
	delete exp0;
	delete res_exp0;
	delete res_exp1;
}

int main() {
	// should_expand_polynomials();
	should_expand_products();

	return 0;
}
