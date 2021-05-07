#include <assert.h>
#include <cstdio>

#include "Core/Simplification/Division.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_divisions() {
	AST* exp0 = div(integer(4), integer(2));
	AST* exp1 = div(symbol("x"), symbol("x"));
	AST* exp2 = div(mul({integer(2), symbol("x")}), symbol("x"));
	AST* exp3 = div(power(symbol("x"), integer(2)), symbol("x"));

	AST* res_exp0 = reduceDivisionAST(exp0);
	AST* res_exp1 = reduceDivisionAST(exp1);
	AST* res_exp2 = reduceDivisionAST(exp2);
	AST* res_exp3 = reduceDivisionAST(exp3);

	assert(res_exp0->kind() == Kind::Integer);
	assert(res_exp0->value() == 2);

	assert(res_exp1->kind() == Kind::Integer);
	assert(res_exp1->value() == 1);

	assert(res_exp2->kind() == Kind::Integer);
	assert(res_exp2->value() == 2);

	assert(res_exp3->kind() == Kind::Symbol);
	assert(res_exp3->identifier() == "x");

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
	should_simplify_divisions();
	return 0;
}
