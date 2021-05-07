#include <assert.h>
#include <cstdio>

#include "Core/Simplification/Subtraction.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_subtractions() {
	AST* exp0 = sub({integer(3), integer(2)});
	AST* exp1 = sub({integer(2), integer(3)});
	AST* exp2 = sub({symbol("x"), symbol("x")});

	AST* res_exp0 = reduceSubtractionAST(exp0);
	AST* res_exp1 = reduceSubtractionAST(exp1);
	AST* res_exp2 = reduceSubtractionAST(exp2);

	assert(res_exp0->kind() == Kind::Integer);
	assert(res_exp0->value() == 1);
	assert(res_exp1->kind() == Kind::Integer);
	assert(res_exp1->value() == -1);
	assert(res_exp2->kind() == Kind::Integer);
	assert(res_exp2->value() == 0);

	delete exp0;
	delete res_exp0;
	delete exp1;
	delete res_exp1;
	delete exp2;
	delete res_exp2;
}

int main() {
	should_simplify_subtractions();
	return 0;
}
