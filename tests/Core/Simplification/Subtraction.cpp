#include <assert.h>
#include <cstdio>

#include "Core/Simplification/Subtraction.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_subtractions() {
	AST* exp0 = sub({
		integer(3),
		integer(2)
	});
	AST* exp1 = sub({
		integer(2),
		integer(3)
	});
	AST* exp2 = sub({
		symbol("x"),
		symbol("x")
	});
	AST* exp3 = sub({
		add({symbol("a"), symbol("b"), symbol("c")}),
		sub({symbol("d"), symbol("e")}),
	});
	AST* exp4 = sub({
		add({symbol("a"), symbol("b"), symbol("c")}),
		sub({symbol("d"), symbol("e")}),
		sub({symbol("f"), symbol("g")}),
	});
	AST* exp5 = sub({
		sub({symbol("a"), symbol("b"), symbol("c")}),
		sub({symbol("d"), symbol("e")}),
		sub({symbol("f"), symbol("g")}),
	});

	AST* res_exp0 = reduceSubtractionAST(exp0);
	AST* res_exp1 = reduceSubtractionAST(exp1);
	AST* res_exp2 = reduceSubtractionAST(exp2);
	AST* res_exp3 = reduceSubtractionAST(exp3);
	AST* res_exp4 = reduceSubtractionAST(exp4);
	AST* res_exp5 = reduceSubtractionAST(exp5);

	assert(res_exp0->kind() == Kind::Integer);
	assert(res_exp0->value() == 1);
	assert(res_exp1->kind() == Kind::Integer);
	assert(res_exp1->value() == -1);
	assert(res_exp2->kind() == Kind::Integer);
	assert(res_exp2->value() == 0);

	assert(res_exp3->kind() == Kind::Addition);
	assert(res_exp3->operand(0)->kind() == Kind::Symbol);
	assert(res_exp3->operand(0)->identifier() == "a");
	assert(res_exp3->operand(1)->kind() == Kind::Symbol);
	assert(res_exp3->operand(1)->identifier() == "b");
	assert(res_exp3->operand(2)->kind() == Kind::Symbol);
	assert(res_exp3->operand(2)->identifier() == "c");
	assert(res_exp3->operand(3)->kind() == Kind::Multiplication);
	assert(res_exp3->operand(3)->operand(0)->kind() == Kind::Integer);
	assert(res_exp3->operand(3)->operand(0)->value() == -1);
	assert(res_exp3->operand(3)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp3->operand(3)->operand(1)->identifier() == "d");
	assert(res_exp3->operand(4)->kind() == Kind::Symbol);
	assert(res_exp3->operand(4)->identifier() == "e");

	assert(res_exp4->kind() == Kind::Addition);
	assert(res_exp4->operand(0)->kind() == Kind::Symbol);
	assert(res_exp4->operand(0)->identifier() == "a");
	assert(res_exp4->operand(1)->kind() == Kind::Symbol);
	assert(res_exp4->operand(1)->identifier() == "b");
	assert(res_exp4->operand(2)->kind() == Kind::Symbol);
	assert(res_exp4->operand(2)->identifier() == "c");
	assert(res_exp4->operand(3)->kind() == Kind::Multiplication);
	assert(res_exp4->operand(3)->operand(0)->kind() == Kind::Integer);
	assert(res_exp4->operand(3)->operand(0)->value() == -1);
	assert(res_exp4->operand(3)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp4->operand(3)->operand(1)->identifier() == "d");
	assert(res_exp4->operand(4)->kind() == Kind::Symbol);
	assert(res_exp4->operand(4)->identifier() == "e");
	assert(res_exp4->operand(5)->kind() == Kind::Multiplication);
	assert(res_exp4->operand(5)->operand(0)->kind() == Kind::Integer);
	assert(res_exp4->operand(5)->operand(0)->value() == -1);
	assert(res_exp4->operand(5)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp4->operand(5)->operand(1)->identifier() == "f");
	assert(res_exp4->operand(6)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp4->operand(6)->operand(1)->identifier() == "g");

	assert(res_exp5->kind() == Kind::Addition);
	assert(res_exp5->operand(0)->kind() == Kind::Symbol);
	assert(res_exp5->operand(0)->identifier() == "a");
	assert(res_exp5->operand(1)->kind() == Kind::Multiplication);
	assert(res_exp5->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(res_exp5->operand(1)->operand(0)->value() == -1);
	assert(res_exp5->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp5->operand(1)->operand(1)->identifier() == "b");
	assert(res_exp5->operand(2)->kind() == Kind::Multiplication);
	assert(res_exp5->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(res_exp5->operand(2)->operand(0)->value() == -1);
	assert(res_exp5->operand(2)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp5->operand(2)->operand(1)->identifier() == "c");
	assert(res_exp5->operand(3)->kind() == Kind::Multiplication);
	assert(res_exp5->operand(3)->operand(0)->kind() == Kind::Integer);
	assert(res_exp5->operand(3)->operand(0)->value() == -1);
	assert(res_exp5->operand(3)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp5->operand(3)->operand(1)->identifier() == "d");
	assert(res_exp5->operand(4)->kind() == Kind::Symbol);
	assert(res_exp5->operand(4)->identifier() == "e");
	assert(res_exp5->operand(5)->kind() == Kind::Multiplication);
	assert(res_exp5->operand(5)->operand(0)->kind() == Kind::Integer);
	assert(res_exp5->operand(5)->operand(0)->value() == -1);
	assert(res_exp5->operand(5)->operand(1)->kind() == Kind::Symbol);
	assert(res_exp5->operand(5)->operand(1)->identifier() == "f");
	assert(res_exp5->operand(6)->kind() == Kind::Symbol);
	assert(res_exp5->operand(6)->identifier() == "g");

	delete exp0;
	delete exp1;
	delete exp2;
	delete exp3;
	delete exp4;
	delete exp5;
	delete res_exp0;
	delete res_exp1;
	delete res_exp2;
	delete res_exp3;
	delete res_exp4;
	delete res_exp5;

}

int main() {
	should_simplify_subtractions();
	return 0;
}
