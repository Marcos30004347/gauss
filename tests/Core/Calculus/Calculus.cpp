#include <assert.h>
#include "Core/Calculus/Calculus.hpp"

using namespace ast;
using namespace algebra;
using namespace calculus;

void should_derivate_expressions()
{
	AST* a = add({
		mul({
			integer(4),
			symbol("x")
		}),
		power(
			symbol("x"),
			integer(2)
		),
		mul({
			integer(5),
			power(
				symbol("x"),
				integer(3)
			)
		})
	});

	AST* x = symbol("x");

	AST* a_dx = derivate(a, x);

	assert(a_dx->kind() == Kind::Addition);
	assert(a_dx->operand(0)->kind() == Kind::Integer);
	assert(a_dx->operand(0)->value() == 4);
	assert(a_dx->operand(1)->kind() == Kind::Multiplication);
	assert(a_dx->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(a_dx->operand(1)->operand(0)->value() == 2);
	assert(a_dx->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(a_dx->operand(1)->operand(1)->identifier() == "x");
	assert(a_dx->operand(2)->kind() == Kind::Multiplication);
	assert(a_dx->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(a_dx->operand(2)->operand(0)->value() == 15);
	assert(a_dx->operand(2)->operand(1)->kind() == Kind::Power);
	assert(a_dx->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(a_dx->operand(2)->operand(1)->operand(0)->identifier() == "x");
	assert(a_dx->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(a_dx->operand(2)->operand(1)->operand(1)->value() == 2);
	
	delete a_dx;
	
	AST* ast_sinh = sinh(x);
	AST* ast_sinh_dx = derivate(ast_sinh, x);

	assert(ast_sinh_dx->kind() == Kind::FunctionCall);
	assert(ast_sinh_dx->funName() == "cosh");
	assert(ast_sinh_dx->numberOfOperands() == 1);
	assert(ast_sinh_dx->operand(0)->kind() == Kind::Symbol);
	assert(ast_sinh_dx->operand(0)->identifier() == "x");

	delete ast_sinh;
	delete ast_sinh_dx;

	AST* ast_cosh = cosh(x);
	AST* ast_cosh_dx = derivate(ast_cosh, x);

	assert(ast_cosh_dx->kind() == Kind::FunctionCall);
	assert(ast_cosh_dx->funName() == "sinh");
	assert(ast_cosh_dx->numberOfOperands() == 1);
	assert(ast_cosh_dx->operand(0)->kind() == Kind::Symbol);
	assert(ast_cosh_dx->operand(0)->identifier() == "x");

	delete ast_cosh;
	delete ast_cosh_dx;

	AST* ast_tanh = tanh(x);
	AST* ast_tanh_dx = derivate(ast_tanh, x);

	assert(ast_tanh_dx->kind() == Kind::Power);
	assert(ast_tanh_dx->operand(0)->kind() == Kind::FunctionCall);
	assert(ast_tanh_dx->operand(0)->funName() == "sech");
	assert(ast_tanh_dx->operand(0)->numberOfOperands() == 1);
	assert(ast_tanh_dx->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(ast_tanh_dx->operand(0)->operand(0)->identifier() == "x");
	assert(ast_tanh_dx->operand(1)->kind() == Kind::Integer);
	assert(ast_tanh_dx->operand(1)->value() == 2);

	delete ast_tanh;
	delete ast_tanh_dx;

	AST* ast_sin = sin(x);
	AST* ast_sin_dx = derivate(ast_sin, x);

	assert(ast_sin_dx->kind() == Kind::FunctionCall);
	assert(ast_sin_dx->funName() == "cos");
	assert(ast_sin_dx->numberOfOperands() == 1);
	assert(ast_sin_dx->operand(0)->kind() == Kind::Symbol);
	assert(ast_sin_dx->operand(0)->identifier() == "x");

	delete ast_sin;
	delete ast_sin_dx;

	AST* ast_cos = cos(x);
	AST* ast_cos_dx = derivate(ast_cos, x);

	assert(ast_cos_dx->kind() == Kind::Multiplication);
	assert(ast_cos_dx->operand(0)->kind() == Kind::Integer);
	assert(ast_cos_dx->operand(0)->value() == -1);
	assert(ast_cos_dx->operand(1)->kind() == Kind::FunctionCall);
	assert(ast_cos_dx->operand(1)->funName() == "sin");
	assert(ast_cos_dx->operand(1)->numberOfOperands() == 1);
	assert(ast_cos_dx->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(ast_cos_dx->operand(1)->operand(0)->identifier() == "x");

	delete ast_cos;
	delete ast_cos_dx;

	AST* ast_tan = tan(x);
	AST* ast_tan_dx = derivate(ast_tan, x);

	assert(ast_tan_dx->kind() == Kind::Power);
	assert(ast_tan_dx->operand(0)->kind() == Kind::FunctionCall);
	assert(ast_tan_dx->operand(0)->funName() == "sec");
	assert(ast_tan_dx->operand(0)->numberOfOperands() == 1);
	assert(ast_tan_dx->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(ast_tan_dx->operand(0)->operand(0)->identifier() == "x");
	assert(ast_tan_dx->operand(1)->kind() == Kind::Integer);
	assert(ast_tan_dx->operand(1)->value() == 2);

	delete ast_tan;
	delete ast_tan_dx;

	AST* ast_cot = cot(x);
	AST* ast_cot_dx = derivate(ast_cot, x);

	assert(ast_cot_dx->kind() == Kind::Multiplication);
	assert(ast_cot_dx->operand(0)->kind() == Kind::Integer);
	assert(ast_cot_dx->operand(0)->value() == -1);
	assert(ast_cot_dx->operand(1)->kind() == Kind::Power);
	assert(ast_cot_dx->operand(1)->operand(0)->kind() == Kind::FunctionCall);
	assert(ast_cot_dx->operand(1)->operand(0)->funName() == "csc");
	assert(ast_cot_dx->operand(1)->operand(0)->numberOfOperands() == 1);
	assert(ast_cot_dx->operand(1)->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(ast_cot_dx->operand(1)->operand(0)->operand(0)->identifier() == "x");
	assert(ast_cot_dx->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(ast_cot_dx->operand(1)->operand(1)->value() == 2);

	delete ast_cot;
	delete ast_cot_dx;

	AST* ast_sec = sec(x);
	AST* ast_sec_dx = derivate(ast_sec, x);

	assert(ast_sec_dx->kind() == Kind::Multiplication);
	assert(ast_sec_dx->operand(0)->kind() == Kind::FunctionCall);
	assert(ast_sec_dx->operand(0)->funName() == "sec");
	assert(ast_sec_dx->operand(0)->numberOfOperands() == 1);
	assert(ast_sec_dx->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(ast_sec_dx->operand(0)->operand(0)->identifier() == "x");
	assert(ast_sec_dx->operand(1)->kind() == Kind::FunctionCall);
	assert(ast_sec_dx->operand(1)->funName() == "tan");
	assert(ast_sec_dx->operand(1)->numberOfOperands() == 1);
	assert(ast_sec_dx->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(ast_sec_dx->operand(1)->operand(0)->identifier() == "x");

	delete ast_sec;
	delete ast_sec_dx;

	AST* ast_csc = csc(x);
	AST* ast_csc_dx = derivate(ast_csc, x);

	assert(ast_csc_dx->kind() == Kind::Multiplication);
	assert(ast_csc_dx->operand(0)->kind() == Kind::Integer);
	assert(ast_csc_dx->operand(0)->value() == -1);
	assert(ast_csc_dx->operand(1)->kind() == Kind::FunctionCall);
	assert(ast_csc_dx->operand(1)->funName() == "cot");
	assert(ast_csc_dx->operand(1)->numberOfOperands() == 1);
	assert(ast_csc_dx->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(ast_csc_dx->operand(1)->operand(0)->identifier() == "x");
	assert(ast_csc_dx->operand(2)->kind() == Kind::FunctionCall);
	assert(ast_csc_dx->operand(2)->funName() == "csc");
	assert(ast_csc_dx->operand(2)->numberOfOperands() == 1);
	assert(ast_csc_dx->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(ast_csc_dx->operand(2)->operand(0)->identifier() == "x");

	delete ast_csc;
	delete ast_csc_dx;

	AST* ast_coth = coth(x);
	AST* ast_coth_dx = derivate(ast_coth, x);
	
	assert(ast_coth_dx->kind() == Kind::Multiplication);
	assert(ast_coth_dx->operand(0)->kind() == Kind::Integer);
	assert(ast_coth_dx->operand(0)->value() == -1);
	assert(ast_coth_dx->operand(1)->kind() == Kind::Power);
	assert(ast_coth_dx->operand(1)->operand(0)->kind() == Kind::FunctionCall);
	assert(ast_coth_dx->operand(1)->operand(0)->funName() == "csch");
	assert(ast_coth_dx->operand(1)->operand(0)->numberOfOperands() == 1);
	assert(ast_coth_dx->operand(1)->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(ast_coth_dx->operand(1)->operand(0)->operand(0)->identifier() == "x");
	assert(ast_coth_dx->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(ast_coth_dx->operand(1)->operand(1)->value() == 2);
	
	delete ast_coth;
	delete ast_coth_dx;

	AST* ast_sech = sech(x);
	AST* ast_sech_dx = derivate(ast_sech, x);
	
	assert(ast_sech_dx->kind() == Kind::Multiplication);
	assert(ast_sech_dx->operand(0)->kind() == Kind::Integer);
	assert(ast_sech_dx->operand(0)->value() == -1);
	assert(ast_sech_dx->operand(1)->kind() == Kind::FunctionCall);
	assert(ast_sech_dx->operand(1)->funName() == "sech");
	assert(ast_sech_dx->operand(1)->numberOfOperands() == 1);
	assert(ast_sech_dx->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(ast_sech_dx->operand(1)->operand(0)->identifier() == "x");
	assert(ast_sech_dx->operand(2)->kind() == Kind::FunctionCall);
	assert(ast_sech_dx->operand(2)->funName() == "tanh");
	assert(ast_sech_dx->operand(2)->numberOfOperands() == 1);
	assert(ast_sech_dx->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(ast_sech_dx->operand(2)->operand(0)->identifier() == "x");

	delete ast_sech;
	delete ast_sech_dx;

	AST* ast_csch = csch(x);
	AST* ast_csch_dx = derivate(ast_csch, x);

	assert(ast_csch_dx->kind() == Kind::Multiplication);
	assert(ast_csch_dx->operand(0)->kind() == Kind::Integer);
	assert(ast_csch_dx->operand(0)->value() == -1);
	assert(ast_csch_dx->operand(1)->kind() == Kind::FunctionCall);
	assert(ast_csch_dx->operand(1)->funName() == "coth");
	assert(ast_csch_dx->operand(1)->numberOfOperands() == 1);
	assert(ast_csch_dx->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(ast_csch_dx->operand(1)->operand(0)->identifier() == "x");
	assert(ast_csch_dx->operand(2)->kind() == Kind::FunctionCall);
	assert(ast_csch_dx->operand(2)->funName() == "csch");
	assert(ast_csch_dx->operand(2)->numberOfOperands() == 1);
	assert(ast_csch_dx->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(ast_csch_dx->operand(2)->operand(0)->identifier() == "x");
	
	delete ast_csch;
	delete ast_csch_dx;

	AST* e0 = power(symbol("x"), fraction(integer(1), integer(2)));
	AST* e0_dx = derivate(e0, x);

	assert(e0_dx->kind() == Kind::Multiplication);
	assert(e0_dx->operand(0)->kind() == Kind::Fraction);
	assert(e0_dx->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(e0_dx->operand(0)->operand(0)->value() == 1);
	assert(e0_dx->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(e0_dx->operand(0)->operand(1)->value() == 2);
	assert(e0_dx->operand(1)->kind() == Kind::Power);
	assert(e0_dx->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(e0_dx->operand(1)->operand(0)->identifier() == "x");
	assert(e0_dx->operand(1)->operand(1)->kind() == Kind::Fraction);
	assert(e0_dx->operand(1)->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(e0_dx->operand(1)->operand(1)->operand(0)->value() == -1);
	assert(e0_dx->operand(1)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(e0_dx->operand(1)->operand(1)->operand(1)->value() == 2);

	delete e0;
	delete e0_dx;

}

int main()
{
	should_derivate_expressions();
}
