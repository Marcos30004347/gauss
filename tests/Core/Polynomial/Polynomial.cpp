#include "Core/Polynomial/Polynomial.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;
using namespace polynomial;

void should_get_polynomial_variable() {
	AST* exp0 = add({
		mul({
			inte(4),
			symb("x")
		}),
		pow(
			symb("x"),
			inte(2)
		),
		mul({
			inte(5),
			pow(
				symb("x"),
				inte(3)
			)
		})
	});

	AST* exp1 = add({
		mul({
			inte(4),
			symb("x")
		}),
		pow(
			symb("y"),
			inte(2)),
		mul({
			inte(5),
			funCall(
				"sin",
				{ symb("x")}
			)
		})
	});

	std::vector<AST*> vars_exp0 = variables(exp0);
	std::vector<AST*> vars_exp1 = variables(exp1);

	assert(vars_exp0.size() == 1);
	assert(vars_exp0[0]->kind() == Kind::Symbol);
	assert(vars_exp0[0]->identifier() == "x");

	assert(vars_exp1.size() == 3);
	assert(vars_exp1[0]->kind() == Kind::Symbol);
	assert(vars_exp1[0]->identifier() == "x");
	assert(vars_exp1[1]->kind() == Kind::Symbol);
	assert(vars_exp1[1]->identifier() == "y");
	assert(vars_exp1[2]->kind() == Kind::FunctionCall);
	assert(vars_exp1[2]->funName() == "sin");
	assert(vars_exp1[2]->numberOfOperands() == 1);
	assert(vars_exp1[2]->operand(0)->kind() == Kind::Symbol);
	assert(vars_exp1[2]->operand(0)->identifier() == "x");


	delete exp0;
	delete exp1;

	for(AST* a : vars_exp0)
		delete a;
	for(AST* a : vars_exp1)
		delete a;
}

void should_get_if_is_polynomial_gpe() {
	AST* exp0 = add({
		mul({
			inte(4),
			symb("x")
		}),
		pow(
			symb("x"),
			inte(2)
		),
		mul({
			inte(5),
			pow(
				symb("x"),
				inte(3)
			)
		})
	});

	AST* exp1 = add({
		mul({
			inte(4),
			symb("x")
		}),
		pow(
			symb("y"),
			inte(2)),
		mul({
			inte(5),
			funCall(
				"sin",
				{ symb("x")}
			)
		})
	});

	std::vector<AST*> vars_exp0 = { symb("x") };
	assert(isPolynomialGPE(exp0, vars_exp0));

	std::vector<AST*> vars_exp1 = { symb("x"),  symb("y"), funCall("sin", { symb("x")})};
	assert(isPolynomialGPE(exp1, vars_exp1));

	delete exp0;
	delete exp1;

	for(AST* a : vars_exp0)
		delete a;
	for(AST* a : vars_exp1)
		delete a;
}

void should_get_degree_of_variables() {
	AST* exp0 = add({
		mul({
			inte(4),
			symb("x")
		}),
		pow(
			symb("x"),
			inte(2)
		),
		mul({
			inte(5),
			pow(
				symb("x"),
				inte(3)
			)
		})
	});

	AST* exp1 = add({
		mul({
			inte(4),
			symb("x")
		}),
		pow(
			symb("y"),
			inte(2)
		),
		pow(
			funCall(
				"sin",
				{ symb("x")}
			),
			inte(5)
		)
	});

	AST* x = symb("x");
	AST* y = symb("y");
	AST* sin_x = funCall("sin",{ symb("x") });
	
	AST* degree_exp0_x = degreeGPE(exp0, x );
	AST* degree_exp1_x = degreeGPE(exp1, x);
	AST* degree_exp1_y = degreeGPE(exp1, y );
	AST* degree_exp1_sin_x = degreeGPE(exp1, sin_x);

	assert(degree_exp0_x->kind() == Kind::Integer);
	assert(degree_exp0_x->value() == 3);
	assert(degree_exp1_x->kind() == Kind::Integer);
	assert(degree_exp1_x->value() == 1);
	assert(degree_exp1_y->kind() == Kind::Integer);
	assert(degree_exp1_y->value() == 2);
	assert(degree_exp1_sin_x->kind() == Kind::Integer);
	assert(degree_exp1_sin_x->value() == 5);

	delete x;
	delete y;
	delete sin_x;
	delete exp0;
	delete exp1;
	delete degree_exp0_x;
	delete degree_exp1_x;
	delete degree_exp1_y;
	delete degree_exp1_sin_x;
}

void should_get_coefficient_gpe() {
	AST* exp0 = mul({
		inte(4),
		pow(
			symb("x"),
			inte(2)
		)
	});

	AST* exp1 = add({
		mul({
			symb("a"),
			pow(
				symb("x"),
				inte(2)
			)
		}),

		mul({
			symb("b"),
			pow(
				symb("x"),
				inte(2)
			)
		}),
	});

	AST* exp2 = sub({
		mul({
			symb("a"),
			pow(
				symb("x"),
				inte(2)
			)
		}),

		mul({
			symb("b"),
			pow(
				symb("x"),
				inte(2)
			)
		}),
	});

	AST* x = pow(symb("x"), inte(2));

	AST* coeff_exp0 = coefficientGPE(exp0, x);
	AST* coeff_exp1 = coefficientGPE(exp1, x);
	AST* coeff_exp2 = coefficientGPE(exp2, x);

	assert(coeff_exp0->kind() == Kind::Integer);
	assert(coeff_exp0->value() == 4);

	assert(coeff_exp1->kind() == Kind::Addition);
	assert(coeff_exp1->operand(0)->kind() == Kind::Symbol);
	assert(coeff_exp1->operand(0)->identifier() == "a");
	assert(coeff_exp1->operand(1)->kind() == Kind::Symbol);
	assert(coeff_exp1->operand(1)->identifier() == "b");

	assert(coeff_exp2->kind() == Kind::Subtraction);
	assert(coeff_exp2->operand(0)->kind() == Kind::Symbol);
	assert(coeff_exp2->operand(0)->identifier() == "a");
	assert(coeff_exp2->operand(1)->kind() == Kind::Symbol);
	assert(coeff_exp2->operand(1)->identifier() == "b");


	delete x;
	delete exp0;
	delete exp1;
	delete exp2;
	delete coeff_exp0;
	delete coeff_exp1;
	delete coeff_exp2;
}

void should_get_leading_coefficient_gpe() {
	AST* exp0 = mul({
		inte(4),
		pow(
			symb("x"),
			inte(2)
		)
	});

	AST* exp1 = add({
		mul({
			symb("a"),
			pow(
				symb("x"),
				inte(2)
			)
		}),

		mul({
			symb("b"),
			pow(
				symb("x"),
				inte(2)
			)
		}),
	});

	AST* exp2 = add({
		mul({
			symb("a"),
			pow(
				symb("x"),
				inte(2)
			)
		}),

		mul({
			symb("b"),
			pow(
				symb("x"),
				inte(3)
			)
		}),
	});

	AST* x = symb("x");

	AST* leadcoeff_exp0 = leadingCoefficientGPE(exp0, x);
	AST* leadcoeff_exp1 = leadingCoefficientGPE(exp1, x);
	AST* leadcoeff_exp2 = leadingCoefficientGPE(exp2, x);
	
	assert(leadcoeff_exp0->kind() == Kind::Integer);
	assert(leadcoeff_exp0->value() == 4);

	assert(leadcoeff_exp1->kind() == Kind::Addition);
	assert(leadcoeff_exp1->operand(0)->kind() == Kind::Symbol);
	assert(leadcoeff_exp1->operand(0)->identifier() == "a");
	assert(leadcoeff_exp1->operand(1)->kind() == Kind::Symbol);
	assert(leadcoeff_exp1->operand(1)->identifier() == "b");

	assert(leadcoeff_exp2->kind() == Kind::Symbol);
	assert(leadcoeff_exp2->identifier() == "b");

	delete x;
	delete exp0;
	delete exp1;
	delete exp2;
	delete leadcoeff_exp0;
	delete leadcoeff_exp1;
	delete leadcoeff_exp2;
}

void should_divided_polynomials() {
	AST* exp0 = add({
		mul({inte(5), pow(symb("x"), inte(2))}),
		mul({inte(4), symb("x")}),
		inte(1)
	});

	AST* exp1 = add({
		mul({inte(2), symb("x")}),
		inte(3)
	});

	AST* x = symb("x");
	std::pair<AST*, AST*> res = divideGPE(exp0, exp1, x);

	assert(res.first->kind() == Kind::Addition);
	assert(res.first->operand(0)->kind() == Kind::Fraction);
	assert(res.first->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(res.first->operand(0)->operand(0)->value() == -7);
	assert(res.first->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(res.first->operand(0)->operand(1)->value() == 4);
	assert(res.first->operand(1)->kind() == Kind::Multiplication);
	assert(res.first->operand(1)->operand(0)->kind() == Kind::Fraction);
	assert(res.first->operand(1)->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(res.first->operand(1)->operand(0)->operand(0)->value() == 5);
	assert(res.first->operand(1)->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(res.first->operand(1)->operand(0)->operand(1)->value() == 2);
	assert(res.first->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(res.first->operand(1)->operand(1)->identifier() == "x");

	assert(res.second->kind() == Kind::Fraction);
	assert(res.second->operand(0)->kind() == Kind::Integer);
	assert(res.second->operand(0)->value() == 25);
	assert(res.second->operand(1)->kind() == Kind::Integer);
	assert(res.second->operand(1)->value() == 4);

	delete x;
	delete exp0;
	delete exp1;
	delete res.first;
	delete res.second;
}

void should_expand_polynomials() {
	AST* u = add({
		pow(symb("x"), inte(5)),
		mul({inte(11), pow(symb("x"), inte(4))}),
		mul({inte(51), pow(symb("x"), inte(3))}),
		mul({inte(124), pow(symb("x"), inte(2))}),
		mul({inte(159), symb("x")}),
		inte(86)
	});

	AST* v = add({
		pow(symb("x"), inte(2)),
		mul({inte(4), symb("x")}),
		inte(5)
	});

	AST* x = symb("x");
	AST* t = symb("t");
	
	AST* e = expandGPE(u, v, x, t);

	printf("%s\n", e->toString().c_str());

	delete x;
	delete t;
	delete u;
	delete v;
	delete e;
}

int main() {
	should_get_polynomial_variable();
	should_get_if_is_polynomial_gpe();
	should_get_degree_of_variables();
	should_get_coefficient_gpe();
	should_get_leading_coefficient_gpe();
	should_divided_polynomials();
	should_expand_polynomials();
	return 0;
}
