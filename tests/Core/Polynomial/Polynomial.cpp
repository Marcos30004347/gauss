#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Expand/Expand.hpp"	
#include "Core/Algebra/List.hpp"	

#include <assert.h>

using namespace ast;
using namespace expand;
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

	assert(e->kind() == Kind::Addition);
	assert(e->operand(0)->kind() == Kind::Integer);
	assert(e->operand(0)->value() == 1);

	assert(e->operand(1)->kind() == Kind::Multiplication);
	assert(e->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(e->operand(1)->operand(0)->value() == 2);
	assert(e->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(e->operand(1)->operand(1)->identifier() == "t");

	assert(e->operand(2)->kind() == Kind::Multiplication);
	assert(e->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(e->operand(2)->operand(0)->value() == 3);
	assert(e->operand(2)->operand(1)->kind() == Kind::Power);
	assert(e->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(e->operand(2)->operand(1)->operand(0)->identifier() == "t");
	assert(e->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(e->operand(2)->operand(1)->operand(1)->value() == 2);

	assert(e->operand(3)->kind() == Kind::Symbol);
	assert(e->operand(3)->identifier() == "x");

	assert(e->operand(4)->kind() == Kind::Multiplication);
	assert(e->operand(4)->operand(0)->kind() == Kind::Symbol);
	assert(e->operand(4)->operand(0)->identifier() == "t");
	assert(e->operand(4)->operand(1)->kind() == Kind::Symbol);
	assert(e->operand(4)->operand(1)->identifier() == "x");

	assert(e->operand(5)->kind() == Kind::Multiplication);
	assert(e->operand(5)->operand(0)->kind() == Kind::Power);
	assert(e->operand(5)->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(e->operand(5)->operand(0)->operand(0)->identifier() == "t");
	assert(e->operand(5)->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(e->operand(5)->operand(0)->operand(1)->value() == 2);
	assert(e->operand(5)->operand(1)->kind() == Kind::Symbol);
	assert(e->operand(5)->operand(1)->identifier() == "x");

	delete x;
	delete t;
	delete u;
	delete v;
	delete e;
}

void should_get_gcd_polynomials() {
	AST* u = add({
		pow(symb("x"), inte(7)),
		mul({inte(-4), pow(symb("x"), inte(5))}),
		mul({inte(-1), pow(symb("x"), inte(2))}),
		inte(4)
	});

	AST* v = add({
		pow(symb("x"), inte(5)),
		mul({inte(-4), pow(symb("x"), inte(3))}),
		mul({inte(-1), pow(symb("x"), inte(2))}),
		inte(4)
	});

	AST* x = symb("x");

	AST* res = gcdGPE(u, v, x);

	assert(res->kind() == Kind::Addition);
	assert(res->operand(0)->kind() == Kind::Integer);
	assert(res->operand(0)->value() == 4);
	assert(res->operand(1)->kind() == Kind::Multiplication);
	assert(res->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(res->operand(1)->operand(0)->value() == -4);
	assert(res->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(res->operand(1)->operand(1)->identifier() == "x");

	assert(res->operand(2)->kind() == Kind::Multiplication);
	assert(res->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(res->operand(2)->operand(0)->value() == -1);
	assert(res->operand(2)->operand(1)->kind() == Kind::Power);
	assert(res->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(res->operand(2)->operand(1)->operand(0)->identifier() == "x");
	assert(res->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(res->operand(2)->operand(1)->operand(1)->value() == 2);

	assert(res->operand(3)->kind() == Kind::Power);
	assert(res->operand(3)->operand(0)->kind() == Kind::Symbol);
	assert(res->operand(3)->operand(0)->identifier() == "x");
	assert(res->operand(3)->operand(1)->kind() == Kind::Integer);
	assert(res->operand(3)->operand(1)->value() == 3);

	delete u;
	delete v;
	delete x;
	delete res;
}

void should_get_extanded_gcd_polynomials() {
	AST* u = add({
		pow(symb("x"), inte(7)),
		mul({inte(-4), pow(symb("x"), inte(5))}),
		mul({inte(-1), pow(symb("x"), inte(2))}),
		inte(4)
	});

	AST* v = add({
		pow(symb("x"), inte(5)),
		mul({inte(-4), pow(symb("x"), inte(3))}),
		mul({inte(-1), pow(symb("x"), inte(2))}),
		inte(4)
	});

	AST* x = symb("x");

	std::vector<AST*> res = extendedEuclideanAlgGPE(u, v, x);

	AST* gcd = res[0];
	AST* A = res[1];
	AST* B = res[2];

	assert(gcd->kind() == Kind::Addition);
	assert(gcd->operand(0)->kind() == Kind::Integer);
	assert(gcd->operand(0)->value() == 4);
	assert(gcd->operand(1)->kind() == Kind::Multiplication);
	assert(gcd->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(gcd->operand(1)->operand(0)->value() == -4);
	assert(gcd->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(gcd->operand(1)->operand(1)->identifier() == "x");
	assert(gcd->operand(2)->kind() == Kind::Multiplication);
	assert(gcd->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(gcd->operand(2)->operand(0)->value() == -1);
	assert(gcd->operand(2)->operand(1)->kind() == Kind::Power);
	assert(gcd->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(gcd->operand(2)->operand(1)->operand(0)->identifier() == "x");
	assert(gcd->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(gcd->operand(2)->operand(1)->operand(1)->value() == 2);
	assert(gcd->operand(3)->kind() == Kind::Power);
	assert(gcd->operand(3)->operand(0)->kind() == Kind::Symbol);
	assert(gcd->operand(3)->operand(0)->identifier() == "x");
	assert(gcd->operand(3)->operand(1)->kind() == Kind::Integer);
	assert(gcd->operand(3)->operand(1)->value() == 3);

	assert(A->kind() == Kind::Multiplication);
	assert(A->operand(0)->kind() == Kind::Integer);
	assert(A->operand(0)->value() == -1);
	assert(A->operand(1)->kind() == Kind::Symbol);
	assert(A->operand(1)->identifier() == "x");

	assert(B->kind() == Kind::Addition);
	assert(B->operand(0)->kind() == Kind::Integer);
	assert(B->operand(0)->value() == 1);
	assert(B->operand(1)->kind() == Kind::Power);
	assert(B->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(B->operand(1)->operand(0)->identifier() == "x");
	assert(B->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(1)->value() == 3);


	AST* k_ = add({
		mul({ A->deepCopy(), u->deepCopy() }),
		mul({ B->deepCopy(), v->deepCopy() })
	});

	AST* k = expandAST(k_);

	assert(k->match(gcd));

	delete k;
	delete u;
	delete v;
	delete x;
	delete A;
	delete B;
	delete k_;
	delete gcd;
}

void should_calculate_monomial_division() {
	AST* u = add({
		pow(symb("x"), inte(3)),
		mul({
			inte(3),
			pow(symb("x"), inte(2)),
			symb("y")
		}),
		mul({
			inte(4),
			symb("x"),
			pow(symb("y"), inte(2))
		})
	});

	AST* v = add({
		mul({ symb("x"), symb("y") }),
		mul({ inte(2), symb("y") }),
		mul({ inte(3), pow(symb("y"), inte(2)) }),
	});

	AST* L = list({symb("x"), symb("y")});

	AST* R = monomialPolyDiv(u,v,L);
	
	assert(R->operand(0)->kind() == Kind::Addition);
	assert(R->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(0)->operand(0)->value() == -6);
	assert(R->operand(0)->operand(1)->kind() == Kind::Multiplication);
	assert(R->operand(0)->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(0)->operand(1)->operand(0)->value() == 3);
	assert(R->operand(0)->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(R->operand(0)->operand(1)->operand(1)->identifier() == "x");
	assert(R->operand(0)->operand(2)->kind() == Kind::Multiplication);
	assert(R->operand(0)->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(0)->operand(2)->operand(0)->value() == -5);
	assert(R->operand(0)->operand(2)->operand(1)->kind() == Kind::Symbol);
	assert(R->operand(0)->operand(2)->operand(1)->identifier() == "y");

	assert(R->operand(1)->kind() == Kind::Addition);
	assert(R->operand(1)->operand(0)->kind() == Kind::Power);
	assert(R->operand(1)->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(0)->operand(0)->identifier() == "x");
	assert(R->operand(1)->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(0)->operand(1)->value() == 3);
	assert(R->operand(1)->operand(1)->kind() == Kind::Multiplication);
	assert(R->operand(1)->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(1)->operand(0)->value() == 12);
	assert(R->operand(1)->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(1)->operand(1)->identifier() == "y");
	assert(R->operand(1)->operand(2)->kind() == Kind::Multiplication);
	assert(R->operand(1)->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(2)->operand(0)->value() == 28);
	assert(R->operand(1)->operand(2)->operand(1)->kind() == Kind::Power);
	assert(R->operand(1)->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(2)->operand(1)->operand(0)->identifier() == "y");
	assert(R->operand(1)->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(2)->operand(1)->operand(1)->value() == 2);
	assert(R->operand(1)->operand(3)->kind() == Kind::Multiplication);
	assert(R->operand(1)->operand(3)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(3)->operand(0)->value() == 15);
	assert(R->operand(1)->operand(3)->operand(1)->kind() == Kind::Power);
	assert(R->operand(1)->operand(3)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(3)->operand(1)->operand(0)->identifier() == "y");
	assert(R->operand(1)->operand(3)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(3)->operand(1)->operand(1)->value() == 3);

	delete u;
	delete v;
	delete L;
	delete R;
}

void should_get_leading_monomial() {
	AST* u = add({
		mul({inte(3), pow(symb("x"), inte(2)), symb("y")}),
		mul({inte(4), symb("x"), pow(symb("y"), inte(2))}),
		pow(symb("y"), inte(3)),
		symb("x"),
		inte(3)
	});

	AST* L0 = list({ symb("x"), symb("y") });
	AST* L1 = list({ symb("y"), symb("x") });

	AST* lm0 = leadingMonomial(u, L0);
	AST* lm1 = leadingMonomial(u, L1);

	assert(lm0->kind() == Kind::Multiplication);
	assert(lm0->operand(0)->kind() == Kind::Integer);
	assert(lm0->operand(0)->value() == 3);
	assert(lm0->operand(1)->kind() == Kind::Power);
	assert(lm0->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(lm0->operand(1)->operand(0)->identifier() == "x");
	assert(lm0->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(lm0->operand(1)->operand(1)->value() == 2);
	assert(lm0->operand(2)->kind() == Kind::Symbol);
	assert(lm0->operand(2)->identifier() == "y");

	assert(lm1->kind() == Kind::Power);
	assert(lm1->operand(0)->kind() == Kind::Symbol);
	assert(lm1->operand(0)->identifier() == "y");
	assert(lm1->operand(1)->kind() == Kind::Integer);
	assert(lm1->operand(1)->value() == 3);

	delete u;
	delete L0;
	delete L1;
	delete lm0;
	delete lm1;
}

void should_rec_divide_polynomials() {
	AST* u = add({
		mul({
			pow(symb("x"), inte(2)),
			pow(symb("y"), inte(2))
		}),
		symb("x"),
	});
	AST* v = add({
		mul({
			symb("x"),
			symb("y"),
		}),
		inte(1)
	});

	AST* L0 = list({ symb("x"), symb("y") });
	AST* L1 = list({ symb("y"), symb("x") });

	AST* Q = symb("Q");

	AST* R0 = recPolyDiv(u, v, L0, Q);
	AST* R1 = recPolyDiv(u, v, L1, Q);

	assert(R0->operand(0)->kind() == Kind::Multiplication);
	assert(R0->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(R0->operand(0)->operand(0)->identifier() == "x");
	assert(R0->operand(0)->operand(1)->kind() == Kind::Symbol);
	assert(R0->operand(0)->operand(1)->identifier() == "y");
	assert(R0->operand(1)->kind() == Kind::Addition);
	assert(R0->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(R0->operand(1)->operand(0)->identifier() == "x");
	assert(R0->operand(1)->operand(1)->kind() == Kind::Multiplication);
	assert(R0->operand(1)->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(R0->operand(1)->operand(1)->operand(0)->value() == -1);
	assert(R0->operand(1)->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(R0->operand(1)->operand(1)->operand(1)->identifier() == "x");
	assert(R0->operand(1)->operand(1)->operand(2)->kind() == Kind::Symbol);
	assert(R0->operand(1)->operand(1)->operand(2)->identifier() == "y");

	assert(R1->operand(0)->kind() == Kind::Addition);
	assert(R1->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(R1->operand(0)->operand(0)->value() == -1);
	assert(R1->operand(0)->operand(1)->kind() == Kind::Multiplication);
	assert(R1->operand(0)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(R1->operand(0)->operand(1)->operand(0)->identifier() == "x");
	assert(R1->operand(0)->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(R1->operand(0)->operand(1)->operand(1)->identifier() == "y");
	assert(R1->operand(1)->kind() == Kind::Addition);
	assert(R1->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(R1->operand(1)->operand(0)->value() == 1);
	assert(R1->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(R1->operand(1)->operand(1)->identifier() == "x");

	delete u;
	delete v;
	delete L0;
	delete L1;
	delete Q;
	delete R0;
	delete R1;
}

int main() {
	should_get_polynomial_variable();
	should_get_if_is_polynomial_gpe();
	should_get_degree_of_variables();
	should_get_coefficient_gpe();
	should_get_leading_coefficient_gpe();
	should_divided_polynomials();
	should_expand_polynomials();
	should_get_gcd_polynomials();
	should_get_extanded_gcd_polynomials();
	should_calculate_monomial_division();
	should_get_leading_monomial();
	should_rec_divide_polynomials();
	return 0;
}
