#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Expand/Expand.hpp"	
#include "Core/Algebra/List.hpp"	
#include "Core/Algebra/Set.hpp"	

#include <assert.h>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace polynomial;

void should_get_polynomial_variable() {
	AST* exp0 = add({
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

	AST* exp1 = add({
		mul({
			integer(4),
			symbol("x")
		}),
		power(
			symbol("y"),
			integer(2)),
		mul({
			integer(5),
			funCall(
				"sin",
				{ symbol("x")}
			)
		})
	});

	AST* vars_exp0 = variables(exp0);
	AST* vars_exp1 = variables(exp1);

	assert(vars_exp0->numberOfOperands() == 1);
	assert(vars_exp0->operand(0)->kind() == Kind::Symbol);
	assert(vars_exp0->operand(0)->identifier() == "x");

	assert(vars_exp1->numberOfOperands() == 3);
	assert(vars_exp1->operand(0)->kind() == Kind::Symbol);
	assert(vars_exp1->operand(0)->identifier() == "x");
	assert(vars_exp1->operand(1)->kind() == Kind::Symbol);
	assert(vars_exp1->operand(1)->identifier() == "y");
	assert(vars_exp1->operand(2)->kind() == Kind::FunctionCall);
	assert(vars_exp1->operand(2)->funName() == "sin");
	assert(vars_exp1->operand(2)->numberOfOperands() == 1);
	assert(vars_exp1->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(vars_exp1->operand(2)->operand(0)->identifier() == "x");


	delete exp0;
	delete exp1;

	delete vars_exp0;
	delete vars_exp1;
}

void should_get_if_is_polynomial_gpe() {
	AST* exp0 = add({
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

	AST* exp1 = add({
		mul({
			integer(4),
			symbol("x")
		}),
		power(
			symbol("y"),
			integer(2)),
		mul({
			integer(5),
			funCall(
				"sin",
				{ symbol("x")}
			)
		})
	});

	std::vector<AST*> vars_exp0 = { symbol("x") };
	// assert(isPolynomialGPE(exp0, vars_exp0));

	std::vector<AST*> vars_exp1 = { symbol("x"),  symbol("y"), funCall("sin", { symbol("x")})};
	// assert(isPolynomialGPE(exp1, vars_exp1));

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

	AST* exp1 = add({
		mul({
			integer(4),
			symbol("x")
		}),
		power(
			symbol("y"),
			integer(2)
		),
		power(
			funCall(
				"sin",
				{ symbol("x")}
			),
			integer(5)
		)
	});

	AST* x = symbol("x");
	AST* y = symbol("y");
	AST* sin_x = funCall("sin",{ symbol("x") });
	
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
		integer(4),
		power(
			symbol("x"),
			integer(2)
		)
	});

	AST* exp1 = add({
		mul({
			symbol("a"),
			power(
				symbol("x"),
				integer(2)
			)
		}),

		mul({
			symbol("b"),
			power(
				symbol("x"),
				integer(2)
			)
		}),
	});

	AST* exp2 = sub({
		mul({
			symbol("a"),
			power(
				symbol("x"),
				integer(2)
			)
		}),

		mul({
			symbol("b"),
			power(
				symbol("x"),
				integer(2)
			)
		}),
	});

	AST* x = symbol("x");
	AST* p = integer(2);
	AST* coeff_exp0 = coefficientGPE(exp0, x, p);
	AST* coeff_exp1 = coefficientGPE(exp1, x, p);
	AST* coeff_exp2 = coefficientGPE(exp2, x, p);

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
	delete p;
	delete exp0;
	delete exp1;
	delete exp2;
	delete coeff_exp0;
	delete coeff_exp1;
	delete coeff_exp2;
}

void should_get_leading_coefficient_gpe() {
	AST* exp0 = mul({
		integer(4),
		power(
			symbol("x"),
			integer(2)
		)
	});

	AST* exp1 = add({
		mul({
			symbol("a"),
			power(
				symbol("x"),
				integer(2)
			)
		}),

		mul({
			symbol("b"),
			power(
				symbol("x"),
				integer(2)
			)
		}),
	});

	AST* exp2 = add({
		mul({
			symbol("a"),
			power(
				symbol("x"),
				integer(2)
			)
		}),

		mul({
			symbol("b"),
			power(
				symbol("x"),
				integer(3)
			)
		}),
	});

	AST* x = symbol("x");

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
		mul({integer(5), power(symbol("x"), integer(2))}),
		mul({integer(4), symbol("x")}),
		integer(1)
	});

	AST* exp1 = add({
		mul({integer(2), symbol("x")}),
		integer(3)
	});

	AST* x = symbol("x");
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
		power(symbol("x"), integer(5)),
		mul({integer(11), power(symbol("x"), integer(4))}),
		mul({integer(51), power(symbol("x"), integer(3))}),
		mul({integer(124), power(symbol("x"), integer(2))}),
		mul({integer(159), symbol("x")}),
		integer(86)
	});

	AST* v = add({
		power(symbol("x"), integer(2)),
		mul({integer(4), symbol("x")}),
		integer(5)
	});

	AST* x = symbol("x");
	AST* t = symbol("t");
	
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
		power(symbol("x"), integer(7)),
		mul({integer(-4), power(symbol("x"), integer(5))}),
		mul({integer(-1), power(symbol("x"), integer(2))}),
		integer(4)
	});

	AST* v = add({
		power(symbol("x"), integer(5)),
		mul({integer(-4), power(symbol("x"), integer(3))}),
		mul({integer(-1), power(symbol("x"), integer(2))}),
		integer(4)
	});

	AST* x = symbol("x");

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

void should_get_extended_gcd_polynomials() {
	AST* u = add({
		power(symbol("x"), integer(7)),
		mul({integer(-4), power(symbol("x"), integer(5))}),
		mul({integer(-1), power(symbol("x"), integer(2))}),
		integer(4)
	});

	AST* v = add({
		power(symbol("x"), integer(5)),
		mul({integer(-4), power(symbol("x"), integer(3))}),
		mul({integer(-1), power(symbol("x"), integer(2))}),
		integer(4)
	});

	AST* x = symbol("x");

	AST* res = extendedEuclideanAlgGPE(u, v, x);
	printf("%s\n", res->toString().c_str());
	AST* gcd = res->operand(0);
	AST* A = res->operand(1);
	AST* B = res->operand(2);
	// AST* l = mul({A, u});
	// AST* t = mul({B, v});
	// AST* p = add({l, t});

	// AST* c = algebraicExpand(p);

	// printf("%s\n", c->toString().c_str());
	// assert(gcd->kind() == Kind::Addition);
	// assert(gcd->operand(0)->kind() == Kind::Integer);
	// assert(gcd->operand(0)->value() == 4);
	// assert(gcd->operand(1)->kind() == Kind::Multiplication);
	// assert(gcd->operand(1)->operand(0)->kind() == Kind::Integer);
	// assert(gcd->operand(1)->operand(0)->value() == -4);
	// assert(gcd->operand(1)->operand(1)->kind() == Kind::Symbol);
	// assert(gcd->operand(1)->operand(1)->identifier() == "x");
	// assert(gcd->operand(2)->kind() == Kind::Multiplication);
	// assert(gcd->operand(2)->operand(0)->kind() == Kind::Integer);
	// assert(gcd->operand(2)->operand(0)->value() == -1);
	// assert(gcd->operand(2)->operand(1)->kind() == Kind::Power);
	// assert(gcd->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	// assert(gcd->operand(2)->operand(1)->operand(0)->identifier() == "x");
	// assert(gcd->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	// assert(gcd->operand(2)->operand(1)->operand(1)->value() == 2);
	// assert(gcd->operand(3)->kind() == Kind::Power);
	// assert(gcd->operand(3)->operand(0)->kind() == Kind::Symbol);
	// assert(gcd->operand(3)->operand(0)->identifier() == "x");
	// assert(gcd->operand(3)->operand(1)->kind() == Kind::Integer);
	// assert(gcd->operand(3)->operand(1)->value() == 3);

	// assert(A->kind() == Kind::Multiplication);
	// assert(A->operand(0)->kind() == Kind::Integer);
	// assert(A->operand(0)->value() == -1);
	// assert(A->operand(1)->kind() == Kind::Symbol);
	// assert(A->operand(1)->identifier() == "x");

	// assert(B->kind() == Kind::Addition);
	// assert(B->operand(0)->kind() == Kind::Integer);
	// assert(B->operand(0)->value() == 1);
	// assert(B->operand(1)->kind() == Kind::Power);
	// assert(B->operand(1)->operand(0)->kind() == Kind::Symbol);
	// assert(B->operand(1)->operand(0)->identifier() == "x");
	// assert(B->operand(1)->operand(1)->kind() == Kind::Integer);
	// assert(B->operand(1)->operand(1)->value() == 3);


	AST* k_ = add({
		mul({ A->deepCopy(), u->deepCopy() }),
		mul({ B->deepCopy(), v->deepCopy() })
	});

	AST* k = expandAST(k_);

	assert(k->match(gcd));

	delete res;
	delete k;
	delete u;
	delete v;
	delete x;
	delete k_;
}

void should_calculate_monomial_division() {
	AST* u = add({
		power(symbol("x"), integer(3)),
		mul({
			integer(3),
			power(symbol("x"), integer(2)),
			symbol("y")
		}),
		mul({
			integer(4),
			symbol("x"),
			power(symbol("y"), integer(2))
		})
	});

	AST* v = add({
		mul({ symbol("x"), symbol("y") }),
		mul({ integer(2), symbol("y") }),
		mul({ integer(3), power(symbol("y"), integer(2)) }),
	});

	AST* L = list({symbol("x"), symbol("y")});

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
		mul({integer(3), power(symbol("x"), integer(2)), symbol("y")}),
		mul({integer(4), symbol("x"), power(symbol("y"), integer(2))}),
		power(symbol("y"), integer(3)),
		symbol("x"),
		integer(3)
	});

	AST* L0 = list({ symbol("x"), symbol("y") });
	AST* L1 = list({ symbol("y"), symbol("x") });

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
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(2))
		}),
		symbol("x"),
	});
	AST* v = add({
		mul({
			symbol("x"),
			symbol("y"),
		}),
		integer(1)
	});

	AST* L0 = list({ symbol("x"), symbol("y") });
	AST* L1 = list({ symbol("y"), symbol("x") });

	AST* Q = symbol("Q");

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

void should_pseudo_divide_polynomials() {
	AST* u = add({
		mul({integer(5), power(symbol("x"), integer(4)),  power(symbol("y"), integer(3))}),
		mul({integer(3), symbol("x"),  symbol("y")}),
		integer(2)
	});

	AST* v = add({
		mul({integer(2), power(symbol("x"), integer(3)),  symbol("y")}),
		mul({integer(2), symbol("x")}),
		integer(3)
	});

	AST* x = symbol("x");
	
	AST* R = pseudoDivision(u, v, x);

	assert(R->kind() == Kind::List);

	assert(R->operand(0)->kind() == Kind::Multiplication);
	assert(R->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(0)->operand(0)->value() == 10);
	assert(R->operand(0)->operand(1)->kind() == Kind::Symbol);
	assert(R->operand(0)->operand(1)->identifier() == "x");
	assert(R->operand(0)->operand(2)->kind() == Kind::Power);
	assert(R->operand(0)->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(0)->operand(2)->operand(0)->identifier() == "y");
	assert(R->operand(0)->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(0)->operand(2)->operand(1)->value() == 4);

	assert(R->operand(1)->kind() == Kind::Addition);
	assert(R->operand(1)->operand(0)->kind() == Kind::Multiplication);
	assert(R->operand(1)->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(0)->operand(0)->value() == 8);
	assert(R->operand(1)->operand(0)->operand(1)->kind() == Kind::Power);
	assert(R->operand(1)->operand(0)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(0)->operand(1)->operand(0)->identifier() == "y");
	assert(R->operand(1)->operand(0)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(0)->operand(1)->operand(1)->value() == 2);
	assert(R->operand(1)->operand(1)->kind() == Kind::Multiplication);
	assert(R->operand(1)->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(1)->operand(0)->value() == 12);
	assert(R->operand(1)->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(1)->operand(1)->identifier() == "x");
	assert(R->operand(1)->operand(1)->operand(2)->kind() == Kind::Power);
	assert(R->operand(1)->operand(1)->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(1)->operand(2)->operand(0)->identifier() == "y");
	assert(R->operand(1)->operand(1)->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(1)->operand(2)->operand(1)->value() == 3);

	assert(R->operand(1)->operand(2)->kind() == Kind::Multiplication);
	assert(R->operand(1)->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(2)->operand(0)->value() == -30);
	assert(R->operand(1)->operand(2)->operand(1)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(2)->operand(1)->identifier() == "x");
	assert(R->operand(1)->operand(2)->operand(2)->kind() == Kind::Power);
	assert(R->operand(1)->operand(2)->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(2)->operand(2)->operand(0)->identifier() == "y");
	assert(R->operand(1)->operand(2)->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(2)->operand(2)->operand(1)->value() == 4);

	assert(R->operand(1)->operand(3)->kind() == Kind::Multiplication);
	assert(R->operand(1)->operand(3)->operand(0)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(3)->operand(0)->value() == -20);
	assert(R->operand(1)->operand(3)->operand(1)->kind() == Kind::Power);
	assert(R->operand(1)->operand(3)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(3)->operand(1)->operand(0)->identifier() == "x");
	assert(R->operand(1)->operand(3)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(3)->operand(1)->operand(1)->value() == 2);
	assert(R->operand(1)->operand(3)->operand(2)->kind() == Kind::Power);
	assert(R->operand(1)->operand(3)->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(R->operand(1)->operand(3)->operand(2)->operand(0)->identifier() == "y");
	assert(R->operand(1)->operand(3)->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(R->operand(1)->operand(3)->operand(2)->operand(1)->value() == 4);

	delete R;
	delete x;
	delete u;
	delete v;
}

void should_normalize_polynomial() {
	AST* u = add({
		mul({
			add({
				mul({integer(2), symbol("y")}),
				integer(1)
			}),
			symbol("x")
		}),
		add({
			mul({
				integer(6),
				symbol("y")
			}),
			integer(3)
		})
	});
	AST* L = list({ symbol("x"), symbol("y") });
	AST* Q = symbol("Q");
	
	AST* u_ = normalizePoly(u, L, Q);

	assert(u_->kind() == Kind::Addition);
	assert(u_->operand(0)->kind() == Kind::Fraction);
	assert(u_->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(u_->operand(0)->operand(0)->value() == 3);
	assert(u_->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(u_->operand(0)->operand(1)->value() == 2);

	assert(u_->operand(1)->kind() == Kind::Multiplication);
	assert(u_->operand(1)->operand(0)->kind() == Kind::Fraction);
	assert(u_->operand(1)->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(u_->operand(1)->operand(0)->operand(0)->value() == 1);
	assert(u_->operand(1)->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(u_->operand(1)->operand(0)->operand(1)->value() == 2);
	assert(u_->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(u_->operand(1)->operand(1)->identifier() == "x");
	assert(u_->operand(2)->kind() == Kind::Multiplication);
	assert(u_->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(u_->operand(2)->operand(0)->value() == 3);
	assert(u_->operand(2)->operand(1)->kind() == Kind::Symbol);
	assert(u_->operand(2)->operand(1)->identifier() == "y");
	assert(u_->operand(3)->kind() == Kind::Multiplication);
	assert(u_->operand(3)->operand(0)->kind() == Kind::Symbol);
	assert(u_->operand(3)->operand(0)->identifier() == "x");
	assert(u_->operand(3)->operand(1)->kind() == Kind::Symbol);
	assert(u_->operand(3)->operand(1)->identifier() == "y");

	
	delete u;
	delete u_;
	delete L;
	delete Q;
}

void should_mv_poly_gcd() {
	AST* u = add({
		mul({integer(-1), symbol("y"), power(symbol("x"), integer(2))}),
		power(symbol("y"), integer(3))
	});

	AST* v = add({
		mul({ symbol("y"), power(symbol("x"), integer(2)) }),
		mul({ integer(2), power(symbol("y"), integer(2)), symbol("x")}),
		power(symbol("y"), integer(3))
	});

	AST* L = list({ symbol("x"), symbol("y") });
	
	AST* Z = symbol("Z");
	
	AST* gcd = mvPolyGCD(u, v, L, Z);

	assert(gcd->kind() == Kind::Addition);
	assert(gcd->operand(0)->kind() == Kind::Multiplication);
	assert(gcd->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(gcd->operand(0)->operand(0)->identifier() == "x");
	assert(gcd->operand(0)->operand(1)->kind() == Kind::Symbol);
	assert(gcd->operand(0)->operand(1)->identifier() == "y");
	assert(gcd->operand(1)->kind() == Kind::Power);
	assert(gcd->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(gcd->operand(1)->operand(0)->identifier() == "y");
	assert(gcd->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(gcd->operand(1)->operand(1)->value() == 2);

	delete u;
	delete v;
	delete L;
	delete Z;
	delete gcd;
}
void should_get_coeff_var_parts_of_monomial() {
	AST* u = mul({
		integer(4),
		integer(5),
		fraction(1,2),
		symbol("x"),
		power(symbol("x"), integer(2)),
		power(symbol("y"), integer(3)),
	});

	AST* S = set({symbol("x"), symbol("y")});
	
	AST* L = coeffVarMonomial(u, S);

	assert(L->numberOfOperands() == 2);
	assert(L->operand(0)->kind() == Kind::Multiplication);
	assert(L->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(L->operand(0)->operand(0)->value() == 4);
	assert(L->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(L->operand(0)->operand(1)->value() == 5);
	assert(L->operand(0)->operand(2)->kind() == Kind::Fraction);
	assert(L->operand(0)->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(L->operand(0)->operand(2)->operand(0)->value() == 1);
	assert(L->operand(0)->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(L->operand(0)->operand(2)->operand(1)->value() == 2);
	
	assert(L->operand(1)->kind() == Kind::Multiplication);
	assert(L->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(L->operand(1)->operand(0)->identifier() == "x");
	assert(L->operand(1)->operand(1)->kind() == Kind::Power);
	assert(L->operand(1)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(L->operand(1)->operand(1)->operand(0)->identifier() == "x");
	assert(L->operand(1)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(L->operand(1)->operand(1)->operand(1)->value() == 2);
	assert(L->operand(1)->operand(2)->kind() == Kind::Power);
	assert(L->operand(1)->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(L->operand(1)->operand(2)->operand(0)->identifier() == "y");
	assert(L->operand(1)->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(L->operand(1)->operand(2)->operand(1)->value() == 3);

	delete u;
	delete S;
	delete L;
}

void should_collect_terms() {
	AST* u = add({
		mul({symbol("a"), symbol("x")}),
		mul({symbol("b"), symbol("x")}),
		symbol("c"),
		symbol("d"),
	});

	AST* S = set({symbol("x")});

	AST* c = collectTerms(u, S);

	assert(c->kind() == Kind::Addition);
	assert(c->operand(0)->kind() == Kind::Multiplication);
	assert(c->operand(0)->operand(0)->kind() == Kind::Addition);
	assert(c->operand(0)->operand(0)->operand(0)->kind() == Kind::Symbol);
	assert(c->operand(0)->operand(0)->operand(0)->identifier() == "a");
	assert(c->operand(0)->operand(0)->operand(1)->kind() == Kind::Symbol);
	assert(c->operand(0)->operand(0)->operand(1)->identifier() == "b");
	assert(c->operand(0)->operand(1)->kind() == Kind::Symbol);
	assert(c->operand(0)->operand(1)->identifier() == "x");
	assert(c->operand(1)->kind() == Kind::Addition);
	assert(c->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(c->operand(1)->operand(0)->identifier() == "c");
	assert(c->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(c->operand(1)->operand(1)->identifier() == "d");


	delete u;
	delete S;
	delete c;
}

void should_algebraic_expand_expressions() {
	AST* u0 = mul({
		add({
			mul({
				symbol("x"),
				power(
					add({symbol("y"), integer(1)}),
					fraction(3, 2)
				)
			}),
			integer(1)
		}),
		sub({
			mul({
				symbol("x"),
				power(
					add({
						symbol("y"),
						integer(1)
					}),
					fraction(3, 2)
				)
			}),
			integer(1)
		})
	});

	AST* u0_	= algebraicExpand(u0);

	AST* result_u0 = add({
		integer(-1),
		power(symbol("x"), integer(2)),
		mul({integer(3), power(symbol("x"), integer(2)), symbol("y") }),
		mul({integer(3), power(symbol("x"), integer(2)), power(symbol("y"), integer(2)) }),
		mul({ power(symbol("x"), integer(2)), power(symbol("y"), integer(3)) }),
	});


	assert(u0_->match(result_u0));

	AST* u1 = power(
		add({
			mul({
				symbol("x"),
				power(
					add({
						symbol("y"),
						integer(1)
					}),
					fraction(1,2)
				)
			}),
			integer(1)
		}),
		integer(4)
	);

	AST* u1_	= algebraicExpand(u1);

	AST* result_u1 = add({
		integer(1),
		mul({ integer(6), power(symbol("x"), integer(2)) }),
		power(symbol("x"), integer(4)),
		mul({ integer(6), power(symbol("x"), integer(2)), symbol("y") }),
		mul({ integer(2), power(symbol("x"), integer(4)), symbol("y") }),
		mul({ power(symbol("x"), integer(4)), power(symbol("y"), integer(2)) }),
		mul({ integer(4), symbol("x"), power(add({integer(1), symbol("y")}), fraction(1, 2)) }),
		mul({ integer(4), power(symbol("x"), integer(3)), power(add({integer(1), symbol("y")}), fraction(3, 2)) })
	});
	assert(u1_->match(result_u1));

	AST* u2 = mul({
		add({symbol("x"), integer(2)}),
		add({symbol("x"), integer(3)}),
		add({symbol("x"), integer(4)}),
	});

	AST* u2_	= algebraicExpand(u2);

	AST* result_u2 = add({
		integer(24),
		mul({ integer(26), symbol("x") }),
		mul({ integer(9), power(symbol("x"), integer(2)) }),
		power(symbol("x"), integer(3))
	});
	assert(u2_->match(result_u2));

	AST* u3 = power(
		add({symbol("x"), symbol("y"), symbol("z")}),
		integer(3)
	);

	AST* u3_	= algebraicExpand(u3);
	
	AST* result_u3 = add({
		power(symbol("x"), integer(3)),
		mul({ integer(3), power(symbol("x"), integer(2)), symbol("y") }),
		mul({ integer(3), symbol("x"), power(symbol("y"), integer(2)) }),
		power(symbol("y"), integer(3)),
		mul({ integer(3), power(symbol("x"), integer(2)), symbol("z") }),
		mul({ integer(6), symbol("x"), symbol("y"), symbol("z") }),
		mul({ integer(3), power(symbol("y"), integer(2)), symbol("z") }),
		mul({ integer(3), symbol("x"), power(symbol("z"), integer(2)) }),
		mul({ integer(3), symbol("y"), power(symbol("z"), integer(2)) }),
		power(symbol("z"), integer(3)),
	});

	assert(u3_->match(result_u3));

	AST* u4 = add({
		power(
			add({
				symbol("x"),
				integer(1)
			}),
			integer(2)
		),
		power(
			add({
				symbol("y"),
				integer(1)
			}),
			integer(2)
		),
	});

	AST* u4_	= algebraicExpand(u4);

	AST* result_u4 = add({
		integer(2),
		mul({integer(2), symbol("x")}),
		power(symbol("x"), integer(2)),
		mul({ integer(2), symbol("y") }),
		power(symbol("y"), integer(2))
	});

	assert(u4_->match(result_u4));

	AST* u5 = power(
		add({
			power(add({ symbol("x"), integer(2) }), integer(2)), 
			integer(3)
		}),
		integer(2)
	);

	AST* u5_	= algebraicExpand(u5);

	AST* result_u5 = add({
		integer(49),
		mul({integer(56), symbol("x")}),
		mul({integer(30), power(symbol("x"), integer(2))}),
		mul({integer(8), power(symbol("x"), integer(3))}),
		power(symbol("x"), integer(4)),
	});

	assert(u5_->match(result_u5));

	delete u0;
	delete u1;
	delete u2;
	delete u3;
	delete u4;
	delete u5;
	delete u0_;
	delete u1_;
	delete u2_;
	delete u3_;
	delete u4_;
	delete u5_;
	delete result_u0;
	delete result_u1;
	delete result_u2;
	delete result_u3;
	delete result_u4;
	delete result_u5;
}

int main() {

	should_get_polynomial_variable();
	should_get_if_is_polynomial_gpe();
	should_get_degree_of_variables();
	should_get_coefficient_gpe();
	should_get_leading_coefficient_gpe();
	should_divided_polynomials();
	// should_expand_polynomials();
	should_get_gcd_polynomials();
	should_get_extended_gcd_polynomials();
	// should_calculate_monomial_division();
	// should_get_leading_monomial();
	// should_rec_divide_polynomials();
	// should_pseudo_divide_polynomials();
	// should_normalize_polynomial();
	// should_mv_poly_gcd();
	// should_get_coeff_var_parts_of_monomial();
	// should_collect_terms();
	should_algebraic_expand_expressions();

	return 0;
}
