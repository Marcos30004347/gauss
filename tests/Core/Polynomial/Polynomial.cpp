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
	
	AST* degree_exp0_x = degree(exp0, x );
	AST* degree_exp1_x = degree(exp1, x);
	AST* degree_exp1_y = degree(exp1, y );
	AST* degree_exp1_sin_x = degree(exp1, sin_x);

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
	AST* coeff_exp0 = coeff(exp0, x, p);
	AST* coeff_exp1 = coeff(exp1, x, p);
	AST* coeff_exp2 = coeff(exp2, x, p);

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

	AST* leadcoeff_exp0 = leadCoeff(exp0, x);
	AST* leadcoeff_exp1 = leadCoeff(exp1, x);
	AST* leadcoeff_exp2 = leadCoeff(exp2, x);
	
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
	AST* res = divideGPE(exp0, exp1, x);

	assert(res->operand(0)->kind() == Kind::Addition);
	assert(res->operand(0)->operand(0)->kind() == Kind::Fraction);
	assert(res->operand(0)->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(res->operand(0)->operand(0)->operand(0)->value() == -7);
	assert(res->operand(0)->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(res->operand(0)->operand(0)->operand(1)->value() == 4);
	assert(res->operand(0)->operand(1)->kind() == Kind::Multiplication);
	assert(res->operand(0)->operand(1)->operand(0)->kind() == Kind::Fraction);
	assert(res->operand(0)->operand(1)->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(res->operand(0)->operand(1)->operand(0)->operand(0)->value() == 5);
	assert(res->operand(0)->operand(1)->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(res->operand(0)->operand(1)->operand(0)->operand(1)->value() == 2);
	assert(res->operand(0)->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(res->operand(0)->operand(1)->operand(1)->identifier() == "x");
	assert(res->operand(1)->kind() == Kind::Fraction);
	assert(res->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(res->operand(1)->operand(0)->value() == 25);
	assert(res->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(res->operand(1)->operand(1)->value() == 4);

	delete x;
	delete exp0;
	delete exp1;
	delete res;
}

// void should_expand_polynomials() {
// 	AST* u = add({
// 		power(symbol("x"), integer(5)),
// 		mul({integer(11), power(symbol("x"), integer(4))}),
// 		mul({integer(51), power(symbol("x"), integer(3))}),
// 		mul({integer(124), power(symbol("x"), integer(2))}),
// 		mul({integer(159), symbol("x")}),
// 		integer(86)
// 	});

// 	AST* v = add({
// 		power(symbol("x"), integer(2)),
// 		mul({integer(4), symbol("x")}),
// 		integer(5)
// 	});

// 	AST* x = symbol("x");
// 	AST* t = symbol("t");
	
// 	AST* e = expandGPE(u, v, x, t);

// 	assert(e->kind() == Kind::Addition);
// 	assert(e->operand(0)->kind() == Kind::Integer);
// 	assert(e->operand(0)->value() == 1);

// 	assert(e->operand(1)->kind() == Kind::Multiplication);
// 	assert(e->operand(1)->operand(0)->kind() == Kind::Integer);
// 	assert(e->operand(1)->operand(0)->value() == 2);
// 	assert(e->operand(1)->operand(1)->kind() == Kind::Symbol);
// 	assert(e->operand(1)->operand(1)->identifier() == "t");

// 	assert(e->operand(2)->kind() == Kind::Multiplication);
// 	assert(e->operand(2)->operand(0)->kind() == Kind::Integer);
// 	assert(e->operand(2)->operand(0)->value() == 3);
// 	assert(e->operand(2)->operand(1)->kind() == Kind::Power);
// 	assert(e->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
// 	assert(e->operand(2)->operand(1)->operand(0)->identifier() == "t");
// 	assert(e->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
// 	assert(e->operand(2)->operand(1)->operand(1)->value() == 2);

// 	assert(e->operand(3)->kind() == Kind::Symbol);
// 	assert(e->operand(3)->identifier() == "x");

// 	assert(e->operand(4)->kind() == Kind::Multiplication);
// 	assert(e->operand(4)->operand(0)->kind() == Kind::Symbol);
// 	assert(e->operand(4)->operand(0)->identifier() == "t");
// 	assert(e->operand(4)->operand(1)->kind() == Kind::Symbol);
// 	assert(e->operand(4)->operand(1)->identifier() == "x");

// 	assert(e->operand(5)->kind() == Kind::Multiplication);
// 	assert(e->operand(5)->operand(0)->kind() == Kind::Power);
// 	assert(e->operand(5)->operand(0)->operand(0)->kind() == Kind::Symbol);
// 	assert(e->operand(5)->operand(0)->operand(0)->identifier() == "t");
// 	assert(e->operand(5)->operand(0)->operand(1)->kind() == Kind::Integer);
// 	assert(e->operand(5)->operand(0)->operand(1)->value() == 2);
// 	assert(e->operand(5)->operand(1)->kind() == Kind::Symbol);
// 	assert(e->operand(5)->operand(1)->identifier() == "x");

// 	delete x;
// 	delete t;
// 	delete u;
// 	delete v;
// 	delete e;
// }

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
		power(
			symbol("x"),
			integer(7)
		),
		mul({
			integer(-4),
			power(
				symbol("x"),
				integer(5)
			)
		}),
		mul({
			integer(-1),
			power(
				symbol("x"),
				integer(2)
			)
		}),
		integer(4)
	});

	AST* v = add({
		power(
			symbol("x"),
			integer(5)
		),
		mul({
			integer(-4),
			power(
				symbol("x"),
				integer(3)
			)}
		),
		mul({
			integer(-1),
			power(
				symbol("x"),
				integer(2)
			)
		}),
		integer(4)
	});

	AST* x = symbol("x");

	AST* res = extendedEuclideanAlgGPE(u, v, x);

	AST* gcd = res->operand(0);
	AST* A = res->operand(1);
	AST* B = res->operand(2);

	AST* a = mul({ integer(-1), symbol("x") });
	AST* b = add({ integer(1), power(symbol("x"), integer(3)) });

	assert(A->match(a));
	assert(B->match(b));

	AST* k_ = add({
		mul({ A->copy(), u->copy() }),
		mul({ B->copy(), v->copy() })
	});

	AST* k = expandAST(k_);

	assert(k->match(gcd));

	delete a;
	delete b;
	delete k;
	delete u;
	delete v;
	delete x;
	delete k_;
	delete res;
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

	AST* lm0 = leadMonomial(u, L0);
	AST* lm1 = leadMonomial(u, L1);

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

	AST* R0_0 = mul({
		symbol("x"),
		symbol("y")
	});

	AST* R0_1 = add({
		symbol("x"),
		mul({
			integer(-1),
			symbol("x"),
			symbol("y")
		})
	});

	assert(R0->operand(0)->match(R0_0));
	assert(R0->operand(1)->match(R0_1));

	AST* R1_0 = add({
		integer(-1),
		mul({
			symbol("x"),
			symbol("y")
		})
	});

	AST* R1_1 = add({
		integer(1),
		symbol("x")
	});

	assert(R1->operand(0)->match(R1_0));
	assert(R1->operand(1)->match(R1_1));

	delete u;
	delete v;
	delete L0;
	delete L1;
	delete Q;
	delete R0;
	delete R0_0;
	delete R0_1;
	delete R1;
	delete R1_0;
	delete R1_1;
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
	assert(R->numberOfOperands() == 2);

	AST* R0 = mul({
		integer(10),
		symbol("x"),
		power(
			symbol("y"),
			integer(4)
		)
	});

	assert(R->operand(0)->match(R0));

	AST* R1 = add({
		mul({
			integer(8),
			power(
				symbol("y"),
				integer(2)
			)
		}),
		mul({
			integer(12),
			symbol("x"),
			power(
				symbol("y"),
				integer(3)
			)
		}),
		mul({
			integer(-30),
			symbol("x"),
			power(
				symbol("y"),
				integer(4)
			)
		}),
		mul({
			integer(-20),
			power(
				symbol("x"),
				integer(2)
			),
			power(
				symbol("y"),
				integer(4)
			)
		})
	});

	assert(R->operand(1)->match(R1));

	delete R;
	delete x;
	delete u;
	delete v;
	delete R0;
	delete R1;
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


	AST* u_res = add({
		fraction(3,2),
		mul({
			fraction(1,2),
			symbol("x")
		}),
		mul({
			integer(3),
			symbol("y"),
		}),
		mul({
			symbol("x"),
			symbol("y"),
		})
	});

	assert(u_->match(u_res));

	delete u;
	delete u_;
	delete u_res;
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

	AST* gcd_res = add({
		mul({
			symbol("x"),
			symbol("y")
		}),
		power(
			symbol("y"),
			integer(2)
		)
	});

	assert(gcd->match(gcd_res));

	delete u;
	delete v;
	delete L;
	delete Z;
	delete gcd;
	delete gcd_res;
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

	AST* L_0 = mul({
		integer(4),
		integer(5),
		fraction(1,2)
	});
	AST* L_1 = mul({
		symbol("x"),
		power(
			symbol("x"),
			integer(2)
		),
		power(
			symbol("y"),
			integer(3)
		)
	});


	assert(L->kind() == Kind::List);
	assert(L->numberOfOperands() == 2);

	assert(L->operand(0)->match(L_0));
	assert(L->operand(1)->match(L_1));

	delete u;
	delete S;
	delete L;
	delete L_0;
	delete L_1;
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

	AST* c_res = add({
		mul({
			add({
				symbol("a"),
				symbol("b"),
			}),
			symbol("x")
		}),
		add({
			symbol("c"),
			symbol("d"),
		})
	});

	assert(c->match(c_res));

	delete u;
	delete S;
	delete c;
	delete c_res;
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

	AST* u6 = div(
		add({
			mul({integer(-32), power(symbol("z"), integer(3))}),
			mul({integer(32), power(symbol("z"), integer(4))}),
			mul({integer(48), power(symbol("z"), integer(5))}),
			mul({integer(-24), power(symbol("z"), integer(6))}),
			mul({integer(-48), power(symbol("z"), integer(7))}),
			mul({integer(-36), power(symbol("z"), integer(8))}),
			mul({integer(-40), power(symbol("z"), integer(9))}),
			mul({integer(-8), power(symbol("z"), integer(10))}),
			mul({integer(-8), power(symbol("z"), integer(11))}),
		}),
		mul({integer(4), power(symbol("z"), integer(2))})
	);

	AST* u6_ = algebraicExpand(u6);

	AST* result_u6 = add({
		mul({integer(-8), symbol("z")}),
		mul({integer(8), power(symbol("z"), integer(2))}),
		mul({integer(12), power(symbol("z"), integer(3))}),
		mul({integer(-6), power(symbol("z"), integer(4))}),
		mul({integer(-12), power(symbol("z"), integer(5))}),
		mul({integer(-9), power(symbol("z"), integer(6))}),
		mul({integer(-10), power(symbol("z"), integer(7))}),
		mul({integer(-2), power(symbol("z"), integer(8))}),
		mul({integer(-2), power(symbol("z"), integer(9))}),
	});

	assert(u6_->match(result_u6));

	delete u0;
	delete u1;
	delete u2;
	delete u3;
	delete u4;
	delete u5;
	delete u6;
	delete u0_;
	delete u1_;
	delete u2_;
	delete u3_;
	delete u4_;
	delete u5_;
	delete u6_;
	delete result_u0;
	delete result_u1;
	delete result_u2;
	delete result_u3;
	delete result_u4;
	delete result_u5;
	delete result_u6;
}

void should_expand_main_operator() {
	AST* u0 = mul({
		symbol("x"),
		add({
			integer(2),
			power(
				add({
					integer(1),
					symbol("x")
				}),
				integer(2)
			)
		})
	});
	AST* r0 = algebraicExpandRoot(u0);
	AST* k0 = add({
		mul({integer(2), symbol("x")}),
		mul({
			symbol("x"),
			power(
				add({integer(1), symbol("x")}),
				integer(2)
			)
		}),
	});
	assert(r0->match(k0));

	AST* u1 = power(
		add({
			symbol("x"),
			power(
				add({
					integer(1),
					symbol("x")
				}),
				integer(2)
			)
		}),
		integer(2)
	);

	AST* r1 = algebraicExpandRoot(u1);

	AST* k1 = add({
		power(symbol("x"), integer(2)),
		mul({
			integer(2),
			symbol("x"),
			power(
				add({
					integer(1),
					symbol("x")
				}),
				integer(2)
			)
		}),
		power(
			add({
				integer(1),
				symbol("x")
			}),
			integer(4)
		)
	});

	assert(r1->match(k1));

	delete u0;
	delete u1;
	delete r0;
	delete r1;
	delete k0;
	delete k1;
}

void should_get_polynomial_content()
{
	AST* u = add({
		mul({integer(4), power(symbol("x"), integer(2))}),
		mul({integer(-1), integer(6), symbol("x")}),
	});

	AST* x = symbol("x");

	AST* L = list({symbol("x")});

	AST* Z = symbol("Z");
	AST* Q = symbol("Q");

	AST* u_cont = cont(u, L, Z);

	assert(u_cont->kind() == Kind::Integer);
	assert(u_cont->value() == 2);

	AST* t = mul({integer(2), symbol("x")});
	AST* t_cont = cont(t, L, Z);

	assert(t_cont->kind() == Kind::Integer);
	assert(t_cont->value() == 2);

	AST* p = mul({integer(-1), symbol("x")});
	AST* p_cont = cont(p, L, Z);

	assert(p_cont->kind() == Kind::Integer);
	assert(p_cont->value() == 1);

	AST* T = list({ symbol("x"), symbol("y") });

	AST* a = add({
		mul({fraction(1, 2), symbol("x"), symbol("y")}),
		mul({integer(6), symbol("y")}),
	});

	AST* a_cont = cont(a, T, Q);

	assert(a_cont->kind() == Kind::Symbol);
	assert(a_cont->identifier() == "y");

	AST* b = add({
		mul({
			add({
				power(symbol("y"), integer(2)),
				mul({integer(2), symbol("y")}),
				integer(1)
			}),
			power(symbol("x"), integer(2))
		}),
		mul({
			sub({
				mul({integer(2), power(symbol("y"), integer(2))}),
				integer(2),
			}),
			symbol("x")
		}),
		add({
			mul({integer(3), symbol("y")}),
			integer(3)
		}),
	});

	AST* b_cont = cont(b, T, Q);

	assert(b_cont->kind() == Kind::Addition);
	assert(b_cont->operand(0)->kind() == Kind::Integer);
	assert(b_cont->operand(0)->value() == 1);
	assert(b_cont->operand(1)->kind() == Kind::Symbol);
	assert(b_cont->operand(1)->identifier() == "y");

	delete T;
	delete L;
	delete Z;
	delete Q;
	delete x;
	delete u;
	delete t;
	delete p;
	delete a;
	delete b;
	delete a_cont;
	delete b_cont;
	delete u_cont;
	delete t_cont;
	delete p_cont;
}



void should_get_polynomial_content_sub_resultant()
{
	AST* u = add({
		mul({ integer(4), power(symbol("x"), integer(2)) }),
		mul({ integer(-1), integer(6), symbol("x") }),
	});

	AST* x = symbol("x");

	AST* L = list({symbol("x")});

	AST* Z = symbol("Z");
	AST* Q = symbol("Q");

	AST* u_cont = cont(u, L, Z);

	assert(u_cont->kind() == Kind::Integer);
	assert(u_cont->value() == 2);

	AST* t = mul({integer(2), symbol("x")});
	AST* t_cont = cont(t, L, Z);

	assert(t_cont->kind() == Kind::Integer);
	assert(t_cont->value() == 2);

	AST* p = mul({integer(-1), symbol("x")});
	AST* p_cont = cont(p, L, Z);

	assert(p_cont->kind() == Kind::Integer);
	assert(p_cont->value() == 1);

	AST* T = list({ symbol("x"), symbol("y") });

	AST* a = add({
		mul({fraction(1, 2), symbol("x"), symbol("y")}),
		mul({integer(6), symbol("y")}),
	});

	AST* a_cont = cont(a, T, Q);

	assert(a_cont->kind() == Kind::Symbol);
	assert(a_cont->identifier() == "y");

	AST* b = add({
		mul({
			add({
				power(symbol("y"), integer(2)),
				mul({integer(2), symbol("y")}),
				integer(1)
			}),
			power(symbol("x"), integer(2))
		}),
		mul({
			sub({
				mul({integer(2), power(symbol("y"), integer(2))}),
				integer(2),
			}),
			symbol("x")
		}),
		add({
			mul({integer(3), symbol("y")}),
			integer(3)
		}),
	});

	AST* b_cont = cont(b, T, Q);

	assert(b_cont->kind() == Kind::Addition);
	assert(b_cont->operand(0)->kind() == Kind::Integer);
	assert(b_cont->operand(0)->value() == 1);
	assert(b_cont->operand(1)->kind() == Kind::Symbol);
	assert(b_cont->operand(1)->identifier() == "y");

	delete T;
	delete L;
	delete Z;
	delete Q;
	delete x;
	delete u;
	delete t;
	delete p;
	delete a;
	delete b;
	delete a_cont;
	delete b_cont;
	delete u_cont;
	delete t_cont;
	delete p_cont;
}

void should_monomial_base_expand_polynomials()
{
	AST* u = add({
		mul({
			power(symbol("a"), integer(2)),
			symbol("b")
		}),
		mul({
			integer(2),
			symbol("a"),
			power(symbol("b"), integer(2)),
		}),
		power(symbol("b"), integer(3)),
		mul({integer(2), symbol("a")}),
		mul({integer(2), symbol("b")}),
		integer(3),
	});

	AST* v = add({symbol("a"), symbol("b")});

	AST* L = list({symbol("a"), symbol("b")});

	AST* t = symbol("t");

	AST* r = monomialBasedPolyExpansion(u, v, L, t);

	assert(r->kind() == Kind::Addition);
	assert(r->operand(0)->kind() == Kind::Integer);
	assert(r->operand(0)->value() == 3);
	assert(r->operand(1)->kind() == Kind::Multiplication);
	assert(r->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(r->operand(1)->operand(0)->value() == 2);
	assert(r->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(r->operand(1)->operand(1)->identifier() == "t");
	assert(r->operand(2)->kind() == Kind::Multiplication);
	assert(r->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(r->operand(2)->operand(0)->identifier() == "b");
	assert(r->operand(2)->operand(1)->kind() == Kind::Power);
	assert(r->operand(2)->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(r->operand(2)->operand(1)->operand(0)->identifier() == "t");
	assert(r->operand(2)->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(r->operand(2)->operand(1)->operand(1)->value() == 2);

	delete u;
	delete v;
	delete L;
	delete t;
	delete r;
}


int main() {

	should_get_polynomial_variable();
	should_get_if_is_polynomial_gpe();
	should_algebraic_expand_expressions();
	should_expand_main_operator();
	should_get_degree_of_variables();
	should_get_coefficient_gpe();
	should_get_leading_coefficient_gpe();
	should_divided_polynomials();
	should_get_gcd_polynomials();
	should_get_extended_gcd_polynomials();
	should_calculate_monomial_division();
	should_get_leading_monomial();
	should_rec_divide_polynomials();
	should_pseudo_divide_polynomials();
	should_normalize_polynomial();
	should_mv_poly_gcd();
	should_get_coeff_var_parts_of_monomial();
	should_collect_terms();
	should_get_polynomial_content();
	should_get_polynomial_content_sub_resultant();
	should_monomial_base_expand_polynomials();

	return 0;
}
