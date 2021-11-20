#include <assert.h>

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Factorization/Wang.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;
using namespace factorization;

void should_get_nondivisors()
{
	AST* F = list({
		integer(-14),
		integer(3),
		integer(-11),
		integer(-17),
	});

	AST* L = list({
		symbol("y"),
		symbol("z"),
	});

	AST* K = symbol("Z");

	AST* d = nondivisors(4, F, 1, L, K);
	// assert(success == 1);

	assert(d->operand(0)->value() == 7);
	assert(d->operand(1)->value() == 3);
	assert(d->operand(2)->value() == 11);
	assert(d->operand(3)->value() == 17);

	delete F;
	delete L;
	delete K;
	delete d;
}

void shoud_get_ground_lead_coeff()
{
	AST* t = add({
		mul({
			add({
				mul({integer(2), power(symbol("y"), integer(2))}),
				mul({integer(3), power(symbol("y"), integer(1))}),
				integer(4)
			}),
			power(symbol("x"), integer(2))
		}),
		integer(5)
	});
	
	AST* L = list({ symbol("x"), symbol("y") });
	
	AST* lc = groundLeadCoeff(t, L);

	assert(lc->is(2));

	delete t;
	delete L;
	delete lc;
}

void should_factor_poly()
{
	AST* t = add({
		mul({
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		integer(-9)
	});
	
	AST* L = list({symbol("x"), symbol("y"), symbol("z")});

	AST* K = symbol("Z");

	AST* F0 = factors(t, L, K);

	delete t;
	delete L;
	delete K;
	delete F0;
}


void should_solve_diophant()
{
	AST* H1 = list({
		add({
			mul({integer(44), power(symbol("x"), integer(2))}),
			mul({integer(42), symbol("x")}),
			integer(1)
		}),
		add({
			mul({integer(126), power(symbol("x"), integer(2))}),
			mul({integer(-9), symbol("x")}),
			integer(28)
		}),
		add({
			mul({integer(187), power(symbol("x"), integer(2))}),
			integer(-23)
		}),
	});


	AST* H2 = list({
		add({
			mul({integer(-4), power(symbol("x"), integer(2)), symbol("y")}),
			mul({integer(-12), power(symbol("x"), integer(2))}),
			mul({integer(-3), symbol("x"), symbol("y")}),
			integer(1)
		}),
		add({
			mul({integer(-9), power(symbol("x"), integer(2)), symbol("y")}),
			mul({integer(-9), symbol("x")}),
			mul({integer(-2), symbol("y")}),
		}),
		add({
			mul({power(symbol("x"), integer(2)), power(symbol("y"), integer(2))}),
			mul({integer(-9), power(symbol("x"), integer(2))}),
			symbol("y"),
			integer(-9)
		}),
	}); 


	AST* H3 = list({
		add({
			mul({integer(-4), power(symbol("x"), integer(2)), symbol("y")}),
			mul({integer(-12), power(symbol("x"), integer(2))}),
			mul({integer(-3), symbol("x"), symbol("y")}),
			integer(1)
		}),
		add({
			mul({integer(-9), power(symbol("x"), integer(2)), symbol("y")}),
			mul({integer(-9), symbol("x")}),
			mul({integer(-2), symbol("y")}),
		}),
		add({
			mul({power(symbol("x"), integer(2)), power(symbol("y"), integer(2))}),
			mul({integer(-9), power(symbol("x"), integer(2))}),
			symbol("y"),
			integer(-9)
		}),
	}); 

	AST* c1 = add({
		mul({integer(-70686), power(symbol("x"), integer(5))}),
		mul({integer(-5863), power(symbol("x"), integer(4))}),
		mul({integer(-17826), power(symbol("x"), integer(3))}),
		mul({integer(2009), power(symbol("x"), integer(2))}),
		mul({integer(5031), symbol("x")}),
		integer(74)
	});

	AST* c2 = add({
		mul({integer(9), power(symbol("x"), integer(5)), power(symbol("y"), integer(4))}),
		mul({integer(12), power(symbol("x"), integer(5)), power(symbol("y"), integer(3))}),
		mul({integer(-45), power(symbol("x"), integer(5)), power(symbol("y"), integer(2))}),
		mul({integer(-108), power(symbol("x"), integer(5)), symbol("y")}),
		mul({integer(-324), power(symbol("x"), integer(5))}),
		mul({integer(18), power(symbol("x"), integer(4)), power(symbol("y"), integer(3))}),
		mul({integer(-216), power(symbol("x"), integer(4)), power(symbol("y"), integer(2))}),
		mul({integer(-810), power(symbol("x"), integer(4)), symbol("y")}),
		mul({integer(2), power(symbol("x"), integer(3)), power(symbol("y"), integer(4))}),
		mul({integer(9), power(symbol("x"), integer(3)), power(symbol("y"), integer(3))}),
		mul({integer(-252), power(symbol("x"), integer(3)), power(symbol("y"), integer(2))}),
		mul({integer(-288), power(symbol("x"), integer(3)), symbol("y")}),
		mul({integer(-945), power(symbol("x"), integer(3))}),
		mul({integer(-30), power(symbol("x"), integer(2)), power(symbol("y"), integer(2)) }),
		mul({integer(-414), power(symbol("x"), integer(2)), symbol("y") }),
		mul({integer(2), symbol("x"), power(symbol("y"), integer(3)) }),
		mul({integer(-54), symbol("x"), power(symbol("y"), integer(2)) }),
		mul({integer(-3), symbol("x"), symbol("y") }),
		mul({integer(81), symbol("x") }),
		mul({integer(12), symbol("y") }),
	});

	AST* c3 = add({
		mul({integer(-36), power(symbol("x"), integer(4)), power(symbol("y"), integer(2))}),
		mul({integer(-108), power(symbol("x"), integer(4)), symbol("y")}),
		mul({integer(-27), power(symbol("x"), integer(3)),  power(symbol("y"), integer(2))}),
		mul({integer(-36), power(symbol("x"), integer(3)),  symbol("y")}),
		mul({integer(-108), power(symbol("x"), integer(3))}),
		mul({integer(-8), power(symbol("x"), integer(2)), power(symbol("y"), integer(2))}),
		mul({integer(-42), power(symbol("x"), integer(2)), symbol("y")}),
		mul({integer(-6), symbol("x"), power(symbol("y"), integer(2))}),
		mul({integer(9), symbol("x")}),
		mul({integer(2), symbol("y")}),
	});

	AST* L1 = list({ symbol("x") });
	AST* I1 = list({});
 	AST* D1 = multivariateDiophant(H1, c1, L1, I1, 5, 6291469, 1);

	AST* R1 = list({
		mul({integer(-3), symbol("x")}),
		integer(-2),
		integer(1),
	});

	assert(D1->match(R1));

	AST* L2 = list({ symbol("x"), symbol("y") });
	AST* I2 = list({integer(-14)});
 	AST* D2 = multivariateDiophant(H2, c2, L2, I2, 5, 6291469, 1);
	
	AST* R2 = list({
		mul({integer(-1), symbol("x"), symbol("y")}),
		mul({integer(-3), symbol("x")}),
		integer(-6),
	});

	assert(D2->match(R2));

 	AST* D3 = multivariateDiophant(H3, c3, L2, I2, 5, 6291469, 1);
	
	AST* R3 = list({
		integer(0),
		integer(0),
		integer(-1),
	});

	assert(D3->match(R3));

	delete H1;
	delete H2;
	delete H3;
	delete c1;
	delete c2;
	delete c3;
	delete L1;
	delete L2;
	delete I1;
	delete I2;
	delete R1;
	delete R2;
	delete R3;
	delete D1;
	delete D2;
	delete D3;
}

void should_get_lead_coeffs()
{
	AST* f = add({
		mul({
			integer(4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(4)), 
			power(symbol("z"), integer(2))
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(3)), 
			power(symbol("z"), integer(3))
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(4))
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(6)), 
			symbol("y"),
			power(symbol("z"), integer(5))
		}),
		mul({
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(4)), 
			power(symbol("z"), integer(3))
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(3)),
			symbol("z") 
		}),
		mul({
			integer(-1),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(5)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(5)), 
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(4)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(3)), 
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(3)), 
			symbol("z")
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(4)), 
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(3)), 
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(2)), 
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(5)), 
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(4)), 
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(3)), 
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(4)),
			symbol("z") 
		}),
		mul({
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(3)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-1),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(9),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			symbol("z")
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(3)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(3)), 
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(3)), 
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(3),
			power(symbol("x"), integer(3)), 
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(3)),
		}),
		mul({
			integer(-6),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(2)),
			symbol("z")
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("x"),
			power(symbol("y"), integer(3)),
			symbol("z"),
		}),
		mul({
			integer(-2),
			symbol("x"),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-3),
			symbol("x"),
			symbol("y"),
			symbol("z"),
		}),
		mul({
			integer(3),
			symbol("x"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-2),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
	});

	AST* K = symbol("Z");

	AST* L = list({symbol("x"), symbol("y"), symbol("z")});
	AST* a = list({integer(-14), integer(3)});
	AST* p = deepReplace(f, L->operand(1), a->operand(0));
	AST* t = deepReplace(p, L->operand(2), a->operand(1));
	AST* k = reduceAST(t);

	AST* d = cont(k, L->operand(0), K);
	AST* s = pp(k, d, L->operand(0), K);

	AST* F = list({
		list({symbol("y"), integer(1)}),
		list({symbol("z"), integer(2) }),
		list({add({symbol("y"), symbol("z")}), integer(2) }),
		list({add({symbol("y"), mul({integer(-1), symbol("z")})}), integer(1) }),
	});

	AST* sF = list({
		integer(-14),
		integer(3),
		integer(-11),
		integer(-17),
	});

	AST* sqf = sqfFactors(s, L->operand(0), K);

	AST* wlc = wangLeadingCoeff(f, d, sqf->operand(1), F, sF, a, L, K);

	assert(wlc->operand(0)->match(f));

	AST* S = set({
		add({
			mul({integer(187), power(symbol("x"), integer(2))}),
			integer(-23)
		}),
		add({
			mul({integer(44), power(symbol("x"), integer(2))}),
			mul({integer(42), symbol("x")}),
			integer(1)
		}),
		add({
			mul({integer(126), power(symbol("x"), integer(2))}),
			mul({integer(-9), symbol("x")}),
			integer(28)
		}),
	});

	AST* Q = set({
		wlc->operand(1)->operand(0)->copy(),
		wlc->operand(1)->operand(1)->copy(),
		wlc->operand(1)->operand(2)->copy(),
	});

	assert(S->match(Q));

	AST* q = set({
		add({
			mul({integer(-4), symbol("y")}),
			mul({integer(-4), symbol("z")}),
		}),
		mul({
			add({symbol("y"), symbol("z")}),
			add({symbol("y"), mul({integer(-1), symbol("z")})}),
		}),
		mul({
			integer(-1),
			symbol("y"),
			power(symbol("z"), integer(2)),
		})
	});

	AST* E = set({
		wlc->operand(2)->operand(0)->copy(),
		wlc->operand(2)->operand(1)->copy(),
		wlc->operand(2)->operand(2)->copy(),
	});

	//	printf("%s\n", wlc->operand(2)->toString().c_str());
	//	printf("%s\n", q->toString().c_str());

	assert(E->match(q));

	delete wlc;
	delete sqf;
	delete f;
	delete S;
	delete Q;
	delete F;
	delete E;
	delete sF;
	delete K;
	delete L;
	delete a;
	delete p;
	delete t;
	delete k;
	delete d;
	delete s;
	delete q;
}

void should_factorize_multivariate_polynomials()
{
	AST* L = list({ symbol("x") });
	AST* K = symbol("Z");

	// AST* u0 = integer(0);
	// AST* U0 = factors(u0, L, K);

	// printf("U0 = %s\n", U0->toString().c_str());

	// AST* u1 = integer(3);
	// AST* U1 = factors(u1, L, K);

	// printf("U1 = %s\n", U1->toString().c_str());
	
	// AST* u2 = integer(-8);
	// AST* U2 = factors(u2, L, K);
	
	// printf("U2 = %s\n", U2->toString().c_str());

	// AST* u3 = add({
	// 	power(symbol("x"), integer(2)),
	// 	integer(-9)
	// });
	// AST* U3 = factors(u3, L, K);
	// printf("U3 = %s\n", U3->toString().c_str());

	AST* u4 = add({
		mul({
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(2)),
			symbol("y")
		}),
		mul({
			integer(9),
			power(symbol("x"), integer(2)),
		}),
		integer(-1)
	});
	// AST* u3 = integer(-8);

	AST* R = list({ symbol("x"), symbol("y") });
	// AST* U4 = factors(u4, R, K);

	// printf("U4 = %s\n", U4->toString().c_str());


	AST* u5 = add({
		mul({ 
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)) 
		}),
		integer(-9)
	});

	AST* T = list({ symbol("x"), symbol("y"), symbol("z") });
	AST* U5 = factors(u5, T, K);
	printf("U5 = %s\n", U5->toString().c_str());

}

void should_get_evaluation_points()
{
	AST* U = add({
		mul({
			integer(4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(4)), 
			power(symbol("z"), integer(2))
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(3)), 
			power(symbol("z"), integer(3))
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(4))
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(6)), 
			symbol("y"),
			power(symbol("z"), integer(5))
		}),
		mul({
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(4)), 
			power(symbol("z"), integer(3))
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(3)),
			symbol("z") 
		}),
		mul({
			integer(-1),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(5)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(5)), 
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(4)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(3)), 
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(3)), 
			symbol("z")
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(4)), 
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(3)), 
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(2)), 
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(5)), 
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(4)), 
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(3)), 
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(4)),
			symbol("z") 
		}),
		mul({
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(3)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-1),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(9),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			symbol("z")
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(3)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(3)), 
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(3)), 
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(3),
			power(symbol("x"), integer(3)), 
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(3)),
		}),
		mul({
			integer(-6),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(2)),
			symbol("z")
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("x"),
			power(symbol("y"), integer(3)),
			symbol("z"),
		}),
		mul({
			integer(-2),
			symbol("x"),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-3),
			symbol("x"),
			symbol("y"),
			symbol("z"),
		}),
		mul({
			integer(3),
			symbol("x"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-2),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
	});

	AST* V = list({
		symbol("y"),
		symbol("z"),
		add({symbol("y"), symbol("z")}),
		add({symbol("y"), mul({ integer(-1), symbol("z")}) }),
	});

	AST* G = integer(4);
	AST* L = list({ symbol("x"), symbol("y"), symbol("z") });
	AST* K = symbol("Z");

	AST* E = getEvaluationPoints(U, G, V, L, K, 3);
	
	printf("%s\n", E->toString().c_str());


	delete U;
	delete G;
	delete L;
	delete K;
	delete E;
}

void should_lift_factors()
{
	AST* f = add({
		mul({
			integer(4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(4)), 
			power(symbol("z"), integer(2))
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(3)), 
			power(symbol("z"), integer(3))
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(6)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(4))
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(6)), 
			symbol("y"),
			power(symbol("z"), integer(5))
		}),
		mul({
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(4)), 
			power(symbol("z"), integer(3))
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(3)),
			symbol("z") 
		}),
		mul({
			integer(-1),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(5)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(5)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(5)), 
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(4)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(3)), 
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(3)), 
			symbol("z")
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(4)), 
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(3)), 
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(4)), 
			power(symbol("y"), integer(2)), 
			power(symbol("z"), integer(2)), 
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(5)), 
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(4)), 
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(4)), 
			symbol("y"),
			power(symbol("z"), integer(3)), 
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(4)),
			symbol("z") 
		}),
		mul({
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(3)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-1),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(9),
			power(symbol("x"), integer(3)), 
			power(symbol("y"), integer(2)),
			symbol("z")
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(3)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(3)), 
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(3)), 
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(3),
			power(symbol("x"), integer(3)), 
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(3)),
		}),
		mul({
			integer(-6),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(2)), 
			power(symbol("y"), integer(2)),
			symbol("z")
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(2)), 
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("x"),
			power(symbol("y"), integer(3)),
			symbol("z"),
		}),
		mul({
			integer(-2),
			symbol("x"),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-3),
			symbol("x"),
			symbol("y"),
			symbol("z"),
		}),
		mul({
			integer(3),
			symbol("x"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-2),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
	});


	AST* U = list({
		add({
			mul({integer(44), power(symbol("x"), integer(2))}),
			mul({integer(42), symbol("x")}),
			integer(1)
		}),
		add({
			mul({integer(126), power(symbol("x"), integer(2))}),
			mul({integer(-9), symbol("x")}),
			integer(28)
		}),
		add({
			mul({integer(187), power(symbol("x"), integer(2))}),
			integer(-23)
		}),
	});

	AST* LC = list({
		add({
			mul({integer(-4), symbol("y")}),
			mul({integer(-4), symbol("z")}),
		}),
		mul({integer(-1), symbol("y"), power(symbol("z"), integer(2))}),
		add({
			power(symbol("y"), integer(2)),
			mul({integer(-1), power(symbol("z"), integer(2))}),
		}),
	});

	AST* a = list({
		integer(-14),
		integer(3),
	});
	
	AST* L = list({symbol("x"), symbol("y"), symbol("z")});
	
	AST* K = symbol("Z");
	
	AST* r = wangEEZ(f, U, LC, a, 6291469, L, K);
	printf("%s\n", r->toString().c_str());
}

int main()
{

	// should_get_nondivisors();
	// should_solve_diophant();
	
	// should_get_lead_coeffs();

	// should_get_evaluation_points();
	// should_factorize_multivariate_polynomials();
	should_lift_factors();
	return 0;
}
