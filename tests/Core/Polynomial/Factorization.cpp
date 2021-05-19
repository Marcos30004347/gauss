#include "Core/Algebra/Set.hpp"	
#include "Core/Algebra/List.hpp"	
#include "Core/Algebra/Algebra.hpp"	
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Factorization.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;
using namespace polynomial;

void should_get_r_matrix() {

	AST* u = add({
		power(symbol("x"), integer(6)),
		power(symbol("x"), integer(5)),
		symbol("x"),
		integer(4)
	});

	AST* x = symbol("x");

	AST* n = integer(6);

	RMatrix(u, x, n, 5);

	// 0 0 1 0 1 0
	// 0 4 3 3 1 1
	// 0 0 1 1 2 0
	// 0 0 3 4 2 0
	// 0 0 2 4 1 0
	// 0 1 3 3 1 4

	assert(getRMatrixValue(0,0) == 0);
	assert(getRMatrixValue(0,1) == 0);
	assert(getRMatrixValue(0,2) == 1);
	assert(getRMatrixValue(0,3) == 0);
	assert(getRMatrixValue(0,4) == 1);
	assert(getRMatrixValue(0,5) == 0);

	assert(getRMatrixValue(1,0) == 0);
	assert(getRMatrixValue(1,1) == 4);
	assert(getRMatrixValue(1,2) == 3);
	assert(getRMatrixValue(1,3) == 3);
	assert(getRMatrixValue(1,4) == 1);
	assert(getRMatrixValue(1,5) == 1);

	assert(getRMatrixValue(2,0) == 0);
	assert(getRMatrixValue(2,1) == 0);
	assert(getRMatrixValue(2,2) == 1);
	assert(getRMatrixValue(2,3) == 1);
	assert(getRMatrixValue(2,4) == 2);
	assert(getRMatrixValue(2,5) == 0);

	assert(getRMatrixValue(3,0) == 0);
	assert(getRMatrixValue(3,1) == 0);
	assert(getRMatrixValue(3,2) == 3);
	assert(getRMatrixValue(3,3) == 4);
	assert(getRMatrixValue(3,4) == 2);
	assert(getRMatrixValue(3,5) == 0);

	assert(getRMatrixValue(4,0) == 0);
	assert(getRMatrixValue(4,1) == 0);
	assert(getRMatrixValue(4,2) == 2);
	assert(getRMatrixValue(4,3) == 4);
	assert(getRMatrixValue(4,4) == 1);
	assert(getRMatrixValue(4,5) == 0);

	assert(getRMatrixValue(5,0) == 0);
	assert(getRMatrixValue(5,1) == 1);
	assert(getRMatrixValue(5,2) == 3);
	assert(getRMatrixValue(5,3) == 3);
	assert(getRMatrixValue(5,4) == 1);
	assert(getRMatrixValue(5,5) == 4);

	destroyRMatrix(6);

	delete u;
	delete x;
	delete n;

}

void should_get_berlekamp_factors() {
	AST* u = add({
		power(symbol("x"), integer(6)),
		power(symbol("x"), integer(5)),
		symbol("x"),
		integer(4)
	});

	AST* x = symbol("x");
	AST* factors = berlekampFactor(u, x, 5);

	AST* F = set({
		add({
			integer(3), 
			mul({
				integer(2),
				symbol("x"),
			}),
			power(
				symbol("x"),
				integer(2)
			)
		}),
		add({
			integer(4), 
			mul({
				integer(3),
				symbol("x"),
			}),
			power(
				symbol("x"),
				integer(2)
			)
		}),
		add({
			integer(2), 
			symbol("x"),
			power(
				symbol("x"),
				integer(2)
			)
		}),
	});

	assert(factors->match(F));

	delete u;
	delete x;
	delete F;
	delete factors;
}

void should_gen_extended_sigma_p() {
	AST* F0 = add({
		symbol("x"),
		integer(1)
	});

	AST* V0 = list({
		symbol("x"),
		sub({ symbol("x"), integer(1) }),
	});

	AST* x = symbol("x");

	AST* S0 = genExtendSigmaP(V0, x, 3);

	AST* k0 = list({
		integer(-1),
		integer(1),
	});

	AST* R0 = genExtendRP(V0, S0, F0, x, 3);

	AST* t0 = list({
		integer(-1),
		integer(-1),
	});

	assert(S0->match(k0));
	assert(R0->match(t0));

	AST* F1 = add({
		mul({integer(-4), power(symbol("x"), integer(4))}),
		power(symbol("x"), integer(3)),
		mul({integer(-4), power(symbol("x"), integer(2))}),
		mul({integer(-3), symbol("x")}),
		integer(4)
	});

	AST* V1 = list({
		add({
			power(symbol("x"), integer(2)),
			mul({integer(-2), symbol("x")}),
			integer(4)
		}),

		add({
			symbol("x"),
			integer(5)
		}),

		add({
			symbol("x"),
			integer(-5)
		}),

		add({
			symbol("x"),
			integer(-2)
		})
	});

	AST* S1 = genExtendSigmaP(V1, x, 11);
	AST* R1 = genExtendRP(V1, S1, F1, x, 11);
	
	AST* t1 = list({
		add({
			integer(4),
			mul({ integer(-2), symbol("x") })
		}),
		integer(0),
		integer(0),
		integer(-2),
	});

	AST* k1 = list({
		add({
			integer(2),
			mul({ integer(3), symbol("x") }),
			mul({ integer(-5), power(symbol("x"), integer(2)) }),
			mul({ integer(-4), power(symbol("x"), integer(3)) }),
			mul({ integer(-2), power(symbol("x"), integer(4)) }),
			mul({ integer(-2), power(symbol("x"), integer(5)) }),
			mul({ integer(-2), power(symbol("x"), integer(6)) }),
		}),
		add({
			integer(5),
			mul({ integer(-2), symbol("x") }),
			mul({ integer(-1), power(symbol("x"), integer(2)) }),
			mul({ integer(4), power(symbol("x"), integer(3)) }),
			mul({ integer(5), power(symbol("x"), integer(4)) }),
			mul({ integer(2), power(symbol("x"), integer(5)) }),
		}),
		add({
			integer(3),
			mul({ integer(-1), symbol("x") }),
			mul({ integer(-1), power(symbol("x"), integer(3)) }),
		}),
		integer(3),
	});

	assert(S1->match(k1));
	assert(R1->match(t1));
	
	delete x;
	delete F0;
	delete F1;
	delete V0;
	delete V1;
	delete k0;
	delete k1;
	delete t0;
	delete t1;
	delete S0;
	delete S1;
	delete R0;
	delete R1;
}

int main() {
	
	should_get_r_matrix();
	should_get_berlekamp_factors();
	should_gen_extended_sigma_p();
	
	return 0;
}
