#include <assert.h>

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Factorization.hpp"
#include "Core/Factorization/Wang.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
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

	long d[4];

	int success = nondivisors(4, F, 1, L, K, d);

	assert(success == 1);
	
	assert(d[0] == 7);
	assert(d[1] == 3);
	assert(d[2] == 11);
	assert(d[3] == 17);

	delete F;
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

int main()
{
	should_get_nondivisors();
	shoud_get_ground_lead_coeff();
	should_factor_poly();
	return 0;
}
