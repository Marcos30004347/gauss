#include <assert.h>

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Factorization.hpp"
#include "Core/Factorization/Berlekamp.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_get_berlekamp_factors() 
{
	AST* u = add({
		power(symbol("x"), integer(6)),
		power(symbol("x"), integer(5)),
		symbol("x"),
		integer(4)
	});

	AST* x = symbol("x");
	AST* p = integer(5);

	AST* factors = berlekampFactors(u, x, 5);

	// AST* factors = berlekamp(u, x, p);

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
	delete p;
	delete x;
	delete F;
	delete factors;
}

void should_factorize_with_berlekamp()
{
	AST* ax = add({
		power(
			symbol("x"),
			integer(6)
		),
		mul({
			integer(-3),
			power(
				symbol("x"),
				integer(5)
			)
		}),
		power(
			symbol("x"),
			integer(4)
		),
		mul({
			integer(-3),
			power(
				symbol("x"),
				integer(3)
			)
		}),
		mul({
			integer(-1),
			power(
				symbol("x"),
				integer(2)
			)
		}),
		mul({
			integer(-3),
			symbol("x"),
		}),
		integer(1)
	});

	AST* x = symbol("x");
	AST* q = integer(11);

	AST* f = berlekampFactors(ax, x, 11);

	AST* F = set({
		add({
			symbol("x"),
			integer(1),
		}),
		add({
			power(symbol("x"), integer(2)),
			mul({
				integer(5),
				symbol("x")
			}),
			integer(3)
		}),
		add({
			power(symbol("x"), integer(3)),
			mul({
				integer(2),
				power(symbol("x"), integer(2))
			}),
			mul({
				integer(3),
				symbol("x")
			}),
			integer(4)
		}),
	});

	assert(f->match(F));

	delete ax;
	delete f;
	delete x;
	delete q;
	delete F;
}

int main()
{
	should_get_berlekamp_factors();	
	should_factorize_with_berlekamp();
}
