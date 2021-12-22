#include <assert.h>

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Factorization/Zassenhaus.hpp"
#include "Core/Factorization/Berlekamp.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_factorize_zassenhaus()
{
	Expr f = add({
		power(symbol("x"), integer(4)),
		integer(-1)
	});

	Expr g = add({
		mul({integer(6), power(symbol("x"), integer(4))}),
		mul({integer(5), power(symbol("x"), integer(3))}),
		mul({integer(15), power(symbol("x"), integer(2))}),
		mul({integer(5), symbol("x")}),
		integer(4)
	});

	Expr x = symbol("x");

	// Expr z = zassenhaus(f, L, K);
	Expr r = zassenhaus(g, x);

	// printf("%s\n", z->toString().c_str());
	printf("%s\n", r.toString().c_str());
}

void should_distinct_degree_factorize()
{
	Expr f = add({
		power(symbol("x"), integer(8)),
		power(symbol("x"), integer(7)),
		mul({integer(-1), power(symbol("x"), integer(6))}),
		power(symbol("x"), integer(5)),
		mul({integer(-1), power(symbol("x"), integer(3))}),
		mul({integer(-1), power(symbol("x"), integer(2))}),
		mul({integer(-1), symbol("x") }),
	});

	Expr x = symbol("x");

	Expr d = cantorZassenhausDDF(f, x, 3);

	Expr F = list({
		list({symbol("x"), integer(1)}),
		list({
			add({
				power(symbol("x"), integer(4)),
				power(symbol("x"), integer(3)),
				symbol("x"),
				integer(-1)
			}),
			integer(2)
		}),
		list({
			add({
				power(symbol("x"), integer(3)),
				mul({integer(-1), symbol("x")}),
				integer(1)
			}),
			integer(3)
		})
	});

	// assert(d->match(F));

	printf("----> %s\n", d.toString().c_str());

	Expr g = add({
		power(symbol("x"), integer(63)),
		integer(1)
	});

	Expr k = cantorZassenhausDDF(g, x, 2);

	printf("----> %s\n", k.toString().c_str());







}


void should_equal_degree_factorize()
{
	Expr a = add({
		power(symbol("x"), integer(15)),
		integer(-1)
	});
	Expr x = symbol("x");

	Expr u = cantorZassenhausDDF(a, x, 11);

	printf("%s\n", u.toString().c_str());

	Expr a1 = add({
		power(symbol("x"), integer(5)),
		integer(-1)
	});

	Expr t = cantorZassenhausEDF(a1, x, 1, 11);

	printf("%s\n", t.toString().c_str());
}

void should_factorize_cantor_zassenhaus()
{
	Expr a = add({
		power(symbol("x"), integer(15)),
		integer(-1)
	});

	Expr x = symbol("x");

	Expr f = cantorZassenhaus(a, x, 11);
	printf("essa = %s\n", f.toString().c_str());
}

int main()
{
	// should_distinct_degree_factorize();
	// should_equal_degree_factorize();
	// should_factorize_cantor_zassenhaus();
	should_factorize_zassenhaus();
	return 0;
}
