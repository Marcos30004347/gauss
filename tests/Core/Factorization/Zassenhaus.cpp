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
	Expr x = Expr("x");

	Expr f = power(x, 8) + power(x, 7) + -1*power(x, 6) + power(x, 5) + -1*power(x, 3) + -1*power(x, 2) + -1*x;

	Expr d = cantorZassenhausDDF(f, x, 3);

	assert(d == list({
				list({x, 1}),
				list({x + power(x, 3) + power(x, 4) + -1, 2}),
				list({-1*x + power(x, 3) + 1, 3}),
			}));
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
	//should_distinct_degree_factorize();
	should_equal_degree_factorize();
	// should_factorize_cantor_zassenhaus();
	// should_factorize_zassenhaus();
	return 0;
}
