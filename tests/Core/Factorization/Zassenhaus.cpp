#include <assert.h>

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Factorization.hpp"
#include "Core/Factorization/Zassenhaus.hpp"
#include "Core/Factorization/Berlekamp.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_factorize_zassenhaus()
{
	AST* f = add({
		power(symbol("x"), integer(4)),
		integer(-1)
	});
	
	AST* g = add({
		mul({integer(6), power(symbol("x"), integer(4))}),
		mul({integer(5), power(symbol("x"), integer(3))}),
		mul({integer(15), power(symbol("x"), integer(2))}),
		mul({integer(5), symbol("x")}),
		integer(4)
	});
	
	AST* x = symbol("x");

	// AST* z = zassenhaus(f, L, K);
	AST* r = zassenhaus(g, x);

	// printf("%s\n", z->toString().c_str());
	printf("%s\n", r->toString().c_str());

	delete f;
	delete g;
	delete x;
	delete r;
}

void should_distinct_degree_factorize()
{
	AST* f = add({
		power(symbol("x"), integer(8)),
		power(symbol("x"), integer(7)),
		mul({integer(-1), power(symbol("x"), integer(6))}),
		power(symbol("x"), integer(5)),
		mul({integer(-1), power(symbol("x"), integer(3))}),
		mul({integer(-1), power(symbol("x"), integer(2))}),
		mul({integer(-1), symbol("x") }),
	});

	AST* x = symbol("x");

	AST* d = cantorZassenhausDDF(f, x, 3);

	AST* F = list({
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

	printf("----> %s\n", d->toString().c_str());

	AST* g = add({
		power(symbol("x"), integer(63)),
		integer(1)
	});

	AST* k = cantorZassenhausDDF(g, x, 2);

	printf("----> %s\n", k->toString().c_str());
	
	delete F;
	delete f;
	delete g;
	delete x;
	delete k;
	delete d;
}


void should_equal_degree_factorize()
{
	AST* a = add({
		power(symbol("x"), integer(15)),
		integer(-1)
	});
	AST* x = symbol("x");

	AST* u = cantorZassenhausDDF(a, x, 11);

	printf("%s\n", u->toString().c_str());
	
	AST* a1 = add({
		power(symbol("x"), integer(5)),
		integer(-1)
	});

	AST* t = cantorZassenhausEDF(a1, x, 1, 11);
	
	printf("%s\n", t->toString().c_str());

	delete a;
	delete x;
	delete u;
	delete a1;
	delete t;
}

void should_factorize_cantor_zassenhaus()
{
	AST* a = add({
		power(symbol("x"), integer(15)),
		integer(-1)
	});

	AST* x = symbol("x");

	AST* f = cantorZassenhaus(a, x, 11);
	printf("essa = %s\n", f->toString().c_str());

	delete a;
	delete x;
	delete f;
}

int main()
{
	// should_distinct_degree_factorize();
	// should_equal_degree_factorize();
	// should_factorize_cantor_zassenhaus();
	should_factorize_zassenhaus();
	return 0;	
}
