#include "Core/Polynomial/Resultant.hpp"

#include <assert.h>


using namespace ast;
using namespace algebra;
using namespace polynomial;

void should_get_univariate_resultant()
{
	AST* ux = add({
		mul({integer(2), power(symbol("x"), integer(3))}),
		mul({integer(-3), symbol("x")}),
		integer(1)
	});

	AST* vx = add({
		mul({integer(3), power(symbol("x"), integer(2))}),
		mul({integer(-4), symbol("x")}),
		integer(3)
	});

	AST* x = symbol("x");

	AST* r0 = univariateResultant(ux, vx, x);

	assert(r0->kind() == Kind::Integer);
	assert(r0->value() == 218);

	delete ux;
	delete vx;
	delete r0;
	delete x;
}

void should_get_multivariate_resultants()
{
	AST* u = add({
		mul({
			power(symbol("x"), integer(3)),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(4)),
			symbol("y"),
		}),
		mul({
			integer(9),
			power(symbol("x"), integer(5)),
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(24),
			power(symbol("x"), integer(3)),
			symbol("y"),
		}),
		mul({
			integer(36),
			power(symbol("x"), integer(3)),
		}),
		mul({
			integer(5),
			symbol("x"),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(45),
			power(symbol("x"), integer(3)),
		}),
		mul({
			integer(2),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(12),
			symbol("y"),
			symbol("x"),
		}),
		mul({
			integer(18),
			power(symbol("x"), integer(2)),
		}),
	});

	AST* v = add({
		mul({
			power(symbol("x"), integer(5)),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)),
			symbol("y"),
		}),
		mul({
			integer(16),
			power(symbol("x"), integer(3)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(4)),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(96),
			power(symbol("x"), integer(3)),
			symbol("y"),
		}),
		mul({
			integer(192),
			power(symbol("x"), integer(2)),
		}),
		mul({
			integer(45),
			power(symbol("x"), integer(3)),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(360),
			symbol("y"),
			power(symbol("x"), integer(2)),
		}),
		mul({
			integer(720),
			symbol("x")
		}),
		mul({
			integer(50),
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(400),
			symbol("y"),
			symbol("x"),
		}),
		integer(800)
	});

	AST* K = symbol("Z");
	AST* L = list({ symbol("x"), symbol("y") });
	
	AST* r0 = srPolynomialResultant(u, v, L, K);

	printf("%s\n", r0->toString().c_str());

	delete L;
	delete K;
	delete r0;
	delete u;
	delete v;
}

int main()
{

	should_get_univariate_resultant();
	should_get_multivariate_resultants();

	return 0;
}
