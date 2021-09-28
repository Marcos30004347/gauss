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

void should_get_multivariate_resultants0()
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

	assert(r0->kind() == Kind::Integer);
	assert(r0->value() == 0);

	delete L;
	delete K;
	delete r0;
	delete u;
	delete v;
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
			power(symbol("x"), integer(2)),
			symbol("y"),
		}),
		mul({
			integer(5),
			symbol("x"),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(2),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(1),
			symbol("y"),
			symbol("x"),
		}),
		mul({
			integer(3),
			power(symbol("x"), integer(2)),
		}),
	});

	AST* v = add({
		mul({
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(5),
			power(symbol("x"), integer(3)),
		}),
		mul({
			integer(3),
			power(symbol("x"), integer(3)),
			symbol("y"),
		}),
		mul({
			integer(4),
			symbol("y"),
			symbol("x"),
		}),
		integer(8)
	});

	AST* K = symbol("Z");
	AST* L = list({ symbol("x"), symbol("y") });
	
	AST* r0 = multivariateResultant(u, v, L, K);
	AST* r1 = srPolynomialResultant(u, v, L, K);

	AST* r = add({
		mul({integer(4), power(symbol("y"), integer(12))}),
		mul({integer(-40), power(symbol("y"), integer(11))}),
		mul({integer(156), power(symbol("y"), integer(10))}),
		mul({integer(288), power(symbol("y"), integer(9))}),
		mul({integer(-1056), power(symbol("y"), integer(8))}),
		mul({integer(1948), power(symbol("y"), integer(7))}),
		mul({integer(20440), power(symbol("y"), integer(6))}),
		mul({integer(42512), power(symbol("y"), integer(5))}),
		mul({integer(-3072), power(symbol("y"), integer(4))}),
		mul({integer(-118024), power(symbol("y"), integer(3))}),
		mul({integer(-133344), power(symbol("y"), integer(2))}),
		mul({integer(-57024), symbol("y")}),
		integer(-8640)
	});

	assert(r0->match(r));
	assert(r1->match(r));
	
	delete L;
	delete K;
	delete r0;
	delete r1;
	delete r;
	delete u;
	delete v;
}

void should_get_remainder_sequence()
{
	AST* u = add({
		power(symbol("x"), integer(8)),
		power(symbol("x"), integer(6)),
		mul({
			integer(-3),
			power(symbol("x"), integer(4)),
		}),
		mul({
			integer(-3),
			power(symbol("x"), integer(3)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("x")
		}),
		integer(-5)
	});

	AST* v = add({
		mul({
			integer(3),
			power(symbol("x"), integer(6)),
		}),
		mul({
			integer(5),
			power(symbol("x"), integer(4)),
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(2)),
		}),
		mul({
			integer(-9),
			symbol("x")
		}),
		integer(21)
	});

	AST* L = list({symbol("x")});
	AST* Q = symbol("Q");

	// PPRS(u, v, x);
	AST* r = polyRemSeq(u, v, L, Q);

	assert(r->operand(0)->kind() == Kind::Integer);
	assert(r->operand(0)->value() == 1);

	assert(r->operand(1)->kind() == Kind::Integer);
	assert(r->operand(1)->value() == 260708);

	delete r;
	delete u;
	delete v;
	delete L;
	delete Q;
}



void should_get_remainder_sequence_mv()
{
	AST* f = add({
		mul({integer(3), symbol("y"), power(symbol("x"), integer(2))}),
		mul({integer(-1), add({power(symbol("y"), integer(3)), integer(4)})}),
	});

	AST* g = add({
		power(symbol("x"), integer(2)),
		mul({power(symbol("y"), integer(3)), symbol("x")}),
		integer(-9)
	});


	AST* L = list({symbol("x"), symbol("y")});
	AST* Z = symbol("Z");

	AST* r = polyRemSeq(f, g, L, Z);

	assert(r->operand(0)->kind() == Kind::Integer);
	assert(r->operand(0)->value() == 1);

	AST* res = add({
		mul({integer(-3), power(symbol("y"), integer(10))}),
		mul({integer(-12), power(symbol("y"), integer(7))}),
		power(symbol("y"), integer(6)),
		mul({integer(-54), power(symbol("y"), integer(4))}),
		mul({integer(8), power(symbol("y"), integer(3))}),
		mul({integer(729), power(symbol("y"), integer(2))}),
		mul({integer(-216), symbol("y")}),
		integer(16)
	});

	assert(r->operand(1)->match(res));

	delete r;
	delete f;
	delete g;
	delete L;
	delete Z;
	delete res;
}


void should_get_remainder_sequence_mv1()
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

	AST* L = list({symbol("x"), symbol("y")});

	AST* Z = symbol("Z");

	AST* r = polyRemSeq(u, v, L, Z);
	
	assert(r->operand(0)->kind() == Kind::Integer);
	assert(r->operand(0)->value() == 1);
	
	assert(r->operand(1)->kind() == Kind::Integer);
	assert(r->operand(1)->value() == 0);

	delete r;
	delete u;
	delete v;
	delete L;
	delete Z;
}

int main()
{
	// should_get_univariate_resultant();
	// should_get_multivariate_resultants();
	// should_get_multivariate_resultants0();
	should_get_remainder_sequence();
	should_get_remainder_sequence_mv();
	should_get_remainder_sequence_mv1();
	return 0;
}
