#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Factorization/Hensel.hpp"

#include "test.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

// void should_replace_leading_coefficients()
// {
// 	Expr ax = add({
// 		mul({integer(3), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(3))}),
// 		mul({integer(1), power(symbol("x"), integer(2))}),
// 		mul({integer(1), power(symbol("x"), integer(1))}),
// 		integer(12)
// 	});

// 	Expr x = symbol("x");
// 	Expr lc = integer(15);

// 	Expr kx = leadCoeffReplace(ax, x, lc);

// 	Expr rx = add({
// 		mul({integer(15), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(3))}),
// 		power(symbol("x"), integer(2)),
// 	  symbol("x"),
// 		integer(12)
// 	});

// 	assert(kx->match(rx));

// 	Expr ux = add({
// 		mul({integer(3), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(3))}),
// 		mul({integer(1), power(symbol("x"), integer(2))}),
// 		mul({integer(1), power(symbol("x"), integer(1))}),
// 		integer(10)
// 	});

// 	Expr px = leadCoeffReplace(ux, x, lc);

// 	Expr zx = add({
// 		mul({integer(15), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(3))}),
// 		power(symbol("x"), integer(2)),
// 	  symbol("x"),
// 		integer(10)
// 	});

// 	assert(px->match(zx));

//
//
//
//
//
//
//
//
// }

// void should_hensel_lift_polynomials()
// {
// 	Expr ax = add({
// 		power(symbol("x"), integer(3)),
// 		mul({integer(10), power(symbol("x"), integer(2))}),
// 		mul({integer(-432), symbol("x")}),
// 		integer(5040)
// 	});

// 	Expr x = symbol("x");

// 	Expr ux_1 = symbol("x");
// 	Expr wx_1 = add({power(symbol("x"), integer(2)), integer(-2)});
// 	Expr p = integer(5);
// 	Expr B = integer(5040);
// 	Expr gamma = undefined();

// 	Expr factors = univariateHensel(ax, x, p, ux_1, wx_1, B, gamma);

// 	Expr ax0 = add({
// 		symbol("x"),
// 		integer(30)
// 	});

// 	Expr ax1 = add({
// 		power(symbol("x"), integer(2)),
// 		mul({ integer(-20), symbol("x") }),
// 		integer(168)
// 	});

// 	assert(factors->kind() == Kind::List);
// 	assert(factors.size() == 2);
// 	assert(factors[0]->match(ax0));
// 	assert(factors[1]->match(ax1));


// 	Expr bx = add({
// 		power(symbol("x"), integer(4)),
// 		integer(1)
// 	});

// 	Expr kx_1 = add({
// 		power(symbol("x"), integer(2)),
// 		integer(2)
// 	});

// 	Expr qx_1 = add({
// 		power(symbol("x"), integer(2)),
// 		integer(-2)
// 	});

// 	Expr m = integer(5);

// 	Expr factors1 = univariateHensel(bx, x, m, kx_1, qx_1, B, gamma);

// 	assert(factors1->kind() == Kind::Fail);

// 	Expr cx = add({
// 		mul({integer(12), power(symbol("x"), integer(3))}),
// 		mul({integer(10), power(symbol("x"), integer(2))}),
// 		mul({integer(-36), symbol("x")}),
// 		integer(35)
// 	});

// 	Expr fx_1 = add({
// 		power(symbol("x"), integer(2)),
// 		integer(2)
// 	});

// 	Expr tx_1 = mul({
// 		integer(2),
// 		symbol("x")
// 	});

// 	Expr r = integer(5);
// 	Expr C = integer(36);

// 	Expr factors2 = univariateHensel(cx, x, r, tx_1, fx_1, C, gamma);

// 	Expr cx0 = add({
// 		integer(5),
// 		mul({ integer(2), symbol("x") }),
// 	});

// 	Expr cx1 = add({
// 		mul({ integer(6), power(symbol("x"), integer(2)) }),
// 		mul({ integer(-10), symbol("x") }),
// 		integer(7)
// 	});

// 	assert(factors2->kind() == Kind::List);
// 	assert(factors2.size() == 2);
// 	assert(factors2[0]->match(cx0));
// 	assert(factors2[1]->match(cx1));

//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
// }

void should_perform_hensel_step()
{
	Expr f = add({
		power(symbol("x"), integer(4)),
		integer(-1),
	});

	Expr p = integer(5);

	Expr g = add({
		power(symbol("x"), integer(3)),
		mul({integer(2), power(symbol("x"), integer(2))}),
		mul({integer(-1), symbol("x")}),
		integer(-2)
	});

	Expr h = add({
		symbol("x"),
		integer(-2)
	});

	Expr s = integer(-2);

	Expr t = add({
		mul({integer(2), power(symbol("x"), integer(2))}),
		mul({integer(-2), symbol("x")}),
		integer(-1)
	});

	Expr x = symbol("x");

	Expr L = henselSep(f, g, h, s, t, x, 5);

	Expr G = add({
		power(symbol("x"), integer(3)),
		mul({
			integer(7),
			power(symbol("x"), integer(2)),
		}),
		mul({integer(-1), symbol("x")}),
		integer(-7)
	});

	Expr H = add({
		symbol("x"),
		integer(-7)
	});

	Expr S = integer(8);

	Expr T = add({
		mul({
			integer(-8),
			power(symbol("x"), integer(2)),
		}),
		mul({integer(-12), symbol("x")}),
		integer(-1)
	});

	assert(L.kind() == Kind::List);
	assert(L.size() == 4);

	assert(L[0] == G);
	assert(L[1] == H);
	assert(L[2] == S);
	assert(L[3] == T);
}

void should_multifactor_hensel_lift()
{
	Expr f = add({
		power(symbol("x"), integer(4)),
		integer(-1)
	});

	Expr H = list({
		add({
			symbol("x"),
			integer(-1)
		}),
		add({
			symbol("x"),
			integer(-2)
		}),
		add({
			symbol("x"),
			integer(2)
		}),
		add({
			symbol("x"),
			integer(1)
		}),
	});

	Expr x = symbol("x");

	Expr F = multifactorHenselLifting(f, H, x, 5, 4);

	printf("%s\n", F.toString().c_str());

}

int main()
{
	// should_replace_leading_coefficients();
	// should_hensel_lift_polynomials();
	TEST(should_perform_hensel_step)
	TEST(should_multifactor_hensel_lift)

	return 0;
}
