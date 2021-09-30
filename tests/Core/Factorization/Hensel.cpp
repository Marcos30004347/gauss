#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Factorization/Hensel.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

// void should_replace_leading_coefficients()
// {
// 	AST* ax = add({
// 		mul({integer(3), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(3))}),
// 		mul({integer(1), power(symbol("x"), integer(2))}),
// 		mul({integer(1), power(symbol("x"), integer(1))}),
// 		integer(12)
// 	});

// 	AST* x = symbol("x");
// 	AST* lc = integer(15);

// 	AST* kx = leadCoeffReplace(ax, x, lc);

// 	AST* rx = add({
// 		mul({integer(15), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(3))}),
// 		power(symbol("x"), integer(2)),
// 	  symbol("x"),
// 		integer(12)
// 	});

// 	assert(kx->match(rx));

// 	AST* ux = add({
// 		mul({integer(3), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(3))}),
// 		mul({integer(1), power(symbol("x"), integer(2))}),
// 		mul({integer(1), power(symbol("x"), integer(1))}),
// 		integer(10)
// 	});

// 	AST* px = leadCoeffReplace(ux, x, lc);

// 	AST* zx = add({
// 		mul({integer(15), power(symbol("x"), integer(4))}),
// 		mul({integer(2), power(symbol("x"), integer(3))}),
// 		power(symbol("x"), integer(2)),
// 	  symbol("x"),
// 		integer(10)
// 	});

// 	assert(px->match(zx));

// 	delete ax;
// 	delete x;
// 	delete lc;
// 	delete kx;
// 	delete rx;
// 	delete ux;
// 	delete px;
// 	delete zx;
// }

void should_hensel_lift_polynomials()
{
	AST* ax = add({
		power(symbol("x"), integer(3)),
		mul({integer(10), power(symbol("x"), integer(2))}),
		mul({integer(-432), symbol("x")}),
		integer(5040)
	});

	AST* x = symbol("x");

	AST* ux_1 = symbol("x");
	AST* wx_1 = add({power(symbol("x"), integer(2)), integer(-2)});
	AST* p = integer(5);
	AST* B = integer(5040);
	AST* gamma = undefined();

	AST* factors = univariateHensel(ax, x, p, ux_1, wx_1, B, gamma);

	AST* ax0 = add({
		symbol("x"),
		integer(30)
	});

	AST* ax1 = add({
		power(symbol("x"), integer(2)),
		mul({ integer(-20), symbol("x") }),
		integer(168)
	});

	assert(factors->kind() == Kind::List);	
	assert(factors->numberOfOperands() == 2);	
	assert(factors->operand(0)->match(ax0));
	assert(factors->operand(1)->match(ax1));


	AST* bx = add({
		power(symbol("x"), integer(4)),
		integer(1)
	});

	AST* kx_1 = add({
		power(symbol("x"), integer(2)),
		integer(2)
	});

	AST* qx_1 = add({
		power(symbol("x"), integer(2)),
		integer(-2)
	});

	AST* m = integer(5);

	AST* factors1 = univariateHensel(bx, x, m, kx_1, qx_1, B, gamma);

	assert(factors1->kind() == Kind::Fail);

	AST* cx = add({
		mul({integer(12), power(symbol("x"), integer(3))}),
		mul({integer(10), power(symbol("x"), integer(2))}),
		mul({integer(-36), symbol("x")}),
		integer(35)
	});

	AST* fx_1 = add({
		power(symbol("x"), integer(2)),
		integer(2)
	});

	AST* tx_1 = mul({
		integer(2),
		symbol("x")
	});
	
	AST* r = integer(5);
	AST* C = integer(36);
	
	AST* factors2 = univariateHensel(cx, x, r, tx_1, fx_1, C, gamma);

	AST* cx0 = add({
		integer(5),
		mul({ integer(2), symbol("x") }),
	});

	AST* cx1 = add({
		mul({ integer(6), power(symbol("x"), integer(2)) }),
		mul({ integer(-10), symbol("x") }),
		integer(7)
	});

	assert(factors2->kind() == Kind::List);	
	assert(factors2->numberOfOperands() == 2);	
	assert(factors2->operand(0)->match(cx0));
	assert(factors2->operand(1)->match(cx1));

	delete fx_1;
	delete tx_1;
	delete r;
	delete cx;
	delete cx0;
	delete cx1;
	delete m;
	delete bx;
	delete kx_1;
	delete qx_1;
	delete ax;
	delete ax0;
	delete ax1;
	delete x;
	delete ux_1;
	delete wx_1;
	delete p;
	delete B;
	delete C;
	delete gamma;
	delete factors;
	delete factors1;
	delete factors2;
}


int main()
{
	// should_replace_leading_coefficients();

	should_hensel_lift_polynomials();

	return 0;
}
