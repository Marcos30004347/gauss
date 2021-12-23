#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Factorization/Hensel.hpp"

#include "test.hpp"
#include <climits>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_perform_hensel_step()
{
	Expr x = Expr("x");
	Expr f = power(x, 4) + -1;
	Expr g = power(x, 3) + 2*power(x, 2) + -1*x + -2;
	Expr h = x + -2;
	Expr s = -2;
	Expr t = 2*power(x, 2) + -2*x + -1;
	Expr L = henselSep(f, g, h, s, t, x, 5);
	Expr G = power(x, 3) + 7*power(x, 2) + -1*x + -7;
	Expr H = x + -7;
	Expr S = 8;
	Expr T = -8*power(x, 2) + -12*x + -1;

	assert(L.kind() == Kind::List);
	assert(L.size() == 4);
	assert(L[0] == G);
	assert(L[1] == H);
	assert(L[2] == S);
	assert(L[3] == T);
}

void should_perform_hensel_step_poly_expr()
{
	Expr x = Expr("x");
	Expr L = list({x});

	Expr f = polyExpr(power(x, 4) + -1, L);
	Expr g = polyExpr(power(x, 3) + 2*power(x, 2) + -1*x + -2, L);
	Expr h = polyExpr(x + -2, L);
	Expr s = polyExpr(-2, L);
	Expr t = polyExpr(2*power(x, 2) + -2*x + -1, L);

	Expr Y = henselSepPolyExpr(f, g, h, s, t, L, 5);

	Expr G = polyExpr(power(x, 3) + 7*power(x, 2) + -1*x + -7, L);
	Expr H = polyExpr(x + -7, L);
	Expr S = polyExpr(8, L);
	Expr T = polyExpr(-8*power(x, 2) + -12*x + -1, L);

	assert(Y.kind() == Kind::List);
	assert(Y.size() == 4);

	assert(Y[0] == G);
	assert(Y[1] == H);
	assert(Y[2] == S);
	assert(Y[3] == T);
}

void should_multifactor_hensel_lift()
{
	Expr x = Expr("x");

	Expr f = power(x, 4) + -1;
	Expr H = list({ x + -1, x + -2, x + 2, x + 1 });

	Expr F = multifactorHenselLifting(f, H, x, 5, 4);

	assert(F == list({ x + -1, x + -182, x + 182, x + 1 }));
}

void should_multifactor_hensel_lift_poly_expr()
{
	Expr x = Expr("x");
	Expr L = list({ x });
	Expr f = polyExpr(power(x, 4) + -1, L);
	Expr H = list({ polyExpr(x + -1, L), polyExpr(x + -2, L), polyExpr(x + 2, L), polyExpr(x + 1, L) });

	Expr F = multifactorHenselLiftingPolyExpr(f, H, L, 5, 4);

	assert(F == list({ polyExpr(x + -1, L), polyExpr(x + -182, L), polyExpr(x + 182, L), polyExpr(x + 1, L) }));
}

int main()
{
	TEST(should_perform_hensel_step)
	TEST(should_perform_hensel_step_poly_expr)
	TEST(should_multifactor_hensel_lift)
	TEST(should_multifactor_hensel_lift_poly_expr)

	return 0;
}
