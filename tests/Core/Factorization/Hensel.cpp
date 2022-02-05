#include "Core/Factorization/Hensel.hpp"

#include "test.hpp"
#include <climits>

using namespace alg;
using namespace polynomial;
using namespace factorization;

void should_perform_hensel_step()
{
	expr x = expr("x");
	expr f = pow(x, 4) + -1;
	expr g = pow(x, 3) + 2*pow(x, 2) + -1*x + -2;
	expr h = x + -2;
	expr s = -2;
	expr t = 2*pow(x, 2) + -2*x + -1;
	expr L = henselSep(f, g, h, s, t, x, 5);
	expr G = pow(x, 3) + 7*pow(x, 2) + -1*x + -7;
	expr H = x + -7;
	expr S = 8;
	expr T = -8*pow(x, 2) + -12*x + -1;

	assert(L.kind() == kind::LIST);
	assert(L.size() == 4);
	assert(L[0] == G);
	assert(L[1] == H);
	assert(L[2] == S);
	assert(L[3] == T);
}

void should_perform_hensel_step_poly_expr()
{
	expr x = expr("x");
	expr L = list({x});

	expr f = polyExpr(pow(x, 4) + -1, L);
	expr g = polyExpr(pow(x, 3) + 2*pow(x, 2) + -1*x + -2, L);
	expr h = polyExpr(x + -2, L);
	expr s = polyExpr(-2, L);
	expr t = polyExpr(2*pow(x, 2) + -2*x + -1, L);

	expr Y = henselSepPolyExpr(f, g, h, s, t, L, 5);

	expr G = polyExpr(pow(x, 3) + 7*pow(x, 2) + -1*x + -7, L);
	expr H = polyExpr(x + -7, L);
	expr S = polyExpr(8, L);
	expr T = polyExpr(-8*pow(x, 2) + -12*x + -1, L);

	assert(Y.kind() == kind::LIST);
	assert(Y.size() == 4);

	assert(Y[0] == G);
	assert(Y[1] == H);
	assert(Y[2] == S);
	assert(Y[3] == T);
}

void should_multifactor_hensel_lift()
{
	expr x = expr("x");

	expr f = pow(x, 4) + -1;
	expr H = list({ x + -1, x + -2, x + 2, x + 1 });

	expr F = multifactorHenselLifting(f, H, x, 5, 4);

	assert(F == list({ x + -1, x + -182, x + 182, x + 1 }));
}

void should_multifactor_hensel_lift_poly_expr()
{
	expr x = expr("x");
	expr L = list({ x });
	expr f = polyExpr(pow(x, 4) + -1, L);
	expr H = list({ polyExpr(x + -1, L), polyExpr(x + -2, L), polyExpr(x + 2, L), polyExpr(x + 1, L) });

	expr F = multifactorHenselLiftingPolyExpr(f, H, L, 5, 4);

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
