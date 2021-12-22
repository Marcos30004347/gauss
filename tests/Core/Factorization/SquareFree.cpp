#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "test.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Factorization/SquareFree.hpp"
#include "Core/Polynomial/Polynomial.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_get_square_free_factorization() {
  Expr x = Expr("x");

  Expr ax = power(x, 8) + -2 * power(x, 6) + 2 * power(x, 2) + -1;
  Expr bx = power(x, 11) + 2 * power(x, 9) + 2 * power(x, 8) + power(x, 6) +
		power(x, 5) + 2 * power(x, 3) + 2 * power(x, 2) + 1;

  assert(squareFreeFactorization(ax, x) ==
         power(-1 + power(x, 2), 3) * (1 + power(x, 2)));
  assert(squareFreeFactorization2(ax, x) ==
         power(-1 + power(x, 2), 3) * (1 + power(x, 2)));
  assert(squareFreeFactorizationFiniteField(bx, x, 3, false) ==
         (x + 1) * power(power(x, 2) + 1, 3) * power(x + 2, 4));
}

void should_get_square_free_factorization_poly_expr() {
  Expr x = Expr("x");

  Expr L = list({x});
	Expr Z = Expr("Z");

  Expr ax = polyExpr(power(x, 8) + -2 * power(x, 6) + 2 * power(x, 2) + -1, L);
  Expr bx = polyExpr(power(x, 11) + 2 * power(x, 9) + 2 * power(x, 8) + power(x, 6) +
										 power(x, 5) + 2 * power(x, 3) + 2 * power(x, 2) + 1, L);

  assert(squareFreeFactorizationPolyExpr(ax, L, Z) ==
         power(1 * power(x, 0) + 1 * power(x, 2), 1) *
				 power(-1 * power(x, 0) + 1 * power(x, 2), 3));

  assert(squareFreeFactorizationPolyExpr2(ax, L, Z) ==
         power(1 * power(x, 0) + 1 * power(x, 2), 1) *
				 power(-1 * power(x, 0) + 1 * power(x, 2), 3));

	assert(squareFreeFactorizationFiniteFieldPolyExpr(bx, L, Z, 3, false) ==
				 power(1*power(x, 0) + 1*power(x, 1), 1)*
				 power(2*power(x, 0) + 1*power(x, 1), 4)*
				 power(1*power(x, 0) + 1*power(x, 2), 3)
				 );
}

void should_compute_square_free_part() {
	Expr x = Expr("x");
	Expr y = Expr("y");

  Expr t = power(x, 3) + 2*power(x, 2)*y + power(y, 2)*x;

	Expr L = list({ x, y });

  Expr K = Expr("Z");

  assert(squareFreePart(t, L, K)[0] == power(x, 2) + x*y);
}

void should_compute_square_free_part_poly_expr() {
	Expr x = Expr("x");
	Expr y = Expr("y");

	Expr L = list({ x, y });

  Expr t = polyExpr(power(x, 3) + 2*power(x, 2)*y + power(y, 2)*x, L);

  Expr K = Expr("Z");

	assert(squareFreePartPolyExpr(t, L, K)[0] ==
				 add({1*power(y, 1)})*power(x, 1) + add({1*power(y, 0)})*power(x, 2));
}


int main() {
  TEST(should_get_square_free_factorization)
	TEST(should_get_square_free_factorization_poly_expr)
	TEST(should_compute_square_free_part)
	TEST(should_compute_square_free_part_poly_expr)

	return 0;
}
