#include "test.hpp"

#include "gauss/Polynomial/Polynomial.hpp"
#include "gauss/Factorization/SquareFree.hpp"

using namespace alg;
using namespace poly;
using namespace factorization;

void should_get_square_free_factorization_poly_expr() {
  expr x = expr("x");

  expr L = list({x});
	expr Z = expr("Z");

  expr ax = polyExpr(pow(x, 8) + -2 * pow(x, 6) + 2 * pow(x, 2) + -1, L);
  expr bx = polyExpr(pow(x, 11) + 2 * pow(x, 9) + 2 * pow(x, 8) + pow(x, 6) +
										 pow(x, 5) + 2 * pow(x, 3) + 2 * pow(x, 2) + 1, L);
  assert(squareFreeFactorizationPolyExpr(ax, L, Z) ==
         pow(1 * pow(x, 0) + 1 * pow(x, 2), 1) *
				 pow(-1 * pow(x, 0) + 1 * pow(x, 2), 3));

  assert(squareFreeFactorizationPolyExpr2(ax, L, Z) ==
         pow(1 * pow(x, 0) + 1 * pow(x, 2), 1) *
				 pow(-1 * pow(x, 0) + 1 * pow(x, 2), 3));

	assert(squareFreeFactorizationFiniteFieldPolyExpr(bx, L, Z, 3, false) ==
				 pow(1*pow(x, 0) + 1*pow(x, 1), 1)*
				 pow(2*pow(x, 0) + 1*pow(x, 1), 4)*
				 pow(1*pow(x, 0) + 1*pow(x, 2), 3)
				 );
}

// void should_compute_square_free_part() {
// 	expr x = expr("x");
// 	expr y = expr("y");

//   expr t = pow(x, 3) + 2*pow(x, 2)*y + pow(y, 2)*x;

// 	expr L = list({ x, y });

//   expr K = expr("Z");

//   assert(squareFreePart(t, L, K)[0] == pow(x, 2) + x*y);
// }

void should_compute_square_free_part_poly_expr() {
	expr x = expr("x");
	expr y = expr("y");

	expr L = list({ x, y });

  expr t = polyExpr(pow(x, 3) + 2*pow(x, 2)*y + pow(y, 2)*x, L);

  expr K = expr("Z");

	assert(squareFreePartPolyExpr(t, L, K)[0] ==
				 create(kind::ADD, {1*pow(y, 1)})*pow(x, 1) + create(kind::ADD, {1*pow(y, 0)})*pow(x, 2));
}


int main() {
  // TEST(should_get_square_free_factorization)
	TEST(should_get_square_free_factorization_poly_expr)
	// TEST(should_compute_square_free_part)
	TEST(should_compute_square_free_part_poly_expr)

	return 0;
}
