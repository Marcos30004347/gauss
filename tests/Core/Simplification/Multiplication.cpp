#include "test.hpp"
#include "Core/Simplification/Multiplication.hpp"
#include <cstdio>

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_products() {
	Expr x = Expr("x");
	Expr y = Expr("y");

	Expr exp0 = Expr(2) * Expr(2);
  Expr exp1 = Expr(Kind::Multiplication, {Expr(2) * Expr(3), Expr(4) * Expr(5)});
  Expr exp2 = x*x;
  Expr exp3 = Expr(Kind::Multiplication, {2*x, 2*x});
	Expr exp4 = 2*x * 4*y * 6*x * 14*y * 10 * x;

  assert(reduceMultiplicationExpr(exp0) == 4);
  assert(reduceMultiplicationExpr(exp1) == 120);
  assert(reduceMultiplicationExpr(exp2) == power(Expr("x"), 2));
	assert(reduceMultiplicationExpr(exp3) == 4*power(Expr("x"), 2));
	assert(reduceMultiplicationExpr(exp4) == 6720 * power(x, 3) * power(y, 2));
}

int main() {
  TEST(should_simplify_products)
  return 0;
}
