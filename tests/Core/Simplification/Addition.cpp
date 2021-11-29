#include "Core/Simplification/Addition.hpp"
#include "Core/AST/AST.hpp"
#include "test.hpp"
#include <cassert>

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_additions() {
  Expr x = Expr("x");
  Expr y = Expr("y");
  Expr z = Expr("z");

  Expr exp0 = Expr(2) + Expr(2);
  Expr exp1 = Expr(2) + Expr(4) + Expr(5) + Expr(6);
  Expr exp2 = Expr(3) + Expr(5) + Expr(6);
  Expr exp3 = Expr(5) + Expr(6) + Expr(3);
  Expr exp4 = Expr(2) + Expr(3) + Expr(4) + Expr(5);
  Expr exp5 = x + x;
  Expr exp6 = 2 + x + 2 + x;
  Expr exp7 = x + 3 + y + 5;
  Expr exp8 = x + 1 + 2 + x + 3 + y + 4 + z + 5;
  Expr exp9 = Expr(1) + Expr(2) + Expr(3) + Expr(4) + Expr(5) + Expr(6) +
              Expr(7) + Expr(8) + Expr(9);
  Expr exp10 =
      Expr(Kind::Addition, {2 * (x + 1) + 3 * x + 4 * (x + 1) + 4,
                            Expr(Kind::Addition, {3 + x + y, 3 * x + 4})});
	Expr exp11 = x + y + z + -x;
	Expr exp12 = -x + 2*x + x + -x + y + -y;
	Expr exp13 = inf() + x + y + 10;
	Expr exp14 = -inf() + x + y + z + 14;
	Expr exp15 = x + y + inf() + -14 + -inf();

	assert(reduceAdditionExpr(exp0) == 4);
  assert(reduceAdditionExpr(exp1) == 17);
  assert(reduceAdditionExpr(exp2) == 14);
  assert(reduceAdditionExpr(exp3) == 14);
  assert(reduceAdditionExpr(exp4) == 14);
  assert(reduceAdditionExpr(exp5) == 2 * x);
  assert(reduceAdditionExpr(exp6) == 2 * x + 4);
  assert(reduceAdditionExpr(exp7) == x + y + 8);
  assert(reduceAdditionExpr(exp8) == 2 * x + y + z + 15);
  assert(reduceAdditionExpr(exp9) == 45);
	assert(reduceAdditionExpr(exp10) == 7*x + y + 6*(x + 1) + 11);
	assert(reduceAdditionExpr(exp11) == y + z);
	assert(reduceAdditionExpr(exp12) == x);
	assert(reduceAdditionExpr(exp13) == inf());
	assert(reduceAdditionExpr(exp14) == -inf());
	assert(reduceAdditionExpr(exp15) == undefined());
}

int main() {
  TEST(should_simplify_additions)
  return 0;
}
