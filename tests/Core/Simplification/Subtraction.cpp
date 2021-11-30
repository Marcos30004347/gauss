#include "Core/Simplification/Subtraction.hpp"
#include "Core/AST/AST.hpp"
#include "test.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_subtractions() {
  Expr x = Expr("x");
  Expr a = Expr("a");
  Expr b = Expr("b");
  Expr c = Expr("c");
  Expr d = Expr("d");
  Expr e = Expr("e");
  Expr f = Expr("f");
  Expr g = Expr("g");

  Expr exp0 = Expr(3) - Expr(2);
  Expr exp1 = Expr(2) - Expr(3);
  Expr exp2 = x - x;
  Expr exp3 = (a + b + c) - (d - e);

  Expr exp4 = Expr(Kind::Subtraction, {
                                          a + b + c,
                                          d - e,
                                          f - g,
                                      });

  Expr exp5 = Expr(Kind::Subtraction, {
                                          a - b - c,
                                          d - e,
                                          f - g,
                                      });

  Expr exp6 = Expr(1) - Expr(2) - Expr(3) - Expr(5) - Expr(7) - x - Expr(4) -
              Expr(6) - a;

	assert(reduceSubtractionExpr(exp0) == 1);
  assert(reduceSubtractionExpr(exp1) == -1);
  assert(reduceSubtractionExpr(exp2) == 0);
  assert(reduceSubtractionExpr(exp3) == a + b + c + -d + e);
  assert(reduceSubtractionExpr(exp4) == a + b + c + -d + e + -f + g);
  assert(reduceSubtractionExpr(exp5) == a + -b + -c + -d + e + -f + g);
	assert(reduceSubtractionExpr(exp6) == -a + -x + -26);
}

int main() {
  TEST(should_simplify_subtractions)
  return 0;
}
