#include "test.hpp"
#include "Core/Simplification/Multiplication.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_products() {
  Expr exp0 = Expr(2) + Expr(2);
  Expr exp1 =
      Expr(Kind::Multiplication, {Expr(2) * Expr(3), Expr(4) * Expr(5)});
  Expr exp2 = Expr("x") * Expr("x");
  Expr exp3 =
      Expr(Kind::Multiplication, {Expr(2) * Expr("x"), Expr(2) * Expr("x")});

  Expr res_exp0 = reduceMultiplicationAST(exp0);
  Expr res_exp1 = reduceMultiplicationAST(exp1);
  Expr res_exp2 = reduceMultiplicationAST(exp2);
  Expr res_exp3 = reduceMultiplicationAST(exp3);

  assert(res_exp0 == 4);
  assert(res_exp1 == 120);
  assert(res_exp2 == power(Expr("x"), 2));
	assert(res_exp3 == 4*power(Expr("x"), 2));
}

int main() {
  TEST(should_simplify_products)
  return 0;
}
