#include "MathSystem/Simplification/Division.hpp"
#include "test.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_divisions() {
  Expr exp0 = Expr(4) / Expr(2);
  Expr exp1 = Expr("x") / Expr("x");
  Expr exp2 = (Expr(2) * Expr("x")) / Expr("x");
  Expr exp3 = power(Expr("x"), 2) / Expr("x");

  Expr res_exp0 = reduceDivisionAST(exp0);
  Expr res_exp1 = reduceDivisionAST(exp1);
  Expr res_exp2 = reduceDivisionAST(exp2);
  Expr res_exp3 = reduceDivisionAST(exp3);

  assert(res_exp0.kind() == Kind::Integer);
  assert(res_exp0.value() == 2);

	assert(res_exp1.kind() == Kind::Integer);
  assert(res_exp1.value() == 1);

  assert(res_exp2.kind() == Kind::Integer);
  assert(res_exp2.value() == 2);

  assert(res_exp3.kind() == Kind::Symbol);
  assert(res_exp3.identifier() == "x");
}

int main() {
  TEST(should_simplify_divisions)
  return 0;
}
