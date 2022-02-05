#include "test.hpp"
#include "MathSystem/Algebra/Algebra.hpp"
#include "MathSystem/Simplification/Rationals.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_rational_number_expressions() {
  Expr exp0 = Expr(1) + Expr(2);
  Expr exp1 = Expr(1) - Expr(2);
  Expr exp2 = Expr(2) - Expr(1);
  Expr exp3 = Expr(2) * Expr(2);
  Expr exp4 = Expr(4) / Expr(2);
  Expr exp5 = power(Expr(4), Expr(2));
  Expr exp6 = fraction(1, 2) + fraction(1, 2);
  Expr exp7 = fraction(1, 2) - fraction(1, 2);
  Expr exp8 = fraction(1, 2) * fraction(1, 2);
  Expr exp9 = fraction(1, 2) / fraction(1, 2);
  Expr exp10 = power(fraction(1, 2), Expr(2));
  Expr exp11 = Expr(3) + fraction(1, 2);
  Expr exp12 = fraction(1, 2) + Expr(3);
  Expr exp13 = Expr(3) - fraction(1, 2);
  Expr exp14 = fraction(1, 2) - Expr(3);

  Expr res_exp0 = reduceRNEAST(exp0);
  assert(res_exp0 == 3);
	Expr res_exp1 = reduceRNEAST(exp1);
  assert(res_exp1 == -1);
	Expr res_exp2 = reduceRNEAST(exp2);
  assert(res_exp2 == 1);
	Expr res_exp3 = reduceRNEAST(exp3);
  assert(res_exp3 == 4);
	Expr res_exp4 = reduceRNEAST(exp4);
  assert(res_exp4 == 2);
	Expr res_exp5 = reduceRNEAST(exp5);
  assert(res_exp5 == 16);
	Expr res_exp6 = reduceRNEAST(exp6);
  assert(res_exp6 == 1);
	Expr res_exp7 = reduceRNEAST(exp7);
  assert(res_exp7 == 0);
	Expr res_exp8 = reduceRNEAST(exp8);
  assert(res_exp8 == fraction(1, 4));
	Expr res_exp9 = reduceRNEAST(exp9);
  assert(res_exp9 == 1);
	Expr res_exp10 = reduceRNEAST(exp10);
  assert(res_exp10 == fraction(1, 4));
	Expr res_exp11 = reduceRNEAST(exp11);
  assert(res_exp11 == fraction(7, 2));
	Expr res_exp12 = reduceRNEAST(exp12);
  assert(res_exp12 == fraction(7, 2));
	Expr res_exp13 = reduceRNEAST(exp13);
  assert(res_exp13 == fraction(5, 2));
	Expr res_exp14 = reduceRNEAST(exp14);
  assert(res_exp14 == fraction(-5, 2));
}

int main() {
  TEST(should_simplify_rational_number_expressions)
  return 0;
}
