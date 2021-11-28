#include "Core/AST/AST.hpp"
#include "test.hpp"
#include "Core/Simplification/Addition.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_additions() {
  Expr exp0 = Expr(2) + Expr(2);
  Expr exp1 = Expr(2) + Expr(4) + Expr(5) + Expr(6);
  Expr exp2 = Expr(3) + Expr(5) + Expr(6);
  Expr exp3 = Expr(5) + Expr(6) + Expr(3);
  Expr exp4 = Expr(2) + Expr(3) + Expr(4) + Expr(5);
  Expr exp5 = Expr("x") + Expr("x");
  Expr exp6 = Expr(2) + Expr("x") + Expr(2) + Expr("x");
	/*
  Expr res_exp0 = reduceAdditionAST(exp0);
  Expr res_exp1 = reduceAdditionAST(exp1);
  Expr res_exp2 = reduceAdditionAST(exp2);
  Expr res_exp3 = reduceAdditionAST(exp3);
  Expr res_exp4 = reduceAdditionAST(exp4);
  Expr res_exp5 = reduceAdditionAST(exp5);
  Expr res_exp6 = reduceAdditionAST(exp6);

  assert(res_exp0 == 4);
  assert(res_exp1 == 17);
  assert(res_exp2 == 14);
  assert(res_exp3 == 14);
  assert(res_exp4 == 14);
  assert(res_exp5 == 2 * Expr("x"));
  assert(res_exp6 == 4 + 2 * Expr("x"));
	*/
	printf("%s\n", reduceAdditionExpr(Expr("x") + Expr(3) + Expr("y") + Expr(5)).toString().c_str());
	//printf("%s\n", reduceAdditionAST(Expr("x") + Expr(3) + Expr("y") + Expr(5)).toString().c_str());


	//printf("%s\n", reduceAdditionExpr(Expr("w") + Expr(1) + Expr(2) + Expr("x") + Expr(3) + Expr("y") + Expr(4) + Expr("z") + Expr(5)).toString().c_str());
	//printf("%s\n", reduceAdditionAST(Expr("w") + Expr(1) + Expr(2) + Expr("x") + Expr(3) + Expr("y") + Expr(4) + Expr("z") + Expr(5)).toString().c_str());

	printf("%s\n", reduceAdditionExpr(Expr("x") + Expr(1) + Expr(2) + Expr("x") + Expr(3) + Expr("y") + Expr(4) + Expr("z") + Expr(5)).toString().c_str());
	//printf("%s\n", reduceAdditionAST(Expr("x") + Expr(1) + Expr(2) + Expr("x") + Expr(3) + Expr("y") + Expr(4) + Expr("z") + Expr(5)).toString().c_str());

	printf("%s\n", reduceAdditionExpr(Expr(1) + Expr(2) + Expr(3) + Expr(4) + Expr(5) + Expr(6) + Expr(7) + Expr(8) + Expr(9)).toString().c_str());
	//printf("%s\n", reduceAdditionAST(Expr(1) + Expr(2) + Expr(3) + Expr(4) + Expr(5) + Expr(6) + Expr(7) + Expr(8) + Expr(9)).toString().c_str());
}

int main() {
  TEST(should_simplify_additions)
  return 0;
}
