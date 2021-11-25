#include <assert.h>
#include <cstdio>

#include "Core/AST/AST.hpp"
#include "Core/Simplification/Subtraction.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

void should_simplify_subtractions() {
  Expr exp0 = Expr(3) - Expr(2);
  Expr exp1 = Expr(2) - Expr(3);
  Expr exp2 = Expr("x") - Expr("x");
  Expr exp3 = Expr(Kind::Subtraction, {
                                          (Expr("a") + Expr("b") + Expr("c")),
                                          (Expr("d") - Expr("e")),
                                      });

  Expr exp4 = Expr(Kind::Subtraction, {
                                          Expr("a") + Expr("b") + Expr("c"),
                                          Expr("d") - Expr("e"),
                                          Expr("f") - Expr("g"),
                                      });

  Expr exp5 = Expr(Kind::Subtraction, {
                                          Expr("a") - Expr("b") - Expr("c"),
                                          Expr("d") - Expr("e"),
                                          Expr("f") - Expr("g"),
                                      });

  Expr res_exp0 = reduceSubtractionAST(exp0);
  Expr res_exp1 = reduceSubtractionAST(exp1);
  Expr res_exp2 = reduceSubtractionAST(exp2);
  Expr res_exp3 = reduceSubtractionAST(exp3);
  Expr res_exp4 = reduceSubtractionAST(exp4);
  Expr res_exp5 = reduceSubtractionAST(exp5);

  assert(res_exp0 == 1);
  assert(res_exp1 == -1);
  assert(res_exp2 == 0);
  assert(res_exp3 ==
         Expr("a") + Expr("b") + Expr("c") + -1 * Expr("d") + Expr("e"));
  assert(res_exp4 == Expr("a") + Expr("b") + Expr("c") + -1 * Expr("d") +
                         Expr("e") + -1 * Expr("f") + Expr("g"));
  assert(res_exp5 == Expr("a") + -1 * Expr("b") + -1 * Expr("c") +
                         -1 * Expr("d") + Expr("e") + -1 * Expr("f") +
                         Expr("g"));
}

int main() {
  should_simplify_subtractions();
  return 0;
}
