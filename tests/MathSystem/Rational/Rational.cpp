#include "MathSystem/Rational/Rational.hpp"
#include "MathSystem/Algebra/Set.hpp"
#include "MathSystem/Polynomial/Polynomial.hpp"
#include "MathSystem/Simplification/Simplification.hpp"
#include <assert.h>

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace poly;
using namespace simplification;

void should_get_numerator() {
  Expr u0 = (symbol("a") + symbol("b") + symbol("c"));

  Expr u1 = (symbol("a") + symbol("b") + symbol("c")) /
            (symbol("d") + symbol("e") + symbol("f"));

  assert(numerator(u0) == Expr("a") + Expr("b") + Expr("c"));
  assert(numerator(u1) == Expr("a") + Expr("b") + Expr("c"));
}

void should_get_denominators() {
  Expr u0 = (symbol("a") + symbol("b") + symbol("c"));

  Expr u1 = (symbol("a") + symbol("b") + symbol("c")) /
            (symbol("d") + symbol("e") + symbol("f"));

  assert(denominator(u0) == 1);
  assert(denominator(u1) == Expr("d") + Expr("e") + Expr("f"));
}

void should_expand_rational_expressions() {
  Expr u = symbol("a") / symbol("b") + symbol("c") / symbol("d") +
           symbol("e") / symbol("f");

  Expr u_ = rationalize(u);
  Expr v = expandRational(u_);

  Expr k = (symbol("b") * symbol("d") * symbol("e") +
            symbol("b") * symbol("c") * symbol("f") +
            symbol("a") * symbol("d") * symbol("f")) /
           (symbol("b") * symbol("d") * symbol("f"));

  assert(v == k);
}

int main() {
  should_get_numerator();
  should_expand_rational_expressions();
  should_get_denominators();
}
