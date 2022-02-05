#include "Core/Calculus/Derivative.hpp"
#include "test.hpp"

#include "Core/AST/AST3.hpp"
#include "Core/Calculus/Calculus.hpp"

using namespace alg;
using namespace calculus;

void should_derivate_expressions() {
  expr x = expr("x");
  assert(derivate(4 * x + pow(x, 2) + 5 * pow(x, 3), x) ==
         4 + 2 * x + 15 * pow(x, 2));

  assert(derivate(sinh(x), x) == cosh(x));

  assert(derivate(cosh(x), x) == sinh(x));

  assert(derivate(tanh(x), x) == pow(sech(x), 2));

  assert(derivate(sin(x), x) == cos(x));

  assert(derivate(cos(x), x) == -sin(x));

  assert(derivate(tan(x), x) == pow(sec(x), 2));

  assert(derivate(cot(x), x) == -pow(csc(x), 2));

  assert(derivate(sec(x), x) == sec(x) * tan(x));

  assert(derivate(csc(x), x) == -cot(x) * csc(x));

  assert(derivate(coth(x), x) == -pow(csch(x), 2));

  assert(derivate(sech(x), x) == -sech(x) * tanh(x));

  assert(derivate(csch(x), x) == -coth(x) * csch(x));

  assert(derivate(pow(x, fraction(1, 2)), x) ==
         fraction(1, 2) * pow(x, fraction(-1, 2)));
}

int main() { TEST(should_derivate_expressions) }
