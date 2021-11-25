#include "Division.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"

using namespace ast;
using namespace algebra;

namespace simplification {

Expr reduceDivisionAST(Expr u) {
  Expr p = power(u[1], -1);
  Expr m = u[0] * reducePowerAST(p);

  return reduceMultiplicationAST(m);
}

} // namespace simplification
