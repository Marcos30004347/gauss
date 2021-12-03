#include "Division.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"
#include <utility>

using namespace ast;
using namespace algebra;

namespace simplification {

Expr reduceDivisionAST(Expr&& u) {
  Expr p = power(u[1], -1);
  Expr m = u[0] * reducePowerExpr(p);

  return reduceMultiplicationAST(m);
}

Expr reduceDivisionAST(Expr& u) {
		return reduceDivisionAST(std::forward<Expr>(u));
}


} // namespace simplification
