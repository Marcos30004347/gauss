#include "Division.hpp"
#include "gauss/Algebra/Algebra.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"
#include <utility>

using namespace ast;
using namespace algebra;

namespace simplification {

Expr reduceDivisionAST(Expr&& u) {
  return reduceMultiplicationExpr(u[0] * reducePowerExpr(power(u[1], -1)));
}

Expr reduceDivisionAST(Expr& u) {
		return reduceDivisionAST(std::forward<Expr>(u));
}


} // namespace simplification
