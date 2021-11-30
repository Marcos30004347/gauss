#include "Simplification.hpp"
#include "Addition.hpp"
#include "Division.hpp"
#include "Factorial.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"
#include "Rationals.hpp"
#include "Subtraction.hpp"
#include "Trigonometry.hpp"

#include <cstdio>

using namespace ast;
using namespace algebra;

namespace simplification {

Expr reduceAST(Expr u) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Symbol ||
      u.kind() == Kind::Infinity || u.kind() == Kind::MinusInfinity ||
      u.kind() == Kind::Undefined || u.kind() == Kind::List ||
      u.kind() == Kind::Set)
    return u;

  if (u.kind() == Kind::Fraction) {
    return reduceRNEAST(u);
  }

  Expr v = mapUnaryAST(u, reduceAST);

  if (v.kind() == Kind::Addition) {
    return reduceAdditionAST(v);
  }
  if (v.kind() == Kind::Subtraction) {
    return reduceSubtractionAST(v);
  }
  if (v.kind() == Kind::Division) {
    return reduceDivisionAST(v);
  }
  if (v.kind() == Kind::Multiplication) {
    return reduceMultiplicationAST(v);
  }
  if (v.kind() == Kind::Power) {
    return reducePowerExpr(v);
  }
  if (v.kind() == Kind::Factorial) {
    return reduceFactorialAST(v);
  }

  return v;
}

} // namespace simplification
