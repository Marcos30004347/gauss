
#include <assert.h>

#include "gauss/Simplification/Simplification.hpp"
#include "Division.hpp"
#include "Expand.hpp"
#include "Multiplication.hpp"
#include "Polynomial.hpp"

using namespace ast;
using namespace simplification;
using namespace algebra;

namespace expand {

/**
 * Tink aboud adding the following expansions:
 *
 * sin(x+y) = sin(x)cos(y)+cos(x)sin(y)
 * cos(2x) = 2cos(x)^2 - 1
 * e^(a+ln(b)) = (e^a)b
 */

Expr expandAST(Expr u) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Fraction ||
      u.kind() == Kind::FunctionCall || u.kind() == Kind::Undefined ||
      u.kind() == Kind::Symbol || u.kind() == Kind::Infinity ||
      u.kind() == Kind::MinusInfinity || u.kind() == Kind::List ||
      u.kind() == Kind::Set) {
    return u;
  }

  Expr k = reduceAST(mapUnaryAST(u, expandAST));


  if (k.kind() == Kind::Power) {
    k = expandMultinomialAST(k);
  } else if (k.kind() == Kind::Multiplication) {
    k = expandMultiplicationAST(k);
  } else if (k.kind() == Kind::Division) {
    k = expandDivisionAST(k);
  }

  return reduceAST(k);
}

} // namespace expand
