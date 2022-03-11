#ifdef WASM_BUILD

#include <emscripten/bind.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "gauss/Algebra/Expression.hpp"

namespace gaussjs {

struct Scope {
  std::vector<alg::expr> ctx;
};

/**
 * @brief Creates a empty scope for expressions.
 *
 * @return Return a new empty scope
 */
Scope createScope() { return Scope(); }

/**
 * @brief Dealocate all resources used by the expressions inside the scope.
 *
 * @param[inout] scope - The scope being destroyed.
 */
void destroyScope(Scope &scope) {
  for (alg::expr e : scope.ctx) {
    e.~expr();
  }
  scope.ctx.~vector<alg::expr>();
}

	alg::expr &scopeAddExpr(Scope &scope, alg::expr e) {
  scope.ctx.push_back(e);

  return scope.ctx[scope.ctx.size() - 1];
}

alg::expr &ExpressionNumber(Scope &scope, double i) {
  double integral, fractional;

  fractional = std::modf(i, &integral);

  unsigned long long n, d;

  alg::toFraction(fractional, 1000, n, d);

  alg::expr r = Int(integral) + alg::fraction(n, d);

  reduce(&r);

  return scopeAddExpr(scope, r);
}

alg::expr &ExpressionSymbol(Scope &scope, std::string s) {
  return scopeAddExpr(scope, alg::expr(s.c_str()));
}

alg::expr &ExpressionAdd(Scope &scope, alg::expr a, alg::expr b) {
  return scopeAddExpr(scope, a + b);
}

alg::expr &ExpressionSub(Scope &scope, alg::expr a, alg::expr b) {
  return scopeAddExpr(scope, a - b);
}

alg::expr &ExpressionMul(Scope &scope, alg::expr a, alg::expr b) {
  return scopeAddExpr(scope, a * b);
}

alg::expr &ExpressionDiv(Scope &scope, alg::expr a, alg::expr b) {
  return scopeAddExpr(scope, a / b);
}

alg::expr &ExpressionPow(Scope &scope, alg::expr a, alg::expr b) {
  return scopeAddExpr(scope, pow(a, b));
}

alg::expr &ExpressionFact(Scope &scope, alg::expr a) {
  return scopeAddExpr(scope, fact(a));
}

alg::expr &ExpressionInfinity(Scope &scope) {
  return scopeAddExpr(scope, alg::inf());
}

alg::expr &ExpressionReduce(Scope &scope, alg::expr &a) {
  alg::expr b = a;

  alg::reduce(&b);

  return scopeAddExpr(scope, b);
}

alg::expr &ExpressionExpand(Scope &scope, alg::expr a) {
  alg::expr b = a;
  alg::expand(&b);
  return scopeAddExpr(scope, b);
}

alg::expr &ExpressionOperand(Scope &scope, alg::expr a, unsigned i) {
  return scopeAddExpr(scope, a[i]);
}

std::string ExpressionToString(alg::expr a) {
  return std::string(to_string(a).c_str());
}

} // namespace gaussjs

EMSCRIPTEN_BINDINGS(gauss_module) {
  emscripten::class_<alg::expr>("expr");

  emscripten::function("num", &gaussjs::ExpressionNumber);
  emscripten::function("sym", &gaussjs::ExpressionSymbol);
  emscripten::function("add", &gaussjs::ExpressionAdd);
  emscripten::function("sub", &gaussjs::ExpressionSub);
  emscripten::function("mul", &gaussjs::ExpressionMul);
  emscripten::function("div", &gaussjs::ExpressionDiv);
  emscripten::function("pow", &gaussjs::ExpressionPow);
  emscripten::function("fact", &gaussjs::ExpressionFact);
  emscripten::function("infinity", &gaussjs::ExpressionInfinity);
  emscripten::function("reduce", &gaussjs::ExpressionReduce);
  emscripten::function("expand", &gaussjs::ExpressionExpand);
  emscripten::function("operand", &gaussjs::ExpressionOperand);
  emscripten::function("toString", &gaussjs::ExpressionToString);
}

#endif
