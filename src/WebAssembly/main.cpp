#ifdef WASM_BUILD

#include <emscripten/bind.h>

#include <cmath>
#include <limits>
#include <string>

#include "MathSystem/Algebra/Expression.hpp"

alg::expr ExpressionNumber(double i) {
		double integral;

		double fractional = std::modf(i, &integral);

		unsigned long long n, d;

		alg::toFraction(fractional, 1000, n, d);

		alg::expr r = alg::Int(integral) + alg::fraction(n, d);

		reduce(&r);

		return r;
}

alg::expr ExpressionSymbol(std::string s) {
  return alg::expr(string(s.c_str()));
}

alg::expr ExpressionAdd(alg::expr a, alg::expr b) { return a + b; }

alg::expr ExpressionSub(alg::expr a, alg::expr b) { return a - b; }

alg::expr ExpressionMul(alg::expr a, alg::expr b) { return a * b; }

alg::expr ExpressionDiv(alg::expr a, alg::expr b) { return a / b; }

alg::expr ExpressionPow(alg::expr a, alg::expr b) { return pow(a, b); }

alg::expr ExpressionFact(alg::expr a) { return fact(a); }

alg::expr ExpressionInfinity() { return alg::inf(); }

alg::expr ExpressionReduce(alg::expr a) {
	alg::expr b = a;
	alg::reduce(&b);

	return b;
}

alg::expr ExpressionExpand(alg::expr a) {
	alg::expr b = a;
	alg::expand(&b);
	return b;
}

alg::expr ExpressionOperand(alg::expr a, unsigned i) { return a[i]; }

std::string ExpressionToString(alg::expr a){
	return std::string(to_string(a).c_str());
}

EMSCRIPTEN_BINDINGS(gauss_module) {
	emscripten::class_<alg::expr>("expr");

  emscripten::function("num", &ExpressionNumber);
  emscripten::function("sym", &ExpressionSymbol);
  emscripten::function("add", &ExpressionAdd);
  emscripten::function("sub", &ExpressionSub);
  emscripten::function("mul", &ExpressionMul);
  emscripten::function("div", &ExpressionDiv);
  emscripten::function("pow", &ExpressionPow);
  emscripten::function("fact", &ExpressionFact);
  emscripten::function("infinity", &ExpressionInfinity);
  emscripten::function("reduce", &ExpressionReduce);
  emscripten::function("expand", &ExpressionExpand);
  emscripten::function("operand", &ExpressionOperand);
  emscripten::function("toString", &ExpressionToString);
}
#endif
