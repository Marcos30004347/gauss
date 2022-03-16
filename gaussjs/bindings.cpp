#ifdef WASM_BUILD

#include <emscripten/bind.h>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "gauss/Gauss.hpp"

namespace gaussjs {

struct Scope {
  std::vector<alg::expr> ctx;
};

Scope *createScope() { return new Scope(); }

void destroyScope(Scope *scope) {
  for (alg::expr e : scope.ctx) {
    e.~expr();
  }

  scope.ctx.~vector<alg::expr>();

  delete scope;
}

alg::expr &scopeAddExpr(Scope *scope, alg::expr &&e) {
  scope->ctx.push_back(e);

  return scope->ctx[scope.ctx.size() - 1];
}

namespace algebra {

typedef gauss::expr expr;
typedef gauss::expr kind;

expr &intFromString(Scope *scope, const char *v) {
        return scopeAddExpr(scope, gauss::algebra::intFromString(v);
}

expr &numberFromDouble(Scope *scope, double v) {
        return scopeAddExpr(scope, gauss::algebra::numberFromDouble(v);
}

expr &intFromLong(Scope *scope, long v) {
  return scopeAddExpr(scope, gauss::algebra::intFromLong(v);
}

expr &symbol(Scope *scope, const char *v) {
  return scopeAddExpr(scope, gauss::algebra::symbol(v));
}

expr pow(Scope *scope, expr a, expr e) {
  return scopeAddExpr(scope, gauss::algebra::pow(a, e));
}

expr &getOperand(expr a, size_t i) { return gauss::algebra::getOperand(a, i); }

void setOperand(expr &a, size_t i, expr b) {
  gauss::algebra::setOperand(a, i, b);
}

expr sqrt(Scope *scope, expr a) {
  return scopeAddExpr(scope, gauss::algebra::sqrt(a));
}

expr root(Scope *scope, expr a, expr b) {
  return scopeAddExpr(scope, gauss::algebra::root(a, b));
}

kind kindOf(expr a) { return gauss::algebra::kindOf(a); }

bool is(expr a, int k) { return gauss::algebra::is(a, k); }

expr &add(scope *scope, expr a, expr b) {
        return scopeAddExpr(scope, gauss::algebra::add(a, b);
}

expr &sub(scope *scope, expr a, expr b) {
        return scopeAddExpr(scope, gauss::algebra::sub(a, b);
}

expr &mul(scope *scope, expr a, expr b) {
        return scopeAddExpr(scope, gauss::algebra::mul(a, b);
}

expr &div(scope *scope, expr a, expr b) {
        return scopeAddExpr(scope, gauss::algebra::div(a, b);
}

expr &expand(scope *scope, expr a) {
        return scopeAddExpr(scope, gauss::algebra::expand(a);
}

expr &reduce(scope *scope, expr a) {
        return scopeAddExpr(scope, gauss::algebra::reduce(a);
}

expr &log(scope *scope, expr x, expr b) {
        return scopeAddExpr(scope, gauss::algebra::log(x, b);
}

expr &exp(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::exp(x);
}

expr &ln(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::ln(x);
}

expr &replace(scope *scope, expr u, expr x, expr v) {
        return scopeAddExpr(scope, gauss::algebra::replace(u, x, v);
}

expr &eval(scope *scope, expr u, expr x, expr v) {
        return scopeAddExpr(scope, gauss::algebra::eval(u, x, v);
}

expr &freeVariables(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::freeVariables(x);
}

namespace trigonometry {

expr &sinh(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::sinh(x);
}

expr &cosh(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::cosh(x);
}

expr &tanh(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::tanh(x);
}

expr &cos(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::cos(x);
}

expr &sin(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::sin(x);
}

expr &tan(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::tan(x);
}

expr &csc(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::csc(x);
}

expr &cot(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::cot(x);
}

expr &sec(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::sec(x);
}

expr &coth(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::coth(x);
}

expr &sech(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::sech(x);
}

expr &csch(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::csch(x);
}

expr &arccos(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::arccos(x);
}

expr &arcsin(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::arcsin(x);
}

expr &arctan(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::arctan(x);
}

expr &arccot(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::arccot(x);
}

expr &arcsec(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::arcsec(x);
}

expr &arccsc(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::arccsc(x);
}

expr &arccosh(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::arccosh(x);
}

expr &arctanh(scope *scope, expr x) {
        return scopeAddExpr(scope, gauss::algebra::trigonometry::arctanh(x);
}

}; // namespace trigonometry

namespace linear {

expr &matrix(scope *scope, unsigned l, unsigned c) {
  return scopeAddExpr(scope, gauss::algebra::linear::matrix(l, c));
}

expr &identity(scope *scope, unsigned l, unsigned c) {
  return scopeAddExpr(scope, gauss::algebra::linear::identity(l, c));
}

expr &matrixGet(scope *scope, expr &A, unsigned i, unsigned j) {
  return scopeAddExpr(scope, gauss::algebra::linear::matrixGet(A, i, j));
}

void matrixSet(scope *scope, expr &A, unsigned i, unsigned j, double a) {
  return scopeAddExpr(scope, gauss::algebra::linear::matrixSet(A, i, j, a));
}

expr &svd(scope *scope, expr &A) {
  return scopeAddExpr(scope, gauss::algebra::linear::svd(A));
}

expr &inverse(scope *scope, expr &A) {
  return scopeAddExpr(scope, gauss::algebra::linear::inverse(A));
}

expr &det(scope *scope, expr &A) {
  return scopeAddExpr(scope, gauss::algebra::linear::det(A));
}

expr &transpose(scope *scope, expr &A) {
  return scopeAddExpr(scope, gauss::algebra::linear::transpose(A));
}

expr &solveLinear(scope *scope, expr &A, expr &b) {
  return scopeAddExpr(scope, gauss::algebra::linear::solveLinear(A, b));
}

} // namespace linear

} // namespace algebra

namespace polynomial {
algebra::expr degreePoly(scope *scope, algebra::expr &f, algebra::expr &x) {
  return scopeAddExpr(scope, gauss::polynomial::degreePoly(f, x));
}

algebra::expr coefficientPoly(scope *scope, algebra::expr &f, algebra::expr &x,
                              algebra::expr &d) {
  return scopeAddExpr(scope, gauss::polynomial::coefficientPoly(f, x, d));
}

algebra::expr leadingCoefficientPoly(scope *scope, algebra::expr &f,
                                     algebra::expr &x) {
  return scopeAddExpr(scope, gauss::polynomial::leadingCoefficientPoly(f, x));
}

algebra::expr rootsOfPoly(scope *scope, algebra::expr &a) {
  return scopeAddExpr(scope, gauss::polynomial::rootsOfPoly(a));
}

algebra::expr factorPoly(scope *scope, algebra::expr f) {
  return scopeAddExpr(scope, gauss::polynomial::factorPoly(f));
}

algebra::expr resultantOfPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::resultantOfPoly(a, b));
}

algebra::expr addPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::addPoly(a, b));
}

algebra::expr subPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::subPoly(a, b));
}

algebra::expr mulPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::mulPoly(a, b));
}

algebra::expr divPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::divPoly(a, b));
}

algebra::expr quoPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::quoPoly(a, b));
}

algebra::expr remPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::remPoly(a, b));
}

algebra::expr gcdPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::gcdPoly(a, b));
}

algebra::expr lcmPoly(scope *scope, algebra::expr a, algebra::expr b) {
  return scopeAddExpr(scope, gauss::polynomial::lcmPoly(a, b));
}

namespace finiteField {

algebra::expr projectPolyFiniteField(scope *scope, algebra::expr a,
                                     long long p) {
  return scopeAddExpr(
      scope, gauss::polynomial::finiteField::projectPolyFiniteField(a, p));
}

algebra::expr addPolyFiniteField(scope *scope, algebra::expr a, algebra::expr b,
                                 long long p) {
  return scopeAddExpr(
      scope, gauss::polynomial::finiteField::addPolyFiniteField(a, b, p));
}

algebra::expr subPolyFiniteField(scope *scope, algebra::expr a, algebra::expr b,
                                 long long p) {
  return scopeAddExpr(
      scope, gauss::polynomial::finiteField::subPolyFiniteField(a, b, p));
}

algebra::expr mulPolyFiniteField(scope *scope, algebra::expr a, algebra::expr b,
                                 long long p) {
  return scopeAddExpr(
      scope, gauss::polynomial::finiteField::mulPolyFiniteField(a, b, p));
}

algebra::expr divPolyFiniteField(scope *scope, algebra::expr a, algebra::expr b,
                                 long long p) {
  return scopeAddExpr(
      scope, gauss::polynomial::finiteField::divPolyFiniteField(a, b, p));
}

algebra::expr quoPolyFiniteField(scope *scope, algebra::expr a, algebra::expr b,
                                 long long p) {
  return scopeAddExpr(
      scope, gauss::polynomial::finiteField::quoPolyFiniteField(a, b, p));
}

algebra::expr remPolyFiniteField(scope *scope, algebra::expr a, algebra::expr b,
                                 long long p) {
  return scopeAddExpr(
      scope, gauss::polynomial::finiteField::remPolyFiniteField(a, b, p));
}

} // namespace finiteField

} // namespace polynomial

namespace calculus {

algebra::expr derivative(scope *scope, algebra::expr a, algebra::expr x) {
  return scopeAddExpr(scope, gauss::calculus::dderivative(scope, a, x));
}

} // namespace calculus

std::string toString(algebra::expr a) { return gauss::toString(a); }

std::string toLatex(algebra::expr a, bool print_as_fractions,
                    unsigned long max_den) {
  return gauss::toLatex(a, print_as_fractions, max_den);
}

}
}

EMSCRIPTEN_BINDINGS(gauss_module) {
  emscripten::class_<gaussjs::Scope>("Scope");

  emscripten::class_<gaussjs::algebra::expr>("expr");

  emscripten::enum_<gaussjs::algebra::kind>("kind")
      .value("ERROR", gaussjs::algebra::kind::ERROR)
      .value("FACT", gaussjs::algebra::kind::FACT)
      .value("POW", gaussjs::algebra::kind::POW)
      .value("MUL", gaussjs::algebra::kind::MUL)
      .value("ADD", gaussjs::algebra::kind::ADD)
      .value("SUB", gaussjs::algebra::kind::SUB)
      .value("DIV", gaussjs::algebra::kind::DIV)
      .value("SQRT", gaussjs::algebra::kind::SQRT)
      .value("INF", gaussjs::algebra::kind::INF)
      .value("UNDEF", gaussjs::algebra::kind::UNDEF)
      .value("SYM", gaussjs::algebra::kind::SYM)
      .value("INT", gaussjs::algebra::kind::INT)
      .value("FRAC", gaussjs::algebra::kind::FRAC)
      .value("FAIL", gaussjs::algebra::kind::FAIL)
      .value("FUNC", gaussjs::algebra::kind::FUNC);

  emscripten::function("createScope", &gaussjs::createScope);
  emscripten::function("destroyScope", &gaussjs::destroyScope);

  emscripten::function("intFromString", &gaussjs::algebra::intFromString);
  emscripten::function("numberFromDouble", &gaussjs::algebra::numberFromDouble);
  emscripten::function("intFromDouble", &gaussjs::algebra::intFromLong);

  emscripten::function("symbol", &gaussjs::algebra::symbol);

  emscripten::function("add", &gaussjs::algebra::add);
  emscripten::function("sub", &gaussjs::algebra::sub);
  emscripten::function("mul", &gaussjs::algebra::mul);
  emscripten::function("div", &gaussjs::algebra::div);

  emscripten::function("pow", &gaussjs::algebra::pow);
  emscripten::function("root", &gaussjs::algebra::root);
  emscripten::function("sqrt", &gaussjs::algebra::sqrt);

  emscripten::function("is", &gaussjs::algebra::is);
  emscripten::function("kindOf", &gaussjs::algebra::kindOf);

  emscripten::function("expand", &gaussjs::algebra::expand);
  emscripten::function("reduce", &gaussjs::algebra::reduce);

  emscripten::function("ln", &gaussjs::algebra::ln);
  emscripten::function("log", &gaussjs::algebra::log);
  emscripten::function("exp", &gaussjs::algebra::exp);

  emscripten::function("eval", &gaussjs::algebra::eval);
  emscripten::function("replace", &gaussjs::algebra::replace);
  emscripten::function("freeVariables", &gaussjs::algebra::freeVariables);

  emscripten::function("sinh", &gaussjs::algebra::trigonometry::sinh);
  emscripten::function("cosh", &gaussjs::algebra::trigonometry::cosh);
  emscripten::function("tanh", &gaussjs::algebra::trigonometry::tanh);
  emscripten::function("cos", &gaussjs::algebra::trigonometry::cos);
  emscripten::function("sin", &gaussjs::algebra::trigonometry::sin);
  emscripten::function("tan", &gaussjs::algebra::trigonometry::tan);
  emscripten::function("csc", &gaussjs::algebra::trigonometry::csc);
  emscripten::function("cot", &gaussjs::algebra::trigonometry::cot);
  emscripten::function("sec", &gaussjs::algebra::trigonometry::sec);
  emscripten::function("coth", &gaussjs::algebra::trigonometry::coth);
  emscripten::function("sech", &gaussjs::algebra::trigonometry::sech);
  emscripten::function("csch", &gaussjs::algebra::trigonometry::csch);
  emscripten::function("arccos", &gaussjs::algebra::trigonometry::arccos);
  emscripten::function("arcsin", &gaussjs::algebra::trigonometry::arcsin);
  emscripten::function("arctan", &gaussjs::algebra::trigonometry::arctan);
  emscripten::function("arccot", &gaussjs::algebra::trigonometry::arccot);
  emscripten::function("arcsec", &gaussjs::algebra::trigonometry::arcsec);
  emscripten::function("arccsc", &gaussjs::algebra::trigonometry::arccsc);
  emscripten::function("arccosh", &gaussjs::algebra::trigonometry::arccosh);
  emscripten::function("arctanh", &gaussjs::algebra::trigonometry::arctanh);

  emscripten::function("svd", &gaussjs::algebra::linear::svd);
  emscripten::function("det", &gaussjs::algebra::linear::det);
  emscripten::function("matrix", &gaussjs::algebra::linear::matrix);
  emscripten::function("transpose", &gaussjs::algebra::linear::transpose);
  emscripten::function("solveLinear", &gaussjs::algebra::linear::solveLinear);
  emscripten::function("identity", &gaussjs::algebra::linear::identity);
  emscripten::function("matrixGet", &gaussjs::algebra::linear::matrixGet);
  emscripten::function("matrixSet", &gaussjs::algebra::linear::matrixSet);

  emscripten::function("degreePoly", &gaussjs::polynomial::degreePoly);
  emscripten::function("coefficientPoly",
                       &gaussjs::polynomial::coefficientPoly);
  emscripten::function("leadingCoefficientPoly",
                       &gaussjs::polynomial::leadingCoefficientPoly);
  emscripten::function("rootsOfPoly", &gaussjs::polynomial::rootsOfPoly);
  emscripten::function("factorPoly", &gaussjs::polynomial::factorPoly);
  emscripten::function("resultantOfPoly",
                       &gaussjs::polynomial::resultantOfPoly);
  emscripten::function("addPoly", &gaussjs::polynomial::addPoly);
  emscripten::function("subPoly", &gaussjs::polynomial::subPoly);
  emscripten::function("mulPoly", &gaussjs::polynomial::mulPoly);
  emscripten::function("divPoly", &gaussjs::polynomial::divPoly);
  emscripten::function("remPoly", &gaussjs::polynomial::remPoly);
  emscripten::function("gcdPoly", &gaussjs::polynomial::gcdPoly);
  emscripten::function("lcmPoly", &gaussjs::polynomial::lcmPoly);

  emscripten::function("addPolyFiniteField",
                       &gaussjs::polynomial::finiteField::addPolyFiniteField);
  emscripten::function("subPolyFiniteField",
                       &gaussjs::polynomial::finiteField::subPolyFiniteField);
  emscripten::function("mulPolyFiniteField",
                       &gaussjs::polynomial::finiteField::mulPolyFiniteField);
  emscripten::function("divPolyFiniteField",
                       &gaussjs::polynomial::finiteField::divPolyFiniteField);
  emscripten::function("quoPolyFiniteField",
                       &gaussjs::polynomial::finiteField::quoPolyFiniteField);
  emscripten::function("remPolyFiniteField",
                       &gaussjs::polynomial::finiteField::remPolyFiniteField);
  emscripten::function(
      "projectPolyFiniteField",
      &gaussjs::polynomial::finiteField::projectPolyFiniteField);

  emscripten::function("derivative", &gaussjs::calculus::derivative);

  emscripten::function("toString", &gaussjs::calculus::toString);
  emscripten::function("toLatex", &gaussjs::calculus::toLatex);
}

#endif
