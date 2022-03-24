#include <cstddef>
#include <cstdint>

#include <emscripten/bind.h>

// TODO: only include on debug builds
#include <sanitizer/lsan_interface.h>

#include "gauss/Gauss.hpp"
#include "gauss/Error/error.hpp"

#include <cmath>
#include <inttypes.h>
#include <limits>
#include <string>
#include <vector>

// namespace gaussjs {

// typedef gauss::expr expr;
// typedef gauss::kind kind;

// struct Scope {
//   std::vector<gauss::expr> ctx;
// };

// Scope createScope() { return Scope(); }

// void destroyScope(Scope &scope) {
// 	//delete scope;
// }

// expr *scopeAddExpr(Scope &scope, gauss::expr &&e) {
//   // scope->ctx.push_back(e);

//   // return &scope->ctx[scope->ctx.size() - 1];
// }

// expr &scopeAAddExpr(Scope &scope, gauss::expr &&e) {
//   scope.ctx.push_back(e);

//   return scope.ctx[scope.ctx.size() - 1];
// }

// namespace algebra {

// expr *intFromString(Scope &scope, const char *v) {
//   return scopeAddExpr(scope, gauss::algebra::intFromString(v));
// }

// expr& numberFromDouble(Scope &scope, double v) {
// 	return scopeAAddExpr(scope, gauss::algebra::numberFromDouble(v)) ;
// }

// expr *intFromLong(Scope &scope, long v) {
//   return scopeAddExpr(scope, gauss::algebra::intFromLong(v));
// }

// expr *symbol(Scope &scope, const char *v) {
//   return scopeAddExpr(scope, gauss::algebra::symbol(v));
// }

// expr *pow(Scope &scope, expr *a, expr *e) {
//   return scopeAddExpr(scope, gauss::algebra::pow(*a, *e));
// }

// expr *getOperand(expr *a, size_t i) {
//   return &gauss::algebra::getOperand(*a, i);
// }

// void setOperand(Scope &scope, expr *a, size_t i, expr *b) {
//   gauss::algebra::setOperand(*a, i, *b);
// }

// expr *sqrt(Scope &scope, expr *a) {
//   return scopeAddExpr(scope, gauss::algebra::sqrt(*a));
// }

// expr *root(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::algebra::root(*a, *b));
// }

// kind kindOf(expr *a) { return gauss::algebra::kindOf(*a); }

// bool is(expr *a, int k) { return gauss::algebra::is(*a, k); }

// expr *add(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::algebra::add(*a, *b));
// }

// expr *sub(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::algebra::sub(*a, *b));
// }

// expr *mul(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::algebra::mul(*a, *b));
// }

// expr *div(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::algebra::div(*a, *b));
// }

// expr *expand(Scope &scope, expr *a) {
//   return scopeAddExpr(scope, gauss::algebra::expand(*a));
// }

// expr *reduce(Scope &scope, expr *a) {
//   return scopeAddExpr(scope, gauss::algebra::reduce(*a));
// }

// expr *log(Scope &scope, expr *x, expr *b) {
//   return scopeAddExpr(scope, gauss::algebra::log(*x, *b));
// }

// expr *exp(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::exp(*x));
// }

// expr *ln(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::ln(*x));
// }

// expr *replace(Scope &scope, expr *u, expr *x, expr *v) {
//   return scopeAddExpr(scope, gauss::algebra::replace(*u, *x, *v));
// }

// expr *eval(Scope &scope, expr *u, expr *x, expr *v) {
//   return scopeAddExpr(scope, gauss::algebra::eval(*u, *x, *v));
// }

// expr *freeVariables(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::freeVariables(*x));
// }

// namespace trigonometry {

// expr *sinh(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::sinh(*x));
// }

// expr *cosh(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::cosh(*x));
// }

// expr *tanh(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::tanh(*x));
// }

// expr *cos(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::cos(*x));
// }

// expr *sin(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::sin(*x));
// }

// expr *tan(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::tan(*x));
// }

// expr *csc(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::csc(*x));
// }

// expr *cot(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::cot(*x));
// }

// expr *sec(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::sec(*x));
// }

// expr *coth(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::coth(*x));
// }

// expr *sech(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::sech(*x));
// }

// expr *csch(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::csch(*x));
// }

// expr *arccos(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::arccos(*x));
// }

// expr *arcsin(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::arcsin(*x));
// }

// expr *arctan(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::arctan(*x));
// }

// expr *arccot(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::arccot(*x));
// }

// expr *arcsec(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::arcsec(*x));
// }

// expr *arccsc(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::arccsc(*x));
// }

// expr *arccosh(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::arccosh(*x));
// }

// expr *arctanh(Scope &scope, expr *x) {
//   return scopeAddExpr(scope, gauss::algebra::trigonometry::arctanh(*x));
// }

// }; // namespace trigonometry

// namespace linear {

// expr *matrix(Scope &scope, unsigned l, unsigned c) {
//   return scopeAddExpr(scope, gauss::algebra::linear::matrix(l, c));
// }

// expr *identity(Scope &scope, unsigned l, unsigned c) {
//   return scopeAddExpr(scope, gauss::algebra::linear::identity(l, c));
// }

// expr *matrixGet(Scope &scope, expr *A, unsigned i, unsigned j) {
//   return scopeAddExpr(scope, gauss::algebra::linear::matrixGet(*A, i, j));
// }

// void matrixSet(expr *A, unsigned i, unsigned j, double a) {
//   return gauss::algebra::linear::matrixSet(*A, i, j, a);
// }

// expr *svd(Scope &scope, expr *A) {
//   return scopeAddExpr(scope, gauss::algebra::linear::svd(*A));
// }

// expr *inverse(Scope &scope, expr *A) {
//   return scopeAddExpr(scope, gauss::algebra::linear::inverse(*A));
// }

// expr *det(Scope &scope, expr *A) {
//   return scopeAddExpr(scope, gauss::algebra::linear::det(*A));
// }

// expr *transpose(Scope &scope, expr *A) {
//   return scopeAddExpr(scope, gauss::algebra::linear::transpose(*A));
// }

// expr *solveLinear(Scope &scope, expr *A, expr *b) {
//   return scopeAddExpr(scope, gauss::algebra::linear::solveLinear(*A, *b));
// }

// } // namespace linear

// } // namespace algebra

// namespace polynomial {
// expr *degreePoly(Scope &scope, expr *f, expr *x) {
//   return scopeAddExpr(scope, gauss::polynomial::degreePoly(*f, *x));
// }

// expr *coefficientPoly(Scope &scope, expr *f, expr *x, expr *d) {
//   return scopeAddExpr(scope, gauss::polynomial::coefficientPoly(*f, *x, *d));
// }

// expr *leadingCoefficientPoly(Scope &scope, expr *f, expr *x) {
//   return scopeAddExpr(scope,
//                        gauss::polynomial::leadingCoefficientPoly(*f, *x));
// }

// expr *rootsOfPoly(Scope &scope, expr *a) {
//   return scopeAddExpr(scope, gauss::polynomial::rootsOfPoly(*a));
// }

// expr *factorPoly(Scope &scope, expr *f) {
//   return scopeAddExpr(scope, gauss::polynomial::factorPoly(*f));
// }

// expr *resultantOfPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::resultantOfPoly(*a, *b));
// }

// expr *addPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::addPoly(*a, *b));
// }

// expr *subPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::subPoly(*a, *b));
// }

// expr *mulPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::mulPoly(*a, *b));
// }

// expr *divPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::divPoly(*a, *b));
// }

// expr *quoPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::quoPoly(*a, *b));
// }

// expr *remPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::remPoly(*a, *b));
// }

// expr *gcdPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::gcdPoly(*a, *b));
// }

// expr *lcmPoly(Scope &scope, expr *a, expr *b) {
//   return scopeAddExpr(scope, gauss::polynomial::lcmPoly(*a, *b));
// }

// namespace finiteField {

// expr *projectPolyFiniteField(Scope &scope, expr *a, long long p) {
//   return scopeAddExpr(
//       scope, gauss::polynomial::finiteField::projectPolyFiniteField(*a, p));
// }

// expr *addPolyFiniteField(Scope &scope, expr *a, expr *b, long long p) {
//   return scopeAddExpr(
//       scope, gauss::polynomial::finiteField::addPolyFiniteField(*a, *b, p));
// }

// expr *subPolyFiniteField(Scope &scope, expr *a, expr *b, long long p) {
//   return scopeAddExpr(
//       scope, gauss::polynomial::finiteField::subPolyFiniteField(*a, *b, p));
// }

// expr *mulPolyFiniteField(Scope &scope, expr *a, expr *b, long long p) {
//   return scopeAddExpr(
//       scope, gauss::polynomial::finiteField::mulPolyFiniteField(*a, *b, p));
// }

// expr *divPolyFiniteField(Scope &scope, expr *a, expr *b, long long p) {
//   return scopeAddExpr(
//       scope, gauss::polynomial::finiteField::divPolyFiniteField(*a, *b, p));
// }

// expr *quoPolyFiniteField(Scope &scope, expr *a, expr *b, long long p) {
//   return scopeAddExpr(
//       scope, gauss::polynomial::finiteField::quoPolyFiniteField(*a, *b, p));
// }

// expr *remPolyFiniteField(Scope &scope, expr *a, expr *b, long long p) {
//   return scopeAddExpr(
//       scope, gauss::polynomial::finiteField::remPolyFiniteField(*a, *b, p));
// }

// } // namespace finiteField

// } // namespace polynomial

// namespace calculus {

// expr *derivative(Scope &scope, expr *a, expr *x) {
//   return scopeAddExpr(scope, gauss::calculus::derivative(*a, *x));
// }

// } // namespace calculus

// std::string toString(expr& a) {
// 	printf("aaaaaaaa\n");
// 	return gauss::toString(a);
// }

// std::string toLatex(expr *a, bool print_as_fractions, unsigned long max_den) {
//   return gauss::toLatex(*a, print_as_fractions, max_den);
// }

// } // namespace gaussjs

EMSCRIPTEN_BINDINGS(gauss) {
  emscripten::class_<gauss::expr>("expr");

  emscripten::enum_<gauss::kind>("kind")
      .value("FACT", gauss::kind::FACT)
      .value("POW", gauss::kind::POW)
      .value("MUL", gauss::kind::MUL)
      .value("ADD", gauss::kind::ADD)
      .value("SUB", gauss::kind::SUB)
      .value("DIV", gauss::kind::DIV)
		  .value("ROOT", gauss::kind::ROOT)
      .value("INF", gauss::kind::INF)
      .value("UNDEF", gauss::kind::UNDEF)
      .value("SYM", gauss::kind::SYM)
      .value("INT", gauss::kind::INT)
      .value("FRAC", gauss::kind::FRAC)
      .value("FAIL", gauss::kind::FAIL)
      .value("FUNC", gauss::kind::FUNC);

	emscripten::enum_<ErrorCode>("ErrorCode")
		.value("INT_BIGGET_THAN_MAX_ULL", ErrorCode::INT_BIGGER_THAN_MAX_ULL)
		.value("DIVISION_BY_ZERO", ErrorCode::DIVISION_BY_ZERO)
		.value("INT_HAVE_NO_MODULAR_INVERSE", ErrorCode::NUMBER_HAVE_NO_MODULAR_INVERSE)
		.value("ARG_IS_INVALID", ErrorCode::ARG_IS_INVALID)
		.value("POLY_HAVE_NON_INTEGER_DEGREE", ErrorCode::POLY_HAVE_NON_INTEGER_DEGREE)
		.value("POLY_HAVE_NON_CONSTANT_COEFFICIENT", ErrorCode::POLY_HAVE_NON_CONSTANT_COEFFICIENT)
		.value("ARG_IS_NOT_SYM_EXPR", ErrorCode::ARG_IS_NOT_SYM_EXPR)
		.value("MATRIX_INVALID_ADDITION", ErrorCode::MATRIX_INVALID_ADDITION)
		.value("MATRIX_INVALID_MUTIPLICATION", ErrorCode::MATRIX_INVALID_MULTIPLICATION)
		.value("ARG_IS_NOT_POW_EXPR", ErrorCode::ARG_IS_NOT_POW_EXPR)
		.value("ARG_IS_NOT_INT_EXPR", ErrorCode::ARG_IS_NOT_INT_EXPR)
		.value("ARG_IS_NOT_FRA_EXPR", ErrorCode::ARG_IS_NOT_FRA_EXPR)
		.value("ARG_IS_NOT_ADD_EXPR", ErrorCode::ARG_IS_NOT_ADD_EXPR)
		.value("ARG_IS_NOT_SUB_EXPR", ErrorCode::ARG_IS_NOT_SUB_EXPR)
		.value("ARG_IS_NOT_ROOT_EXPR", ErrorCode::ARG_IS_NOT_ROOT_EXPR)
		.value("ARG_IS_NOT_POLY_EXPR", ErrorCode::ARG_IS_NOT_POLY_EXPR)
		.value("ARG_IS_NOT_LIST_EXPR", ErrorCode::ARG_IS_NOT_LIST_EXPR)
		.value("ARG_IS_IMAGINARY", ErrorCode::ARG_IS_IMAGINARY)
		.value("ARG_IS_NOT_UNIVARIATE_POLY", ErrorCode::ARG_IS_NOT_UNIVARIATE_POLY);

	emscripten::function("errorArg", &errorArg);
	emscripten::function("errorCode", &errorCode);

  emscripten::function("intFromString", &gauss::algebra::intFromString,
                       emscripten::allow_raw_pointers());

  emscripten::function("numberFromDouble", &gauss::algebra::numberFromDouble,
                       emscripten::allow_raw_pointers());

  emscripten::function("intFromLong", &gauss::algebra::intFromLong,
                       emscripten::allow_raw_pointers());

  emscripten::function("symbol", &gauss::algebra::symbol,
                       emscripten::allow_raw_pointers());

  emscripten::function("add", &gauss::algebra::add,
                       emscripten::allow_raw_pointers());
  emscripten::function("sub", &gauss::algebra::sub,
                       emscripten::allow_raw_pointers());
  emscripten::function("mul", &gauss::algebra::mul,
                       emscripten::allow_raw_pointers());
  emscripten::function("div", &gauss::algebra::div,
                       emscripten::allow_raw_pointers());

  emscripten::function("pow", &gauss::algebra::pow,
                       emscripten::allow_raw_pointers());
  emscripten::function("root", &gauss::algebra::root,
                       emscripten::allow_raw_pointers());
  emscripten::function("sqrt", &gauss::algebra::sqrt,
                       emscripten::allow_raw_pointers());

  emscripten::function("getOperand", &gauss::algebra::getOperand,
                       emscripten::allow_raw_pointers());
  emscripten::function("setOperand", &gauss::algebra::setOperand,
                       emscripten::allow_raw_pointers());

  emscripten::function("is", &gauss::algebra::is,
                       emscripten::allow_raw_pointers());
  emscripten::function("kindOf", &gauss::algebra::kindOf,
                       emscripten::allow_raw_pointers());

  emscripten::function("expand", &gauss::algebra::expand,
                       emscripten::allow_raw_pointers());
  emscripten::function("reduce", &gauss::algebra::reduce,
                       emscripten::allow_raw_pointers());

  emscripten::function("ln", &gauss::algebra::ln,
                       emscripten::allow_raw_pointers());
  emscripten::function("log", &gauss::algebra::log,
                       emscripten::allow_raw_pointers());
  emscripten::function("exp", &gauss::algebra::exp,
                       emscripten::allow_raw_pointers());

  emscripten::function("eval", &gauss::algebra::eval,
                       emscripten::allow_raw_pointers());
  emscripten::function("replace", &gauss::algebra::replace,
                       emscripten::allow_raw_pointers());
  emscripten::function("freeVariables", &gauss::algebra::freeVariables,
                       emscripten::allow_raw_pointers());

  emscripten::function("sinh", &gauss::algebra::trigonometry::sinh,
                       emscripten::allow_raw_pointers());
  emscripten::function("cosh", &gauss::algebra::trigonometry::cosh,
                       emscripten::allow_raw_pointers());
  emscripten::function("tanh", &gauss::algebra::trigonometry::tanh,
                       emscripten::allow_raw_pointers());
  emscripten::function("cos", &gauss::algebra::trigonometry::cos,
                       emscripten::allow_raw_pointers());
  emscripten::function("sin", &gauss::algebra::trigonometry::sin,
                       emscripten::allow_raw_pointers());
  emscripten::function("tan", &gauss::algebra::trigonometry::tan,
                       emscripten::allow_raw_pointers());
  emscripten::function("csc", &gauss::algebra::trigonometry::csc,
                       emscripten::allow_raw_pointers());
  emscripten::function("cot", &gauss::algebra::trigonometry::cot,
                       emscripten::allow_raw_pointers());
  emscripten::function("sec", &gauss::algebra::trigonometry::sec,
                       emscripten::allow_raw_pointers());
  emscripten::function("coth", &gauss::algebra::trigonometry::coth,
                       emscripten::allow_raw_pointers());
  emscripten::function("sech", &gauss::algebra::trigonometry::sech,
                       emscripten::allow_raw_pointers());
  emscripten::function("csch", &gauss::algebra::trigonometry::csch,
                       emscripten::allow_raw_pointers());
  emscripten::function("arccos", &gauss::algebra::trigonometry::arccos,
                       emscripten::allow_raw_pointers());
  emscripten::function("arcsin", &gauss::algebra::trigonometry::arcsin,
                       emscripten::allow_raw_pointers());
  emscripten::function("arctan", &gauss::algebra::trigonometry::arctan,
                       emscripten::allow_raw_pointers());
  emscripten::function("arccot", &gauss::algebra::trigonometry::arccot,
                       emscripten::allow_raw_pointers());
  emscripten::function("arcsec", &gauss::algebra::trigonometry::arcsec,
                       emscripten::allow_raw_pointers());
  emscripten::function("arccsc", &gauss::algebra::trigonometry::arccsc,
                       emscripten::allow_raw_pointers());
  emscripten::function("arccosh", &gauss::algebra::trigonometry::arccosh,
                       emscripten::allow_raw_pointers());
  emscripten::function("arctanh", &gauss::algebra::trigonometry::arctanh,
                       emscripten::allow_raw_pointers());

  emscripten::function("svd", &gauss::algebra::linear::svd,
                       emscripten::allow_raw_pointers());
  emscripten::function("det", &gauss::algebra::linear::det,
                       emscripten::allow_raw_pointers());
  emscripten::function("matrix", &gauss::algebra::linear::matrix,
                       emscripten::allow_raw_pointers());
  emscripten::function("transpose", &gauss::algebra::linear::transpose,
                       emscripten::allow_raw_pointers());
  emscripten::function("solveLinear", &gauss::algebra::linear::solveLinear,
                       emscripten::allow_raw_pointers());
  emscripten::function("identity", &gauss::algebra::linear::identity,
                       emscripten::allow_raw_pointers());
  emscripten::function("matrixGet", &gauss::algebra::linear::matrixGet,
                       emscripten::allow_raw_pointers());
  emscripten::function("matrixSet", &gauss::algebra::linear::matrixSet,
                       emscripten::allow_raw_pointers());

  emscripten::function("degreePoly", &gauss::polynomial::degreePoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("coefficientPoly", &gauss::polynomial::coefficientPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("leadingCoefficientPoly",
                       &gauss::polynomial::leadingCoefficientPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("rootsOfPoly", &gauss::polynomial::rootsOfPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("factorPoly", &gauss::polynomial::factorPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("resultantOfPoly", &gauss::polynomial::resultantOfPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("addPoly", &gauss::polynomial::addPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("subPoly", &gauss::polynomial::subPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("mulPoly", &gauss::polynomial::mulPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("divPoly", &gauss::polynomial::divPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("remPoly", &gauss::polynomial::remPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("gcdPoly", &gauss::polynomial::gcdPoly,
                       emscripten::allow_raw_pointers());
  emscripten::function("lcmPoly", &gauss::polynomial::lcmPoly,
                       emscripten::allow_raw_pointers());

  emscripten::function("addPolyFiniteField",
                       &gauss::polynomial::finiteField::addPolyFiniteField,
                       emscripten::allow_raw_pointers());
  emscripten::function("subPolyFiniteField",
                       &gauss::polynomial::finiteField::subPolyFiniteField,
                       emscripten::allow_raw_pointers());
  emscripten::function("mulPolyFiniteField",
                       &gauss::polynomial::finiteField::mulPolyFiniteField,
                       emscripten::allow_raw_pointers());
  emscripten::function("divPolyFiniteField",
                       &gauss::polynomial::finiteField::divPolyFiniteField,
                       emscripten::allow_raw_pointers());
  emscripten::function("quoPolyFiniteField",
                       &gauss::polynomial::finiteField::quoPolyFiniteField,
                       emscripten::allow_raw_pointers());
  emscripten::function("remPolyFiniteField",
                       &gauss::polynomial::finiteField::remPolyFiniteField,
                       emscripten::allow_raw_pointers());
  emscripten::function(
      "projectPolyFiniteField",
      &gauss::polynomial::finiteField::projectPolyFiniteField,
      emscripten::allow_raw_pointers());

  emscripten::function("derivative", &gauss::calculus::derivative,
                       emscripten::allow_raw_pointers());

  emscripten::function("toString", &gauss::toString,
                       emscripten::allow_raw_pointers());
  emscripten::function("toLatex", &gauss::toLatex,
                       emscripten::allow_raw_pointers());

  emscripten::function("doLeakCheck", &__lsan_do_recoverable_leak_check);
};
