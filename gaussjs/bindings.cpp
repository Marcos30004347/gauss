#include <cstddef>
#include <cstdint>

#include <emscripten/bind.h>

// // TODO: only include on debug builds
// #include <sanitizer/lsan_interface.h>

#include "gauss/Gauss.hpp"
#include "gauss/Error/error.hpp"

#include <cmath>
#include <inttypes.h>
#include <limits>
#include <string>
#include <vector>

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

  emscripten::function("numberFromString", &gauss::algebra::numberFromString,
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
  emscripten::function("arcsech", &gauss::algebra::trigonometry::arcsech,
                       emscripten::allow_raw_pointers());
  emscripten::function("arccsch", &gauss::algebra::trigonometry::arccsch,
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

  // emscripten::function("doLeakCheck", &__lsan_do_recoverable_leak_check);
};
