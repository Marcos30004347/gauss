#include <climits>
#include <cstddef>
#include <limits>
#include <vector>

#include "Algebra/Expression.hpp"
#include "Error/error.hpp"
#include "Gauss.hpp"

#include "gauss/Error/error.hpp"
#include "gauss/Algebra/Expression.hpp"
#include "gauss/Algebra/Reduction.hpp"
#include "gauss/Algebra/Trigonometry.hpp"
#include "gauss/Calculus/Derivative.hpp"
#include "gauss/Factorization/SquareFree.hpp"
#include "gauss/GaloisField/GaloisField.hpp"
#include "gauss/Polynomial/Polynomial.hpp"
#include "gauss/Polynomial/Resultant.hpp"
#include "gauss/Polynomial/Roots.hpp"
#include "gauss/Primes/Primes.hpp"

namespace gauss {

namespace algebra {

expr numberFromDouble(double v) {
  Int n = 0, d = 1;

  double integral = 0;
  double fractional = 0;

  fractional = std::modf(v, &integral);

  alg::decimalToFraction(fractional, 99999999999999, n, d);

  alg::expr r = Int(integral) + alg::fraction(n, d);

  return alg::reduce(r);
}

expr intFromString(const char *v) {
  // TODO: accept rational numbers with '.'
  return expr(Int::fromString(v));
}

expr intFromLong(long v) { return expr(v); }

expr symbol(std::string s) { return expr(s.c_str()); }

expr pow(expr a, expr b) { return alg::pow(a, b); }

expr &getOperand(expr a, size_t i) { return a[i]; }

void setOperand(expr &a, size_t i, expr b) { a[i] = b; }

kind kindOf(expr a) { return a.kind(); }

expr root(expr a, expr b) { return alg::sqrt(a, b); }

bool is(expr a, int k) { return alg::is(&a, k); }

expr sqrt(expr a) { return alg::sqrt(a, 2); }

expr add(expr a, expr b) { return a + b; }

expr sub(expr a, expr b) { return a - b; }

expr mul(expr a, expr b) { return a * b; }

expr div(expr a, expr b) { return a / b; }

expr expand(expr a) {
  expand(&a);

  return a;
}

expr reduce(expr a) { return alg::reduce(a); }

expr replace(expr u, expr x, expr v) {
  if (x.kind() != alg::kind::SYM) {
		raise(error(ErrorCode::ARG_IS_NOT_SYM_EXPR, 1));
  }

  return alg::replace(u, x, v);
}

expr eval(expr u, expr x, expr v) {
  if (x.kind() != alg::kind::SYM) {
		raise(error(ErrorCode::ARG_IS_NOT_SYM_EXPR, 1));
  }

  return algebra::expand(algebra::replace(u, x, v));
}

expr freeVariables(expr u) { return alg::freeVariables(u); }

bool isEqual(expr a, expr b) {
	alg::expand(&a);
	alg::expand(&b);

	return alg::reduce(a - b) == 0;
}

expr numerator(expr a) {
	if(is(&a, kind::FRAC | kind::DIV)) {
		return a[0];
	}

	return a;
}

expr denominator(expr a) {
	if(is(&a, kind::FRAC | kind::DIV)) {
		return a[1];
	}

	return 1;
}

expr powDegree(expr a) {
	if(!is(&a, kind::POW)) {
		raise(error(ErrorCode::ARG_IS_NOT_POW_EXPR, 0));
	}

	return a[1];
}

expr powBase(expr a) {
	if(!is(&a, kind::POW)) {
		raise(error(ErrorCode::ARG_IS_NOT_POW_EXPR, 0));
	}

  return a[0];
}

expr rootIndex(expr a) {
	if(!is(&a, kind::SQRT)) {
		raise(error(ErrorCode::ARG_IS_NOT_ROOT_EXPR, 0));
	}

	return a[1];
}

expr rootRadicand(expr a) {
	if(!is(&a, kind::SQRT)) {
		raise(error(ErrorCode::ARG_IS_NOT_ROOT_EXPR, 0));
	}

  return a[0];
}

expr log(expr a, expr b) { return alg::log(a, b); }

expr ln(expr a) { return alg::ln(a); }

expr exp(expr a) { return alg::exp(a); }

expr abs(expr a) { return alg::abs(a); }

expr trigonometry::sinh(expr x) { return alg::trig::sinh(x); }

expr trigonometry::cosh(expr x) { return alg::trig::cosh(x); }

expr trigonometry::tanh(expr x) { return alg::trig::tanh(x); }

expr trigonometry::cos(expr x) { return alg::trig::cos(x); }

expr trigonometry::sin(expr x) { return alg::trig::sin(x); }

expr trigonometry::tan(expr x) { return alg::trig::tan(x); }

expr trigonometry::csc(expr x) { return alg::trig::csc(x); }

expr trigonometry::cot(expr x) { return alg::trig::cot(x); }

expr trigonometry::sec(expr x) { return alg::trig::sec(x); }

expr trigonometry::coth(expr x) { return alg::trig::coth(x); }

expr trigonometry::sech(expr x) { return alg::trig::sech(x); }

expr trigonometry::csch(expr x) { return alg::trig::csch(x); }

expr trigonometry::arccos(expr x) { return alg::trig::arccos(x); }

expr trigonometry::arcsin(expr x) { return alg::trig::arcsin(x); }

expr trigonometry::arctan(expr x) { return alg::trig::arctan(x); }

expr trigonometry::arccot(expr x) { return alg::trig::arccot(x); }

expr trigonometry::arcsec(expr x) { return alg::trig::arcsec(x); }

expr trigonometry::arccsc(expr x) { return alg::trig::arccsc(x); }

expr trigonometry::arccosh(expr x) { return alg::trig::arccosh(x); }

expr trigonometry::arctanh(expr x) { return alg::trig::arctanh(x); }

expr linear::matrix(unsigned int l, unsigned int c) { return alg::mat(l, c); }

expr linear::identity(unsigned int l, unsigned int c) {
  return alg::identity_matrix(l, c);
}

expr linear::matrixGet(expr A, unsigned int i, unsigned int j) {
  return alg::mat_get(A, i, j);
}

void linear::matrixSet(expr A, unsigned int i, unsigned int j, double v) {
  return alg::mat_set(A, i, j, numberFromDouble(v));
}

expr linear::svd(expr A) { return alg::svd_matrix(A); }

expr linear::inverse(expr A) { return alg::inverse_matrix(A); }

expr linear::det(expr A) { return alg::determinant_matrix(A); }

expr linear::transpose(expr A) { return alg::transpose_matrix(A); }

expr linear::solveLinear(expr A, expr b) {
  return alg::solve_linear_system(A, b);
}

} // namespace algebra

expr polynomial::factorPoly(expr poly) {
  expr L = poly::getVariableListForPolyExpr(poly);

  expr p = poly::polyExpr(poly, L);

  return poly::factorPolyExprAndExpand(p, L, expr("Q"));
}

expr polynomial::degreePoly(expr f, expr x) { return poly::degree(f, x); }

expr polynomial::coefficientPoly(expr f, expr x, expr d) {
  return poly::coeff(f, x, d);
}

expr polynomial::leadingCoefficientPoly(expr f, expr x) {
  return poly::coeff(f, x, poly::degree(f, x));
}

expr polynomial::resultantOfPoly(expr a, expr b) {
  expr T = poly::normalizeToPolyExprs(a, b);
  return poly::resultantPolyExpr(T[1], T[2], T[0], expr("Q"));
}

expr polynomial::rootsOfPoly(expr a) { return poly::realPolyRoots(a); }

expr polynomial::addPoly(expr a, expr b) {
  alg::expand(&a);
  alg::expand(&b);

  expr c = a + b;

  alg::reduce(&c);

  return c;
}

expr polynomial::subPoly(expr a, expr b) {
  alg::expand(&a);
  alg::expand(&b);

  expr c = a - b;

  alg::reduce(&c);

  return c;
}

expr polynomial::mulPoly(expr a, expr b) {
  alg::expand(&a);
  alg::expand(&b);

  expr c = a * b;

  alg::reduce(&c);

  return c;
}

expr polynomial::divPoly(expr a, expr b) {
  expr T = poly::normalizeToPolyExprs(a, b);
  expr L = T[0];

  expr K = expr("Z");

  expr D = poly::divPolyExpr(T[1], T[2], L, K);

  return D[0] + D[1];
}

expr polynomial::quoPoly(expr a, expr b) {
  expr T = poly::normalizeToPolyExprs(a, b);
  expr L = T[0];

  expr K = expr("Z");

  expr D = poly::divPolyExpr(T[1], T[2], L, K);

  return D[0];
}

expr polynomial::remPoly(expr a, expr b) {
  expr T = poly::normalizeToPolyExprs(a, b);
  expr L = T[0];

  expr K = expr("Z");

  expr D = poly::divPolyExpr(T[1], T[2], L, K);

  return D[1];
}

expr polynomial::gcdPoly(expr a, expr b) {
  expr T = poly::normalizeToPolyExprs(a, b);

  expr L = T[0];

  expr K = expr("Q");

  expr D = poly::gcdPolyExpr(T[1], T[2], L, K);

  return D[1];
}

expr polynomial::lcmPoly(expr a, expr b) {
  expr T = poly::normalizeToPolyExprs(a, b);

  expr L = T[0];

  expr K = expr("Q");

  expr D = poly::lcmPolyExpr(T[1], T[2], L, K);

  return D[1];
}

expr polynomial::finiteField::projectPolyFiniteField(expr a, long long p) {
  expr L = poly::getVariableListForPolyExpr(a);

  expr u = poly::polyExpr(a, L);

  return galoisField::gfPolyExpr(u, p, false);
}

expr polynomial::finiteField::addPolyFiniteField(expr a, expr b, long long p) {
  expr T = poly::normalizeToPolyExprs(a, b);

  expr D = galoisField::addPolyExprGf(T[1], T[2], p);

  return D;
}

expr polynomial::finiteField::subPolyFiniteField(expr a, expr b, long long p) {
  expr T = poly::normalizeToPolyExprs(a, b);

  expr D = galoisField::subPolyExprGf(T[1], T[2], p);

  return D;
}

expr polynomial::finiteField::mulPolyFiniteField(expr a, expr b, long long p) {
  expr T = poly::normalizeToPolyExprs(a, b);

  expr D = galoisField::mulPolyExprGf(T[1], T[2], p);

  return D;
}

expr polynomial::finiteField::divPolyFiniteField(expr a, expr b, long long p) {
  expr T = poly::normalizeToPolyExprs(a, b);

  expr L = T[0];

  expr D = galoisField::divPolyExprGf(T[1], T[2], L, p);

  return D[0] + D[1];
}

expr polynomial::finiteField::quoPolyFiniteField(expr a, expr b, long long p) {
  expr T = poly::normalizeToPolyExprs(a, b);

  expr L = T[0];

  expr D = galoisField::quoPolyExprGf(T[1], T[2], L, p);

  return D;
}

expr polynomial::finiteField::remPolyFiniteField(expr a, expr b, long long p) {
  expr T = poly::normalizeToPolyExprs(a, b);

  expr L = T[0];

  expr D = galoisField::remPolyExprGf(T[1], T[2], L, p);

  return D;
}

expr calculus::derivative(expr a, expr x) { return calc::derivate(a, x); }

std::string toString(expr a) { return alg::to_string(&a); }

std::string toLatex(expr a, bool p, unsigned long k) {
  return alg::to_latex(&a, p, k);
}

expr algebra::prime(size_t i) { return intFromLong(primes[i]); }

expr algebra::primeFactors(expr a) {
  if (!is(&a, kind::INT)) {
		raise(error(ErrorCode::ARG_IS_NOT_INT_EXPR, 0));
  }

  if (abs(a.value()) > std::numeric_limits<unsigned long long>::max()) {
		raise(error(ErrorCode::INT_BIGGER_THAN_MAX_ULL, 0));
  }

  char s = 1;

  Int v = a.value();

  if (v < 0) {
    s = -1;
    v = -v;
  }

  std::vector<unsigned long long> f = primes.factorsOf(v.longValue());

  expr F = 1;
  size_t start = 0;
  if (s == -1) {
    F = -1;
  } else {
    F = Int(f[0]);
    start = 1;
  }

  for (size_t i = start; i < f.size(); i++) {
    F = F * Int(f[i]);
  }

  return F;
}

} // namespace gauss
