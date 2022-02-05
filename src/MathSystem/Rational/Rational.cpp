#include "Rational.hpp"
#include "MathSystem/Algebra/Set.hpp"
#include "MathSystem/Polynomial/Polynomial.hpp"
#include "MathSystem/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace rational {

bool isRationalExpression(Expr u, Expr S) {
  Expr n = numerator(u);
  Expr d = denominator(u);

  bool r = isGerenalPolynomial(n, S) && isGerenalPolynomial(d, S);

  return r;
}

Expr rationalVariables(Expr u) {
  Expr n = numerator(u);
  Expr d = denominator(u);

  Expr K = variables(n);
  Expr J = variables(d);

  Expr R = unification(K, J);

  return R;
}

Expr rationalizeSum(Expr u, Expr v) {
  Expr m = numerator(u);
  Expr r = denominator(u);
  Expr n = numerator(v);
  Expr s = denominator(v);

  if (r == 1 && s == 1) {

    Expr t = u + v;
    Expr k = reduceAST(t);

    return k;
  }

  Expr num_a = m * s;
  Expr num_b = n * r;
  Expr den = r * s;

  Expr num = rationalizeSum(num_a, num_b);

  return div(num, den);
}

Expr rationalize(Expr u) {

  if (u.kind() == Kind::Power) {
    return power(rationalize(u[0]), reduceAST(u[1]));
  }

  if (u.kind() == Kind::Multiplication) {
    Expr f = u[0];

    Expr k = reduceAST(u / f);

    return rationalize(f) * rationalize(k);
  }

  if (u.kind() == Kind::Addition) {
    Expr f = u[0];

    Expr k = reduceAST(u - f);

    Expr g = rationalize(f);
    Expr r = rationalize(k);

    Expr t = rationalizeSum(g, r);

    return t;
  }

  return u;
}

Expr numerator(Expr u) {
  if (u.kind() == Kind::Fraction || u.kind() == Kind::Division)
    return u[0];

  if (u.kind() == Kind::Power) {
    if (u[1].kind() == Kind::Integer && u[1].value() < 0) {
      return integer(1);
    }
    return u;
  }

  if (u.kind() == Kind::Multiplication) {
    if (u.size() == 1) {
      return numerator(u[0]);
    }
    Expr v = u[0];

    Expr h_ = div(u, v);

    Expr h = reduceAST(h_);

    Expr r_ = mul({numerator(v), numerator(h)});

    Expr r = reduceAST(r_);

    return r;
  }

  return u;
}

Expr denominator(Expr u) {
  if (u.kind() == Kind::Fraction || u.kind() == Kind::Division)
    return u[1];

  if (u.kind() == Kind::Power) {
    if (u[1] < 0) {
      return reduceAST(power(u, -1));
    }

    return 1;
  }

  if (u.kind() == Kind::Multiplication) {
    if (u.size() == 1) {
      return denominator(u[0]);
    }

    Expr v = u[0];

    Expr h = reduceAST(u / v);
    Expr r = reduceAST(denominator(v) * denominator(h));

    return r;
  }

  return 1;
}

Expr expandRational(Expr u) {
  Expr n_ = numerator(u);
  Expr d_ = denominator(u);

  Expr n = algebraicExpand(n_);
  Expr d = algebraicExpand(d_);

  Expr k = div(n, d);

  return k;
}

} // namespace rational
