/**
 * @file GaloisField.cpp
 * @author Marcos Vincius Moreira Santos (marcos30004347@gmail.com)
 * @brief This file implement some of the Finite field methods used in the
 * library
 * @version 0.1
 * @date 2021-10-11
 *
 * @ref Michael Monagan - In-place Arithmetic for Polynomials over Zn
 *
 * @copyright Copyright (c) 2021
 */

#include "GaloisField.hpp"

#include "Core/AST/AST.hpp"
#include "Core/AST/Integer.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Debug/Assert.hpp"

#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"

#include <climits>
#include <limits>
#include <random>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace galoisField {

Int mod(Int a, Int b, bool symmetric) {

  Int n = (b + (a % b)) % b;

  if (symmetric) {
    if (0 <= n && n <= b / 2) {
      return n;
    }

    return n - b;
  }

  return n;
}

Int randomGf(Int p, bool symmetric) {
  std::random_device dev;

  std::mt19937 rng(dev());

  std::uniform_int_distribution<std::mt19937::result_type> dist(
      std::numeric_limits<long long>::min(),
      std::numeric_limits<long long>::max());

  return mod((long long)dist(rng), p, symmetric);
}

Int inverseGf(Int a, Int b, bool symmetric) {
  Int t, nt, r, nr, q, tmp;

  if (b < 0)
    b = -b;
  if (a < 0)
    a = b - (-a % b);

  t = 0;
  nt = 1;
  r = b;
  nr = a % b;

  while (nr != 0) {
    q = r / nr;
    tmp = nt;
    nt = t - q * nt;
    t = tmp;
    tmp = nr;
    nr = r - q * nr;
    r = tmp;
  }

  if (r > 1) {
    printf("%s have no inverse mod %s\n", a.to_string().c_str(),
           b.to_string().c_str());
    exit(1);
  }

  if (t < 0)
    t += b;

  return mod(t, b, symmetric);
}

Int quoGf(Int s, Int t, Int p, bool symmetric) {
  return mod((s * inverseGf(t, p, symmetric)), p, symmetric);
}

Expr gf(Expr u, Expr x, Int s, bool symmetric) {
  Expr k = algebraicExpand(u);

  if (k.kind() == Kind::Fail || k.kind() == Kind::Undefined) {
    return k;
  }

  if (k.kind() == Kind::MinusInfinity || k.kind() == Kind::Infinity) {
    return undefined();
  }

  if (k.kind() == Kind::Integer) {
    Int p = k.value();

    return integer(mod(p, s, symmetric));
  }

  if (k.kind() == Kind::Symbol) {
    return k;
  }

  if (k.kind() == Kind::Fraction) {
    assert(k[0].kind() == Kind::Integer,
           "numerator of a fraction needs to be a integer");

    assert(k[1].kind() == Kind::Integer,
           "denominator of a fraction needs to be a integer");

    Int n = k[0].value();
    Int d = k[1].value();

    return mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric);
  }

  if (k.kind() == Kind::Derivative) {
    Expr p = Expr(Kind::Derivative, {gf(k[0], x, s, symmetric), k[1]});

    return p;
  }

  if (k.kind() == Kind::Integral) {
    Expr p = Expr(Kind::Integral, {gf(k[0], x, s, symmetric), k[1]});

    return p;
  }

  if (k.kind() == Kind::Factorial) {
    if (k[0].kind() == Kind::Integer) {
      Expr f = reduceAST(k);
      return gf(f, x, s, symmetric);
    }

    return k;
  }

  if (k.kind() == Kind::Division) {
    Expr p = div(gf(k[0], x, s, symmetric), gf(k[1], x, s, symmetric));

    return gf(reduceAST(p), x, s, symmetric);
  }

  if (k.kind() == Kind::Power) {
    Expr p = power(gf(k[0], x, s, symmetric), k[1]);

    return p;
  }

  if (k.kind() == Kind::FunctionCall) {
    return k;
  }

  if (k.kind() == Kind::Multiplication) {
    Expr p = Expr(Kind::Multiplication);

    for (long i = 0; i < k.size(); i++) {
      p.insert(gf(k[i], x, s, symmetric));
    }

    return p;
  }

  if (k.kind() == Kind::Addition || k.kind() == Kind::Subtraction) {
    Expr p = Expr(k.kind());

    Expr d = degree(k, x);

    for (Int n = 0; n <= d.value(); n++) {
      Expr c = coeff(k, x, n);
      p.insert(gf(c, x, s, symmetric) * power(x, n));
    }

    return reduceAST(p);
  }

  return k;
}

Expr groundGf(Expr u, Int s, bool symmetric) {
  Expr k = u;

  if (k.kind() == Kind::Fail || k.kind() == Kind::Undefined) {
    return k;
  }

  if (k.kind() == Kind::MinusInfinity || k.kind() == Kind::Infinity) {
    return undefined();
  }

  if (k.kind() == Kind::Integer) {
    return integer(mod(k.value(), s, symmetric));
  }

  if (k.kind() == Kind::Symbol) {
    return k;
  }

  if (k.kind() == Kind::Fraction) {
    assert(k[0].kind() == Kind::Integer,
           "numerator of a fraction needs to be a integer");

    assert(k[1].kind() == Kind::Integer,
           "denominator of a fraction needs to be a integer");

    Int n = k[0].value();
    Int d = k[1].value();

    return mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric);
  }

  if (k.kind() == Kind::Derivative) {
    Expr p = Expr(Kind::Derivative, {groundGf(k[0], s, symmetric), k[1]});

    return p;
  }

  if (k.kind() == Kind::Integral) {
    return Expr(Kind::Integral, {groundGf(k[0], s, symmetric), k[1]});
  }

  if (k.kind() == Kind::Factorial) {
    if (k[0].kind() == Kind::Integer) {
      return groundGf(fact(k.value()), s, symmetric);
    }

    return k;
  }

  if (k.kind() == Kind::Division) {
    Expr p = div(groundGf(k[0], s, symmetric), groundGf(k[1], s, symmetric));
    Expr t = reduceAST(p);

    return groundGf(t, s, symmetric);
  }

  if (k.kind() == Kind::Power) {
    return power(groundGf(k[0], s, symmetric), k[1]);
  }

  if (k.kind() == Kind::FunctionCall) {
    return k;
  }

  if (k.kind() == Kind::Multiplication) {
    Expr p = Expr(Kind::Multiplication);

    for (long i = 0; i < k.size(); i++) {
      p.insert(groundGf(k[i], s, symmetric));
    }

    return reduceAST(p);
  }

  if (k.kind() == Kind::Addition || k.kind() == Kind::Subtraction) {
    Expr p = Expr(k.kind());

    for (long i = 0; i < k.size(); i++) {
      p.insert(groundGf(k[i], s, symmetric));
    }

    Expr r = reduceAST(p);

    return r;
  }

  return k;
}

Expr gf(Expr u, Int s, bool symmetric) {
  Expr k = algebraicExpand(u);

  if (k.kind() == Kind::Fail || k.kind() == Kind::Undefined) {
    return k;
  }

  if (k.kind() == Kind::MinusInfinity || k.kind() == Kind::Infinity) {
    return undefined();
  }

  if (k.kind() == Kind::Integer) {
    Int p = k.value();

    return integer(mod(p, s, symmetric));
  }

  if (k.kind() == Kind::Symbol) {
    return k;
  }

  if (k.kind() == Kind::Fraction) {
    assert(k[0].kind() == Kind::Integer,
           "numerator of a fraction needs to be a integer");

    assert(k[1].kind() == Kind::Integer,
           "denominator of a fraction needs to be a integer");

    Int n = k[0].value();
    Int d = k[1].value();

    return integer(
        mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric));
  }

  if (k.kind() == Kind::Derivative) {
    Expr p = Expr(Kind::Derivative, {gf(k[0], s, symmetric), k[1]});

    return p;
  }

  if (k.kind() == Kind::Integral) {
    Expr p = Expr(Kind::Integral, {gf(k[0], s, symmetric), k[1]});

    return p;
  }

  if (k.kind() == Kind::Factorial) {
    if (k[0].kind() == Kind::Integer) {
      Expr f = reduceAST(k);

      Expr p = gf(f, s, symmetric);

      return p;
    }

    return k;
  }

  if (k.kind() == Kind::Division) {
    Expr p = div(gf(k[0], s, symmetric), gf(k[1], s, symmetric));
    Expr t = reduceAST(p);
    Expr r = gf(t, s, symmetric);

    return r;
  }

  if (k.kind() == Kind::Power) {
    Expr p = power(gf(k[0], s, symmetric), k[1]);

    return p;
  }

  if (k.kind() == Kind::FunctionCall) {
    return k;
  }

  if (k.kind() == Kind::Multiplication) {
    Expr p = Expr(Kind::Multiplication);

    for (long i = 0; i < k.size(); i++) {
      p.insert(gf(k[i], s, symmetric));
    }

    return p;
  }

  if (k.kind() == Kind::Addition || k.kind() == Kind::Subtraction) {
    Expr p = Expr(k.kind());

    for (long i = 0; i < k.size(); i++) {
      p.insert(gf(k[i], s, symmetric));
    }

    Expr r = reduceAST(p);

    return r;
  }

  return k;
}

Expr divPolyGf(Expr a, Expr b, Expr x, Int p, bool symmetric) {
  Expr da = degree(a, x);
  Expr db = degree(b, x);

  assert(da.kind() == Kind::Integer,
         "degree of polynomial should be an integer\n");
  assert(db.kind() == Kind::Integer,
         "degree of polynomial should be an integer\n");

  if (da.value() < db.value()) {
    return list({0, a});
  }

  long long k, j, lb, d;
  Int s, e;

  Expr dq, dr, q, r;
  Expr t1, t2, t3, ex;

  std::vector<Expr> A = std::vector<Expr>(da.value().longValue() + 1);
  std::vector<Expr> B = std::vector<Expr>(db.value().longValue() + 1);

  for (k = da.value().longValue(); k >= 0; k--) {
    A[k] = coeff(a, x, k);
  }

  for (k = db.value().longValue(); k >= 0; k--) {
    B[k] = coeff(b, x, k);
  }

  dq = da.value() - db.value();

  dr = db.value() - 1;

  t1 = leadCoeff(b, x);

  lb = inverseGf(t1.value(), p, symmetric).longValue();

  for (k = da.value().longValue(); k >= 0; k--) {
    t1 = A[k];

    s = max(0, k - dq.value());
    e = min(dr.value(), k);

    for (j = s.longValue(); j <= e; j++) {
      t2 = mulPoly(B[j], A[k - j + db.value().longValue()]);
      t3 = subPoly(t1, t2);
      t1 = t3;
    }

    t3 = reduceAST(t1);
    t1 = t3;

    t2 = mod(t1.value(), p, symmetric);
    t1 = t2;

    if (t1.value() < 0) {
      t2 = integer(t1.value() + p);
      t1 = t2;
    }

    if (da.value() - k <= dq.value()) {
      t3 = integer(lb);

      t2 = mulPoly(t1, t3);
      t1 = reduceAST(t2);
      t2 = integer(mod(t1.value(), p, symmetric));
      t1 = t2;
    }

    A[k] = t1;
  }

  q = add({0});
  r = add({0});

  d = 0;

  for (k = da.value().longValue() - dq.value().longValue(); k <= da.value();
       k++) {
    q.insert(A[k] * power(x, integer(d)));
    d = d + 1;
  }

  d = 0;

  for (k = 0; k <= dr.value(); k++) {
    r.insert(A[k]*power(x, integer(d)));

    d = d + 1;
  }

  t1 = gf(q, x, p, symmetric);
  t2 = gf(r, x, p, symmetric);

  return list({t1, t2});
}

Expr remPolyGf(Expr a, Expr b, Expr x, Int p, bool symmetric) {
  return divPolyGf(a, b, x, p, symmetric)[1];
}

Expr quoPolyGf(Expr a, Expr b, Expr x, Int p, bool symmetric) {
  return divPolyGf(a, b, x, p, symmetric)[0];
}

Expr monicPolyGf(Expr f, Expr x, Int p, bool symmetric) {
  if (f == 0) {
    return 0;
  }

  Expr lc = leadCoeff(f, x);

  Expr F = quoPolyGf(f, lc, x, p, symmetric);

  return list({lc, F});
}

Expr gcdPolyGf(Expr a, Expr b, Expr x, Int p, bool symmetric) {
  Expr da = degree(a, x);
  Expr db = degree(b, x);

  if (da.kind() == Kind::MinusInfinity || db.value() > da.value()) {
    return gcdPolyGf(b, a, x, p, symmetric);
  }

  Expr t;

  while (b != 0 && db.kind() != Kind::MinusInfinity && db.value() >= 0) {
    t = a;
    a = b;
    b = remPolyGf(t, b, x, p, symmetric);

    db = degree(b, x);
  }

  b = monicPolyGf(a, x, p, symmetric);

  return b[1];
}

Expr addPolyGf(Expr f, Expr g, Expr x, Int p, bool symmetric) {
  Expr u = addPoly(f, g);

  return gf(u, x, p, symmetric);
}

Expr subPolyGf(Expr f, Expr g, Expr x, Int p, bool symmetric) {
  Expr t, u;

  u = subPoly(f, g);

  t = gf(u, x, p, symmetric);

  return t;
}

Expr mulPolyGf(Expr f, Expr g, Expr x, Int p, bool symmetric) {
  Expr t, u;

  u = mulPoly(f, g);

  t = gf(u, x, p, symmetric);

  return t;
}

Expr powModPolyGf(Expr f, Expr g, Expr x, Int n, Int p, bool symmetric) {
  if (n == 0)
    return 1;

  Expr a = f;
  Expr b = 1;
  Expr t = 0;

  while (n > 1) {
    if (n % 2 == 0) {
      t = mulPolyGf(a, a, x, p, symmetric);
      a = remPolyGf(t, g, x, p, symmetric);

      n = n / 2;
    } else {
      t = mulPolyGf(a, b, x, p, symmetric);
      b = remPolyGf(t, g, x, p, symmetric);
      t = mulPolyGf(a, a, x, p, symmetric);
      a = remPolyGf(t, g, x, p, symmetric);

      n = (n - 1) / 2;
    }
  }

  t = mulPolyGf(a, b, x, p, symmetric);

  return remPolyGf(t, g, x, p, symmetric);
}

Expr randPolyGf(Int d, Expr x, Int p, bool symmetric) {
  Int k = 0;

  if (d == 0) {
    return randomGf(p, symmetric);
  }

  if (d == 1) {
    k = randomGf(p, symmetric);

    if (k == 0)
      return x;

    return x + k;
  }

  Expr r = Expr(Kind::Addition, {power(x, d)});

  for (Int i = d - 1; i >= 2; i--) {
    k = randomGf(p, symmetric);

    if (k != 0) {
      r = r + k * power(x, i);
    }
  }

  k = randomGf(p, symmetric);

  if (k != 0) {
    r = r + k * x;
  }

  k = randomGf(p, symmetric);

  if (k != 0) {
    r = r + k;
  }

  return r;
}

Expr extendedEuclidGf(Expr f, Expr g, Expr x, Int p, bool sym) {
  if (f == 0 || g == 0) {
    return list({1, 0, 0});
  }

  Expr t, s, i, lc, k1, t0, t3, s0, s1, Q, R, T;

  Expr t1 = monicPolyGf(f, x, p, sym);
  Expr t2 = monicPolyGf(g, x, p, sym);

  Expr p0 = t1[0];
  Expr r0 = t1[1];

  Expr p1 = t2[0];
  Expr r1 = t2[1];

  if (f == 0) {
    return list({0, inverseGf(p1.value(), p, sym), r1});
  }

  if (g == 0) {
    return list({inverseGf(p0.value(), p, sym), 0, r0});
  }

  s0 = inverseGf(p0.value(), p, sym);

  s1 = 0;
  t0 = 0;

  t1 = inverseGf(p1.value(), p, sym);

  while (true) {
    T = divPolyGf(r0, r1, x, p, sym);

    Q = T[0];
    R = T[1];

    if (R == 0) {
      break;
    }

    T = monicPolyGf(R, x, p, sym);

    r0 = r1;

    lc = T[0];
    r1 = T[1];

    i = inverseGf(lc.value(), p, sym);

    k1 = mulPolyGf(s1, Q, x, p, sym);
    s = subPolyGf(s0, k1, x, p, sym);

    k1 = mulPolyGf(t1, Q, x, p, sym);
    t = subPolyGf(t0, k1, x, p, sym);

    s0 = s1;
    t0 = t1;

    s1 = mulPolyGf(s, i, x, p, sym);
    t1 = mulPolyGf(t, i, x, p, sym);
  }

  return list({r1, s1, t1});
}

} // namespace galoisField
