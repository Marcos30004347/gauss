#include "SquareFree.hpp"

#include "Core/AST/AST.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/GaloisField/GaloisField.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace galoisField;
using namespace simplification;

namespace factorization {

Expr squareFreeFactorization(Expr ax, Expr x) {
  Expr ox, bx, cx, wx, yx, zx, qx, tx;

  long i = 1;

  ox = 1;

  bx = derivate(ax, x);
  cx = gcdGPE(ax, bx, x);
  wx = quotientGPE(ax, cx, x);

  while (cx != 1) {
    yx = gcdGPE(wx, cx, x);
    zx = quotientGPE(wx, yx, x);

    ox = mul({ox, power(zx, integer(i))});

    i = i + 1;

    wx = yx;

    qx = quotientGPE(cx, yx, x);

    cx = qx;
  }

  ox = mul({ox, power(wx, integer(i))});

  tx = reduceAST(ox);

  return tx;
}

Expr squareFreeFactorization2(Expr ax, Expr x) {
  Expr ox, bx, cx, wx, yx, kx, zx, gx, tx, rx, ux;

  unsigned int i = 1;

  ox = integer(1);

  bx = derivate(ax, x);
  cx = gcdGPE(ax, bx, x);

  if (cx == 1) {
    wx = ax;
  } else {
    wx = quotientGPE(ax, cx, x);
    yx = quotientGPE(bx, cx, x);

    kx = derivate(wx, x);
    zx = subPoly(yx, kx);

    while (zx != 0) {
      gx = gcdGPE(wx, zx, x);

      ox = mul({ox, power(gx, integer(i))});

      i = i + 1;

      tx = quotientGPE(wx, gx, x);

      wx = tx;

      yx = quotientGPE(zx, gx, x);

      rx = derivate(wx, x);

      zx = subPoly(yx, rx);
    }
  }

  ox = mul({ox, power(wx, integer(i))});

  ux = reduceAST(ox);

  return ux;
}

Expr squareFreeFactorizationFiniteField(Expr ax, Expr x, Int p,
                                        bool symmetric) {
  unsigned int i = 1;

  Expr ox = 1;
  Expr ux = derivate(ax, x);
  Expr bx = gf(ux, p, symmetric);

  if (bx != 0) {
    Expr cx = gcdPolyGf(ax, bx, x, p, symmetric);
    Expr wx = quoPolyGf(ax, cx, x, p, symmetric);

    while (wx != 1) {
      Expr yx = gcdPolyGf(wx, cx, x, p, symmetric);
      Expr zx = quoPolyGf(wx, yx, x, p, symmetric);

      ox = mul({ox, power(zx, integer(i))});

      i = i + 1;

      wx = yx;

      Expr kx = quoPolyGf(cx, yx, x, p, symmetric);

      cx = kx;
    }


		if (cx != 1) {
      Expr kx = add({});
      Expr deg = degree(cx, x);


      for (Int i = 0; i <= deg.value(); i++) {
        kx.insert(coeff(cx, x, i) * power(x, i / p));
      }

      cx = reduceAST(kx);

      Expr sx = squareFreeFactorizationFiniteField(cx, x, p);

      cx = sx;

      ox = mul({ox, power(cx, p)});
    }

  } else {
    Expr deg = degree(ax, x);
    Expr kx = add({});

    for (Int i = 0; i <= deg.value(); i++) {
      kx.insert(coeff(ax, x, i) * power(x, i / p));
    }

    ax = reduceAST(kx);

    Expr sx = squareFreeFactorizationFiniteField(ax, x, p);

    ox = power(sx, p);
  }

  return reduceAST(ox);
}

bool isSquareFreeInZp(Expr f, Expr x, long p, bool symmetric) {
  bool r = false;

  Expr lc, t, k, v, g;

  if (f == 0) {
    return true;
  }

  lc = leadCoeff(f, x);

  v = quoPolyGf(f, lc, x, p, symmetric);

  k = derivate(v, x);

  t = gf(k, p, symmetric);

  g = gcdPolyGf(v, t, x, p, symmetric);

  r = g == 1;

  return r;
}

Expr squareFreePart(Expr f, Expr L, Expr K) {
  Expr g, u, v, s;

  long i;

  g = f;

  for (i = 0; i < L.size(); i++) {
    u = reduceAST(derivate(f, L[i]));
    v = mvPolyGCD(g, u, L, K);
    g = v;
  }

  s = recQuotient(f, g, L, K);
  g = pp(s, L, K);

  Expr R = list({});

  for (i = 0; i < L.size(); i++) {
    if (!g.freeOf(L[i])) {
      R.insert(L[i]);
    }
  }

  return list({g, R});
}

bool isSquareFree(ast::Expr f, ast::Expr x, ast::Expr K) {
  long e = 1;

  Expr k, n, g;

  if (f == 0) {
    return true;
  }

  k = derivate(f, x);
  g = gcdGPE(f, k, x);
  n = degree(g, x);

  e = n.value().longValue();

  return !e;
}

} // namespace factorization
