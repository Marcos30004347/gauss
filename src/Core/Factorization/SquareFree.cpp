#include "SquareFree.hpp"

#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/GaloisField/GaloisField.hpp"
#include "Core/Polynomial/Polynomial.hpp"
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

Expr squareFreeFactorizationPolyExpr(Expr ax, Expr L, Expr Z) {
	assert(L.size() <= 1, "only univariate poly expressions allowed");
	assert(Z.identifier() == "Z", "only the integer field allowed");

  Int i = 1;

  Expr ox = Expr(Kind::Multiplication);

  Expr bx = diffPolyExpr(ax, L[0]);

  Expr cx = gcdPolyExpr(ax, bx, L, Z);
  Expr wx = quoPolyExpr(ax, cx, L, Z);

  while (!isConstantPolyExpr(cx, 1)) {
    Expr yx = gcdPolyExpr(wx, cx, L, Z);
    Expr zx = quoPolyExpr(wx, yx, L, Z);

		if(!isConstantPolyExpr(zx, 1)) {
			ox = ox*power(zx, i);
		}

		i = i + 1;

    wx = yx;

    cx = quoPolyExpr(cx, yx, L, Z);
  }

  if(isConstantPolyExpr(wx, 1)) {
		return ox;
	}

	return ox*power(wx, i);
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

Expr squareFreeFactorizationPolyExpr2(Expr ax, Expr L, Expr Z) {
 	assert(L.size() <= 1, "only univariate poly expressions allowed");
	assert(Z.identifier() == "Z", "only the integer field allowed");

	Expr ox, bx, cx, wx, yx, kx, zx, gx, tx, rx, ux;

  Int i = 1;

  ox = Expr(Kind::Multiplication);

  bx = diffPolyExpr(ax, L[0]);
  cx = gcdPolyExpr(ax, bx, L, Z);

  if (isConstantPolyExpr(cx, 1)) {
    wx = ax;
  } else {
    wx = quoPolyExpr(ax, cx, L, Z);
    yx = quoPolyExpr(bx, cx, L, Z);

    kx = diffPolyExpr(wx, L[0]);
    zx = subPolyExpr(yx, kx);

    while (!isConstantPolyExpr(zx, 0)) {
      gx = gcdPolyExpr(wx, zx, L, Z);

			if(!isConstantPolyExpr(gx, 1)) {
				ox = ox*power(gx, i);
			}

      i = i + 1;

      tx = quoPolyExpr(wx, gx, L, Z);

      wx = tx;

      yx = quoPolyExpr(zx, gx, L, Z);

      rx = diffPolyExpr(wx, L[0]);

      zx = subPolyExpr(yx, rx);
    }
  }

	if(!isConstantPolyExpr(wx, 1)) {
		ox = ox*power(wx, i);
	}

  return ox;
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



Expr squareFreeFactorizationFiniteFieldPolyExpr(Expr ax, Expr L, Expr Z, Int p,
                                        bool symmetric) {
 	assert(L.size() <= 1, "only univariate poly expressions allowed");
	assert(Z.identifier() == "Z", "only the integer field allowed");

  Int i = 1;

  Expr ox = Expr(Kind::Multiplication);
  Expr ux = diffPolyExpr(ax, L[0]);
  Expr bx = gfPolyExpr(ux, p, symmetric);

  if (!isConstantPolyExpr(bx, 0)) {
    Expr cx = gcdPolyExprGf(ax, bx, L, p, symmetric);
    Expr wx = quoPolyExprGf(ax, cx, L, p, symmetric);

    while (!isConstantPolyExpr(wx, 1)) {
      Expr yx = gcdPolyExprGf(wx, cx, L, p, symmetric);
      Expr zx = quoPolyExprGf(wx, yx, L, p, symmetric);

			if(!isConstantPolyExpr(zx, 1)) {
				ox = ox * power(zx, i);
			}

      i = i + 1;

      wx = yx;

      Expr kx = quoPolyExprGf(cx, yx, L, p, symmetric);

      cx = kx;
    }


		if (!isConstantPolyExpr(cx, 1)) {
      Expr kx = Expr(Kind::Addition);

      for (Int i = 0; i < cx.size(); i++) {
				kx.insert(cx[i][0]*power(cx[i][1][0], cx[i][1][1].value()/p));
      }

      cx = kx;

      Expr sx = squareFreeFactorizationFiniteFieldPolyExpr(cx, L, Z, p, symmetric);

      cx = sx;

			for(Int i =0; i < cx.size(); i++) {
				if(!isConstantPolyExpr(cx[i], 1)) {
					ox = ox * power(cx[i][0], cx[i][1].value() * p);
				}
			}
    }

  } else {
		Expr kx = Expr(Kind::Addition);

		for (Int i = 0; i < ax.size(); i++) {
			kx.insert(ax[i][0]*power(ax[i][1][0], ax[i][1][1].value()/p));
		}

    ax = kx;

    Expr sx = squareFreeFactorizationFiniteFieldPolyExpr(ax, L, Z, p, symmetric);

		for(Int i =0; i < sx.size(); i++) {
			if(!isConstantPolyExpr(sx[i], 1)) {
				ox = ox * power(sx[i][0], sx[i][1].value() * p);
			}
		}
  }

  return ox;
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

    //printf("u = %s\n", algebraicExpand(u).toString().c_str());
    //printf("g = %s\n", algebraicExpand(g).toString().c_str());
    //printf("--> L = %s\n", L.toString().c_str());

    g = gcdPoly(g, u, L, K);

		//printf("--> G = %s\n", g.toString().c_str());
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

Expr squareFreePartPolyExpr(Expr f, Expr L, Expr K) {
  Expr g, u, v, s;

  long i;

  g = f;
  for (i = 0; i < L.size(); i++) {
		u = diffPolyExpr(f, L[i]);

		//printf("--> u = %s\n", algebraicExpand(u).toString().c_str());
    //printf("--> g = %s\n", algebraicExpand(g).toString().c_str());
    //printf("--> L = %s\n", L.toString().c_str());

		g = gcdPolyExpr(g, u, L, K);

		//printf("--> G = %s\n", g.toString().c_str());
  }
  s = quoPolyExpr(f, g, L, K);
  g = ppPolyExpr(s, L, K);
	// TODO: remove following lines and
	// return just the square free part

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

bool isSquareFreePolyExpr(Expr f, Expr L, Expr K) {
  if (isZeroPolyExpr(f)) return true;

  Expr k = diffPolyExpr(f, L[0]);
  Expr g = gcdPolyExpr(f, k, L, K);
  Expr n = degreePolyExpr(g);

  return n == 0;
}



} // namespace factorization
