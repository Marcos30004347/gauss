#include "SquareFree.hpp"

#include "Core/AST/AST3.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/GaloisField/GaloisField.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include <cstddef>

using namespace alg;
using namespace calculus;
using namespace polynomial;
using namespace galoisField;

namespace factorization {

expr squareFreeFactorization(expr ax, expr x) {
  expr ox, bx, cx, wx, yx, zx, qx, tx;

  long i = 1;

  ox = 1;

	bx = derivate(ax, x);
	// TODO: use gcdPoly instead of GPE
	cx = gcdGPE(ax, bx, x);
	// TODO: use quoPoly instead of GPE
  wx = quotientGPE(ax, cx, x);

  while (cx != 1) {
    yx = gcdGPE(wx, cx, x);
    zx = quotientGPE(wx, yx, x);

    ox = (ox * pow(zx, integer(i)));

    i = i + 1;

    wx = yx;

    qx = quotientGPE(cx, yx, x);

    cx = qx;
  }

  ox = (ox * pow(wx, integer(i)));

  tx = reduce(ox);

  return tx;
}

expr squareFreeFactorizationPolyExpr(expr ax, expr L, expr Z) {
	assert(L.size() <= 1);
	assert(Z.identifier() == "Z");

  Int i = 1;

  expr ox = expr(kind::MUL);

  expr bx = diffPolyExpr(ax, L[0]);

  expr cx = gcdPolyExpr(ax, bx, L, Z);
  expr wx = quoPolyExpr(ax, cx, L, Z);

  while (!isConstantPolyExpr(cx, 1)) {
    expr yx = gcdPolyExpr(wx, cx, L, Z);
    expr zx = quoPolyExpr(wx, yx, L, Z);

		if(!isConstantPolyExpr(zx, 1)) {
			ox = ox*pow(zx, i);
		}

		i = i + 1;

    wx = yx;

    cx = quoPolyExpr(cx, yx, L, Z);
  }

  if(isConstantPolyExpr(wx, 1)) {
		return ox;
	}

	return ox*pow(wx, i);
}


expr squareFreeFactorization2(expr ax, expr x) {
  expr ox, bx, cx, wx, yx, kx, zx, gx, tx, rx, ux;

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

      ox = (ox * pow(gx, integer(i)));

      i = i + 1;

      tx = quotientGPE(wx, gx, x);

      wx = tx;

      yx = quotientGPE(zx, gx, x);

      rx = derivate(wx, x);

      zx = subPoly(yx, rx);
    }
  }

  ox = (ox * pow(wx, integer(i)));

  ux = reduce(ox);

  return ux;
}

expr squareFreeFactorizationPolyExpr2(expr ax, expr L, expr Z) {
 	assert(L.size() <= 1);
	assert(Z.identifier() == "Z");

	expr ox, bx, cx, wx, yx, kx, zx, gx, tx, rx, ux;

  Int i = 1;

  ox = expr(kind::MUL);

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
				ox = ox*pow(gx, i);
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
		ox = ox*pow(wx, i);
	}

  return ox;
}


expr squareFreeFactorizationFiniteField(expr ax, expr x, Int p,
                                        bool symmetric) {
  unsigned int i = 1;

  expr ox = 1;
  expr ux = derivate(ax, x);
  expr bx = gf(ux, p, symmetric);

  if (bx != 0) {
    expr cx = gcdPolyGf(ax, bx, x, p, symmetric);
    expr wx = quoPolyGf(ax, cx, x, p, symmetric);

    while (wx != 1) {
      expr yx = gcdPolyGf(wx, cx, x, p, symmetric);
      expr zx = quoPolyGf(wx, yx, x, p, symmetric);

      ox = (ox * pow(zx, integer(i)));

      i = i + 1;

      wx = yx;

      expr kx = quoPolyGf(cx, yx, x, p, symmetric);

      cx = kx;
    }


		if (cx != 1) {
      expr kx = create(kind::ADD);
      expr deg = degree(cx, x);


      for (Int i = 0; i <= deg.value(); i++) {
        kx.insert(coeff(cx, x, i) * pow(x, i / p));
      }

      cx = reduce(kx);

      expr sx = squareFreeFactorizationFiniteField(cx, x, p);

      cx = sx;

      ox = (ox * pow(cx, p));
    }

  } else {
    expr deg = degree(ax, x);
    expr kx = create(kind::ADD);

    for (Int i = 0; i <= deg.value(); i++) {
      kx.insert(coeff(ax, x, i) * pow(x, i / p));
    }

    ax = reduce(kx);

    expr sx = squareFreeFactorizationFiniteField(ax, x, p);

    ox = pow(sx, p);
  }

  return reduce(ox);
}



expr squareFreeFactorizationFiniteFieldPolyExpr(expr ax, expr L, expr Z, Int p,
                                        bool symmetric) {
 	assert(L.size() <= 1);
	assert(Z.identifier() == "Z");

  Int i = 1;

  expr ox = expr(kind::MUL);
  expr ux = diffPolyExpr(ax, L[0]);
  expr bx = gfPolyExpr(ux, p, symmetric);

  if (!isConstantPolyExpr(bx, 0)) {
    expr cx = gcdPolyExprGf(ax, bx, L, p, symmetric);
    expr wx = quoPolyExprGf(ax, cx, L, p, symmetric);

    while (!isConstantPolyExpr(wx, 1)) {
      expr yx = gcdPolyExprGf(wx, cx, L, p, symmetric);
      expr zx = quoPolyExprGf(wx, yx, L, p, symmetric);

			if(!isConstantPolyExpr(zx, 1)) {
				ox = ox * pow(zx, i);
			}

      i = i + 1;

      wx = yx;

      expr kx = quoPolyExprGf(cx, yx, L, p, symmetric);

      cx = kx;
    }


		if (!isConstantPolyExpr(cx, 1)) {
      expr kx = expr(kind::ADD);

      for (Int i = 0; i < cx.size(); i++) {
				kx.insert(cx[i][0]*pow(cx[i][1][0], cx[i][1][1].value()/p));
      }

      cx = kx;

      expr sx = squareFreeFactorizationFiniteFieldPolyExpr(cx, L, Z, p, symmetric);

      cx = sx;

			for(Int i =0; i < cx.size(); i++) {
				if(!isConstantPolyExpr(cx[i], 1)) {
					ox = ox * pow(cx[i][0], cx[i][1].value() * p);
				}
			}
    }

  } else {
		expr kx = expr(kind::ADD);

		for (Int i = 0; i < ax.size(); i++) {
			kx.insert(ax[i][0]*pow(ax[i][1][0], ax[i][1][1].value()/p));
		}

    ax = kx;

    expr sx = squareFreeFactorizationFiniteFieldPolyExpr(ax, L, Z, p, symmetric);

		for(Int i =0; i < sx.size(); i++) {
			if(!isConstantPolyExpr(sx[i], 1)) {
				ox = ox * pow(sx[i][0], sx[i][1].value() * p);
			}
		}
  }

  return ox;
}


bool isSquareFreeInZp(expr f, expr x, long p, bool symmetric) {
  bool r = false;

  expr lc, t, k, v, g;

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

expr squareFreePart(expr f, expr L, expr K) {
  expr g, u, v, s;

  g = f;

  for (size_t i = 0; i < L.size(); i++) {
    u = reduce(derivate(f, L[i]));

    //printf("u = %s\n", algebraicExpand(u).toString().c_str());
    //printf("g = %s\n", algebraicExpand(g).toString().c_str());
    //printf("--> L = %s\n", L.toString().c_str());

    g = gcdPoly(g, u, L, K);

		//printf("--> G = %s\n", g.toString().c_str());
	}

  s = recQuotient(f, g, L, K);
  g = pp(s, L, K);

  expr R = list({});

  for (size_t i = 0; i < L.size(); i++) {
    if (!g.freeOf(L[i])) {
      R.insert(L[i]);
    }
  }

  return list({g, R});
}

expr squareFreePartPolyExpr(expr f, expr L, expr K) {
  expr g, u, v, s;

  g = f;

  for (size_t i = 0; i < L.size(); i++) {
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

	expr R = list({});

  for (size_t i = 0; i < L.size(); i++) {
    if (!g.freeOf(L[i])) {
      R.insert(L[i]);
    }
  }

  return list({g, R});
}


bool isSquareFree(expr f, expr x, expr K) {
  long e = 1;

  expr k, n, g;

  if (f == 0) {
    return true;
  }

  k = derivate(f, x);
  g = gcdGPE(f, k, x);
  n = degree(g, x);

  e = n.value().longValue();

  return !e;
}

bool isSquareFreePolyExpr(expr f, expr L, expr K) {
  if (isZeroPolyExpr(f)) return true;

  expr k = diffPolyExpr(f, L[0]);
  expr g = gcdPolyExpr(f, k, L, K);
  expr n = degreePolyExpr(g);

  return n == 0;
}



} // namespace factorization
