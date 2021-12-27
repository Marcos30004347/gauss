#include "Hensel.hpp"

#include "Core/AST/AST.hpp"
#include "Core/AST/Integer.hpp"
#include "Core/GaloisField/GaloisField.hpp"
#include "Core/Polynomial/Polynomial.hpp"

#include <bits/types/FILE.h>
#include <cmath>

using namespace ast;
using namespace algebra;
using namespace galoisField;
using namespace polynomial;

namespace factorization {

Expr leadCoeffReplace(Expr ux, Expr x, Expr c) {
  Expr lc = leadCoeff(ux, x);

  Expr de = degree(ux, x);

  Expr px = mul({lc, power(x, de)});

  Expr rx = sub({ux, px});

  Expr kx = algebraicExpand(rx);

  Expr ox = add({kx, mul({c, power(x, de)})});

  Expr zx = algebraicExpand(ox);

  return zx;
}

Expr leadCoeffReplacePolyExpr(Expr ux, Expr c) {
  ux[ux.size() - 1][0] = c;
  return ux;
}

Expr normalize(Expr ux, Expr x) {
  Expr lc = leadCoeff(ux, x);
  Expr px = power(lc, integer(-1)) * ux;

  return algebraicExpand(px);
}

Expr normalizePolyExpr(Expr ux, Expr L, Expr K) {
  Expr lc = raisePolyExpr(leadCoeffPolyExpr(ux), 0, L[0]);
  return quoPolyExpr(ux, lc, L, K);
}

Expr henselSep(Expr f, Expr g, Expr h, Expr s, Expr t, Expr x, Int m,
               bool symmetric) {
  Expr e, q, r, G, H, b, c, d, S, T, t1, t2, t3, t4;
  Expr Z = Expr("Z");
  Expr L = list({x});

  Int m2 = m * m;

  t1 = mulPoly(g, h);

  e = subPoly(f, t1);
  e = gf(e, m2, symmetric);

  t2 = mulPoly(s, e);
  t1 = recPolyDiv(t2, h, L, Z);

  q = gf(t1[0], m2, symmetric);
  r = gf(t1[1], m2, symmetric);

  t2 = mulPoly(t, e);
  t3 = mulPoly(q, g);

  t4 = addPoly(t2, t3);

  G = addPoly(g, t4);
  G = gf(G, m2, symmetric);

  H = addPoly(h, r);
  H = gf(H, m2, symmetric);

  t2 = mulPoly(s, G);
  t3 = mulPoly(t, H);
  t4 = addPoly(t2, t3);

  b = subPoly(t4, 1);
  b = gf(b, m2, symmetric);

  t2 = mulPoly(s, b);
  t3 = recPolyDiv(t2, H, L, Z);

  c = gf(t3[0], m2, symmetric);
  d = gf(t3[1], m2, symmetric);

  S = subPoly(s, d);
  S = gf(S, m2, symmetric);

  t2 = mulPoly(t, b);
  t3 = mulPoly(c, G);
  t1 = addPoly(t2, t3);

  T = subPoly(t, t1);
  T = gf(T, m2, symmetric);

  return list({G, H, S, T});
}

// Expr henselSepPolyExpr(Expr f, Expr g, Expr h, Expr s, Expr t, Expr L, Int m,
//                        bool symmetric) {
//   assert(L.kind() == Kind::List && L.size() <= 1);

//   Expr e, q, r, G, H, b, c, d, S, T, t1, t2, t3, t4;

//   Int m2 = m * m;

//   t2 = mulPolyExprGf(g, h, m2, symmetric);
//   t3 = subPolyExpr(f, t2);

//   e = subPolyExprGf(f, t2, m2, symmetric);

//   t2 = mulPolyExprGf(s, e, m2, symmetric);
//   t1 = divPolyExprGf(t2, h, L, m2, symmetric);

//   q = t1[0];
//   r = t1[1];
//   t2 = mulPolyExpr(t, e);
//   t3 = mulPolyExpr(q, g);

//   t4 = addPolyExpr(t2, t3);

//   G = addPolyExprGf(g, t4, m2, symmetric);
//   H = addPolyExprGf(h, r, m2, symmetric);

//   t2 = mulPolyExprGf(s, G, m2, symmetric);
//   t3 = mulPolyExprGf(t, H, m2, symmetric);
//   t4 = addPolyExprGf(t2, t3, m2, symmetric);

//   b = subPolyExprGf(t4, raisePolyExpr(1, 0, L[0]), m2, symmetric);

//   t2 = mulPolyExprGf(s, b, m2, symmetric);
//   t3 = divPolyExprGf(t2, H, L, m2, symmetric);

//   c = t3[0];
//   d = t3[1];

//   S = subPolyExprGf(s, d, m2, symmetric);

//   t2 = mulPolyExprGf(t, b, m2, symmetric);

//   t3 = mulPolyExprGf(c, G, m2, symmetric);
//   t1 = addPolyExprGf(t2, t3, m2, symmetric);

//   T = subPolyExprGf(t, t1, m2, symmetric);

//   return list({G, H, S, T});
// }



Expr henselSepPolyExpr(Expr f, Expr g, Expr h, Expr s, Expr t, Expr L, Int m,
               bool symmetric) {
  Expr e, q, r, G, H, b, c, d, S, T, t1, t2, t3, t4;

  Expr Z = Expr("Z");

  Expr o = raisePolyExpr(1, 0, L[0]);

  Int m2 = m * m;

  t1 = mulPolyExpr(g, h);
  e = subPolyExpr(f, t1);

  e = gfPolyExpr(e, m2, symmetric);

  t2 = mulPolyExpr(s, e);
  t1 = divPolyExpr(t2, h, L, Z);

  q = gfPolyExpr(t1[0], m2, symmetric);
  r = gfPolyExpr(t1[1], m2, symmetric);

  t2 = mulPolyExpr(t, e);
  t3 = mulPolyExpr(q, g);

  t4 = addPolyExpr(t2, t3);

  G = addPolyExpr(g, t4);
  G = gfPolyExpr(G, m2, symmetric);

  H = addPolyExpr(h, r);
  H = gfPolyExpr(H, m2, symmetric);

  t2 = mulPolyExpr(s, G);
  t3 = mulPolyExpr(t, H);
  t4 = addPolyExpr(t2, t3);

  b = subPolyExpr(t4, o);
  b = gfPolyExpr(b, m2, symmetric);

  t2 = mulPolyExpr(s, b);
  t3 = divPolyExpr(t2, H, L, Z);

  c = gfPolyExpr(t3[0], m2, symmetric);
  d = gfPolyExpr(t3[1], m2, symmetric);

  S = subPolyExpr(s, d);
  S = gfPolyExpr(S, m2, symmetric);

  t2 = mulPolyExpr(t, b);
  t3 = mulPolyExpr(c, G);
  t1 = addPolyExpr(t2, t3);

  T = subPolyExpr(t, t1);
  T = gfPolyExpr(T, m2, symmetric);

  return list({G, H, S, T});
}


Expr multifactorHenselLifting(Expr v, Expr H, Expr x, Int p, Int l,
                              bool symmetric) {
  Int i, j, r, k, d;

  Int a;

  Expr f, fi, lc, t1, g, h, s;
  Expr t, e, T, H0, H1, F0, F1, F;

  lc = leadCoeff(v, x);

  f = v;

  r = H.size();

  if (r == 1) {
    Int pl = pow(p, l);

    a = inverseGf(lc.value(), pl, symmetric);

    t1 = mod(a, pl, symmetric);

    fi = mulPolyGf(f, t1, x, pl, symmetric);

    return list({fi});
  }

  k = r / 2;
  d = std::ceil(log2(l.longValue()));
  g = lc;

  h = 1;

  for (i = 0; i < k; i++) {
    g = mulPolyGf(g, H[i], x, p, symmetric);
  }

  for (i = k; i < r; i++) {
    h = mulPolyGf(h, H[i], x, p, symmetric);
  }

  e = extendedEuclidGf(g, h, x, p, symmetric);

  s = gf(e[1], p, symmetric);
  t = gf(e[2], p, symmetric);

  for (j = 1; j <= d; j++) {
    T = henselSep(f, g, h, s, t, x, pow(p, pow(2, j - 1)), symmetric);
    g = T[0];
    h = T[1];
    s = T[2];
    t = T[3];
  }

  H0 = list({});
  H1 = list({});

  for (i = 0; i < k; i++) {
    H0.insert(H[i]);
  }

  for (i = k; i < r; i++) {
    H1.insert(H[i]);
  }

  F0 = multifactorHenselLifting(g, H0, x, p, l, symmetric);
  F1 = multifactorHenselLifting(h, H1, x, p, l, symmetric);

  F = list({});

  while (F0.size() > 0) {
    F.insert(F0[0]);
    F0.remove(0L);
  }

  while (F1.size() > 0) {
    F.insert(F1[0]);
    F1.remove(0L);
  }

  return F;
}

Expr multifactorHenselLiftingPolyExpr(Expr v, Expr H, Expr L, Int p, Int l,
                                      bool symmetric) {
  assert(L.kind() == Kind::List && L.size() == 1);

  Int i, j, r, k, d;

  Int a;

  Expr f, fi, lc, t1, g, h, s;
  Expr t, e, T, H0, H1, F0, F1, F;

  Int pl = pow(p, l);

  lc = leadCoeffPolyExpr(v);

  f = v;

  r = H.size();

  if (r == 1) {
    a = inverseGf(lc.value(), pl, symmetric);

    fi = mulPolyExprGf(f, mod(a, pl, symmetric), pl, symmetric);

    return list({fi});
  }

  k = r / 2;
  d = std::ceil(log2(l.longValue()));

  g = lc;

  h = gfPolyExpr(H[k], p, symmetric); // integer(1);

  // TODO: use pow
  for (i = 0; i < k; i++) {
    g = mulPolyExprGf(g, H[i], p, symmetric);
  }
  // TODO: use pow
  for (i = k + 1; i < r; i++) {
    h = mulPolyExprGf(h, H[i], p, symmetric);
  }

  e = extendedEuclidPolyExprGf(g, h, L, p, symmetric);

  s = e[1];
  t = e[2];

  for (j = 1; j <= d; j++) {
    T = henselSepPolyExpr(f, g, h, s, t, L, pow(p, pow(2, j - 1)), symmetric);

    g = T[0];
    h = T[1];
    s = T[2];
    t = T[3];
  }

  H0 = list({});
  H1 = list({});

  for (i = 0; i < k; i++) {
    H0.insert(H[i]);
  }

  for (i = k; i < r; i++) {
    H1.insert(H[i]);
  }

  F0 = multifactorHenselLiftingPolyExpr(g, H0, L, p, l, symmetric);
  F1 = multifactorHenselLiftingPolyExpr(h, H1, L, p, l, symmetric);

  F = list({});

  for (Int i = 0; i < F0.size(); i++) {
    F.insert(F0[i]);
  }

  for (Int i = 0; i < F1.size(); i++) {
    F.insert(F1[i]);
  }

  return F;
}

// Expr univariateHensel(Expr ax, Expr x, Expr p, Expr ux_1, Expr wx_1, Expr B,
// Expr zeta, bool symmetric)
// {

// 	Expr tmp = nullptr;
// 	Expr gam = nullptr;
// 	Expr tal = nullptr;
// 	Expr qx  = nullptr;
// 	Expr rx  = nullptr;
// 	Expr cx  = nullptr;

// 	Expr Z = symbol("Z");

// 	Expr L = list({ x });

// 	// 1. Define polynomial and its modulo p factors
// 	Expr alph = leadCoeff(ax, x);

// 	if(zeta->kind() == Kind::Undefined)
// 	{
// 		zeta = alph;
// 	}
// 	else
// 	{
// 		zeta = zeta;
// 	}

// 	tmp = mul({ zeta, ax });

// 	ax = algebraicExpand(tmp);

//

// 	// Normalization maybe wrong
// 	Expr kx = normalize(ux_1, x);
// 	Expr zx = normalize(wx_1, x);

// 	ux_1 = mulPolyGf(zeta, kx, x, p.value(), symmetric);
// 	wx_1 = mulPolyGf(alph, zx, x, p.value(), symmetric);

//
//

// // 	// printf("u[1](x) = %s\n", ux_1->toString().c_str());
// // 	// printf("w[1](x) = %s\n", wx_1->toString().c_str());

// 	// 2. Apply extended Euclidean algorithm to ux_1, wx_2 defined in Zp[x]
// 	Expr l = extendedEuclidGf(ux_1, wx_1, x, p.value(), symmetric);
// //extendedEuclideanAlgGPE_sZp(ux_1, wx_1, x, p.value());

// 	Expr sx = l[1];
// 	Expr tx = l[2];

// // 	// printf("s(x) = %s\n", sx->toString().c_str());
// // 	// printf("t(x) = %s\n", tx->toString().c_str());

//

// 	// 3. Initialization for iteration
// 	Expr ux = leadCoeffReplace(ux_1, x, zeta);
// 	Expr wx = leadCoeffReplace(wx_1, x, alph);

// 	Expr fx = sub({ ax, mul({ ux, wx }) });

// 	Expr ex = algebraicExpand(fx);

//

// 	long modulus = p.value();

// 	// 4. Iterate until either the factorization in Z[x] is obtained
// 	// or else the bound on modulus is reached
// 	while(ex->isNot(0) && modulus < 2 * B.value() * zeta.value())
// 	{
// 		// 4.1 Solve in the domain Zp[x] the polynomial equation
// 		tmp = div(ex, integer(modulus)); // MAYBE DIV IN sZp[x]

// 		cx = algebraicExpand(tmp);

//

// 		// tmp = mul({ sx, cx });

//

// 		gam = mulPolyGf(sx, cx, x, p.value(), symmetric);// sZp(tmp, x,
// p.value());

// 		//

// 		// tmp = mul({ tx, cx });

//

// 		tal = mulPolyGf(tx, cx, x, p.value(), symmetric);// sZp(tmp, x,
// p.value());

//

// 		tmp = quoPolyGf(gam, wx_1, x, p.value(), symmetric);

// 		qx = tmp[0];
// 		rx = tmp[1];

//

//

// 		gam = rx;

// 		tmp = mulPolyGf(qx, ux_1, x, p.value(), symmetric); //add({ tal,
// mul({ qx, ux_1 }) });

//

// 		tal = addPolyGf(tal, tmp, x, p.value(), symmetric); //sZp(tmp,
// x, p.value());

//

// 		// 4.2 Update factors and compute the error
// 		ux = add({ ux, mul({ tal, integer(modulus) })});
// 		wx = add({ wx, mul({ gam, integer(modulus) })});

// 		tmp = sub({ ax, mul({ ux, wx }) });

//

// 		ex = algebraicExpand(tmp);

//

// 		tmp = algebraicExpand(ux);

//

// 		ux = tmp;

// 		tmp = algebraicExpand(wx);

//

// 		wx = tmp;

// 		modulus = modulus * p.value();

//
// 	}

// 	Expr lf = nullptr;

// 	// 5. Check termination status
// 	if(ex->is(0))
// 	{
// 		// Factorization obtained - remove contents
// 		Expr d = cont(ux, x);

// 		Expr px = quotientGPE(ux, d, x);

// 		Expr k = integer(zeta.value() / d.value());

// 		Expr kx = quotientGPE(wx, k, x);

//
//

// 		// Note: a(x) <- a(x)/y would restore a(x) to its input value
// 		lf = list({ px, kx });
// 	}
// 	else
// 	{
// 		lf = new AST(Kind::Fail);
// 	}
//
//

//
//
//
//

//
//
//
//
//
//
//
//

// 	return lf;
// }
} // namespace factorization
