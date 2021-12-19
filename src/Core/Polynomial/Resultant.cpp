#include "Resultant.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"

// TODO: optimize all the methods in this file when new polynomial multiplication methods gets added

using namespace ast;
using namespace algebra;
using namespace simplification;

namespace polynomial {

// Expr univariateResultant(Expr ux, Expr vx, Expr x) {
//   Expr m = degree(ux, x);
//   Expr n = degree(vx, x);

//   if (n == 0) {

//     return power(vx, m);
//   }

//   Expr r = remainderGPE(ux, vx, x);

//   if (r == 0) {

//     return integer(0);
//   }

//   Expr s = degree(r, x);
//   Expr l = coeff(vx, x, n);

//   Expr k = mul({power(integer(-1), mul({m, n})), power(l, sub({m, s})),
//                 univariateResultant(vx, r, x)});

//   Expr e = algebraicExpand(k);

//   return e;
// }

// Expr multivariateResultant(Expr u, Expr v, Expr L, Expr K) {
//   Expr x = L[0];

//   Expr m = degree(u, x);
//   Expr n = degree(v, x);

//   // if(m->isLessThan(n))
//   if (m.value() < n.value()) {
//     Expr k = mul(
//         {power(integer(-1), mul({m, n})), multivariateResultant(v, u, L, K)});

//     Expr t = algebraicExpand(k);

//     return t;
//   }

//   if (n == 0) {
//     Expr k = power(v, m);

//     return k;
//   }

//   Expr delta = add({m, mul({integer(-1), n}), integer(1)});

//   Expr r = pseudoRemainder(u, v, x);

//   if (r == 0) {

//     return integer(0);
//   }

//   Expr s = degree(r, x);

//   Expr e =
//       mul({power(integer(-1), mul({m, n})), multivariateResultant(v, r, L, K)});

//   Expr w = algebraicExpand(e);

//   Expr l = leadCoeff(v, x);

//   Expr k = add({mul({delta, n}), mul({integer(-1), m}), s});

//   Expr z = power(l, k);

//   Expr f = algebraicExpand(z);

//   Expr g = recQuotient(w, f, L, K);

//   return g;
// }

Expr polynomialResultantRec(Expr u, Expr v, Expr L, Expr K, Expr i,
                              Expr delta_prev, Expr gamma_prev) {
  assert(u != 0, "Polynomial should be non-zero");
  assert(v != 0, "Polynomial should be non-zero");

  Expr x = L[0];
  Expr m = degree(u, x);
  Expr n = degree(v, x);

  if (m.value() < n.value()) {
    Expr e =
        mul({power(integer(-1), mul({m, n})),
             polynomialResultantRec(v, u, L, K, i, delta_prev, gamma_prev)});

    Expr k = reduceAST(e);

    return k;
  }

  if (n == 0) {
    Expr e = power(v, m);

    Expr k = reduceAST(e);

    return k;
  }

  Expr r = pseudoRemainder(u, v, x);

  if (r == 0) {

    return integer(0);
  }

  Expr delta = add({m, mul({integer(-1), n}), integer(1)});

  Expr R = rest(L);

  Expr gamma = undefined();
  Expr beta = undefined();

  if (i == 1) {
    gamma = integer(-1);

    Expr tmp = power(integer(-1), delta);

    beta = reduceAST(tmp);

  } else {
    Expr f = leadCoeff(u, x);

    Expr tmp1 = power(mul({integer(-1), f}), sub({delta_prev, integer(1)}));
    Expr tmp2 = power(gamma_prev, sub({delta_prev, integer(2)}));

    Expr tmp3 = algebraicExpand(tmp1);
    Expr tmp4 = algebraicExpand(tmp2);

    gamma = recQuotient(tmp3, tmp4, R, K);

    Expr tmp5 = mul({integer(-1), f, power(gamma, sub({delta, integer(1)}))});

    beta = algebraicExpand(tmp5);
  }

  Expr t = recQuotient(r, beta, L, K);

  r = t;

  Expr tmp1 = add({i, integer(1)});
  Expr tmp2 = reduceAST(tmp1);

  Expr tmp3 = mul({power(integer(-1), mul({m, n})), power(beta, n),
                   polynomialResultantRec(v, r, L, K, tmp2, delta, gamma)});

  Expr w = algebraicExpand(tmp3);

  Expr l = coeff(v, x, n);

  Expr s = degree(r, x);

  Expr k = add({mul({delta, n}), mul({integer(-1), m}), s});

  Expr tmp4 = power(l, k);

  Expr f = algebraicExpand(tmp4);

  Expr o = recQuotient(w, f, L, K);

  return o;
}

Expr polynomialResultant(Expr u, Expr v, Expr L, Expr K) {
  Expr x = L[0];

  Expr m = degree(u, x);
  Expr n = degree(v, x);

  Expr cont_u = cont(u, L, K);
  Expr pp_u = recQuotient(u, cont_u, L, K);
  Expr cont_v = cont(v, L, K);
  Expr pp_v = recQuotient(v, cont_v, L, K);

  Expr i = integer(1);
  Expr delta = integer(0);
  Expr g = integer(0);

  Expr s = polynomialResultantRec(pp_u, pp_v, L, K, i, delta, g);

  return reduceAST(power(cont_u, n) * power(cont_v, m) * s);
}

Expr polyRemSeqRec(Expr Gi2, Expr Gi1, Expr L, Expr hi2, Expr K) {
  Expr Gi, hi1, d, t1, t2, t3, t4, t5, t6, nk, cnt, ppk, r, x;

  x = L[0];

  if (Gi1 == 0) {
    return list({integer(1), integer(0)});
  }

  t4 = pseudoRemainder(Gi2, Gi1, x);

  //printf("r = \n%s\n", t4.toString().c_str());

  if (t4 == 0) {

    nk = degree(Gi1, x);

    if (nk.value() > 0) {

      cnt = cont(Gi1, L, K);

      ppk = recQuotient(Gi1, cnt, L, K);

      r = list({ppk, integer(0)});

    } else {
      r = list({integer(1), Gi1});
    }

    return r;
  }


  d = degree(Gi2, x) - degree(Gi1, x);
  t1 = d + 1;
  t2 = power(-1, t1);
  t3 = t2 * t4;
  t4 = algebraicExpand(t3);

  t1 = leadCoeff(Gi2, x);
  t2 = power(hi2, d);
  t3 = mul({t1, t2});

  t5 = algebraicExpand(t3);

	//printf("a = \n%s\n\n", t4.toString().c_str());
	//printf("b = \n%s\n\n", t5.toString().c_str());

  Gi = recQuotient(t4, t5, L, K);

	//printf("Gi = \n%s\n\n", Gi.toString().c_str());
	//printf("d = \n%s\n\n", d.toString().c_str());

  hi1 = algebraicExpand(power(leadCoeff(Gi1, x), d) * power(hi2, 1 - d)); //h4

  //printf("c = \n%s\n\n", hi1.toString().c_str());
  return polyRemSeqRec(Gi1, Gi, L, hi1, K);
}

Expr polyRemSeq(Expr F1, Expr F2, Expr L, Expr K) {

  if (F1.kind() == Kind::Integer && F2.kind() == Kind::Integer) {
    return gcd(F1.value(), F2.value());
  }

  Expr x = L[0];


  Expr m = degree(F1, x);
  Expr n = degree(F2, x);
	F1 = algebraicExpand(F1);
	F2 = algebraicExpand(F2);
	//printf("F = \n%s\n", F1.toString().c_str());
	//printf("G = \n%s\n", F2.toString().c_str());
	//printf("%s\n", m.toString().c_str());
	//printf("%s\n", n.toString().c_str());

  if (m.value() < n.value()) {
    return polyRemSeq(F2, F1, L, K);
  }

  Expr d, t1, t2, t3, t4, t5;
  Expr G1, G2, G3, h2, nk, ppk, cnt, r;

  G1 = F1;
  G2 = F2;

  // compute G[3]
  d  = degree(F1, x) - degree(F2, x);
  t3 = power(-1, d + 1);

	//printf("b = \n%s\n", t3.toString().c_str());

	t4 = pseudoRemainder(G1, G2, x);

	//printf("prem = \n%s\n", t4.toString().c_str());

  G3 = reduceAST(mulPoly(power(-1, d + 1), t4));
	//printf("G3 = \n%s\n", G3.toString().c_str());

  if (G3 == 0) {
    nk = degree(G2, x);

    if (nk.value() > 0) {
      cnt = cont(G2, L, K);
      ppk = recQuotient(G2, cnt, L, K);

      r = list({ppk, integer(0)});
    } else {
      r = list({integer(1), G2});
    }

    return r;
  }

  // compute h[2]
  h2 = reduceAST(power(leadCoeff(G2, x), d));

  return polyRemSeqRec(G2, G3, L, h2, K);
}


Expr colPolyResultantRec(Expr u, Expr v, Expr L, Expr K, Expr i,
                              Expr delta_prev, Expr gamma_prev) {
  assert(u != 0, "Polynomial should be non-zero");
  assert(v != 0, "Polynomial should be non-zero");

	Expr m = degreePolyExpr(u);
  Expr n = degreePolyExpr(v);

  if (m.value() < n.value()) {
		Expr t = colPolyResultantRec(v, u, L, K, i, delta_prev, gamma_prev);
		Expr g = polyExpr(pow(-1, m.value() * n.value()), L);

		return mulPolyExpr(g, t);
  }

  if (n == 0) {
    return powPolyExpr(v, m.value());
  }

  Expr r = pseudoRemPolyExpr(u, v, L[0]);

  if (r == 0) {
    return polyExpr(0, L);
  }

	Expr delta = m.value() - n.value() + 1;

  Expr R = rest(L);

  Expr gama = undefined();
  Expr beta = undefined();

  if (i == 1) {
    gama = -1;
    beta = pow(-1, delta.value());
  } else {
		Expr k = - 1;
    Expr f = leadCoeffPolyExpr(u);
		Expr r = mulPolyExpr(k, f);

		Expr tmp1 = powPolyExpr(r, delta_prev.value() - 1);
		Expr tmp2 = pow(gamma_prev.value(), delta_prev.value() - 2);

    gama = recQuoPolyExpr(tmp1, tmp2, R, K);

		tmp1 = polyExpr(pow(gama.value(), delta.value() - 1), R);

		beta = mulPolyExpr(r, tmp1);
  }

  Expr t = recQuoPolyExpr(r, beta, L, K);

  r = t;

	Expr tmp1 = colPolyResultantRec(v, r, L, K, i.value() + 1, delta, gama);
	Expr tmp2 = pow(-1, m.value() * n.value()) * pow(beta.value(), n.value());

	Expr w = mulPolyExpr(tmp2, tmp1);

  Expr l = leadCoeffPolyExpr(v);

  Expr s = degreePolyExpr(r);

	Int k = delta.value() * n.value() + -1*m.value() + s.value();

  Expr f = powPolyExpr(l, k);

  return recQuoPolyExpr(w, f, L, K);
}

Expr colPolyResultant(Expr u, Expr v, Expr L, Expr K) {
  // Expr x = L[0];

  Expr m = degreePolyExpr(u);
  Expr n = degreePolyExpr(v);

	// TODO write colPolyCont
  Expr cont_u = cont(u, L, K);

  Expr pp_u = recQuoPolyExpr(u, cont_u, L, K);

	// TODO write colPolyCont
  Expr cont_v = cont(v, L, K);

  Expr pp_v = recQuoPolyExpr(v, cont_v, L, K);

  Expr i = integer(1);
  Expr delta = integer(0);
  Expr g = integer(0);

  Expr s = colPolyResultantRec(pp_u, pp_v, L, K, i, delta, g);

	Expr k = pow(cont_u.value(), n.value()) * pow(cont_v.value(), m.value());

	return mulPolyExpr(k, s);
}


Expr colPolyRemSeqRec(Expr Gi2, Expr Gi1, Expr L, Expr hi2, Expr K) {
  Expr Gi, hi1, d, t1, t2, t3, t4, t5, t6, nk, cnt, ppk, r;

  if (Gi1 == 0) {
    return list({polyExpr(1, L), polyExpr(0, L)});
  }

  t4 = pseudoRemPolyExpr(Gi2, Gi1, L[0]);

  if (t4 == 0) {

    nk = degreePolyExpr(Gi1);

    if (nk.value() > 0) {
			// TODO: write contColPoly
      cnt = polyExpr(cont(algebraicExpand(Gi1), L, K), L);

      ppk = recQuoPolyExpr(Gi1, cnt, L, K);

      r = list({ ppk, polyExpr(0, L) });
    } else {
      r = list({ polyExpr(1, L), Gi1 });
    }

    return r;
  }

  d = degreePolyExpr(Gi2) - degreePolyExpr(Gi1);

	t2 = pow(-1, d.value() + 1);
  t4 = mulPolyExpr(t2, t4);

  t1 = leadCoeffPolyExpr(Gi2);
  t2 = powPolyExpr(hi2, d.value());
  t5 = mulPolyExpr(t2, t1);

  Gi = recQuoPolyExpr(t4, t5, L, K);

	Expr a = powPolyExpr(leadCoeffPolyExpr(Gi1), d.value());
	Expr b = powPolyExpr(hi2, Int(1) - d.value());

	hi1 = mulPolyExpr(a, b);

  return colPolyRemSeqRec(Gi1, Gi, L, hi1, K);
}


Expr colPolyRemSeq(Expr F1, Expr F2, Expr L, Expr K) {

  if (F1.kind() == Kind::Integer && F2.kind() == Kind::Integer) {
    return gcd(F1.value(), F2.value());
  }

  //Expr x = L[0];

  Expr m = degreePolyExpr(F1);
  Expr n = degreePolyExpr(F2);

  if (m.value() < n.value()) {
    return colPolyRemSeq(F2, F1, L, K);
  }

  Expr t1, t2, t3, t4, t5;
  Expr G1, G2, G3, h2, nk, ppk, cnt, r;

  G1 = F1;
  G2 = F2;

  // compute G[3]
  Int d  = m.value() - n.value();

	t3 = pow(-1, d + 1);

	t4 = pseudoRemPolyExpr(G1, G2, L[0]);

  Expr k = pow(-1, d + 1);

	G3 = mulPolyExpr(k, t4);

  if (G3 == 0) {
    nk = degreePolyExpr(G2);

    if (nk.value() > 0) {
			// TODO: write contColPoly
      cnt = polyExpr(cont(algebraicExpand(G2), L, K), L);
      ppk = recQuoPolyExpr(G2, cnt, L, K);

      r = list({ppk, polyExpr(0, L)});
    } else {
      r = list({polyExpr(1, L), G2});
    }

    return r;
  }

  // compute h[2]
	h2 = powPolyExpr(leadCoeffPolyExpr(G2), d);

  return colPolyRemSeqRec(G2, G3, L, h2, K);
}


} // namespace polynomial
