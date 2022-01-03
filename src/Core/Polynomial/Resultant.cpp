#include "Resultant.hpp"

#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace simplification;

namespace polynomial {

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

  Expr delta = reduceAST(add({m, mul({integer(-1), n}), integer(1)}));

  Expr R = rest(L);

  Expr gamma = undefined();
  Expr beta = undefined();

  if (i == 1) {
    gamma = integer(-1);

    Expr tmp = power(integer(-1), delta);

    beta = reduceAST(tmp);

  } else {
    Expr f = leadCoeff(u, x);
		Expr r = algebraicExpand(-1*f);

		Expr tmp1 = power(r, sub({delta_prev, integer(1)}));
    Expr tmp2 = power(gamma_prev, sub({delta_prev, integer(2)}));

    Expr tmp3 = algebraicExpand(tmp1);
    Expr tmp4 = algebraicExpand(tmp2);

		gamma = recQuotient(tmp3, tmp4, R, K);

		tmp1 = power(gamma, delta.value() - 1);

		Expr tmp5 = mul({r, tmp1});

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

  return recQuotient(w, f, L, K);
}

Expr polynomialResultant(Expr u, Expr v, Expr L, Expr K) {
  Expr x = L[0];

  Expr m = degree(u, x);
  Expr n = degree(v, x);

  Expr cont_u = cont(u, L, K);
  Expr pp_u = recQuotient(u, cont_u, L, K);
  Expr cont_v = cont(v, L, K);
  Expr pp_v = recQuotient(v, cont_v, L, K);

  Expr i = 1;
  Expr d = 0;
  Expr g = 0;

  Expr s = polynomialResultantRec(pp_u, pp_v, L, K, i, d, g);

  return reduceAST(power(cont_u, n) * power(cont_v, m) * s);
}

Expr polyRemSeqRec(Expr Gi2, Expr Gi1, Expr L, Expr hi2, Expr K) {
  Expr Gi, hi1, d, t1, t2, t3, t4, t5, t6, nk, cnt, ppk, r, x;
  x = L[0];

  if (Gi1 == 0) {
    return list({integer(1), integer(0)});
  }

  t4 = pseudoRemainder(Gi2, Gi1, x);

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

  Gi = recQuotient(t4, t5, L, K);

	Expr c = leadCoeff(Gi1, x);
	Expr a = power(c, d);
	Expr b = power(hi2, 1 - d);

  hi1 = algebraicExpand(a * b); //h4

	return polyRemSeqRec(Gi1, Gi, L, hi1, K);
}

Expr polyRemSeq(Expr F1, Expr F2, Expr L, Expr K) {

	Expr U = polyExpr(F1, L);
  Expr V = polyExpr(F2, L);

	Expr G = remSeqPolyExpr(U, V, L, K);

	return list({ expandPolyExpr(G[0]), expandPolyExpr(G[1]) });

  // if (F1.kind() == Kind::Integer && F2.kind() == Kind::Integer) {
  //   return gcd(F1.value(), F2.value());
  // }

  // Expr x = L[0];

  // Expr m = degree(F1, x);
  // Expr n = degree(F2, x);

	// F1 = algebraicExpand(F1);
	// F2 = algebraicExpand(F2);

	// if (m.value() < n.value()) {
  //   return polyRemSeq(F2, F1, L, K);
  // }

  // Expr d, t1, t2, t4, t5;
  // Expr G1, G2, G3, h2, nk, ppk, cnt, r;

  // G1 = F1;
  // G2 = F2;


  // d  = degree(F1, x) - degree(F2, x);

	// t4 = pseudoRemainder(G1, G2, x);

  // G3 = reduceAST(mulPoly(power(-1, d + 1), t4));

  // if (G3 == 0) {
  //   nk = degree(G2, x);

  //   if (nk.value() > 0) {
  //     cnt = cont(G2, L, K);
  //     ppk = recQuotient(G2, cnt, L, K);

  //     r = list({ppk, integer(0)});
  //   } else {
  //     r = list({integer(1), G2});
  //   }

  //   return r;
  // }

  // // compute h[2]
  // h2 = reduceAST(power(leadCoeff(G2, x), d));

  // return polyRemSeqRec(G2, G3, L, h2, K);
}


Expr resultantPolyExprRec(Expr u, Expr v, Expr L, Expr K, Expr i,
                              Expr delta_prev, Expr gamma_prev) {
  assert(u != 0, "Polynomial should be non-zero");
  assert(v != 0, "Polynomial should be non-zero");

	Expr m = degreePolyExpr(u);
  Expr n = degreePolyExpr(v);

	if (m.value() < n.value()) {
		Expr t = resultantPolyExprRec(v, u, L, K, i, delta_prev, gamma_prev);

		Expr g = polyExpr(pow(-1, m.value() * n.value()), L);

		return mulPolyExpr(g, t);
  }

  if (isZeroPolyExpr(n)) {
    return powPolyExpr(v, m.value());
  }

	Expr r = pseudoRemPolyExpr(u, v, L);

  if (isZeroPolyExpr(r)) {
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
		Expr tmp2 = polyExpr(pow(gamma_prev.value(), delta_prev.value() - 2), R);

		gama = quoPolyExpr(tmp1, tmp2, R, K);

		tmp1 = powPolyExpr(gama, delta.value() - 1);

		beta = raisePolyExpr(mulPolyExpr(r, tmp1), 0, L[0]);
  }

  r = quoPolyExpr(r, beta, L, K);

	Expr tmp1 = resultantPolyExprRec(v, r, L, K, i.value() + 1, delta, gama);

	Expr tmp2 = mulPolyExpr(pow(-1, m.value() * n.value()), powPolyExpr(beta, n.value()));

	Expr w = mulPolyExpr(tmp2, tmp1);

  Expr l = raisePolyExpr(leadCoeffPolyExpr(v), 0, L[0]);

  Expr s = degreePolyExpr(r);

	Int k = delta.value() * n.value() + -1*m.value() + s.value();

	Expr f = powPolyExpr(l, k);

  return quoPolyExpr(w, f, L, K);
}

Expr resultantPolyExpr(Expr u, Expr v, Expr L, Expr K) {
  // Expr x = L[0];
	Expr R = rest(L);

  Expr m = degreePolyExpr(u);
  Expr n = degreePolyExpr(v);

	Expr U = contAndPpPolyExpr(u, L, K);

	Expr ct_u = U[0];
  Expr pp_u = U[1];

	Expr V = contAndPpPolyExpr(v, L, K);

	Expr ct_v = V[0];
  Expr pp_v = V[1];

  Expr i = 1;
  Expr d = 0;
  Expr g = 0;

	Expr s = resultantPolyExprRec(pp_u, pp_v, L, K, i, d, g);

	Expr a = powPolyExpr(ct_u, n.value());
	Expr b = powPolyExpr(ct_v, m.value());
	Expr k = mulPolyExpr(a, b);

	k = raisePolyExpr(k, 0, L[0]);

	return mulPolyExpr(k, s);
}


Expr remSeqPolyExprRec(Expr Gi2, Expr Gi1, Expr L, Expr hi2, Expr K) {
  Expr Gi, hi1, d, t1, t2, t3, t4, t5, t6, nk, cnt, ppk, r;

	if (isZeroPolyExpr(Gi1)) {
    return list({ polyExpr(0, L), polyExpr(1, L) });
  }

  t4 = pseudoRemPolyExpr(Gi2, Gi1, L);

  if (isZeroPolyExpr(t4)) {
		return list({ Gi1, polyExpr(1, L) });

		nk = degreePolyExpr(Gi1);

    if (nk.value() > 0) {
      return list({ polyExpr(0, L), ppPolyExpr(Gi1, L, K) });
    }

    return list({ polyExpr(1, L), Gi1 });
  }

  d = degreePolyExpr(Gi2).value() - degreePolyExpr(Gi1).value();

	t2 = pow(-1, d.value() + 1);
  t4 = mulPolyExpr(t2, t4);

	Expr R = rest(L);

	t1 = leadCoeffPolyExpr(Gi2);
	t2 = powPolyExpr(hi2, d.value());
  t5 = mulPolyExpr(t2, t1);
	t5 = raisePolyExpr(t5, 0, L[0]);

	Gi = quoPolyExpr(t4, t5, L, K);

	Expr c = leadCoeffPolyExpr(Gi1);
	Expr a = powPolyExpr(c, d.value());
	Expr b = powPolyExpr(hi2, Int(1) - d.value());

	if(b.kind() == Kind::Division) {
		hi1 = quoPolyExpr(a, b[1], L, K);
	} else {
		hi1 = mulPolyExpr(a, b);
	}

	return remSeqPolyExprRec(Gi1, Gi, L, hi1, K);
}


Expr remSeqPolyExpr(Expr G1, Expr G2, Expr L, Expr K) {

	if(isZeroPolyExpr(G1)) {
		return  list({ polyExpr(0, L), polyExpr(0, L) });
	}

	if(isZeroPolyExpr(G2)) {
		return  list({ G1, polyExpr(1, L) });
	}

  if (G1.kind() == Kind::Integer && G2.kind() == Kind::Integer) {
    return gcd(G1.value(), G2.value());
  }

  Expr m = degreePolyExpr(G1);
  Expr n = degreePolyExpr(G2);

  if (m.value() < n.value()) {
    return remSeqPolyExpr(G2, G1, L, K);
  }

  Expr t1, t2, t4, t5;
  Expr G3, h2, nk, ppk, cnt, r;

  Int d  = m.value() - n.value();

	t4 = pseudoRemPolyExpr(G1, G2, L);

  Expr k = pow(-1, d + 1);

	G3 = mulPolyExpr(k, t4);

  if (isZeroPolyExpr(G3)) {
		//return list({ G1, G2 });
    nk = degreePolyExpr(G2);

    if (nk.value() > 0) {
      return list({ ppPolyExpr(G2, L, K), polyExpr(1, L) });
    }

    return list({ polyExpr(1, L), G2 });
  }

  Expr R = rest(L);

	h2 = powPolyExpr(leadCoeffPolyExpr(G2), d);

	return remSeqPolyExprRec(G2, G3, L, h2, K);
}

} // namespace polynomial
