#include "Resultant.hpp"

#include "MathSystem/Algebra/Expression.hpp"
#include "MathSystem/Debug/Assert.hpp"
#include "MathSystem/Polynomial/Polynomial.hpp"

using namespace alg;

namespace polynomial {

expr polynomialResultantRec(expr u, expr v, expr L, expr K, expr i,
                              expr delta_prev, expr gamma_prev) {
  assert(u != 0, "Polynomial should be non-zero");
  assert(v != 0, "Polynomial should be non-zero");

	expr x = L[0];
	expr m = degree(u, x);
  expr n = degree(v, x);

  if (m.value() < n.value()) {
    expr e =
			create(kind::MUL, { pow(integer(-1), m*n),
             polynomialResultantRec(v, u, L, K, i, delta_prev, gamma_prev)});

    expr k = reduce(e);

    return k;
  }

  if (n == 0) {
    expr e = pow(v, m);

    expr k = reduce(e);

    return k;
  }

  expr r = pseudoRemainder(u, v, x);

  if (r == 0) {
    return integer(0);
  }

  expr delta = reduce(m + -n + 1);

  expr R = rest(L);

  expr gamma = undefined();
  expr beta = undefined();

  if (i == 1) {
    gamma = integer(-1);

    expr tmp = pow(integer(-1), delta);

    beta = reduce(tmp);

  } else {
    expr f = leadCoeff(u, x);
		expr r = expand(-1*f);

		expr tmp1 = pow(r, delta_prev - 1);
    expr tmp2 = pow(gamma_prev, delta_prev - 2);

    expr tmp3 = expand(tmp1);
    expr tmp4 = expand(tmp2);

		gamma = recQuotient(tmp3, tmp4, R, K);

		tmp1 = pow(gamma, delta.value() - 1);

		expr tmp5 = r*tmp1;

    beta = expand(tmp5);
  }
  expr t = recQuotient(r, beta, L, K);

  r = t;

  expr tmp1 = i + 1;
  expr tmp2 = reduce(tmp1);

  expr tmp3 = pow(-1, m*n)*pow(beta, n)*
                   polynomialResultantRec(v, r, L, K, tmp2, delta, gamma);

  expr w = expand(tmp3);

  expr l = coeff(v, x, n);

  expr s = degree(r, x);

  expr k = delta*n + -m + s;

  expr tmp4 = pow(l, k);

  expr f = expand(tmp4);

  return recQuotient(w, f, L, K);
}

expr polynomialResultant(expr u, expr v, expr L, expr K) {
  expr x = L[0];

  expr m = degree(u, x);
  expr n = degree(v, x);

  expr cont_u = cont(u, L, K);
  expr pp_u = recQuotient(u, cont_u, L, K);
  expr cont_v = cont(v, L, K);
  expr pp_v = recQuotient(v, cont_v, L, K);

  expr i = 1;
  expr d = 0;
  expr g = 0;

  expr s = polynomialResultantRec(pp_u, pp_v, L, K, i, d, g);

  return reduce(pow(cont_u, n) * pow(cont_v, m) * s);
}

expr polyRemSeqRec(expr Gi2, expr Gi1, expr L, expr hi2, expr K) {
  expr Gi, hi1, d, t1, t2, t3, t4, t5, t6, nk, cnt, ppk, r, x;
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
  t2 = pow(-1, t1);
  t3 = t2 * t4;
  t4 = expand(t3);

  t1 = leadCoeff(Gi2, x);
  t2 = pow(hi2, d);
  t3 = t1 * t2;

  t5 = expand(t3);

  Gi = recQuotient(t4, t5, L, K);

	expr c = leadCoeff(Gi1, x);
	expr a = pow(c, d);
	expr b = pow(hi2, 1 - d);

  hi1 = expand(a * b); //h4

	return polyRemSeqRec(Gi1, Gi, L, hi1, K);
}

expr polyRemSeq(expr F1, expr F2, expr L, expr K) {

	expr U = polyExpr(F1, L);
  expr V = polyExpr(F2, L);

	expr G = remSeqPolyExpr(U, V, L, K);

	return list({ expandPolyExpr(G[0]), expandPolyExpr(G[1]) });

  // if (F1.kind() == Kind::Integer && F2.kind() == Kind::Integer) {
  //   return gcd(F1.value(), F2.value());
  // }

  // expr x = L[0];

  // expr m = degree(F1, x);
  // expr n = degree(F2, x);

	// F1 = expand(F1);
	// F2 = expand(F2);

	// if (m.value() < n.value()) {
  //   return polyRemSeq(F2, F1, L, K);
  // }

  // expr d, t1, t2, t4, t5;
  // expr G1, G2, G3, h2, nk, ppk, cnt, r;

  // G1 = F1;
  // G2 = F2;


  // d  = degree(F1, x) - degree(F2, x);

	// t4 = pseudoRemainder(G1, G2, x);

  // G3 = reduce(mulPoly(pow(-1, d + 1), t4));

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
  // h2 = reduce(pow(leadCoeff(G2, x), d));

  // return polyRemSeqRec(G2, G3, L, h2, K);
}


expr resultantPolyExprRec(expr u, expr v, expr L, expr K, expr i,
                              expr delta_prev, expr gamma_prev) {
	assert(!isZeroPolyExpr(u), "Polynomial should be non-zero");
  assert(!isZeroPolyExpr(v), "Polynomial should be non-zero");


	expr m = degreePolyExpr(u);
  expr n = degreePolyExpr(v);

  // printf("A\n");

  if (m.value() < n.value()) {
		expr t = resultantPolyExprRec(v, u, L, K, i, delta_prev, gamma_prev);

		expr g = polyExpr(pow(-1, m.value() * n.value()), L);

		return mulPolyExpr(g, t);
  }

  // printf("B\n");
  if (isZeroPolyExpr(n)) {
    return powPolyExpr(v, m.value());
	}

  // printf("D\n");

	expr r = pseudoRemPolyExpr(u, v, L);

	// printf("E\n");

  if (isZeroPolyExpr(r)) {
    return polyExpr(0, L);
  }

	expr delta = m.value() - n.value() + 1;

	expr R = rest(L);

  expr gama = undefined();
  expr beta = undefined();

	// printf("F\n");

 if (i == 1) {
    gama = -1;
    beta = pow(-1, delta.value());
  } else {
		// printf("-> A\n");
		expr k = - 1;

		expr f = leadCoeffPolyExpr(u);

		// printf("-> B\n");
		expr r = mulPolyExpr(k, f);

		// printf("-> C\n");
		expr tmp1 = powPolyExpr(r, delta_prev.value() - 1);
		// printf("-> D\n");

		expr t = powPolyExpr(gamma_prev, delta_prev.value() - 2);

		expr tmp2 = polyExpr(t, R);

		// printf("-> E\n");
		gama = quoPolyExpr(tmp1, tmp2, R, K);

		// printf("-> F\n");
		tmp1 = powPolyExpr(gama, delta.value() - 1);

		// printf("-> G\n");
		beta = raisePolyExpr(mulPolyExpr(r, tmp1), 0, L[0]);
  }

  // printf("CC\n");
  r = quoPolyExpr(r, beta, L, K);
  // printf("DD\n");
	expr tmp1 = resultantPolyExprRec(v, r, L, K, i.value() + 1, delta, gama);

	expr tmp2 = mulPolyExpr(pow(-1, m.value() * n.value()), powPolyExpr(beta, n.value()));

	expr w = mulPolyExpr(tmp2, tmp1);

  expr l = raisePolyExpr(leadCoeffPolyExpr(v), 0, L[0]);

  expr s = degreePolyExpr(r);

	Int k = delta.value() * n.value() + -1*m.value() + s.value();

	expr f = powPolyExpr(l, k);

  return quoPolyExpr(w, f, L, K);
}

expr resultantPolyExpr(expr u, expr v, expr L, expr K) {
	expr R = rest(L);

  expr m = degreePolyExpr(u);
  expr n = degreePolyExpr(v);

	expr U = contAndPpPolyExpr(u, L, K);

	expr ct_u = U[0];
  expr pp_u = U[1];

	expr V = contAndPpPolyExpr(v, L, K);

	expr ct_v = V[0];
  expr pp_v = V[1];

  expr i = 1;
  expr d = 0;
  expr g = 0;
	// printf("aaaa\n");
	expr s = resultantPolyExprRec(pp_u, pp_v, L, K, i, d, g);

	expr a = powPolyExpr(ct_u, n.value());
	expr b = powPolyExpr(ct_v, m.value());
	expr k = mulPolyExpr(a, b);

	k = raisePolyExpr(k, 0, L[0]);

	return mulPolyExpr(k, s);
}


expr remSeqPolyExprRec(expr Gi2, expr Gi1, expr L, expr hi2, expr K) {
  expr Gi, hi1, d, t1, t2, t3, t4, t5, t6, nk, cnt, ppk, r;

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

	expr R = rest(L);

	t1 = leadCoeffPolyExpr(Gi2);
	t2 = powPolyExpr(hi2, d.value());
  t5 = mulPolyExpr(t2, t1);
	t5 = raisePolyExpr(t5, 0, L[0]);

	Gi = quoPolyExpr(t4, t5, L, K);

	expr c = leadCoeffPolyExpr(Gi1);
	expr a = powPolyExpr(c, d.value());
	expr b = powPolyExpr(hi2, Int(1) - d.value());

	if(b.kind() == kind::DIV) {
		hi1 = quoPolyExpr(a, b[1], L, K);
	} else {
		hi1 = mulPolyExpr(a, b);
	}

	return remSeqPolyExprRec(Gi1, Gi, L, hi1, K);
}


expr remSeqPolyExpr(expr G1, expr G2, expr L, expr K) {

	if(isZeroPolyExpr(G1)) {
		return  list({ polyExpr(0, L), polyExpr(0, L) });
	}

	if(isZeroPolyExpr(G2)) {
		return  list({ G1, polyExpr(1, L) });
	}

  if (G1.kind() == kind::INT && G2.kind() == kind::INT) {
    return gcd(G1.value(), G2.value());
  }

  expr m = degreePolyExpr(G1);
  expr n = degreePolyExpr(G2);

  if (m.value() < n.value()) {
    return remSeqPolyExpr(G2, G1, L, K);
  }

  expr t1, t2, t4, t5;
  expr G3, h2, nk, ppk, cnt, r;

  Int d  = m.value() - n.value();

	t4 = pseudoRemPolyExpr(G1, G2, L);

  expr k = pow(-1, d + 1);

	G3 = mulPolyExpr(k, t4);

  if (isZeroPolyExpr(G3)) {
		//return list({ G1, G2 });
    nk = degreePolyExpr(G2);

    if (nk.value() > 0) {
      return list({ ppPolyExpr(G2, L, K), polyExpr(1, L) });
    }

    return list({ polyExpr(1, L), G2 });
  }

  expr R = rest(L);

	h2 = powPolyExpr(leadCoeffPolyExpr(G2), d);

	return remSeqPolyExprRec(G2, G3, L, h2, K);
}

} // namespace polynomial
