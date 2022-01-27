
#include "Wang.hpp"
#include "Berlekamp.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "SquareFree.hpp"
#include "Utils.hpp"
#include "Zassenhaus.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/GaloisField/GaloisField.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Simplification/Simplification.hpp"

#include <cstddef>
#include <cstdio>
#include <limits>
#include <random>

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace galoisField;
using namespace simplification;

namespace factorization {

Expr nondivisors(Int G, Expr F, Int d, Expr L, Expr K) {
	// TODO: remove
  assert(G != 0, "G needs to be different from zero!");
  assert(d != 0, "d needs to be different from zero!");

  long i, j, k;

  Int q, r;

  Expr Fi;

  k = F.size();

  Int *x = new Int[k + 1];

  x[0] = G * d;
  for (i = 1; i <= k; i++) {
    Fi = F[i - 1];

    q = norm(Fi, L, K);

    for (j = i - 1; j >= 0; j--) {
      r = x[j];

      while (abs(r) != 1) {
        r = gcd(r, q);
        q = q / r;
      }

      if (q == 1) {
        return fail();
      }
    }

    x[i] = q;
  }

  Expr p = list({});

  for (i = 1; i <= k; i++) {
    p.insert(x[i]);
  }

  delete[] x;

  return p;
}

Expr expandList(Expr L) {
  Expr K = list({});
  for (long i = 0; i < L.size(); i++) {
    K.insert(algebraicExpand(L[i]));
  }
  return K;
}

Expr nondivisorsPolyExpr(Int G, Expr F, Int d, Expr L, Expr K) {
  assert(G != 0, "G needs to be different from zero!");
  assert(d != 0, "d needs to be different from zero!");

  long i, j, k;

  Int q, r;

  Expr Fi;

  k = F.size();

  Int *x = new Int[k + 1];

  x[0] = G * d;

	for (i = 1; i <= k; i++) {
    Fi = F[i - 1];

    q = normPolyExpr(Fi, L, K);

    for (j = i - 1; j >= 0; j--) {
      r = x[j];

      while (abs(r) != 1) {
        r = gcd(r, q);
        q = q / r;
      }

      if (q == 1) {
        return fail();
      }
    }

    x[i] = q;
  }

  Expr p = list({});

  for (i = 1; i <= k; i++) {
    p.insert(x[i]);
  }

  delete[] x;

  return p;
}

Expr trialDivision(Expr& f, Expr& F, Expr& L, Expr K) {
  Expr v, q, r, d, t = list({});

  bool stop = false;
  long i, k;

  for (i = 0; i < F.size(); i++) {
    k = 0;

    stop = false;

    v = algebraicExpand(F[i]);

    while (!stop) {
      d = recPolyDiv(f, v, L, K);

      q = d[0];
      r = d[1];

      if (r == 0) {
        f = q;

        k = k + 1;
      }

      stop = r != 0;
    }

		if (k > 0) {
      t.insert(list({F[i], integer(k)}));
    }
  }

  return t;
}

Expr trialDivisionPolyExpr(Expr& f, Expr& F, Expr& L, Expr K) {
  Expr t = list({});
  for (Int i = 0; i < F.size(); i++) {
    Int k = 0;

    bool stop = false;

    Expr v = F[i];

    while (!stop) {
      Expr d = divPolyExpr(f, v, L, K);

      Expr q = d[0];
      Expr r = d[1];

      if (isZeroPolyExpr(r)) {
        f = q;

        k = k + 1;
      }

      stop = !isZeroPolyExpr(r);
    }

    if (k > 0) {
      t.insert(list({ F[i], k }));
    }
  }

  return t;
}

Expr sqfFactors(Expr f, Expr x, Expr K) {
  Expr n, cn, pr, lc, L, t1, F;

  L = list({x});
  cn = cont(f, L, K);
  pr = pp(f, cn, L, K);
  lc = leadCoeff(pr, x);

  if (lc.value() < 0) {
    cn = groundInvert(cn);
    pr = groundInvert(pr);
  }

  n = degree(pr, x);

  if (n.value() <= 0) {
    return list({cn, list({})});
  }

  if (n.value() == 1) {
    return list({cn, list({pr})});
  }

  F = zassenhaus(pr, x, K);

  return list({cn, F});
}


Expr sqfFactorsPolyExpr(Expr f, Expr L, Expr K) {
  assert(L.size() == 1, "should be a univaritate polynomial");

  Expr n, t1, F;

	Expr CP = contAndPpPolyExpr(f, L, K);

	Expr cn = CP[0];
  Expr pr = CP[1];

  Expr lc = leadCoeffPolyExpr(pr);

  assert(lc.kind() == Kind::Integer, "lead coefficient should be an integer");

  if (lc.value() < 0) {
    cn = groundInvertPolyExpr(cn);
    pr = groundInvertPolyExpr(pr);
  }

  n = degreePolyExpr(pr);

  assert(n.kind() == Kind::Integer, "degree of polynomial should be an integer");

	if (n.value() <= 0) {
    return list({cn, list({})});
  }

  if (n.value() == 1) {
    return list({cn, list({pr})});
  }

  return list({ cn, zassenhausPolyExpr(pr, L, K) });
}



Expr factors(Expr f, Expr L, Expr K) {
  if (f == 0) {
    return list({0, list({})});
  }

  if (L.size() == 1) {
    Expr cnt = cont(f, L, K);
    Expr ppr = pp(f, cnt, L, K);

    // TODO: maybe unnecessary
    Expr lc = leadCoeff(ppr, L[0]);
    assert(lc.kind() == Kind::Integer, "not integer coefficient");
    if (lc.value() < 0) {
      cnt = mulPoly(cnt, -1);
      ppr = mulPoly(ppr, -1);
    }

    Expr n = degree(ppr, L[0]);

    if (n == 0) {
      return list({cnt, list({})});
    }

    if (n == 1) {
      return list({cnt, list({ppr, 1})});
    }

    Expr g = squareFreePart(ppr, L, K);
    Expr H = zassenhaus(ppr, L[0], K);

    Expr F = trialDivision(ppr, H, L, K);

    return list({cnt, F});
  }

  if (L.size() == 0) {
    return list({f, list({})});
  }

  Expr x = L[0];
  Expr R = rest(L);
  Expr cnt = groundContPoly(f, L, K);
  Expr prp = groundPPPoly(f, L, K);

  Expr lc = groundLeadCoeffPoly(prp, L);

  // TODO: mayne unnecessary
  if (lc.value() < 0) {
    cnt = mulPoly(cnt, -1);
    prp = mulPoly(prp, -1);
  }

  bool is_const = true;

  // TODO: maybe just verify if prp is integer or fraction
  for (Int i = 0; i < L.size() && is_const; i++) {
    Expr d = degree(prp, L[i]);
    if (d != 0)
      is_const = false;
  }
  if (is_const) {
    return list({cnt, list({})});
  }

  Expr G = cont(prp, L, K);
  Expr g = pp(prp, G, L, K);

  Expr F = list({});

  Expr n = degree(g, x);

  if (n.value() > 0) {
		Expr S = squareFreePart(g, L, K);
    Expr H = factorsWang(S[0], S[1], K);

    F = trialDivision(f, H, L, K);
  }

  Expr t1 = factors(G, R, K);

  for (size_t i = 0; i < t1[1].size(); i++) {
    F.insert(list({t1[1][i][0], t1[1][i][1]}));
  }

  return list({cnt, F});
}





Expr factorsPolyExpr(Expr f, Expr L, Expr K) {
	if (isZeroPolyExpr(f)) {
    return list({f, list({})});
  }

  if (L.size() == 1) {
		Expr CP = contAndPpPolyExpr(f, L, K);

		Expr cnt = CP[0];
    Expr ppr = CP[1];

    Expr lc = leadCoeffPolyExpr(ppr);

    assert(lc.kind() == Kind::Integer, "lc should be an integer");

		if (lc.value() < 0) {
			cnt = groundInvertPolyExpr(cnt);
			ppr = groundInvertPolyExpr(ppr);
    }

		cnt = raisePolyExpr(cnt, 0, L[0]);

    Expr n = degreePolyExpr(ppr);

    assert(n.kind() == Kind::Integer, "n should be an integer");

		if (n == 0) {
			return list({cnt, list({})});
    }

    if (n == 1) {
      return list({cnt, list({ppr, 1})});
    }

    Expr g = squareFreePartPolyExpr(ppr, L, K);
    Expr H = zassenhausPolyExpr(ppr, L, K);
    Expr F = trialDivisionPolyExpr(ppr, H, L, K);

    return list({cnt, F});
  }

  Expr R = rest(L);

	Expr cnt = groundContPolyExpr(f);
  Expr prp = groundPPPolyExpr(f);

  Expr lc = groundLeadCoeffPolyExpr(prp);

	assert(lc.kind() == Kind::Integer, "lc should be an integer");

  if (lc.value() < 0) {
		cnt = groundInvertPolyExpr(cnt);
		prp = groundInvertPolyExpr(prp);
  }

  bool is_const = true;

  // TODO: maybe just verify if prp is integer or fraction
  for (Int i = 0; i < L.size() && is_const; i++) {
    Expr d = degreePolyExpr(prp);
    if (d != 0) is_const = false;
  }

  if (is_const) {
    return list({cnt, list({})});
  }

	Expr O = contAndPpPolyExpr(prp, L, K);

	Expr G = O[0];
  Expr g = O[1];
  Expr F = list({});

  Expr n = degreePolyExpr(g);

  if (n.value() > 0) {
		Expr S = squareFreePartPolyExpr(g, L, K);
    Expr H = factorsWangPolyExpr(S[0], S[1], K);

    F = trialDivisionPolyExpr(f, H, L, K);
  }

  Expr t1 = factorsPolyExpr(G, R, K);

  for (size_t i = 0; i < t1[1].size(); i++) {
    F.insert(list({raisePolyExpr(t1[1][i][0], 0, L[0]), t1[1][i][1]}));
  }

  return list({cnt, F});
}

Expr eval(Expr f, Expr L, Expr a, Int j = 0) {
	// TODO: remove
	if (a.size() == 0) {
    return f;
  }

  Expr k;

  Expr p = deepReplace(f, L[j], a[0]);

  for (Int i = j + 1; i < L.size(); i++) {
    k = deepReplace(p, L[i], a[i - j]);

    p = k;
  }

  k = reduceAST(p);

  return k;
}

Expr comp(Expr f, Expr x, Expr a) {
  Expr g = deepReplace(f, x, a);
  Expr t = reduceAST(g);

  return t;
}

Expr testEvaluationPoints(Expr U, Expr G, Expr F, Expr a, Expr L, Expr K) {
  assert(G.kind() == Kind::Integer, "Gamma parameter needs to be an integer");

  long i;

  Expr x, V, U0, g, delta, pr, lc, t1, R, E, d;

  x = L[0];

  R = rest(L);

  assert(R.size() == a.size(), "Wrong numbers of test points");

	// 	Test Wang condition 1: V(a1, ..., ar) = lc(f(x))(a1,...,ar) != 0
  V = leadCoeff(U, x);
  g = eval(V, R, a, 0);

  if (g == 0) {
    return fail();
  }

  // Test Wang condition 3: U0(x) = U(x, a1, ..., at) is square free
  U0 = eval(U, R, a, 0);

  if (!isSquareFree(U0, x, K)) {

    return fail();
  }

  // Test Wang condition 2: For each F[i], E[i] = F[i](a1, ..., ar)
  // has at least one prime division p[i] which does not divide
  // any E[j] j < i, Gamma, or the content of U0
  delta = cont(U0, L, K);
  pr = pp(U0, delta, L, K);

	lc = groundLeadCoeffPoly(pr, L);

  if (lc.value() < 0) {
    t1 = groundInvert(delta);

    delta = t1;

    t1 = groundInvert(pr);

    pr = t1;
  }

  E = list({});

  for (i = 0; i < F.size(); i++) {
    E.insert(eval(F[i], R, a, 0));
  }

  if (delta.value() == 0) {
    return fail();
  }

  d = nondivisors(G.value(), E, delta.value(), R, K);

  if (d.kind() == Kind::Fail) {

    return fail();
  }

  // return the content of U0, the primitive part of U0 and the E[i] = F[i](a1,
  // ... ar) paper defines delta = cn, and u(x) = pp(U0)
  return list({delta, pr, E});
}

Expr testEvaluationPointsPolyExpr(Expr& U, Expr& G, Expr& F, Expr& a, Expr& L, Expr& K) {
  assert(G.kind() == Kind::Integer, "Gamma parameter needs to be an integer");

  long i;

  Expr T = list({ L[0] });

  Expr R = rest(L);

  assert(R.size() == a.size(), "Wrong numbers of test points");
  // 	Test Wang condition 1: V(a1, ..., ar) = lc(f(x))(a1,...,ar) != 0
  Expr V = leadCoeffPolyExpr(U);
  Expr g = evalTailPolyExpr(V, R, a);

  if (isZeroPolyExpr(g)) return fail();

  // Test Wang condition 3: U0(x) = U(x, a1, ..., at) is square free
  Expr U0 = evalTailPolyExpr(U, R, a);

	if (!isSquareFreePolyExpr(U0, T, K)) return fail();

  // Test Wang condition 2: For each F[i], E[i] = F[i](a1, ..., ar)
  // has at least one prime division p[i] which does not divide
  // any E[j] j < i, Gamma, or the content of U0
	Expr CP = contAndPpPolyExpr(U0, T, K);

	Expr dt = CP[0];
  Expr pr = CP[1];

	Expr lc = groundLeadCoeffPolyExpr(pr);

  if (lc.value() < 0) {
    dt = groundInvertPolyExpr(dt);
    pr = groundInvertPolyExpr(pr);
  }

  Expr E = list({});

  for (i = 0; i < F.size(); i++) {
    E.insert(evalTailPolyExpr(F[i], R, a));
  }

	if (isZeroPolyExpr(dt)) return fail();

	Expr d = nondivisorsPolyExpr(G.value(), E, dt.value(), L, K);

  if (d == fail()) return fail();

  // return the content of U0, the primitive part of U0 and the E[i] = F[i](a1,
  // ... ar) paper defines delta = cn, and u(x) = pp(U0)
  return list({dt, pr, E});
}


Int degreeSum(Expr f, Expr L) {
  Expr n;

  Int s = 0;

  for (long i = 0; i < L.size(); i++) {
    n = degree(f, L[i]);

    if (n.kind() == Kind::Integer) {
      s += n.value();
    }
  }

  return s;
}


Int sumDegreesPolyExpr(Expr f, Expr L) {
  Expr n;

  Int s = 0;
  for (long i = 0; i < L.size(); i++) {
    n = degreePolyExpr(f, L[i]);
		s += n.value();
  }

  return s;
}

Int mignotteBound(Expr f, Expr L, Expr K) {
  Expr l = groundLeadCoeffPoly(f, L);

  Int a = norm(f, L, K);
  Int b = abs(l.value());
  Int n = degreeSum(f, L);

  return std::sqrt(n.longValue() + 1) * std::pow(2, n.longValue()) * a.longValue() * b.longValue();
}


Int mignotteBoundPolyExpr(Expr f, Expr L, Expr K) {
	Expr l = groundLeadCoeffPolyExpr(f);
  Int a = normPolyExpr(f, L, K);
  Int b = abs(l.value());
  Int n = sumDegreesPolyExpr(f, L);

  return std::sqrt(n.longValue() + 1) * std::pow(2, n.longValue()) * a.longValue() * b.longValue();
}

// return l, such that p^l is a bound to the coefficients of the factors of f in
// K[L]
Int mignoteExpoent(Expr f, Expr L, Expr K, Int p) {
  return std::ceil(std::log(2 * mignotteBound(f, L, K).longValue() + 1) /
                   std::log(p.longValue()));
}

// return l, such that p^l is a bound to the coefficients of the factors of f in
// K[L]
Int mignoteExpoentPolyExpr(Expr f, Expr L, Expr K, Int p) {
  return std::ceil(std::log(2 * mignotteBoundPolyExpr(f, L, K).longValue() + 1) /
                   std::log(p.longValue()));
}

Expr getEvaluationPoints(Expr& f, Expr& G, Expr& F, Expr& L, Expr K, Int p,
                         Expr c) {

  Expr t1, t2, t3;

  Int o = -1;
  Int r = -1;

  Int t = L.size() - 1;

  while (c.size() < 3) {
    for (t = 0; t < 5; t++) {
      Expr a = list({});

      for (size_t i = 0; i < L.size() - 1; i++) {
        a.insert(mod(random(), p, true));
      }

      Expr s = testEvaluationPoints(f, G, F, a, L, K);

      if (s == fail()) {
        continue;
			}

      Expr delta = s[0];
      Expr pr_u0 = s[1];
      Expr E = s[2];

      Expr pr = sqfFactors(pr_u0, L[0], K)[1];

      // Verify that the sets a[i] are
      // given same low r value
      o = pr.size();

      if (r == -1) {
        // first set
        r = o;

        c = set({list({
            delta, // paper delta
            pr_u0, // paper pr(U0)
            E,     // paper ~F[i]
            pr,    // paper u[i](x)*...*u[r](x)
            a      // paper a[i]
        })});
      } else {
        // reset current configuration if the found one leads to a
        // smaller number of factors
        if (o < r) {
          r = o;

          c = set({});
        }

        // If this config leads to the same number
        // of factors save it
        if (o == r) {
          t2 = list({
              delta, // paper delta
              pr_u0, // paper pr(U0)
              E,     // paper ~F[i]
              pr,    // paper u[i](x)*...*u[r](x)
              a      // paper a[i]
          });


          c = unification(c, set({t2}));
        }
      }

      if (c.size() < 3) {
        break;
      }
    }

    p = p + 3;
  }

  return c;
}

Expr getEvaluationPointsPolyExpr(Expr& f, Expr& G, Expr& F, Expr& L, Expr K, Int p,
                                 Expr c) {
  Int o = -1;
  Int r = -1;

  Int t = L.size() - 1;

	while (c.size() < 3) {
    for (t = 0; t < 5; t++) {
      Expr a = list({});

      for (size_t i = 0; i < L.size() - 1; i++) {
        a.insert(mod(random(), p, true));
      }

			Expr s = testEvaluationPointsPolyExpr(f, G, F, a, L, K);

      if (s == fail()) {
				continue;
			}

      Expr dt = s[0];
      Expr u0 = s[1];

      Expr E = s[2];

      Expr T = list({L[0]});
      Expr pr = sqfFactorsPolyExpr(u0, T, K)[1];

      o = pr.size();

      if (r == -1) {
        r = o;
        c = set({list({dt, u0, E, pr, a})});
      } else {
        if (o < r) {
          r = o;
          c = set({});
        }

        if (o == r) {
          c = unification(c, set({list({dt, u0, E, pr, a})}));
        }
      }

      if (c.size() < 3)
        break;
    }

    p = p + 3;
  }

  return c;
}

Expr wangLeadingCoeff(Expr f, Expr delta, Expr u, Expr F, Expr sF, Expr a,
                      Expr L, Expr K) {
	assert(K.identifier() == "Z", "only the integer Z field can be used");
  /**
   * From the Wang's paper:
   *
   * If none of u[1](x),...,u[r](x) is extraneous, then U factors into r
   * distinct irreductible polynomials U = prod i = 1 to r G[i](x[0], ...,
   * x[t]).
   *
   * Let C[i](x2, ..., x[t]) = lc(G[i]), ~C[i] = C[i](a[1], ..., a[t - 1]),
   * and G[i](x, a[1], ..., a[t - 1]) = delta[i] * u[i] where delta[i]
   * is some divisor of delta.
   *
   * Lemma: If there are no extraneous factors, then for all i and m, F[k]^m
   * divides C[i], then ~C[i] = ~F[1]^s1 * ... * F[k]^s[k]*w where w | G, s[i]
   * >= 0 and s[k] < m. Thus p[k]^m dows not divided ~C[i], which implies that
   * ~F[k]^m does not divide lc(u[i])*delta
   *
   * This lemma enables one to distribute all F[k] first, then all F[k-1], etc.
   * Thus D[i](x[2], ... ,x[t]) can be determined as products of powers of F[i],
   * Now let ~D[i] = D[i](a2, ..., a[t]). If delta == 1, then C[i] =
   * lc(u[i]/~D[i])*D[i]. Otherwise, if delta != 1, the following steps are
   * carried out for all i = 1, ..., r to correctly distribute the factors of
   * delta.
   *
   * 	1. Let d = gcd(lc(u[i]), ~D[i]) and C[i] = D[i] * lc(u[i]) / d
   *  2. Let u[i] = (~D[i] / d) * u[i]
   *  3. Let delta = delta / (D[i] / d)
   *
   * The process ends when delta = 1, Otherwise, let u[i] = delta * u[i], C[i] =
   * delta * C[i], u = delta^(r - 1) * U. In this case, when the true factors
   * over Z of U are found, they may have integer contents which should be
   * removed.
   *
   * In the above process, if any factors of V[n] is not distributed, then there
   * are extraneous factors, and the program goes back for different
   * substitutions that lower r.
   */

  bool extraneous = false;

  Int i, m, k;

  Int d, di, dt;

  Expr x, R, Di, D, ui, sFk, C, ci, Fk, lc, sDi, U, t1, t2;

  x = L[0];
  R = rest(L);

  D = list({});
  C = list({});
  U = list({});

  // 1. Distribute D[i]
  bool *was_set = new bool[sF.size()];

  for (i = 0; i < sF.size(); i++) {
    was_set[i.longValue()] = false;
  }

  for (i = 0; i < u.size(); i++) {
    ui = u[i];

    lc = leadCoeff(ui, x);

    Di = 1;

    /**
     * Aplying lemma: It there are no extraneous factors, then,
     * for all i and m, F[k]^m divides C[i] if and only if ~F[k]^m
     * divides lc(u[i])*delta
     *
     */
    d = lc.value() * delta.value();

    // Distribute F[k], then k - 1, ...

    Int size = sF.size();
    for (k = size - 1; k >= 0; k--) {
      m = 0;
      sFk = sF[k];

      // find expoent m
      while (d % sFk.value() == 0) {
        d = d / sFk.value();
        m = m + 1;
      }

      Fk = F[k];

      if (m != 0) {
        Di = Di * power(Fk, m);
        was_set[k.longValue()] = true;
      }
    }

    lc = reduceAST(Di);

    Di = lc;

    D.insert(Di);
  }

  extraneous = false;

  for (i = 0; i < sF.size(); i++) {
    if (!was_set[i.longValue()]) {
      // Extraneous factors found, stop and try again for another set of a[i]'s
      extraneous = true;
      break;
    }
  }

  if (extraneous) {
    return fail();
  }

  dt = delta.value();

  // otherwise, if delta != 1, the following steps
  // are carried out, for all i = 1, ..., r, to
  // correctly distribute the factors of delta
  for (i = 0; i < D.size(); i++) {
    ui = u[i];
    Di = D[i];

    lc = leadCoeff(ui, x);
    sDi = eval(Di, R, a, 0);

		assert(sDi.kind() == Kind::Integer, "integer expected");
    assert(lc.kind() == Kind::Integer, "integer expected");

    di = sDi.value();

    // if delta == 1, Then Ci = (lc(ui) / sD[i])*D[i]
    if (delta == 1) {
      ci = integer(lc.value() / sDi.value());
    } else {
      // * 	1. Let d = gcd(lc(u[i]), ~D[i]) and C[i] = D[i] * lc(u[i]) / d
      d = gcd(lc.value(), di);

      ci = integer(lc.value() / d);

      di = di / d;

      // *  2. Let u[i] = (~D[i] / d) * u[i]
      t2 = mulPoly(ui, di);

      ui = t2;

      // *  3. Let delta = delta / (D[i] / d)
      dt = dt / di;
    }

    // Ci = (lc(ui)/d) * Di = ci * Di
    t1 = mulPoly(ci, Di);

    ci = t1;

    C.insert(ci);
    U.insert(ui);
  }

  // Now if dt == 1, the process ends
  if (dt == 1) {
    return list({f, U, C});
  }

  // otherwise, let ui = delta*ui, Ci = delta*Ci, U = delta^(r - 1) * U
  for (i = 0; i < C.size(); i++) {
    U[i] = mulPoly(U[i], dt);
    C[i] = mulPoly(C[i], dt);
  }

  dt = pow(dt, u.size() - 1);
  t2 = mulPoly(f, dt);

  return list({t2, U, C});
}

Expr wangLeadingCoeffPolyExpr(Expr f, Expr delta, Expr u, Expr F, Expr sF,
                              Expr a, Expr L, Expr K) {
  assert(K.identifier() == "Z", "only the integer Z field can be used");
  /**
   * From the Wang's paper:
   *
   * If none of u[1](x),...,u[r](x) is extraneous, then U factors into r
   * distinct irreductible polynomials U = prod i = 1 to r G[i](x[0], ...,
   * x[t]).
   *
   * Let C[i](x2, ..., x[t]) = lc(G[i]), ~C[i] = C[i](a[1], ..., a[t - 1]),
   * and G[i](x, a[1], ..., a[t - 1]) = delta[i] * u[i] where delta[i]
   * is some divisor of delta.
   *
   * Lemma: If there are no extraneous factors, then for all i and m, F[k]^m
   * divides C[i], then ~C[i] = ~F[1]^s1 * ... * F[k]^s[k]*w where w | G, s[i]
   * >= 0 and s[k] < m. Thus p[k]^m dows not divided ~C[i], which implies that
   * ~F[k]^m does not divide lc(u[i])*delta
   *
   * This lemma enables one to distribute all F[k] first, then all F[k-1], etc.
   * Thus D[i](x[2], ... ,x[t]) can be determined as products of powers of F[i],
   * Now let ~D[i] = D[i](a2, ..., a[t]). If delta == 1, then C[i] =
   * lc(u[i]/~D[i])*D[i]. Otherwise, if delta != 1, the following steps are
   * carried out for all i = 1, ..., r to correctly distribute the factors of
   * delta.
   *
   * 	1. Let d = gcd(lc(u[i]), ~D[i]) and C[i] = D[i] * lc(u[i]) / d
   *  2. Let u[i] = (~D[i] / d) * u[i]
   *  3. Let delta = delta / (D[i] / d)
   *
   * The process ends when delta = 1, Otherwise, let u[i] = delta * u[i], C[i] =
   * delta * C[i], u = delta^(r - 1) * U. In this case, when the true factors
   * over Z of U are found, they may have integer contents which should be
   * removed.
   *
   * In the above process, if any factors of V[n] is not distributed, then there
   * are extraneous factors, and the program goes back for different
   * substitutions that lower r.
   */

  bool extraneous = false;

  Int i, m, k;

  Int d, di, dt;

  Expr  ci;

  // Expr x = L[0];
  Expr R = rest(L);

  Expr D = list({});
  Expr C = list({});
  Expr U = list({});

  // 1. Distribute D[i]
  bool *was_set = new bool[sF.size()];

  for (i = 0; i < sF.size(); i++) {
    was_set[i.longValue()] = false;
  }

  for (i = 0; i < u.size(); i++) {
    Expr ui = u[i];

    Expr lc = leadCoeffPolyExpr(ui);

		assert(lc.kind() == Kind::Integer, "lc should be and integer");

		Expr Di = polyExpr(1, R);

    /**
     * Aplying lemma: It there are no extraneous factors, then,
     * for all i and m, F[k]^m divides C[i] if and only if ~F[k]^m
     * divides lc(u[i])*delta
     *
     */
    d = lc.value() * delta.value();

    // Distribute F[k], then k - 1, ...

    Int size = sF.size();

		for (k = size - 1; k >= 0; k--) {
      m = 0;

      Expr sFk = sF[k];

      // find expoent m
      while (d % sFk.value() == 0) {
        d = d / sFk.value();
        m = m + 1;
      }

      if (m != 0) {
			  Expr t1 = powPolyExpr(F[k], m);

				Di = mulPolyExpr(Di, t1);

        was_set[k.longValue()] = true;
      }
    }

    D.insert(Di);
  }

  extraneous = false;

  for (i = 0; i < sF.size(); i++) {
    if (!was_set[i.longValue()]) {
      // Extraneous factors found, stop and try again for another set of a[i]'s
      extraneous = true;
      break;
    }
  }

  if (extraneous) {
    return fail();
  }

  dt = delta.value();

  // otherwise, if delta != 1, the following steps
  // are carried out, for all i = 1, ..., r, to
  // correctly distribute the factors of delta
  for (i = 0; i < D.size(); i++) {
    Expr ui = u[i];
    Expr Di = D[i];

    Expr lc = leadCoeffPolyExpr(ui);
    Expr sDi = evalTailPolyExpr(Di, R, a);

    di = sDi.value();

    // if delta == 1, Then Ci = (lc(ui) / sD[i])*D[i]
    if (delta == 1) {
      ci = lc.value() / sDi.value();
    } else {
      // * 	1. Let d = gcd(lc(u[i]), ~D[i]) and C[i] = D[i] * lc(u[i]) / d
      d = gcd(lc.value(), di);

      ci = lc.value() / d;

      di = di / d;

      // *  2. Let u[i] = (~D[i] / d) * u[i]
		  Expr dti = di;
			ui = mulPolyExpr(ui, dti);

      // *  3. Let delta = delta / (D[i] / d)
      dt = dt / di;
    }

		// Ci = (lc(ui)/d) * Di = ci * Di
    ci = mulPolyExpr(ci, Di);

    C.insert(ci);
    U.insert(ui);
  }

  // Now if dt == 1, the process ends
  if (dt == 1) {
    return list({f, U, C});
  }

	delta = dt;

	// otherwise, let ui = delta*ui, Ci = delta*Ci, U = delta^(r - 1) * U
  for (i = 0; i < C.size(); i++) {
    U[i] = mulPolyExpr(U[i], delta);
    C[i] = mulPolyExpr(C[i], delta);
  }

  delta = pow(dt, (Int)u.size() - 1);

	Expr t2 = mulPolyExpr(f, delta);

	return list({t2, U, C});
}

Expr EEAlift(Expr& a, Expr& b, Expr& x, Int p, Int k) {
  Int j, modulus;
  Expr t1, t2, t3, t4, t5, _sig, _tal, tal, sig;

  Expr amodp, bmodp, smodp, tmodp, s, t, G, e, c, mod, q;

  amodp = gf(a, p, true);
  bmodp = gf(b, p, true);

	G = extendedEuclidGf(amodp, bmodp, x, p);

  s = G[1];
  t = G[2];

  smodp = s;
  tmodp = t;

  modulus = p;

  for (j = 1; j <= k - 1; j++) {
    t1 = integer(1);
    t2 = mulPoly(s, a);
    t3 = mulPoly(t, b);
    t4 = subPoly(t1, t2);
    t5 = subPoly(t4, t3);

    e = reduceAST(t5);

    mod = integer(modulus);

    c = quoPolyGf(e, mod, x, p, true);

    _sig = mulPoly(smodp, c);
    _tal = mulPoly(tmodp, c);

    t3 = divideGPE(_sig, bmodp, x);

    q = t3[0];

    sig = t3[1];

    t3.remove(0L);
    t3.remove(0L);

    t3 = mulPoly(q, amodp);
    t5 = addPoly(_tal, t3);

    tal = gf(t5, p, true);

    t3 = mulPoly(sig, mod);
    t5 = addPoly(s, t3);

    s = t5;

    t3 = mulPoly(tal, mod);
    t5 = addPoly(t, t3);

    t = t5;

    modulus = modulus * p;
  }

  return list({s, t});
}


Expr EEAliftPolyExpr(Expr& a, Expr& b, Expr& L, Int p, Int k, Expr K) {
  Int j, modulus;

	Expr t1, t2, t3, t4, t5, _sig, _tal, tal, sig;

  Expr s, t, G, e, c, mod, q;

  Expr amodp = gfPolyExpr(a, p, true);
  Expr bmodp = gfPolyExpr(b, p, true);

  G = extendedEuclidPolyExprGf(amodp, bmodp, L, p);

	s = G[1];
  t = G[2];

  Expr smodp = s;
  Expr tmodp = t;

  modulus = p;

  for (j = 1; j <= k - 1; j++) {
		t1 = polyExpr(1, L);

    t2 = mulPolyExpr(s, a);
    t3 = mulPolyExpr(t, b);
    t4 = subPolyExpr(t1, t2);

		e = subPolyExpr(t4, t3);

    //mod = modulus;
		c = groundQuoPolyExprGf(e, modulus, p, true);
    //c = quoPolyExprGf(e, mod, L, p, true);

    _sig = mulPolyExpr(smodp, c);
    _tal = mulPolyExpr(tmodp, c);

    t3 = divPolyExpr(_sig, bmodp, L, K);

    q   = t3[0];
    sig = t3[1];

    t3 = mulPolyExpr(q, amodp);
    t5 = addPolyExpr(_tal, t3);

    tal = gfPolyExpr(t5, p, true);

    t3 = mulPolyExpr(sig, mod);
    s = addPolyExpr(s, t3);

    t3 = mulPolyExpr(tal, mod);

		t = addPolyExpr(t, t3);

    modulus = modulus * p;
  }

  return list({s, t});
}


Expr multiTermEEAlift(Expr& a, Expr& L, Int p, Int k) {
  long j, r;

  Expr q, t1, t2, s, bet, sig;

  r = a.size();

  q = list({a[r - 1]});

  for (j = r - 2; j >= 1; j--) {
    t1 = mulPoly(a[j], q[0L]);
    q.insert(reduceAST(t1), 0L);
  }

  bet = list({integer(1)});

  t2 = list({});
  s = list({});

  for (j = 0; j < r - 1; j++) {
    t1 = list({q[j], a[j]});

    sig = multivariateDiophant(t1, bet[bet.size() - 1], L, t2, 0, p, k);

    bet.insert(sig[0]);
    s.insert(sig[1]);

    sig.remove(1L);
    sig.remove(0L);
  }

  s.insert(bet[r - 1]);

  return s;
}



Expr multiTermEEAliftPolyExpr(Expr& a, Expr& L, Int p, Int k, Expr K) {
  long j, r;

  Expr q, t1, t2, s, bet, sig;

  r = a.size();

  q = list({a[r - 1]});

	for (j = r - 2; j >= 1; j--) {
    q.insert(mulPolyExpr(a[j], q[0]), 0);
  }

  bet = list({ polyExpr(1, L) });

  t2 = list({});
  s = list({});

  for (j = 0; j < r - 1; j++) {
    t1 = list({q[j], a[j]});

		sig = multivariateDiophantPolyExpr(t1, bet[bet.size() - 1], L, t2, 0, p, k, K);

		bet.insert(sig[0]);
    s.insert(sig[1]);
  }

  s.insert(bet[r - 1]);

  return s;
}


Expr replaceAndReduce(Expr f, Expr x, Expr a) {
  Expr g = deepReplace(f, x, a);
  Expr r = reduceAST(g);

  return r;
}

Expr diff(Expr f, Int j, Expr x) {
  Expr g = f;

  for (int i = 0; i < j; i++) {
    g = derivate(g, x);
  }

  return reduceAST(g);
}


Expr diffNthPolyExpr(Expr f, Int j, Expr x) {
  Expr g = f;

	for (int i = 0; i < j; i++) {
    g = diffPolyExpr(g, x);
  }

  return g;
}


Expr taylorExpansionCoeffAt(Expr f, Int m, Int j, Expr L, Expr a) {
  Expr g = diff(f, m, L[j]);
  Expr t = replaceAndReduce(g, L[j], a);

  Expr n = integer(fact(m));
  Expr K = symbol("Z");
  Expr q = recQuotient(t, n, L, K);

  return q;
}


Expr taylorExpansionCoeffAtPolyExpr(Expr f, Int m, Int j, Expr L, Expr a, Expr K) {
	assert(K.identifier() == "Z", "only Z is allowed");

	Expr g = diffNthPolyExpr(f, m, L[j]);

	Expr t = evalPolyExpr(g, L[j], a);

  Expr n = fact(m);

	return groundDivPolyExpr(t, n);
}

Expr multivariateDiophant(Expr& a, Expr& c, Expr& L, Expr& I, Int d, Int p, Int k) {
  long long i, j;

  Int m, v, r;

  Expr K, x1, ds, monomial, cm, e, sig, R, xv, av, A, t1, t2, t3, b, anew, Inew,
      cnew;

  // 1. Initialization
  r = a.size();
  v = 1 + I.size();

  if (v > 1) {
    K = symbol("Z");

    xv = L[L.size() - 1];
    av = I[I.size() - 1];

    // 2.1. Multivariate case
    A = 1;

    for (i = 0; i < r; i++) {
      t1 = mulPoly(A, a[i]);

      A = t1;
    }

    t1 = reduceAST(A);

    A = t1;

    b = list({});
    anew = list({});

    for (j = 0; j < r; j++) {
      b.insert(recQuotient(A, a[j], L, K));
    }

    for (j = 0; j < a.size(); j++) {
      anew.insert(replaceAndReduce(a[j], xv, av));
    }

		cnew = replaceAndReduce(c, xv, av);

		Inew = I;
    Inew.remove(Inew.size() - 1);

    R = L;
    R.remove(R.size() - 1);

    sig = multivariateDiophant(anew, cnew, R, Inew, d, p, k);

    t1 = integer(0);

    for (j = 0; j < sig.size(); j++) {
      t2 = mulPoly(sig[j], b[j]);
      t3 = addPoly(t1, t2);

      t1 = t3;
    }

    t2 = subPoly(c, t1);

    e = gf(t2, pow(p, k), true);

    monomial = integer(1);

    for (m = 1; m < d; m++) {
      if (e == 0)
        break;

      t1 = subPoly(xv, av);
      t2 = mulPoly(monomial, t1);

      monomial = t2;

      cm = taylorExpansionCoeffAt(e, m, v - 1, L, av);

      if (cm != 0) {
        ds = multivariateDiophant(anew, cm, L, Inew, d, p, k);

        for (j = 0; j < ds.size(); j++) {
          t1 = ds[j];
          ds.remove(j);

          t2 = mulPoly(t1, monomial);

          ds.insert(t2, j);
        }

        for (j = 0; j < ds.size(); j++) {
          t1 = ds[j];

          t2 = sig[j];

          sig.remove(j);

          t3 = addPoly(t1, t2);

          sig.insert(t3, j);
        }

        t1 = integer(0);
        for (j = 0; j < r; j++) {
          t2 = mulPoly(ds[j], b[j]);
          t3 = addPoly(t1, t2);

          t1 = t3;
        }

        t2 = subPoly(e, t1);

        e = gf(t2, pow(p, k), true);
      }
    }

  } else {
    x1 = L[0];

    sig = list({});
    for (j = 0; j < r; j++) {
      sig.insert(integer(0));
    }

    Expr C = c;

    while (C != 0) {
      t1 = degree(C, x1);
      m = t1.value();
      cm = leadCoeff(C, x1);
			ds = univariateDiophant(a, L, m, p, k);

      for (i = 0; i < ds.size(); i++) {
        t2 = mulPoly(ds[i], cm);

        t3 = addPolyGf(sig[i], t2, x1, pow(p, k), true);

        sig.remove(i);
        sig.insert(t3, i);
      }

      t1 = mul({cm, power(x1, m)});

      t2 = subPoly(C, t1);

      C = reduceAST(t2);
    }
  }

  for (j = 0; j < sig.size(); j++) {
    t2 = sig[j];

    sig.remove(j);

    sig.insert(gf(t2, pow(p, k), true), j);
  }

  return sig;
}

Expr multivariateDiophantPolyExpr(Expr& a, Expr& c, Expr& L, Expr& I, Int d, Int p, Int k, Expr K) {
  long long i, j;

  Int m, v, r;

  Expr ds, cm, t1, t2, t3, sig;

	// 1. Initialization
  r = a.size();
  v = 1 + I.size();

  if (v > 1) {
    Expr xv = L[L.size() - 1];
    Expr av = I[I.size() - 1];

    // 2.1. Multivariate case
    Expr A = 1;

    for (i = 0; i < r; i++) {
      A = mulPolyExpr(A, a[i]);
    }

    Expr b = list({});
    Expr anew = list({});

    for (j = 0; j < r; j++) {
			b.insert(quoPolyExpr(A, a[j], L, K));
    }

    for (j = 0; j < a.size(); j++) {
			anew.insert(evalPolyExpr(a[j], xv, av));
    }

    Expr cnew = evalPolyExpr(c, xv, av);

		Expr I_ = I;
    Expr L_ = L;

		I_.remove(I_.size() - 1);
		L_.remove(L_.size() - 1);

    sig = multivariateDiophantPolyExpr(anew, cnew, L_, I_, d, p, k, K);

		t1 = polyExpr(0, L);

    for (j = 0; j < sig.size(); j++) {
			sig[j] = insertSymbolPolyExpr(sig[j], L[L.size() - 1], 0, L.size() - 1);
		}

		for (j = 0; j < sig.size(); j++) {
      t2 = mulPolyExpr(sig[j], b[j]);
      t1 = addPolyExpr(t1, t2);
    }

    t2 = subPolyExpr(c, t1);

		Expr e = gfPolyExpr(t2, pow(p, k), true);
		Expr monomial = 1;

    for (m = 1; m < d; m++) {
      if (isZeroPolyExpr(e)) {
				break;
			}

			t1 = subPolyExpr(polyExpr(xv, L), polyExpr(av, L));

      monomial = mulPolyExpr(monomial, t1);

      cm = taylorExpansionCoeffAtPolyExpr(e, m, v - 1, L, av, K);

			if (!isZeroPolyExpr(cm)) {
        ds = multivariateDiophantPolyExpr(anew, cm, L_, I_, d, p, k, K);

        for (j = 0; j < ds.size(); j++) {
          ds[j] = mulPolyExpr(ds[j], monomial);
        }

        for (j = 0; j < ds.size(); j++) {
					sig[j] = addPolyExpr(ds[j], sig[j]);
        }

        t1 = polyExpr(0, L);

        for (j = 0; j < r; j++) {
          t2 = mulPolyExpr(ds[j], b[j]);
          t1 = addPolyExpr(t1, t2);
        }

        t2 = subPolyExpr(e, t1);

        e = gfPolyExpr(t2, pow(p, k), true);
      }
    }

  } else {
		sig = list({});

    for (j = 0; j < r; j++) {
      sig.insert(polyExpr(0, L));
    }

    Expr C = c;

    while (!isZeroPolyExpr(C)) {
			m = degreePolyExpr(C).value();
      cm = leadCoeffPolyExpr(C);

      ds = univariateDiophantPolyExpr(a, L, m, p, k, K);

      for (i = 0; i < ds.size(); i++) {
        t2 = mulPolyExpr(ds[i], cm);
        sig[i] = addPolyExprGf(sig[i], t2, pow(p, k), true);
      }

      t1 = polyExpr(cm * power(L[0], m), L);
			C = subPolyExpr(C, t1);
    }
  }

  for (j = 0; j < sig.size(); j++) {
    sig[j] =  gfPolyExpr(sig[j], pow(p, k), true);
  }

  return sig;
}


Expr univariateDiophant(Expr a, Expr L, Int m, Int p, Int k) {
  Expr x, s, t1, t2, t3, t4, result, u, v;

  x = L[0];

  long long r, j;

  r = a.size();

  result = list({});
  if (r > 2) {
    s = multiTermEEAlift(a, L, p, k);

    for (j = 0; j < r; j++) {
      t1 = power(x, integer(m));
      t2 = mulPoly(s[j], t1);

      result.insert(remPolyGf(t2, a[j], x, pow(p, k), true));
    }

  } else {
    s = EEAlift(a[1], a[0], x, p, k);

    t1 = power(x, integer(m));

    t2 = mulPoly(s[0], t1);

    t3 = divPolyGf(t2, a[0], x, pow(p, k), true);

    u = t3[0];
    v = t3[1];

    t3.remove(0L);
    t3.remove(0L);

    t1 = power(x, integer(m));
    t2 = mulPoly(s[1], t1);

    t3 = mulPoly(u, a[1]);
    t4 = addPolyGf(t2, t3, x, pow(p, k));

    result.insert(v);
    result.insert(t4);
  }

  return result;
}

Expr univariateDiophantPolyExpr(Expr a, Expr L, Int m, Int p, Int k, Expr K) {
  Int r, j;

  r = a.size();
	Expr result = list({});

  if (r > 2) {
		Expr s = multiTermEEAliftPolyExpr(a, L, p, k, K);

    for (j = 0; j < r; j++) {
			Expr t1 = polyExpr(power(L[0], m), L);
      Expr t2 = mulPolyExpr(s[j], t1);

      result.insert(remPolyExprGf(t2, a[j], L, pow(p, k), true));
    }

  } else {
		Expr s = EEAliftPolyExpr(a[1], a[0], L, p, k, K);

		Expr t1 = polyExpr(power(L[0], m), L);
    Expr t2 = mulPolyExpr(s[0], t1);
    Expr t3 = divPolyExprGf(t2, a[0], L, pow(p, k), true);

    Expr u = t3[0];
    Expr v = t3[1];

		t2 = mulPolyExpr(s[1], t1);

		t3 = mulPolyExpr(u, a[1]);

		t1 = addPolyExprGf(t2, t3, pow(p, k));

    result.insert(v);
    result.insert(t1);
  }

  return result;
}


Expr wangEEZ(Expr U, Expr u, Expr lc, Expr a, Int p, Expr L, Expr K) {
  Expr G, C, S, Ri, Y, s, ni, I, J, h, T, X, M, m, c, ti, Ui, ui, t1, t2, t3,
      lci, rij, xi, t4, t5, t6, t7, t8, t9, ai, si;

  Int r, i, j, k, t, w, z;

  // Compute U[i] where
  // U[i] = U(x,...,x[i], a[3], ..., a[t]);
  S = list({U});

  t = a.size() - 1;

  for (i = t; i >= 1; i--) {
    ai = a[i];
    xi = L[i + 1];
    si = comp(S[0], xi, ai);

    S.insert(gf(si, p, true), 0);
  }

  r = u.size();

  // Construct sequence of polynomials Rij
  for (j = 2; j <= t + 2; j++) {
    // Get list of expanded coefficients
    G = list({});
    for (i = 0; i < u.size(); i++) {
      G.insert(algebraicExpand(u[i]));
    }

    ai = a[j - 2];
    si = S[j - 2];

    I = list({});
    J = list({});

    for (k = 0; k < j - 2; k++)
      I.insert(a[k]);

    for (k = j - 1; k < a.size(); k++)
      J.insert(a[k]);

    for (i = 0; i < r; i++) {
      xi = L[0];

      t2 = eval(lc[i], L, J, j);
      t3 = groundGf(t2, p, true);
      // Replace leading coefficient by
      // pre computed coefficient

      t5 = leadCoeff(u[i], xi);
      t6 = degree(u[i], xi);

      t5 = t5 * power(xi, t6);
      t6 = t3 * power(xi, t6);

      // Compute R1
      t7 = algebraicExpand(u[i] - t5); // reduceAST(t4 - t5);
      t8 = algebraicExpand(t7 + t6);   // reduceAST(t7 + t6);

      u[i] = t8;
    }

    X = 1;

    for (k = 0; k < r; k++)
      X = mulPoly(X, algebraicExpand(u[k]));

    xi = L[j - 1];

    Ri = subPoly(si, X);
    ni = degree(Ri, xi);

    M = 1;

    m = xi + -ai;

    for (k = 0; k < ni.value(); k++) {
      // if R[m] = 0 mod(xk - ak)^(m+1), S[i][m+1] = R[i][m]
      if (Ri == 0) {
        break;
      }

      xi = L[j - 1];

      M = mulPoly(M, m);

      t1 = diff(Ri, k + 1, xi);
      C = comp(t1, xi, ai);
      if (C != 0) {
        Y = list({});

        for (z = 0; z <= j - 2; z++) {
          Y.insert(L[z]);
        }

        C = recQuotient(C, fact(k + 1), Y, K);
        T = multivariateDiophant(G, C, Y, I, ni.value(), p, 1);

        for (i = 0; i < r; i++) {
          ui = u[i];
          ti = T[i];

          t1 = mulPoly(ti, M);

          t2 = ui + t1;

          u[i] = groundGf(reduceAST(t2), p, true);

        }

        X = 1;

        for (i = 0; i < r; i++)
          X = mulPoly(X, algebraicExpand(u[i]));

        t1 = subPoly(si, X);
        Ri = gf(t1, p, true);
      }
    }
  }

  X = 1;

  for (k = 0; k < u.size(); k++) {
    X = mulPoly(X, algebraicExpand(u[k]));
  }

  if (X != U) {
    return fail();
  }

  return u;
}

Expr wangEEZPolyExpr(Expr U, Expr u, Expr lc, Expr a, Int p, Expr L, Expr K) {
  Expr G, C, S, Ri, Y, s, ni, I, J, h, T, X, M, m, c, ti, Ui, ui, t1, t2, t3,
      lci, rij, xi, t4, t5, t6, t7, t8, t9, ai, si;

  Int r, i, j, k, t, w, z;

  // Compute U[i] where
  // U[i] = U(x,...,x[i], a[3], ..., a[t]);
  S = list({ U });

  t = a.size() - 1;

  for (i = t; i >= 1; i--) {
    ai = a[i];
    xi = L[i + 1];

    si = evalPolyExpr(S[0], xi, ai);

    S.insert(gfPolyExpr(si, p, true), 0);
  }

  r = u.size();

  // Construct sequence of polynomials Rij
  for (j = 2; j <= t + 2; j++) {
    // Get list of expanded coefficients
    G = u;
    ai = a[j - 2];
    si = S[j - 2];

    I = list({});
    J = list({});

    for (k = 0; k < j - 2; k++)
      I.insert(a[k]);

    for (k = j - 1; k < a.size(); k++)
      J.insert(a[k]);

    for (i = 0; i < r; i++) {
			t2 = evalTailPolyExpr(lc[i], L, J, j);
			t3 = gfPolyExpr(t2, p, true);

      // Replace leading coefficient by
      // pre computed coefficient
			t6 = degreePolyExpr(u[i]);

			t5 = leadCoeffPolyExpr(u[i]);

			t5 = raisePolyExpr(t5, t6.value(), L[0]);

			u[i] = subPolyExpr(u[i], t5);

			u[i] = insertSymbolPolyExpr(u[i], L[j - 1], 0, j - 1);

			t5 = raisePolyExpr(t3, t6.value(), L[0]);

			u[i] = addPolyExpr(u[i], t5);
    }


		X = 1;

		for (k = 0; k < r; k++) {
			X = mulPolyExpr(X, u[k]);
		}

		xi = L[j - 1];


		Ri = subPolyExpr(si, X);
    ni = degreePolyExpr(Ri);

		Expr W = list({});

		for(Int z = 0; z <= j - 1; z++) {
			W.insert(L[z]);
		}

    M = polyExpr(1, W);
		m = polyExpr(xi + -1*ai.value(), W);

    for (k = 0; k < ni.value(); k++) {
      // if R[m] = 0 mod(xk - ak)^(m+1), S[i][m+1] = R[i][m]
      if (isZeroPolyExpr(Ri)) {
        break;
      }

      xi = L[j - 1];

      M = mulPolyExpr(M, m);
      t1 = diffNthPolyExpr(Ri, k + 1, xi);
      C = evalPolyExpr(t1, xi, ai);

      if (!isZeroPolyExpr(C)) {
        Y = list({});

				for (z = 0; z <= j - 2; z++) {
          Y.insert(L[z]);
        }

				Expr fk = fact(k + 1);

				C = groundDivPolyExpr(C, fk);

				T = multivariateDiophantPolyExpr(G, C, Y, I, ni.value(), p, 1, K);

				for (i = 0; i < r; i++) {
          ui = u[i];

					ti = insertSymbolPolyExpr(T[i], L[j - 1], 0, j - 1);
					t1 = mulPolyExpr(ti, M);
          t2 = addPolyExpr(ui, t1);

          u[i] = gfPolyExpr(t2, p, true);
        }

				Expr U = list({});

				for(i = 0; i < j; i++) {
					U.insert(L[i]);
				}

				X = polyExpr(1, U);

				for (i = 0; i < r; i++) {
          X = mulPolyExpr(X, u[i]);
				}

        t1 = subPolyExpr(si, X);
        Ri = gfPolyExpr(t1, p, true);
      }
    }
  }

  X = 1;

  for (k = 0; k < u.size(); k++) {
    X = mulPolyExpr(X, u[k]);
  }

  if (X != U) return fail();

  return u;
}



Expr factorsWangRec(Expr& f, Expr& L, Expr K, Int mod) {
  long long i = 0;
  long long j = 0;

  Int B = mignotteBound(f, L, K);

  long p = primes[0];

  while (p <= B) {
    p = primes[++i];
  }

  Int nrm1 = std::numeric_limits<long long>::min();
  Int nrm2 = std::numeric_limits<long long>::min();

  // First step: factor lc(f)
  Expr x = L[0];

  Expr lc = leadCoeff(f, x);

  Expr R = rest(L);

  Expr H = factors(lc, R, K);

  Expr G = H[0];

  Expr Vn = list({});

  for (i = 0; i < H[1].size(); i++) {
    Vn.insert(H[1][i][0]);
  }

  Expr a = list({});

  for (i = 0; i < L.size() - 1; i++) {
    a.insert(0);
  }

  Expr S = set({});

	// Test all zeros evaluation points
  Expr Q = testEvaluationPoints(f, G, Vn, a, L, K);

  if (Q != fail()) {

    H = sqfFactors(Q[1], L[0], K)[1];

    if (H.size() == 1)
      return list({f});

    S.insert(list({Q[0], Q[1], Q[2], H, a}));
  }

  // Second step: find integers a1, ..., ar
  S = getEvaluationPoints(f, G, Vn, L, K, mod, S);

  j = 0;

  for (i = 0; i < S.size(); i++) {
    nrm2 = norm(S[i][1], x);

    if (nrm2 > nrm1) {
      nrm1 = nrm2;
      j = i;
    }
  }
  Expr c = S[j];

  Expr delta = c[0];
  Expr pp_u0 = c[1];
  Expr sF = c[2];
  Expr u = c[3];

  a = c[4];

  Expr WLC = wangLeadingCoeff(f, delta, u, Vn, sF, a, L, K);


  if (WLC == fail()) {
    // try again
    return factorsWangRec(f, L, K, mod + 1);
  }

  Expr h = WLC[0];
  Expr U = WLC[1];
  Expr LC = WLC[2];

  Expr E = wangEEZ(f, U, LC, a, p, L, K);

  if (E == fail()) {
    return factorsWangRec(f, L, K, mod + 1);
  }


  Expr F = list({});

  for (i = 0; i < E.size(); i++) {
    Expr w = groundPPPoly(E[i], L, K);
    Expr l = groundLeadCoeffPoly(w, L);

    if (l.value() < 0) {
      w = reduceAST(reduceAST(mul({integer(-1), E[i]})));
    } else {
      w = E[i];
    }

    F.insert(w);
  }

  return E;
}




Expr factorsWangPolyExprRec(Expr& f, Expr& L, Expr K, Int mod) {
  // printf("\n\nWANG START\n\n");

  long long i = 0;
  long long j = 0;

  Int B = mignotteBoundPolyExpr(f, L, K);

  long p = primes[0];

  while (p <= B) {
    p = primes[++i];
  }

  Int nrm1 = std::numeric_limits<long long>::min();
  Int nrm2 = std::numeric_limits<long long>::min();

  // First step: factor lc(f)
  Expr lc = leadCoeffPolyExpr(f);

  Expr R = rest(L);

	Expr H = factorsPolyExpr(lc, R, K);

  Expr G = groundLeadCoeffPolyExpr(H[0]);

  Expr Vn = list({});

  for (i = 0; i < H[1].size(); i++) {
    Vn.insert(H[1][i][0]);
  }

  Expr a = list({});

  for (i = 0; i < L.size() - 1; i++) {
    a.insert(0);
  }

  Expr S = set({});

  // Test all zeros evaluation points
  Expr Q = testEvaluationPointsPolyExpr(f, G, Vn, a, L, K);


	if (Q != fail()) {
		H = sqfFactorsPolyExpr(Q[1], list({L[0]}), K)[1];

    if (H.size() == 1)
      return list({f});

    S.insert(list({Q[0], Q[1], Q[2], H, a}));
  }

	// Second step: find integers a1, ..., ar
  S = getEvaluationPointsPolyExpr(f, G, Vn, L, K, mod, S);


	j = 0;

  for (i = 0; i < S.size(); i++) {
    nrm2 = normPolyExpr(S[i][1]);

    if (nrm2 > nrm1) {
      nrm1 = nrm2;
      j = i;
    }
  }

	Expr c = S[j];

  Expr delta = c[0];
  Expr pp_u0 = c[1];
  Expr sF = c[2];
  Expr u = c[3];

  a = c[4];

	Expr WLC = wangLeadingCoeffPolyExpr(f, delta, u, Vn, sF, a, L, K);

  if (WLC == fail()) {
    return factorsWangPolyExprRec(f, L, K, mod + 1);
  }

  Expr h = WLC[0];
  Expr U = WLC[1];
  Expr LC = WLC[2];

	Expr E = wangEEZPolyExpr(f, U, LC, a, p, L, K);

	if (E == fail()) {
    return factorsWangPolyExprRec(f, L, K, mod + 1);
  }

  Expr F = list({});

  for (i = 0; i < E.size(); i++) {
    Expr w = groundPPPolyExpr(E[i]);
    Expr l = groundLeadCoeffPolyExpr(w);

    if (l.value() < 0) {
			w = groundInvertPolyExpr(E[i]);
    } else {
      w = E[i];
    }

    F.insert(w);
  }

  return E;
}

Expr factorsWang(Expr& f, Expr& L, Expr K) { return factorsWangRec(f, L, K, 3); }
Expr factorsWangPolyExpr(Expr& f, Expr& L, Expr K) { return factorsWangPolyExprRec(f, L, K, 3); }

} // namespace factorization
