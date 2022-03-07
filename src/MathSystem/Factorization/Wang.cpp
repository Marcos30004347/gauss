#include "Wang.hpp"
#include "Berlekamp.hpp"
#include "SquareFree.hpp"
#include "Utils.hpp"
#include "Zassenhaus.hpp"

#include "MathSystem/Primes/Primes.hpp"
#include "MathSystem/Calculus/Calculus.hpp"
#include "MathSystem/Algebra/Expression.hpp"
#include "MathSystem/Polynomial/Polynomial.hpp"
#include "MathSystem/GaloisField/GaloisField.hpp"

#include <cstddef>
#include <limits>

using namespace alg;
using namespace calc;
using namespace poly;
using namespace galoisField;

namespace factorization {

// expr nondivisors(Int G, expr F, Int d, expr L, expr K) {
//   assert(G != 0);
//   assert(d != 0);

//   long i, j, k;

//   Int q, r;

//   expr Fi;

//   k = F.size();

//   Int *x = new Int[k + 1];

//   x[0] = G * d;

//   for (i = 1; i <= k; i++) {
//     Fi = F[i - 1];

//     q = norm(Fi, L, K);

//     for (j = i - 1; j >= 0; j--) {
//       r = x[j];

//       while (abs(r) != 1) {
//         r = gcd(r, q);
//         q = q / r;
//       }

//       if (q == 1) {
//         return fail();
//       }
//     }

//     x[i] = q;
//   }

//   expr p = list({});

//   for (i = 1; i <= k; i++) {
//     p.insert(x[i]);
//   }

//   delete[] x;

//   return p;
// }

expr expandList(expr L) {
  expr K = list({});

  for (size_t i = 0; i < L.size(); i++) {
    K.insert(expand(L[i]));
  }

  return K;
}

expr nondivisorsPolyExpr(Int G, expr F, Int d, expr L, expr K) {
  assert(G != 0);
  assert(d != 0);

  long i, j, k;

	Int q, r;

  expr Fi;

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

				delete[] x;

        return fail();
      }
    }

    x[i] = q;
  }

  expr p = list({});

  for (i = 1; i <= k; i++) {
    p.insert(x[i]);
  }

  delete[] x;

  return p;
}

// expr trialDivision(expr& f, expr& F, expr& L, expr K) {
//   expr v, q, r, d, t = list({});

//   bool stop = false;

//   for (size_t i = 0; i < F.size(); i++) {
//     size_t k = 0;

//     stop = false;

//     v = expand(F[i]);

//     while (!stop) {
//       d = recPolyDiv(f, v, L, K);

//       q = d[0];
//       r = d[1];

//       if (r == 0) {
//         f = q;

//         k = k + 1;
//       }

//       stop = r != 0;
//     }

// 		if (k > 0) {
//       t.insert(list({F[i], integer(k)}));
//     }
//   }

//   return t;
// }

expr trialDivisionPolyExpr(expr& f, expr& F, expr& L, expr K) {
  expr t = list({});
  for (Int i = 0; i < F.size(); i++) {
    Int k = 0;

    bool stop = false;

    expr v = F[i];

    while (!stop) {
      expr d = divPolyExpr(f, v, L, K);

      expr q = d[0];
      expr r = d[1];

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

// expr sqfFactors(expr f, expr x, expr K) {

// 	expr n, cn, pr, lc, L, t1, F;

//   L = list({x});
//   cn = cont(f, L, K);
//   pr = pp(f, cn, L, K);
//   lc = leadCoeff(pr, x);

//   if (lc.value() < 0) {
//     cn = groundInvert(cn);
//     pr = groundInvert(pr);
//   }

//   n = degree(pr, x);

//   if (n.value() <= 0) {
//     return list({cn, list({})});
//   }

//   if (n.value() == 1) {
//     return list({cn, list({pr})});
//   }

//   F = zassenhaus(pr, x, K);

//   return list({cn, F});
// }



expr sqfFactorsPolyExpr(expr f, expr L, expr K) {

	assert(L.size() == 1);

  expr n, t1, F;

	expr CP = contAndPpPolyExpr(f, L, K);

	expr cn = CP[0];
  expr pr = CP[1];

  expr lc = leadCoeffPolyExpr(pr);

  assert(lc.kind() == kind::INT);

  if (lc.value() < 0) {
    cn = groundInvertPolyExpr(cn);
    pr = groundInvertPolyExpr(pr);
  }

  n = degreePolyExpr(pr);

  assert(n.kind() == kind::INT);

	if (n.value() <= 0) {
    return list({cn, list({})});
  }

  if (n.value() == 1) {
    return list({cn, list({pr})});
  }
  return list({ cn, zassenhausPolyExpr(pr, L, K) });
}



// expr factors(expr f, expr L, expr K) {

// 	if (f == 0) {
//     return list({0, list({})});
//   }

//   if (L.size() == 1) {
//     expr cnt = cont(f, L, K);
//     expr ppr = pp(f, cnt, L, K);

//     // TODO: maybe unnecessary
//     expr lc = leadCoeff(ppr, L[0]);
//     assert(lc.kind() == kind::INT);
//     if (lc.value() < 0) {
//       cnt = mulPoly(cnt, -1);
//       ppr = mulPoly(ppr, -1);
//     }

//     expr n = degree(ppr, L[0]);

//     if (n == 0) {
//       return list({cnt, list({})});
//     }

//     if (n == 1) {
//       return list({cnt, list({ppr, 1})});
//     }

//     expr g = squareFreePart(ppr, L, K);
//     expr H = zassenhaus(ppr, L[0], K);

//     expr F = trialDivision(ppr, H, L, K);

//     return list({cnt, F});
//   }

//   if (L.size() == 0) {
//     return list({f, list({})});
//   }

//   expr x = L[0];
//   expr R = rest(L);
//   expr cnt = groundContPoly(f, L, K);
//   expr prp = groundPPPoly(f, L, K);

//   expr lc = groundLeadCoeffPoly(prp, L);

//   // TODO: mayne unnecessary
//   if (lc.value() < 0) {
//     cnt = mulPoly(cnt, -1);
//     prp = mulPoly(prp, -1);
//   }

//   bool is_const = true;

//   // TODO: maybe just verify if prp is integer or fraction
//   for (Int i = 0; i < L.size() && is_const; i++) {
//     expr d = degree(prp, L[i]);
//     if (d != 0)
//       is_const = false;
//   }
//   if (is_const) {
//     return list({cnt, list({})});
//   }

//   expr G = cont(prp, L, K);
//   expr g = pp(prp, G, L, K);

//   expr F = list({});

//   expr n = degree(g, x);

//   if (n.value() > 0) {
// 		expr S = squareFreePart(g, L, K);
//     expr H = factorsWang(S[0], S[1], K);

//     F = trialDivision(f, H, L, K);
//   }

//   expr t1 = factors(G, R, K);

//   for (size_t i = 0; i < t1[1].size(); i++) {
//     F.insert(list({t1[1][i][0], t1[1][i][1]}));
//   }

//   return list({cnt, F});
// }

expr factorsPolyExpr(expr f, expr L, expr K) {
	if (isZeroPolyExpr(f)) {
    return list({f, list({})});
  }

  if (L.size() == 1) {
		expr CP = contAndPpPolyExpr(f, L, K);

		expr cnt = CP[0];
    expr ppr = CP[1];

    expr lc = leadCoeffPolyExpr(ppr);

    assert(lc.kind() == kind::INT);

		if (lc.value() < 0) {
			cnt = groundInvertPolyExpr(cnt);
			ppr = groundInvertPolyExpr(ppr);
    }


		cnt = raisePolyExpr(cnt, 0, L[0]);
    expr n = degreePolyExpr(ppr);

    assert(n.kind() == kind::INT);


		if (n == 0) {
			list k = list({});
			 list r = list({cnt, list({})});
			 return r;
    }

    if (n == 1) {
      return list({cnt, list({ppr, 1})});
    }

    expr g = squareFreePartPolyExpr(ppr, L, K);
    expr H = zassenhausPolyExpr(ppr, L, K);
    expr F = trialDivisionPolyExpr(ppr, H, L, K);

    return list({cnt, F});
  }


	expr R = rest(L);

	expr cnt = groundContPolyExpr(f);
  expr prp = groundPPPolyExpr(f);

  expr lc = groundLeadCoeffPolyExpr(prp);

	assert(lc.kind() == kind::INT);

  if (lc.value() < 0) {
		cnt = groundInvertPolyExpr(cnt);
		prp = groundInvertPolyExpr(prp);
  }

  bool is_const = true;

  // TODO: maybe just verify if prp is integer or fraction
  for (Int i = 0; i < L.size() && is_const; i++) {
    expr d = degreePolyExpr(prp);
    if (d != 0) is_const = false;
  }

  if (is_const) {
    return list({cnt, list({})});
  }

	expr O = contAndPpPolyExpr(prp, L, K);

	expr G = O[0];
  expr g = O[1];
  expr F = list({});

  expr n = degreePolyExpr(g);

  if (n.value() > 0) {
		expr S = squareFreePartPolyExpr(g, L, K);
    expr H = factorsWangPolyExpr(S[0], S[1], K);

    F = trialDivisionPolyExpr(f, H, L, K);
  }

  expr t1 = factorsPolyExpr(G, R, K);
  for (size_t i = 0; i < t1[1].size(); i++) {
    F.insert(list({raisePolyExpr(t1[1][i][0], 0, L[0]), t1[1][i][1]}));
  }

  return list({cnt, F});
}

expr eval(expr f, expr L, expr a, Int j = 0) {
	// TODO: remove

	if (a.size() == 0) {
    return f;
  }

  expr k;

  expr p = replace(f, L[j], a[0]);

  for (Int i = j + 1; i < L.size(); i++) {
    k = replace(p, L[i], a[i - j]);

    p = k;
  }

  k = reduce(p);

  return k;
}

expr comp(expr f, expr x, expr a) {

	expr g = replace(f, x, a);
  expr t = reduce(g);

  return t;
}

// expr testEvaluationPoints(expr U, expr G, expr F, expr a, expr L, expr K) {
//   assert(G.kind() == kind::INT);

//   expr x, V, U0, g, delta, pr, lc, t1, R, E, d;

//   x = L[0];

//   R = rest(L);

//   assert(R.size() == a.size());

// 	// 	Test Wang condition 1: V(a1, ..., ar) = lc(f(x))(a1,...,ar) != 0
//   V = leadCoeff(U, x);
//   g = eval(V, R, a, 0);

//   if (g == 0) {
//     return fail();
//   }

//   // Test Wang condition 3: U0(x) = U(x, a1, ..., at) is square free
//   U0 = eval(U, R, a, 0);

//   if (!isSquareFree(U0, x, K)) {

//     return fail();
//   }

//   // Test Wang condition 2: For each F[i], E[i] = F[i](a1, ..., ar)
//   // has at least one prime division p[i] which does not divide
//   // any E[j] j < i, Gamma, or the content of U0
//   delta = cont(U0, L, K);
//   pr = pp(U0, delta, L, K);

// 	lc = groundLeadCoeffPoly(pr, L);

//   if (lc.value() < 0) {
//     t1 = groundInvert(delta);

//     delta = t1;

//     t1 = groundInvert(pr);

//     pr = t1;
//   }

//   E = list({});

//   for (size_t i = 0; i < F.size(); i++) {
//     E.insert(eval(F[i], R, a, 0));
//   }

//   if (delta.value() == 0) {
//     return fail();
//   }

//   d = nondivisors(G.value(), E, delta.value(), R, K);

//   if (d.kind() == kind::FAIL) {
//     return fail();
//   }

//   // return the content of U0, the primitive part of U0 and the E[i] = F[i](a1,
//   // ... ar) paper defines delta = cn, and u(x) = pp(U0)
//   return list({delta, pr, E});
// }

expr testEvaluationPointsPolyExpr(expr& U, expr& G, expr& F, expr& a, expr& L, expr& K) {
  assert(G.kind() == kind::INT);

  expr T = list({ L[0] });

  expr R = rest(L);

  assert(R.size() == a.size());

  // 	Test Wang condition 1: V(a1, ..., ar) = lc(f(x))(a1,...,ar) != 0
  expr V = leadCoeffPolyExpr(U);
  expr g = evalTailPolyExpr(V, R, a);

  if (isZeroPolyExpr(g)) return fail();

  // Test Wang condition 3: U0(x) = U(x, a1, ..., at) is square free
  expr U0 = evalTailPolyExpr(U, R, a);

	if (!isSquareFreePolyExpr(U0, T, K)) return fail();

  // Test Wang condition 2: For each F[i], E[i] = F[i](a1, ..., ar)
  // has at least one prime division p[i] which does not divide
  // any E[j] j < i, Gamma, or the content of U0
	expr CP = contAndPpPolyExpr(U0, T, K);

	expr dt = CP[0];
  expr pr = CP[1];

	expr lc = groundLeadCoeffPolyExpr(pr);

  if (lc.value() < 0) {
    dt = groundInvertPolyExpr(dt);
    pr = groundInvertPolyExpr(pr);
  }

  expr E = list({});

  for (size_t i = 0; i < F.size(); i++) {
    E.insert(evalTailPolyExpr(F[i], R, a));
  }

	if (isZeroPolyExpr(dt)) return fail();

	expr d = nondivisorsPolyExpr(G.value(), E, dt.value(), L, K);

  if (d == fail()) return fail();

  // return the content of U0, the primitive part of U0 and the E[i] = F[i](a1,
  // ... ar) paper defines delta = cn, and u(x) = pp(U0)
  return list({dt, pr, E});
}


// Int degreeSum(expr f, expr L) {
//   expr n;

//   Int s = 0;

//   for (size_t i = 0; i < L.size(); i++) {
//     n = degree(f, L[i]);

//     if (n.kind() == kind::INT) {
//       s += n.value();
//     }
//   }

//   return s;
// }


Int sumDegreesPolyExpr(expr f, expr L) {
  expr n;

  Int s = 0;

  for (size_t i = 0; i < L.size(); i++) {
    n = degreePolyExpr(f, L[i]);
		s += n.value();
  }

  return s;
}

// Int mignotteBound(expr f, expr L, expr K) {
//   expr l = groundLeadCoeffPoly(f, L);

//   Int a = norm(f, L, K);
//   Int b = abs(l.value());
//   Int n = degreeSum(f, L);

//   return std::sqrt(n.longValue() + 1) * std::pow(2, n.longValue()) * a.longValue() * b.longValue();
// }


Int mignotteBoundPolyExpr(expr f, expr L, expr K) {

	expr l = groundLeadCoeffPolyExpr(f);
  Int a = normPolyExpr(f, L, K);
  Int b = abs(l.value());
  Int n = sumDegreesPolyExpr(f, L);

  return std::sqrt(n.longValue() + 1) * std::pow(2, n.longValue()) * a.longValue() * b.longValue();
}

// return l, such that p^l is a bound to the coefficients of the factors of f in
// K[L]
// Int mignoteExpoent(expr f, expr L, expr K, Int p) {
//   return std::ceil(std::log(2 * mignotteBound(f, L, K).longValue() + 1) /
//                    std::log(p.longValue()));
// }

// return l, such that p^l is a bound to the coefficients of the factors of f in
// K[L]
Int mignoteExpoentPolyExpr(expr f, expr L, expr K, Int p) {
  return std::ceil(std::log(2 * mignotteBoundPolyExpr(f, L, K).longValue() + 1) /
                   std::log(p.longValue()));
}

// expr getEvaluationPoints(expr& f, expr& G, expr& F, expr& L, expr K, Int p,
//                          expr c) {

//   expr t1, t2, t3;

//   Int o = -1;
//   Int r = -1;

//   Int t = L.size() - 1;

//   while (c.size() < 3) {
//     for (t = 0; t < 5; t++) {
//       expr a = list({});

//       for (size_t i = 0; i < L.size() - 1; i++) {
//         a.insert(mod(random(), p, true));
//       }

//       expr s = testEvaluationPoints(f, G, F, a, L, K);

//       if (s == fail()) {
//         continue;
// 			}

//       expr delta = s[0];
//       expr pr_u0 = s[1];
//       expr E = s[2];

//       expr pr = sqfFactors(pr_u0, L[0], K)[1];

//       // Verify that the sets a[i] are
//       // given same low r value
//       o = pr.size();

//       if (r == -1) {
//         // first set
//         r = o;

//         c = set({list({
//             delta, // paper delta
//             pr_u0, // paper pr(U0)
//             E,     // paper ~F[i]
//             pr,    // paper u[i](x)*...*u[r](x)
//             a      // paper a[i]
//         })});
//       } else {
//         // reset current configuration if the found one leads to a
//         // smaller number of factors
//         if (o < r) {
//           r = o;

//           c = set({});
//         }

//         // If this config leads to the same number
//         // of factors save it
//         if (o == r) {
//           t2 = list({
//               delta, // paper delta
//               pr_u0, // paper pr(U0)
//               E,     // paper ~F[i]
//               pr,    // paper u[i](x)*...*u[r](x)
//               a      // paper a[i]
//           });


//           c = unification(c, set({t2}));
//         }
//       }

//       if (c.size() < 3) {
//         break;
//       }
//     }

//     p = p + 3;
//   }

//   return c;
// }

expr getEvaluationPointsPolyExpr(expr& f, expr& G, expr& F, expr& L, expr K, Int p,
                                 expr c) {
  Int o = -1;
  Int r = -1;

  Int t = L.size() - 1;

	set points = {};

	while (c.size() < 3) {
    for (t = 0; t < 3; t++) {
			expr a = list({});

      for (size_t i = 0; i < L.size() - 1; i++) {
        a.insert(mod(random(), p, true));
      }

			if(exists(points, a)) {
				continue;
			}

			points.insert(a);

			expr s = testEvaluationPointsPolyExpr(f, G, F, a, L, K);

      if (s == fail()) {
				continue;
			}

      expr dt = s[0];
      expr u0 = s[1];

      expr E = s[2];

      expr T = list({L[0]});
      expr pr = sqfFactorsPolyExpr(u0, T, K)[1];

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

// expr wangLeadingCoeff(expr f, expr delta, expr u, expr F, expr sF, expr a,
//                       expr L, expr K) {
// 	assert(K.identifier() == "Z");

//   /**
//    * From the Wang's paper:
//    *
//    * If none of u[1](x),...,u[r](x) is extraneous, then U factors into r
//    * distinct irreductible polynomials U = prod i = 1 to r G[i](x[0], ...,
//    * x[t]).
//    *
//    * Let C[i](x2, ..., x[t]) = lc(G[i]), ~C[i] = C[i](a[1], ..., a[t - 1]),
//    * and G[i](x, a[1], ..., a[t - 1]) = delta[i] * u[i] where delta[i]
//    * is some divisor of delta.
//    *
//    * Lemma: If there are no extraneous factors, then for all i and m, F[k]^m
//    * divides C[i], then ~C[i] = ~F[1]^s1 * ... * F[k]^s[k]*w where w | G, s[i]
//    * >= 0 and s[k] < m. Thus p[k]^m dows not divided ~C[i], which implies that
//    * ~F[k]^m does not divide lc(u[i])*delta
//    *
//    * This lemma enables one to distribute all F[k] first, then all F[k-1], etc.
//    * Thus D[i](x[2], ... ,x[t]) can be determined as products of powers of F[i],
//    * Now let ~D[i] = D[i](a2, ..., a[t]). If delta == 1, then C[i] =
//    * lc(u[i]/~D[i])*D[i]. Otherwise, if delta != 1, the following steps are
//    * carried out for all i = 1, ..., r to correctly distribute the factors of
//    * delta.
//    *
//    * 	1. Let d = gcd(lc(u[i]), ~D[i]) and C[i] = D[i] * lc(u[i]) / d
//    *  2. Let u[i] = (~D[i] / d) * u[i]
//    *  3. Let delta = delta / (D[i] / d)
//    *
//    * The process ends when delta = 1, Otherwise, let u[i] = delta * u[i], C[i] =
//    * delta * C[i], u = delta^(r - 1) * U. In this case, when the true factors
//    * over Z of U are found, they may have integer contents which should be
//    * removed.
//    *
//    * In the above process, if any factors of V[n] is not distributed, then there
//    * are extraneous factors, and the program goes back for different
//    * substitutions that lower r.
//    */

//   bool extraneous = false;

//   Int i, m, k;

//   Int d, di, dt;

//   expr x, R, Di, D, ui, sFk, C, ci, Fk, lc, sDi, U, t1, t2;

//   x = L[0];
//   R = rest(L);

//   D = list({});
//   C = list({});
//   U = list({});

//   // 1. Distribute D[i]
//   bool *was_set = new bool[sF.size()];

//   for (i = 0; i < sF.size(); i++) {
//     was_set[i.longValue()] = false;
//   }

//   for (i = 0; i < u.size(); i++) {
//     ui = u[i];

//     lc = leadCoeff(ui, x);

//     Di = 1;

//     /**
//      * Aplying lemma: It there are no extraneous factors, then,
//      * for all i and m, F[k]^m divides C[i] if and only if ~F[k]^m
//      * divides lc(u[i])*delta
//      *
//      */
//     d = lc.value() * delta.value();

//     // Distribute F[k], then k - 1, ...

//     Int size = sF.size();
//     for (k = size - 1; k >= 0; k--) {
//       m = 0;
//       sFk = sF[k];

//       // find expoent m
//       while (d % sFk.value() == 0) {
//         d = d / sFk.value();
//         m = m + 1;
//       }

//       Fk = F[k];

//       if (m != 0) {
//         Di = Di * pow(Fk, m);
//         was_set[k.longValue()] = true;
//       }
//     }

//     lc = reduce(Di);

//     Di = lc;

//     D.insert(Di);
//   }

//   extraneous = false;

//   for (i = 0; i < sF.size(); i++) {
//     if (!was_set[i.longValue()]) {
//       // Extraneous factors found, stop and try again for another set of a[i]'s
//       extraneous = true;
//       break;
//     }
//   }

//   if (extraneous) {
//     return fail();
//   }

//   dt = delta.value();

//   // otherwise, if delta != 1, the following steps
//   // are carried out, for all i = 1, ..., r, to
//   // correctly distribute the factors of delta
//   for (i = 0; i < D.size(); i++) {
//     ui = u[i];
//     Di = D[i];

//     lc = leadCoeff(ui, x);
//     sDi = eval(Di, R, a, 0);

// 		assert(sDi.kind() == kind::INT);
//     assert(lc.kind() == kind::INT);

//     di = sDi.value();

//     // if delta == 1, Then Ci = (lc(ui) / sD[i])*D[i]
//     if (delta == 1) {
//       ci = integer(lc.value() / sDi.value());
//     } else {
//       // * 	1. Let d = gcd(lc(u[i]), ~D[i]) and C[i] = D[i] * lc(u[i]) / d
//       d = gcd(lc.value(), di);

//       ci = integer(lc.value() / d);

//       di = di / d;

//       // *  2. Let u[i] = (~D[i] / d) * u[i]
//       t2 = mulPoly(ui, di);

//       ui = t2;

//       // *  3. Let delta = delta / (D[i] / d)
//       dt = dt / di;
//     }

//     // Ci = (lc(ui)/d) * Di = ci * Di
//     t1 = mulPoly(ci, Di);

//     ci = t1;

//     C.insert(ci);
//     U.insert(ui);
//   }

//   // Now if dt == 1, the process ends
//   if (dt == 1) {
//     return list({f, U, C});
//   }

//   // otherwise, let ui = delta*ui, Ci = delta*Ci, U = delta^(r - 1) * U
//   for (i = 0; i < C.size(); i++) {
//     U[i] = mulPoly(U[i], dt);
//     C[i] = mulPoly(C[i], dt);
//   }

//   dt = pow(dt, u.size() - 1);
//   t2 = mulPoly(f, dt);

//   return list({t2, U, C});
// }

expr wangLeadingCoeffPolyExpr(expr f, expr delta, expr u, expr F, expr sF,
                              expr a, expr L, expr K) {
  assert(K.identifier() == "Z");

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

  expr  ci;

  // expr x = L[0];
  expr R = rest(L);

  expr D = list({});
  expr C = list({});
  expr U = list({});

  // 1. Distribute D[i]
  bool *was_set = new bool[sF.size()];

  for (i = 0; i < sF.size(); i++) {
    was_set[i.longValue()] = false;
  }

  for (i = 0; i < u.size(); i++) {
    expr ui = u[i];

    expr lc = leadCoeffPolyExpr(ui);

		assert(lc.kind() == kind::INT);

		expr Di = polyExpr(1, R);

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

      expr sFk = sF[k];

      // find expoent m
      while (d % sFk.value() == 0) {
        d = d / sFk.value();
        m = m + 1;
      }

      if (m != 0) {
			  expr t1 = powPolyExpr(F[k], m);

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

	delete[] was_set;

  if (extraneous) {
    return fail();
  }

  dt = delta.value();
  // otherwise, if delta != 1, the following steps
  // are carried out, for all i = 1, ..., r, to
  // correctly distribute the factors of delta
  for (i = 0; i < D.size(); i++) {
    expr ui = u[i];
    expr Di = D[i];

    expr lc = leadCoeffPolyExpr(ui);
    expr sDi = evalTailPolyExpr(Di, R, a);

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
		  expr dti = di;
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

	expr t2 = mulPolyExpr(f, delta);

	return list({t2, U, C});
}

// expr EEAlift(expr& a, expr& b, expr& x, Int p, Int k) {
//   Int j, modulus;
//   expr t1, t2, t3, t4, t5, _sig, _tal, tal, sig;

//   expr amodp, bmodp, smodp, tmodp, s, t, G, e, c, mod, q;

//   amodp = gf(a, p, true);
//   bmodp = gf(b, p, true);

// 	G = extendedEuclidGf(amodp, bmodp, x, p);

//   s = G[1];
//   t = G[2];

//   smodp = s;
//   tmodp = t;

//   modulus = p;

//   for (j = 1; j <= k - 1; j++) {
//     t1 = integer(1);
//     t2 = mulPoly(s, a);
//     t3 = mulPoly(t, b);
//     t4 = subPoly(t1, t2);
//     t5 = subPoly(t4, t3);

//     e = reduce(t5);

//     mod = integer(modulus);

//     c = quoPolyGf(e, mod, x, p, true);

//     _sig = mulPoly(smodp, c);
//     _tal = mulPoly(tmodp, c);

//     t3 = divideGPE(_sig, bmodp, x);

//     q = t3[0];

//     sig = t3[1];

//     t3.remove(0L);
//     t3.remove(0L);

//     t3 = mulPoly(q, amodp);
//     t5 = addPoly(_tal, t3);

//     tal = gf(t5, p, true);

//     t3 = mulPoly(sig, mod);
//     t5 = addPoly(s, t3);

//     s = t5;

//     t3 = mulPoly(tal, mod);
//     t5 = addPoly(t, t3);

//     t = t5;

//     modulus = modulus * p;
//   }

//   return list({s, t});
// }


expr EEAliftPolyExpr(expr& a, expr& b, expr& L, Int p, Int k, expr K) {
  Int j, modulus;

	expr t1, t2, t3, t4, t5, _sig, _tal, tal, sig;

  expr s, t, G, e, c, mod, q;

  expr amodp = gfPolyExpr(a, p, true);
  expr bmodp = gfPolyExpr(b, p, true);

  G = extendedEuclidPolyExprGf(amodp, bmodp, L, p);

	s = G[1];
  t = G[2];

  expr smodp = s;
  expr tmodp = t;

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


// expr multiTermEEAlift(expr& a, expr& L, Int p, Int k) {
//   long j, r;

//   expr q, t1, t2, s, bet, sig;

//   r = a.size();

//   q = list({a[r - 1]});

//   for (j = r - 2; j >= 1; j--) {
//     t1 = mulPoly(a[j], q[0L]);
//     q.insert(reduce(t1), 0L);
//   }

//   bet = list({integer(1)});

//   t2 = list({});
//   s = list({});

//   for (j = 0; j < r - 1; j++) {
//     t1 = list({q[j], a[j]});

//     sig = multivariateDiophant(t1, bet[bet.size() - 1], L, t2, 0, p, k);

//     bet.insert(sig[0]);
//     s.insert(sig[1]);

//     sig.remove(1L);
//     sig.remove(0L);
//   }

//   s.insert(bet[r - 1]);

//   return s;
// }



expr multiTermEEAliftPolyExpr(expr& a, expr& L, Int p, Int k, expr K) {
  long j, r;

  expr q, t1, t2, s, bet, sig;

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


expr replaceAndReduce(expr f, expr x, expr a) {
  expr g = replace(f, x, a);
  expr r = reduce(g);

  return r;
}

expr diff(expr f, Int j, expr x) {
  expr g = f;

  for (int i = 0; i < j; i++) {
    g = derivate(g, x);
  }

  return reduce(g);
}


expr diffNthPolyExpr(expr f, Int j, expr x) {
  expr g = f;

	for (int i = 0; i < j; i++) {
    g = diffPolyExpr(g, x);
  }

  return g;
}


// expr taylorExpansionCoeffAt(expr f, Int m, Int j, expr L, expr a) {
//   expr g = diff(f, m, L[j]);
//   expr t = replaceAndReduce(g, L[j], a);

//   expr n = integer(fact(m));
//   expr K = symbol("Z");
//   expr q = recQuotient(t, n, L, K);

//   return q;
// }


expr taylorExpansionCoeffAtPolyExpr(expr f, Int m, Int j, expr L, expr a, expr K) {
	assert(K.identifier() == "Z");

	expr g = diffNthPolyExpr(f, m, L[j]);

	expr t = evalPolyExpr(g, L[j], a);

  expr n = fact(m);

	return groundDivPolyExpr(t, n);
}

// expr multivariateDiophant(expr& a, expr& c, expr& L, expr& I, Int d, Int p, Int k) {
//   Int m, v, r;

//   expr K, x1, ds, monomial, cm, e, sig, R, xv, av, A, t1, t2, t3, b, anew, Inew,
//       cnew;

//   // 1. Initialization
//   r = a.size();
//   v = 1 + I.size();

//   if (v > 1) {
//     K = symbol("Z");

//     xv = L[L.size() - 1];
//     av = I[I.size() - 1];

//     // 2.1. Multivariate case
//     A = 1;

//     for (size_t i = 0; i < r; i++) {
//       t1 = mulPoly(A, a[i]);

//       A = t1;
//     }

//     t1 = reduce(A);

//     A = t1;

//     b = list({});
//     anew = list({});

//     for (size_t j = 0; j < r; j++) {
//       b.insert(recQuotient(A, a[j], L, K));
//     }

//     for (size_t j = 0; j < a.size(); j++) {
//       anew.insert(replaceAndReduce(a[j], xv, av));
//     }

// 		cnew = replaceAndReduce(c, xv, av);

// 		Inew = I;
//     Inew.remove(Inew.size() - 1);

//     R = L;
//     R.remove(R.size() - 1);

//     sig = multivariateDiophant(anew, cnew, R, Inew, d, p, k);

//     t1 = integer(0);

//     for (size_t j = 0; j < sig.size(); j++) {
//       t2 = mulPoly(sig[j], b[j]);
//       t3 = addPoly(t1, t2);

//       t1 = t3;
//     }

//     t2 = subPoly(c, t1);

//     e = gf(t2, pow(p, k), true);

//     monomial = integer(1);

//     for (m = 1; m < d; m++) {
//       if (e == 0)
//         break;

//       t1 = subPoly(xv, av);
//       t2 = mulPoly(monomial, t1);

//       monomial = t2;

//       cm = taylorExpansionCoeffAt(e, m, v - 1, L, av);

//       if (cm != 0) {
//         ds = multivariateDiophant(anew, cm, L, Inew, d, p, k);

//         for (size_t j = 0; j < ds.size(); j++) {
//           t1 = ds[j];
//           ds.remove(j);

//           t2 = mulPoly(t1, monomial);

//           ds.insert(t2, j);
//         }

//         for (size_t j = 0; j < ds.size(); j++) {
//           t1 = ds[j];

//           t2 = sig[j];

//           sig.remove(j);

//           t3 = addPoly(t1, t2);

//           sig.insert(t3, j);
//         }

//         t1 = integer(0);
//         for (size_t j = 0; j < r; j++) {
//           t2 = mulPoly(ds[j], b[j]);
//           t3 = addPoly(t1, t2);

//           t1 = t3;
//         }

//         t2 = subPoly(e, t1);

//         e = gf(t2, pow(p, k), true);
//       }
//     }

//   } else {
//     x1 = L[0];

//     sig = list({});

//     for (size_t j = 0; j < r; j++) {
//       sig.insert(integer(0));
//     }

//     expr C = c;

//     while (C != 0) {
//       t1 = degree(C, x1);
//       m = t1.value();
//       cm = leadCoeff(C, x1);
// 			ds = univariateDiophant(a, L, m, p, k);

//       for (size_t i = 0; i < ds.size(); i++) {
//         t2 = mulPoly(ds[i], cm);

//         t3 = addPolyGf(sig[i], t2, x1, pow(p, k), true);

//         sig.remove(i);
//         sig.insert(t3, i);
//       }

//       t1 = (cm * pow(x1, m));

//       t2 = subPoly(C, t1);

//       C = reduce(t2);
//     }
//   }

//   for (size_t j = 0; j < sig.size(); j++) {
//     t2 = sig[j];

//     sig.remove(j);

//     sig.insert(gf(t2, pow(p, k), true), j);
//   }

//   return sig;
// }

expr multivariateDiophantPolyExpr(expr& a, expr& c, expr& L, expr& I, Int d, Int p, Int k, expr K) {
  //long long i, j;

  Int m, v, r;

  expr ds, cm, t1, t2, t3, sig;

	// 1. Initialization
  r = a.size();
  v = 1 + I.size();

  if (v > 1) {
		// printf("========= A\n");
    expr xv = L[L.size() - 1];
    expr av = I[I.size() - 1];

    // 2.1. Multivariate case
    expr A = 1;

    for (size_t i = 0; i < r; i++) {
      A = mulPolyExpr(A, a[i]);
    }

    expr b = list({});
    expr anew = list({});

    for (size_t j = 0; j < r; j++) {
			b.insert(quoPolyExpr(A, a[j], L, K));
    }

    for (size_t j = 0; j < a.size(); j++) {
			anew.insert(evalPolyExpr(a[j], xv, av));
			// printf("eval %s for %s in %s = %s\n", to_string(a[j]).c_str(), to_string(xv).c_str(), to_string(av).c_str(), to_string(anew[j]).c_str());
    }

    expr cnew = evalPolyExpr(c, xv, av);

		expr I_ = I;
    expr L_ = L;

		I_.remove(I_.size() - 1);
		L_.remove(L_.size() - 1);
    sig = multivariateDiophantPolyExpr(anew, cnew, L_, I_, d, p, k, K);

		t1 = polyExpr(0, L);

    for (size_t j = 0; j < sig.size(); j++) {
			sig[j] = insertSymbolPolyExpr(sig[j], L[L.size() - 1], 0, L.size() - 1);
		}

		for (size_t j = 0; j < sig.size(); j++) {
      t2 = mulPolyExpr(sig[j], b[j]);
      t1 = addPolyExpr(t1, t2);
    }

    t2 = subPolyExpr(c, t1);

		expr e = gfPolyExpr(t2, pow(p, k), true);

		expr monomial = 1;

    for (m = 1; m < d; m++) {
      if (isZeroPolyExpr(e)) {
				break;
			}

			t1 = subPolyExpr(polyExpr(xv, L), polyExpr(av, L));

      monomial = mulPolyExpr(monomial, t1);

      cm = taylorExpansionCoeffAtPolyExpr(e, m, v - 1, L, av, K);

			if (!isZeroPolyExpr(cm)) {
        ds = multivariateDiophantPolyExpr(anew, cm, L_, I_, d, p, k, K);

        for (size_t j = 0; j < ds.size(); j++) {
          ds[j] = mulPolyExpr(ds[j], monomial);
        }

        for (size_t j = 0; j < ds.size(); j++) {
					sig[j] = addPolyExpr(ds[j], sig[j]);
        }

        t1 = polyExpr(0, L);

        for (size_t j = 0; j < r; j++) {
          t2 = mulPolyExpr(ds[j], b[j]);
          t1 = addPolyExpr(t1, t2);
        }

        t2 = subPolyExpr(e, t1);

        e = gfPolyExpr(t2, pow(p, k), true);
      }
    }

  } else {
		sig = list({});

    for (size_t j = 0; j < r; j++) {
      sig.insert(polyExpr(0, L));
    }

    expr C = c;

    while (!isZeroPolyExpr(C)) {
			m = degreePolyExpr(C).value();
      cm = leadCoeffPolyExpr(C);
      ds = univariateDiophantPolyExpr(a, L, m, p, k, K);
      for (size_t i = 0; i < ds.size(); i++) {
        t2 = mulPolyExpr(ds[i], cm);
        sig[i] = addPolyExprGf(sig[i], t2, pow(p, k), true);
      }

      t1 = polyExpr(cm * pow(L[0], m), L);
			C = subPolyExpr(C, t1);
    }
  }

  for (size_t j = 0; j < sig.size(); j++) {
    sig[j] =  gfPolyExpr(sig[j], pow(p, k), true);
  }

  return sig;
}


// expr univariateDiophant(expr a, expr L, Int m, Int p, Int k) {
//   expr x, s, t1, t2, t3, t4, result, u, v;

//   x = L[0];

//   long long r, j;

//   r = a.size();

//   result = list({});
//   if (r > 2) {
//     s = multiTermEEAlift(a, L, p, k);

//     for (j = 0; j < r; j++) {
//       t1 = pow(x, integer(m));
//       t2 = mulPoly(s[j], t1);

//       result.insert(remPolyGf(t2, a[j], x, pow(p, k), true));
//     }

//   } else {
//     s = EEAlift(a[1], a[0], x, p, k);

//     t1 = pow(x, integer(m));

//     t2 = mulPoly(s[0], t1);

//     t3 = divPolyGf(t2, a[0], x, pow(p, k), true);

//     u = t3[0];
//     v = t3[1];

//     t3.remove(0L);
//     t3.remove(0L);

//     t1 = pow(x, integer(m));
//     t2 = mulPoly(s[1], t1);

//     t3 = mulPoly(u, a[1]);
//     t4 = addPolyGf(t2, t3, x, pow(p, k));

//     result.insert(v);
//     result.insert(t4);
//   }

//   return result;
// }

expr univariateDiophantPolyExpr(expr a, expr L, Int m, Int p, Int k, expr K) {
  Int r, j;

  r = a.size();
	expr result = list({});

  if (r > 2) {
		expr s = multiTermEEAliftPolyExpr(a, L, p, k, K);

    for (j = 0; j < r; j++) {
			expr t1 = polyExpr(pow(L[0], m), L);
      expr t2 = mulPolyExpr(s[j], t1);

      result.insert(remPolyExprGf(t2, a[j], L, pow(p, k), true));
    }

  } else {
		expr s = EEAliftPolyExpr(a[1], a[0], L, p, k, K);
		expr t1 = polyExpr(pow(L[0], m), L);
    expr t2 = mulPolyExpr(s[0], t1);
    expr t3 = divPolyExprGf(t2, a[0], L, pow(p, k), true);

    expr u = t3[0];
    expr v = t3[1];

		t2 = mulPolyExpr(s[1], t1);

		t3 = mulPolyExpr(u, a[1]);

		t1 = addPolyExprGf(t2, t3, pow(p, k));

    result.insert(v);
    result.insert(t1);
  }

  return result;
}


// expr wangEEZ(expr U, expr u, expr lc, expr a, Int p, expr L, expr K) {
//   expr G, C, S, Ri, Y, s, ni, I, J, h, T, X, M, m, c, ti, Ui, ui, t1, t2, t3,
//       lci, rij, xi, t4, t5, t6, t7, t8, t9, ai, si;

//   Int r, i, j, k, t, w, z;

//   // Compute U[i] where
//   // U[i] = U(x,...,x[i], a[3], ..., a[t]);
//   S = list({U});

//   t = a.size() - 1;

//   for (i = t; i >= 1; i--) {
//     ai = a[i];
//     xi = L[i + 1];
//     si = comp(S[0], xi, ai);

//     S.insert(gf(si, p, true), 0);
//   }

//   r = u.size();

//   // Construct sequence of polynomials Rij
//   for (j = 2; j <= t + 2; j++) {
//     // Get list of expanded coefficients
//     G = list({});
//     for (i = 0; i < u.size(); i++) {
//       G.insert(expand(u[i]));
//     }

//     ai = a[j - 2];
//     si = S[j - 2];

//     I = list({});
//     J = list({});

//     for (k = 0; k < j - 2; k++)
//       I.insert(a[k]);

//     for (k = j - 1; k < a.size(); k++)
//       J.insert(a[k]);

//     for (i = 0; i < r; i++) {
//       xi = L[0];

//       t2 = eval(lc[i], L, J, j);
//       t3 = groundGf(t2, p, true);
//       // Replace leading coefficient by
//       // pre computed coefficient

//       t5 = leadCoeff(u[i], xi);
//       t6 = degree(u[i], xi);

//       t5 = t5 * pow(xi, t6);
//       t6 = t3 * pow(xi, t6);

//       // Compute R1
//       t7 = expand(u[i] - t5); // reduceAST(t4 - t5);
//       t8 = expand(t7 + t6);   // reduceAST(t7 + t6);

//       u[i] = t8;
//     }

//     X = 1;

//     for (k = 0; k < r; k++)
//       X = mulPoly(X, expand(u[k]));

//     xi = L[j - 1];

//     Ri = subPoly(si, X);
//     ni = degree(Ri, xi);

//     M = 1;

//     m = xi + -ai;

//     for (k = 0; k < ni.value(); k++) {
//       // if R[m] = 0 mod(xk - ak)^(m+1), S[i][m+1] = R[i][m]
//       if (Ri == 0) {
//         break;
//       }

//       xi = L[j - 1];

//       M = mulPoly(M, m);

//       t1 = diff(Ri, k + 1, xi);
//       C = comp(t1, xi, ai);
//       if (C != 0) {
//         Y = list({});

//         for (z = 0; z <= j - 2; z++) {
//           Y.insert(L[z]);
//         }

//         C = recQuotient(C, fact(k + 1), Y, K);
//         T = multivariateDiophant(G, C, Y, I, ni.value(), p, 1);

//         for (i = 0; i < r; i++) {
//           ui = u[i];
//           ti = T[i];

//           t1 = mulPoly(ti, M);

//           t2 = ui + t1;

//           u[i] = groundGf(reduce(t2), p, true);

//         }

//         X = 1;

//         for (i = 0; i < r; i++)
//           X = mulPoly(X, expand(u[i]));

//         t1 = subPoly(si, X);
//         Ri = gf(t1, p, true);
//       }
//     }
//   }

//   X = 1;

//   for (k = 0; k < u.size(); k++) {
//     X = mulPoly(X, expand(u[k]));
//   }

//   if (X != U) {
//     return fail();
//   }

//   return u;
// }

expr wangEEZPolyExpr(expr U, expr u, expr lc, expr a, Int p, expr L, expr K) {
  expr G, C, S, Ri, Y, s, ni, I, J, h, T, X, M, m, c, ti, Ui, ui, t1, t2, t3,
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

		if(isZeroPolyExpr(Ri)) {
			continue;
		}

		ni = degreePolyExpr(Ri);

		expr W = list({});

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

				expr fk = fact(k + 1);

				C = groundDivPolyExpr(C, fk);

				T = multivariateDiophantPolyExpr(G, C, Y, I, ni.value(), p, 1, K);

				for (i = 0; i < r; i++) {
          ui = u[i];

					ti = insertSymbolPolyExpr(T[i], L[j - 1], 0, j - 1);
					t1 = mulPolyExpr(ti, M);
          t2 = addPolyExpr(ui, t1);

          u[i] = gfPolyExpr(t2, p, true);
        }

				expr U = list({});

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



// expr factorsWangRec(expr& f, expr& L, expr K, Int mod) {
//   long long i = 0;
//   long long j = 0;

//   Int B = mignotteBound(f, L, K);

//   long p = primes[0];

//   while (p <= B) {
//     p = primes[++i];
//   }

//   Int nrm1 = std::numeric_limits<long long>::min();
//   Int nrm2 = std::numeric_limits<long long>::min();

//   // First step: factor lc(f)
//   expr x = L[0];

//   expr lc = leadCoeff(f, x);

//   expr R = rest(L);

//   expr H = factors(lc, R, K);

//   expr G = H[0];

//   expr Vn = list({});

//   for (size_t i = 0; i < H[1].size(); i++) {
//     Vn.insert(H[1][i][0]);
//   }

//   expr a = list({});

//   for (size_t i = 0; i < L.size() - 1; i++) {
//     a.insert(0);
//   }

//   expr S = set({});

// 	// Test all zeros evaluation points
//   expr Q = testEvaluationPoints(f, G, Vn, a, L, K);

//   if (Q != fail()) {

//     H = sqfFactors(Q[1], L[0], K)[1];

//     if (H.size() == 1)
//       return list({f});

//     S.insert(list({Q[0], Q[1], Q[2], H, a}));
//   }

//   // Second step: find integers a1, ..., ar
//   S = getEvaluationPoints(f, G, Vn, L, K, mod, S);

//   j = 0;

//   for (size_t i = 0; i < S.size(); i++) {
//     nrm2 = norm(S[i][1], x);

//     if (nrm2 > nrm1) {
//       nrm1 = nrm2;
//       j = i;
//     }
//   }
//   expr c = S[j];

//   expr delta = c[0];
//   expr pp_u0 = c[1];
//   expr sF = c[2];
//   expr u = c[3];

//   a = c[4];

//   expr WLC = wangLeadingCoeff(f, delta, u, Vn, sF, a, L, K);


//   if (WLC == fail()) {
//     // try again
//     return factorsWangRec(f, L, K, mod + 1);
//   }

//   expr h = WLC[0];
//   expr U = WLC[1];
//   expr LC = WLC[2];

//   expr E = wangEEZ(f, U, LC, a, p, L, K);

//   if (E == fail()) {
//     return factorsWangRec(f, L, K, mod + 1);
//   }


//   expr F = list({});

//   for (size_t i = 0; i < E.size(); i++) {
//     expr w = groundPPPoly(E[i], L, K);
//     expr l = groundLeadCoeffPoly(w, L);

//     if (l.value() < 0) {
//       w = reduce(-E[i]);
//     } else {
//       w = E[i];
//     }

//     F.insert(w);
//   }

//   return E;
// }

expr factorsWangPolyExprRec(expr& f, expr& L, expr K, Int mod) {
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
  expr lc = leadCoeffPolyExpr(f);

  expr R = rest(L);

	expr H = factorsPolyExpr(lc, R, K);
  expr G = groundLeadCoeffPolyExpr(H[0]);

  expr Vn = list({});

  for (size_t i = 0; i < H[1].size(); i++) {
    Vn.insert(H[1][i][0]);
  }

  expr a = list({});

  for (size_t i = 0; i < L.size() - 1; i++) {
    a.insert(0);
  }

  expr S = set({});

  // Test all zeros evaluation points
  expr Q = testEvaluationPointsPolyExpr(f, G, Vn, a, L, K);


	if (Q != fail()) {
		H = sqfFactorsPolyExpr(Q[1], list({L[0]}), K)[1];

    if (H.size() == 1)
      return list({f});

    S.insert(list({Q[0], Q[1], Q[2], H, a}));
  }

	// Second step: find integers a1, ..., ar
  S = getEvaluationPointsPolyExpr(f, G, Vn, L, K, mod, S);


	j = 0;

  for (size_t i = 0; i < S.size(); i++) {
    nrm2 = normPolyExpr(S[i][1]);

    if (nrm2 > nrm1) {
      nrm1 = nrm2;
      j = i;
    }
  }

	expr c = S[j];

  expr delta = c[0];
  expr pp_u0 = c[1];
  expr sF = c[2];
  expr u = c[3];

  a = c[4];

	expr WLC = wangLeadingCoeffPolyExpr(f, delta, u, Vn, sF, a, L, K);

  if (WLC == fail()) {
    return factorsWangPolyExprRec(f, L, K, mod + 1);
  }

  expr h = WLC[0];
  expr U = WLC[1];
  expr LC = WLC[2];

	expr E = wangEEZPolyExpr(f, U, LC, a, p, L, K);

	if (E == fail()) {
    return factorsWangPolyExprRec(f, L, K, mod + 1);
  }

  expr F = list({});

  for (size_t i = 0; i < E.size(); i++) {
    expr w = groundPPPolyExpr(E[i]);
    expr l = groundLeadCoeffPolyExpr(w);

    if (l.value() < 0) {
			w = groundInvertPolyExpr(E[i]);
    } else {
      w = E[i];
    }

    F.insert(w);
  }

  return E;
}

// expr factorsWang(expr& f, expr& L, expr K) { return factorsWangRec(f, L, K, 3); }
expr factorsWangPolyExpr(expr& f, expr& L, expr K) { return factorsWangPolyExprRec(f, L, K, 3); }

} // namespace factorization
