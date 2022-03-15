#include "Zassenhaus.hpp"
#include "gauss/Algebra/Expression.hpp"
#include "Hensel.hpp"
#include "SquareFree.hpp"
#include "Utils.hpp"

#include "gauss/GaloisField/GaloisField.hpp"
#include "gauss/Calculus/Derivative.hpp"
#include "gauss/Primes/Primes.hpp"

#include <cmath>
#include <cstddef>

using namespace alg;
using namespace calc;
using namespace poly;
using namespace galoisField;

namespace factorization {

void subsetsRec(expr &arr, expr &data, list &s, Int start, Int end, Int index,
                Int r) {
  Int i, j;

  if (index == r) {
    s.insert(set({}));

    for (j = 0; j < r; j++) {
      s[s.size() - 1].insert(data[j]);
    }
  } else {
    for (i = start; i <= end && end - i + 1 >= r - index; i++) {
      data.insert(arr[i]);
      subsetsRec(arr, data, s, i + 1, end, index + 1, r);
      data.remove(data.size() - 1);
    }
  }
}

list subset(expr s, Int r) {
  long n = s.size();

  expr d;

  d = set({});

  list res = list({});

  subsetsRec(s, d, res, 0, n - 1, 0, r);

  return res;
}

// from Algorithms for Computer Algebra Geddes
// expr cantorZassenhausDDF(expr ax, expr x, Int p) {
//   long i;

//   expr wx, t, G, gx, n;

//   i = 1;

//   wx = x;

//   G = list({});

//   gx = 1;

//   n = degree(ax, x);

//   while (n != -inf() && n.value() >= 2 * i) {
//     wx = powModPolyGf(wx, ax, x, p, p, true);

//     t = subPolyGf(wx, x, x, p, true);
//     gx = gcdPolyGf(ax, t, x, p, true);

//     if (gx != 1) {
//       G.insert(list({gx, i}));

//       ax = quoPolyGf(ax, gx, x, p, true);

//       t = remPolyGf(wx, ax, x, p, true);

//       wx = t;
//     }

//     n = degree(ax, x);

//     i = i + 1;
//   };

//   if (ax != 1) {
//     G.insert(list({ax, degree(ax, x)}));
//   }

//   return G;
// }

// from Algorithms for Computer Algebra Geddes
expr cantorZassenhausDDFPolyExpr(expr ax, expr L, Int p) {
  long i;
  assert(L.kind() == kind::LIST && L.size() <= 1);

  expr wx, x, t, G, gx, n;

  i = 1;

  x = polyExpr(L[0], L);
  wx = polyExpr(L[0], L);

  G = list({});

  gx = 1;


	n = degreePolyExpr(ax);

  while (n != -inf() && n.value() >= 2 * i) {
    wx = powModPolyExprGf(wx, ax, L, p, p, true);

    t = subPolyExprGf(wx, x, p, true);

    gx = gcdPolyExprGf(ax, t, L, p, true);

    if (!isConstantPolyExpr(gx, 1)) {
      G.insert(list({gx, i}));

      ax = quoPolyExprGf(ax, gx, L, p, true);
      t = remPolyExprGf(wx, ax, L, p, true);

      wx = t;
    }

    n = degreePolyExpr(ax);

    i = i + 1;
  };

  if (!isConstantPolyExpr(ax, 1)) {
    G.insert(list({ax, degreePolyExpr(ax)}));
  }

  return G;
}

// from Algorithms for Computer Algebra Geddes
// expr cantorZassenhausEDF(expr a, expr x, Int n, Int p) {
//   Int m, i;

//   expr g, da, F, v, h, k, f1, f2, t;

//   da = degree(a, x);

//   if (da.value() <= n) {
//     return list({a});
//   }

//   m = da.value() / n;

//   F = list({a});


// 	while (F.size() < m) {
// 		v = randPolyGf(2 * n - 1, x, p);
//     if (p == 2) {
//       h = v;
//       for (i = 0; i < pow(2, n * m - 1); i++) {
//         // TODO change all trues to symmetric
//         h = powModPolyGf(h, a, x, 2, p, true);
//         v = addPolyGf(v, h, x, p, true);
//       }
//     } else {

//       v = powModPolyGf(v, a, x, (pow(p, n) - 1) / 2, p, true);
//       v = subPolyGf(v, 1, x, p, true);
//       g = gcdPolyGf(a, v, x, p, true);
//     }

//     g = gcdPolyGf(a, v, x, p, true);

//     if (g != 1 && g != a) {
//       k = quoPolyGf(a, g, x, p, true);

//       f1 = cantorZassenhausEDF(g, x, n, p);
//       f2 = cantorZassenhausEDF(k, x, n, p);

//       F = append(f1, f2);
//     }
//   }
//   return F;
// }

// from Algorithms for Computer Algebra Geddes
expr cantorZassenhausEDFPolyExpr(expr a, expr L, Int n, Int p) {
  assert(L.kind() == kind::LIST && L.size() <= 1);

  Int m, i;

  expr g, da, F, v, h, k, f1, f2, t, o;

  da = degreePolyExpr(a);

  if (da.value() <= n) {
    return list({a});
  }

  m = da.value() / n;

  F = list({a});

  o = polyExpr(1, L);

  while (F.size() < m) {
		v = randPolyExprGf(2 * n - 1, L, p);

    if (p == 2) {
      h = v;
      for (i = 0; i < pow(2, n * m - 1); i++) {
        h = powModPolyExprGf(h, a, 2, p, true);
        v = addPolyExprGf(v, h, p, true);
      }
    } else {
      v = powModPolyExprGf(v, a, L, (pow(p, n) - 1) / 2, p, true);

      v = subPolyExprGf(v, o, p, true);
      g = gcdPolyExprGf(a, v, L, p, true);
    }

    g = gcdPolyExprGf(a, v, L, p, true);

    if (!isConstantPolyExpr(g, 1) && g != a) {
      k = quoPolyExprGf(a, g, L, p, true);

      f1 = cantorZassenhausEDFPolyExpr(g, L, n, p);
      f2 = cantorZassenhausEDFPolyExpr(k, L, n, p);

      F = append(f1, f2);
    }
  }

  return F;
}

// expr cantorZassenhaus(expr f, expr x, Int m) {
//   expr U = monicPolyGf(f, x, m);

//   expr lc = U[0];
//   expr u = U[1];

//   expr n = degree(u, x);

//   if (n.value() == 0) {
//     return list({lc, list({})});
//   }

//   expr F = cantorZassenhausDDF(u, x, m);

//   expr P = list({});

//   for (Int i = 0; i < F.size(); i++) {
//     expr T = cantorZassenhausEDF(F[i][0], x, F[i][1].value(), m);

//     for (Int i = 0; i < T.size(); i++) {
//       P.insert(T[i]);
//     }
//   }

//   P = sortTerms(P);

//   return list({lc, P});
// }

expr cantorZassenhausPolyExpr(expr f, expr L, Int m) {
  assert(L.kind() == kind::LIST && L.size() <= 1);

  expr U = monicPolyExprGf(f, L, m);

  expr lc = U[0];
  expr u = U[1];

  expr n = degreePolyExpr(u);

  if (n.value() == 0) {
    return list({lc, list({})});
  }

  expr F = cantorZassenhausDDFPolyExpr(u, L, m);

  expr P = list({});

  for (Int i = 0; i < F.size(); i++) {
    expr T = cantorZassenhausEDFPolyExpr(F[i][0], L, F[i][1].value(), m);
    for (Int i = 0; i < T.size(); i++) {
      P.insert(T[i]);
    }
  }

  P = sortTerms(P);

  return list({lc, P});
}

// From modern computer algebra by Gathen
// expr zassenhaus(expr f, expr x, expr K) {
//   assert(K.identifier() == "Z");

//   bool stop = false;

//   Int d, s, i, j, l, p, A, B, C, gamma, gcd;

//   expr g, n, b, F, D, E, H, Z, G, T, S, M, u, v, gi, L, I;

//   L = list({x});

//   n = degree(f, x);

//   if (n == 1) {
//     return list({f});
//   }

//   A = norm(f, x);

//   d = n.value();

//   b = leadCoeff(f, x);

//   B = Int(std::abs(std::sqrt(d.doubleValue() + 1))) * pow(2, d) * A * b.value();

//   C = pow(d + 1, 2 * d) * pow(A, 2 * d - 1);

//   gamma = 2 * C.ceil_log2();

//   double y = gamma.doubleValue();
//   // choose a prime number p such that f be square free in Zp[x]
//   // and such that p dont divide lc(f)
//   for (size_t i = 1; primes[i] <= 2 * y * std::log(y); i++) {
//     p = primes[i];

//     if (b.value() % p == 0) {
//       continue;
//     }
//     F = gf(f, p, true);
//     D = derivate(F, x);

//     E = gf(D, p, true);
//     D = gcdPolyGf(F, E, x, p, false);

//     gcd = D.value();

//     if (b.value() % p > 0 && gcd == 1) {
//       break;
//     }
//   }

//   l = std::ceil((2 * B + 1).ceil_log2().doubleValue() /
//                 std::log(p.doubleValue()));

//   I = cantorZassenhaus(f, x, p);

//   Z = I[1];

//   g = multifactorHenselLifting(f, Z, x, p, l);

//   T = set({});

//   for (i = 0; i < g.size(); i++) {
//     T.insert(i);
//   }

//   F = list({});

//   s = 1;

//   while (2 * s <= T.size()) {
//     stop = false;

//     M = subset(T, s);
//     for (j = 0; j < M.size(); j++) {
//       S = M[j];

//       H = b;
//       G = b;

//       Z = difference(T, S);

//       for (i = 0; i < S.size(); i++) {
//         G = G * g[S[i].value()];
//       }

//       for (i = 0; i < Z.size(); i++) {
//         H = H * g[Z[i].value()];
//       }

//       u = gf(G, pow(p, l), true);
//       v = gf(H, pow(p, l), true);

//       G = u;
//       H = v;

//       if (norm(G, x) > pow(p, l) / 2) {
//         continue;
//       }

//       if (norm(H, x) > pow(p, l) / 2) {
//         continue;
//       }

//       if (l1norm(G, x) * l1norm(H, x) <= B) {
//         T = Z;

//         F.insert(pp(G, L, K));

//         f = pp(H, L, K);

//         b = leadCoeff(f, x);

//         stop = true;
//       }

//       if (stop)
//         break;
//     }

//     if (!stop)
//       s = s + 1;
//   }

//   F.insert(f);

//   return F;
// }

// From modern computer algebra by Gathen
expr zassenhausPolyExpr(expr f, expr L, expr K) {
  assert(L.kind() == kind::LIST && L.size() <= 1);
  assert(K.identifier() == "Z");

  bool stop = false;

  Int d, s, l, p, A, B, C, gamma, gcd;

  expr g, n, b, F, D, E, H, Z, G, T, S, u, v, gi, I;

	n = degreePolyExpr(f);

  if (n == 1) {
    return list({f});
  }

  A = normPolyExpr(f);

  d = n.value();

  b = leadCoeffPolyExpr(f);

  assert(b.kind() == kind::INT);

  B = Int(std::abs(std::sqrt(d.doubleValue() + 1))) * pow(Int(2), d) * A * b.value();

  C = pow(n.value() + 1, 2 * n.value()) * pow(A, 2 * n.value() - 1);

  gamma = 2 * C.ceil_log2();

  double y = gamma.doubleValue();
  // choose a prime number p such that f be square free in Zp[x]
  // and such that p dont divide lc(f)
  for (size_t i = 1; primes[i] <= 2 * y * std::log(y); i++) {
    p = primes[i];

    if (b.value() % p == 0) {
      continue;
    }

    F = gfPolyExpr(f, p, true);

    D = diffPolyExpr(F, L[0]);

    E = gfPolyExpr(D, p, true);

    D = gcdPolyExprGf(F, E, L, p, false);
    // gcd = D.value();
    if (b.value() % p > 0 && isConstantPolyExpr(D, 1)) {
      break;
    }
  }

  l = std::ceil((2 * B + 1).ceil_log2().doubleValue() /
                std::log(p.doubleValue()));

  I = cantorZassenhausPolyExpr(f, L, p);

  Z = I[1];

  g = multifactorHenselLiftingPolyExpr(f, Z, L, p, l);

  T = set({});

  for (size_t i = 0; i < g.size(); i++) {
    T.insert(integer(i));
  }

	F = list({});

  s = 1;

  while (2 * s <= T.size()) {
    stop = false;

    list M = subset(T, s);

		for (size_t j = 0; j < M.size(); j++) {
      S = M[j];

      H = polyExpr(b, L); // mul({ b });
      G = polyExpr(b, L); // mul({ b });

      Z = difference(T, S);

      for (size_t i = 0; i < S.size(); i++) {
        G = mulPolyExpr(G, g[S[i].value()]);
      }

      for (size_t i = 0; i < Z.size(); i++) {
        H = mulPolyExpr(H, g[Z[i].value()]);
      }

      G = gfPolyExpr(G, pow(p, l), true);
      H = gfPolyExpr(H, pow(p, l), true);

      if (normPolyExpr(G) > pow(p, l) / 2) {
        continue;
      }

      if (normPolyExpr(H) > pow(p, l) / 2) {
        continue;
      }

      if (l1normPolyExpr(G) * l1normPolyExpr(H) <= B) {
        T = Z;

        F.insert(ppPolyExpr(G, L, K));

        f = ppPolyExpr(H, L, K);

        b = leadCoeffPolyExpr(f);

        stop = true;
      }

      if (stop)
        break;
    }

    if (!stop)
      s = s + 1;
  }

  F.insert(f);

  return F;
}

} // namespace factorization
