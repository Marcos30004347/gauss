#include "Zassenhaus.hpp"
#include "Berlekamp.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Hensel.hpp"
#include "SquareFree.hpp"
#include "Utils.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/GaloisField/GaloisField.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Simplification/Simplification.hpp"

#include <cmath>




// #include <chrono>


// std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
// std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();


using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace galoisField;
using namespace simplification;

namespace factorization {

void subsetsRec(Expr &arr, Expr &data, Expr &s, Int start, Int end, Int index,
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

Expr subset(Expr s, Int r) {
  long n = s.size();

  Expr d, res;

  d = set({});

  res = set({});

  subsetsRec(s, d, res, 0, n - 1, 0, r);

  return res;
}

// from Algorithms for Computer Algebra Geddes
Expr cantorZassenhausDDF(Expr ax, Expr x, Int p) {
  long i;

  Expr wx, t, G, gx, n;

  i = 1;

  wx = x;

  G = list({});

  gx = 1;

  n = degree(ax, x);

  while (n != -inf() && n.value() >= 2 * i) {
    wx = powModPolyGf(wx, ax, x, p, p, true);

    t = subPolyGf(wx, x, x, p, true);
    gx = gcdPolyGf(ax, t, x, p, true);

    if (gx != 1) {
      G.insert(list({gx, i}));

      ax = quoPolyGf(ax, gx, x, p, true);

      t = remPolyGf(wx, ax, x, p, true);

      wx = t;
    }

    n = degree(ax, x);

    i = i + 1;
  };

  if (ax != 1) {
    G.insert(list({ax, degree(ax, x)}));
  }

  return G;
}

// from Algorithms for Computer Algebra Geddes
Expr cantorZassenhausDDFPolyExpr(Expr ax, Expr L, Int p) {
  long i;
  assert(L.kind() == Kind::List && L.size() <= 1,
         "L should be a list with only one element");

  Expr wx, x, t, G, gx, n;

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
Expr cantorZassenhausEDF(Expr a, Expr x, Int n, Int p) {
  Int m, i;

  Expr g, da, F, v, h, k, f1, f2, t;

  da = degree(a, x);

  if (da.value() <= n) {
    return list({a});
  }

  m = da.value() / n;

  F = list({a});

  while (F.size() < m) {
    v = randPolyGf(2 * n - 1, x, p);
    if (p == 2) {
      h = v;
      for (i = 0; i < pow(2, n * m - 1); i++) {
        // TODO change all trues to symmetric
        h = powModPolyGf(h, a, x, 2, p, true);
        v = addPolyGf(v, h, x, p, true);
      }
    } else {
      v = powModPolyGf(v, a, x, (pow(p, n) - 1) / 2, p, true);
      v = subPolyGf(v, 1, x, p, true);
      g = gcdPolyGf(a, v, x, p, true);
    }

    g = gcdPolyGf(a, v, x, p, true);

    if (g != 1 && g != a) {
      k = quoPolyGf(a, g, x, p, true);

      f1 = cantorZassenhausEDF(g, x, n, p);
      f2 = cantorZassenhausEDF(k, x, n, p);

      F = append(f1, f2);
    }
  }

  return F;
}

// from Algorithms for Computer Algebra Geddes
Expr cantorZassenhausEDFPolyExpr(Expr a, Expr L, Int n, Int p) {
  assert(L.kind() == Kind::List && L.size() <= 1,
         "L should be a list with only one element");

  Int m, i;

  Expr g, da, F, v, h, k, f1, f2, t, o;

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

Expr cantorZassenhaus(Expr f, Expr x, Int m) {
  Expr U = monicPolyGf(f, x, m);

  Expr lc = U[0];
  Expr u = U[1];

  Expr n = degree(u, x);

  if (n.value() == 0) {
    return list({lc, list({})});
  }

	Expr F = cantorZassenhausDDF(u, x, m);

  Expr P = list({});

  for (Int i = 0; i < F.size(); i++) {
    Expr T = cantorZassenhausEDF(F[i][0], x, F[i][1].value(), m);

    for (Int i = 0; i < T.size(); i++) {
      P.insert(T[i]);
    }
  }

  P = sortTerms(P);

  return list({lc, P});
}

Expr cantorZassenhausPolyExpr(Expr f, Expr L, Int m) {
  assert(L.kind() == Kind::List && L.size() <= 1,
         "L should be a list with at most one element");

  Expr U = monicPolyExprGf(f, L, m);

  Expr lc = U[0];
  Expr u = U[1];

  Expr n = degreePolyExpr(u);

  if (n.value() == 0) {
    return list({lc, list({})});
  }

  Expr F = cantorZassenhausDDFPolyExpr(u, L, m);

  Expr P = list({});

  for (Int i = 0; i < F.size(); i++) {
    Expr T = cantorZassenhausEDFPolyExpr(F[i][0], L, F[i][1].value(), m);
    for (Int i = 0; i < T.size(); i++) {
      P.insert(T[i]);
    }
  }

  P = sortTerms(P);

  return list({lc, P});
}

// From modern computer algebra by Gathen
Expr zassenhaus(Expr f, Expr x, Expr K) {
  assert(K.identifier() == "Z", "");

  bool stop = false;

  Int s, i, j, l, p, A, B, C, gamma, gcd;

  Expr g, n, b, F, D, E, H, Z, G, T, S, M, u, v, gi, L, I;

  L = list({x});

  n = degree(f, x);

  if (n == 1) {
    return list({f});
  }

  A = norm(f, x);

  b = leadCoeff(f, x);

  B = Int(std::abs(std::sqrt(n.value().longValue() + 1))) *
      Int(pow(2, n.value())) * A * b.value();

  C = pow(n.value() + 1, 2 * n.value()) * pow(A, 2 * n.value() - 1);

  gamma = std::ceil(
      2 * (2 * n.value().longValue() * log2(n.value().longValue() + 1) +
           (2 * n.value().longValue() - 1) * log2(A.longValue())));

  // choose a prime number p such that f be square free in Zp[x]
  // and such that p dont divide lc(f)

  for (i = 1; primes[i.longValue()] <=
              2 * gamma.longValue() * std::log(gamma.longValue());
       i++) {
    p = primes[i.longValue()];

    if (b.value() % p == 0) {
      continue;
    }
    F = gf(f, p, true);

    D = derivate(F, x);

    E = gf(D, p, true);

    D = gcdPolyGf(F, E, x, p, false);

    gcd = D.value();

    if (b.value() % p > 0 && gcd == 1) {
      break;
    }
  }

  l = std::ceil(std::log(2 * B.longValue() + 1) / std::log(p.longValue()));

  I = cantorZassenhaus(f, x, p);

  Z = I[1];

  g = multifactorHenselLifting(f, Z, x, p, l);

  T = set({});

  for (i = 0; i < g.size(); i++) {
    T.insert(integer(i));
  }

  F = list({});

  s = 1;

  while (2 * s <= T.size()) {
    stop = false;

    M = subset(T, s);

    for (j = 0; j < M.size(); j++) {
      S = M[j];

      H = mul({b});
      G = mul({b});

      Z = difference(T, S);

      for (i = 0; i < S.size(); i++) {
        gi = g[S[i].value()];
        G.insert(gi);
      }

      for (i = 0; i < Z.size(); i++) {
        gi = g[Z[i].value()];
        H.insert(gi);
      }

      u = gf(G, pow(p, l), true);
      v = gf(H, pow(p, l), true);

      G = u;
      H = v;

      if (norm(G, x) > pow(p, l) / 2) {
        continue;
      }

      if (norm(H, x) > pow(p, l) / 2) {
        continue;
      }

      if (l1norm(G, x) * l1norm(H, x) <= B) {
        T = Z;

        F.insert(pp(G, L, K));

        f = pp(H, L, K);

        b = leadCoeff(f, x);

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

// From modern computer algebra by Gathen
Expr zassenhausPolyExpr(Expr f, Expr L, Expr K) {
  assert(L.kind() == Kind::List && L.size() <= 1,
         "L should be a list with at most one element");
  assert(K.identifier() == "Z", "");

  bool stop = false;

  Int s, i, j, l, p, A, B, C, gamma, gcd;

  Expr g, n, b, F, D, E, H, Z, G, T, S, M, u, v, gi, I;
  n = degreePolyExpr(f);

  if (n == 1) {
    return list({f});
  }

  A = normPolyExpr(f);

  b = leadCoeffPolyExpr(f);

  assert(b.kind() == Kind::Integer, "not a poly expression");

  B = Int(std::abs(std::sqrt(n.value().longValue() + 1))) *
      Int(pow(2, n.value())) * A * b.value();

  C = pow(n.value() + 1, 2 * n.value()) * pow(A, 2 * n.value() - 1);

  gamma = std::ceil(
      2 * (2 * n.value().longValue() * log2(n.value().longValue() + 1) +
           (2 * n.value().longValue() - 1) * log2(A.longValue())));

  // choose a prime number p such that f be square free in Zp[x]
  // and such that p dont divide lc(f)

  for (i = 1; primes[i.longValue()] <=
              2 * gamma.longValue() * std::log(gamma.longValue());
       i++) {
    p = primes[i.longValue()];

    if (b.value() % p == 0) {
      continue;
    }

    F = gfPolyExpr(f, p, true);

    D = diffPolyExpr(F, L[0]);

    E = gfPolyExpr(D, p, true);

    D = gcdPolyExprGf(F, E, L, p, false);

    // gcd = D.value();
    if (b.value() % p > 0 && gcd == 1) {
      break;
    }
  }
  l = std::ceil(std::log(2 * B.longValue() + 1) / std::log(p.longValue()));

  I = cantorZassenhausPolyExpr(f, L, p);
  Z = I[1];

  g = multifactorHenselLiftingPolyExpr(f, Z, L, p, l);

  T = set({});

  for (i = 0; i < g.size(); i++) {
    T.insert(i);
  }

  F = list({});

  s = 1;

  while (2 * s <= T.size()) {
    stop = false;

    M = subset(T, s);

    for (j = 0; j < M.size(); j++) {
      S = M[j];

      H = polyExpr(b, L); // mul({ b });
      G = polyExpr(b, L); // mul({ b });

      Z = difference(T, S);

      for (i = 0; i < S.size(); i++) {
        G = mulPolyExpr(G, g[S[i].value()]);
      }

      for (i = 0; i < Z.size(); i++) {
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
