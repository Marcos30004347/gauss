#include "Polynomial.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Exponential/Exponential.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Resultant.hpp"

#include <climits>
#include <cstddef>
#include <cstdio>
#include <map>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

using namespace ast;
using namespace expand;
using namespace simplification;
using namespace algebra;
using namespace calculus;

namespace polynomial {

Int collectDegree(Expr &u, Expr &x) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Fraction) {
    return 0;
  }

  if (u.kind() == Kind::Symbol) {
    if (u.identifier() == x.identifier()) {
      return 1;
    }

    return 0;
  }

  if (u.kind() == Kind::Power) {
    if (u[0] == x) {
      return u[1].value().longValue();
    }

    return 0;
  }

  Int d = 0;

  for (Int j = 0; j < u.size(); j++) {
    d = max(d, collectDegree(u[j], x));
  }

  return d;
}

Expr collectCoeff(Expr &u, Expr &x, Int d) {
  if (u == x && d == 1)
    return 1;
  if (u.kind() == x.kind()) {
    if (d == 0 && u != x) {
      return u;
    }

    return 0;
  }

  if (u.kind() == Kind::Power) {
    if (d == 0) {
      if (u[0] == x)
        return u[1] == 0 ? 1 : 0;
      return u;
    }
    return u[0] == x && u[1] == d ? 1 : 0;
  }

  if (u.kind() == Kind::Multiplication) {
    Expr c = Expr(Kind::Multiplication);
    bool f = 0;

    for (Int i = 0; i < u.size(); i++) {
      if (collectCoeff(u[i], x, d) == 1) {
        f = 1;
      } else {
        c.insert(u[i]);
      }
    }

    if (c.size() == 0)
      c = 1;

    if (d == 0 && f == false)
      return c;

    return f ? c : 0;
  }

  return d == 0 ? u : 0;
}

Expr collectRec(Expr &u, Expr &L, Int i) {
  if (i == L.size()) {
    if (i == 0 && u.kind() == Kind::Power) {
      return Expr(Kind::Addition, {1 * u});
    }

    if (i == 0 && u.kind() == Kind::Multiplication) {
      return Expr(Kind::Addition, {u});
    }

    if (i == 0 && u.kind() == Kind::Symbol) {
      return Expr(Kind::Addition, {power(u, 0)});
    }

    if (i == 0 && u.kind() == Kind::FunctionCall) {
      return Expr(Kind::Addition, {power(u, 0)});
    }
    return u;
  }

  if (u.kind() == Kind::Multiplication && u.size() == 2 &&
      u[1].kind() == Kind::Power && u[1][0] == L[i]) {
    return Expr(Kind::Addition, {u});
  }

  Expr c = 1;

  if (u.kind() == Kind::Power && u[0] == L[i]) {
    return Expr(Kind::Addition, {collectRec(c, L, i + 1) * u});
  }

  if (u.kind() == Kind::Symbol || u.kind() == Kind::FunctionCall) {
    if (u == L[i]) {
      return add({collectRec(c, L, i + 1) * power(L[i], 1)});
    }

    return add({collectRec(u, L, i + 1) * power(L[i], 0)});
  }

  Int d = collectDegree(u, L[i]);

  if (u.kind() == Kind::Multiplication) {
    Int k = collectDegree(u, L[i]);
    Expr c = collectCoeff(u, L[i], k);

    return Expr(Kind::Addition,
                {collectRec(c, L, i + 1) * power(L[i], Int(k))});
  }

  if (u.kind() == Kind::Addition) {
    std::map<Int, Expr> coeffs;

    for (long j = 0; j < u.size(); j++) {
      Int k = collectDegree(u[j], L[i]);
      Expr c = collectCoeff(u[j], L[i], k);
      // printf("%s %s %s\n", u[j].toString().c_str(), L[i].toString().c_str(),
      // k.to_string().c_str()); printf("%s\n", c.toString().c_str());
      if (c == 0)
        continue;

      if (coeffs.count(k) == 0)
        coeffs.insert(std::pair<Int, Expr>({k, c}));
      else
        coeffs[k] = coeffs[k] + c;
    }

    Expr g = Expr(Kind::Addition);

    for (std::map<Int, Expr>::iterator it = coeffs.begin(); it != coeffs.end();
         it++) {
      if (it->second != 0) {
        g.insert(collectRec(it->second, L, i + 1) * power(L[i], it->first));
      }
    }

    return g;
  }
  return Expr(Kind::Addition, {collectRec(u, L, i + 1) * power(L[i], d)});
}

Expr polyExpr(Expr &&u, Expr &&L) { return collectRec(u, L, 0); }
Expr polyExpr(Expr &u, Expr &L) { return collectRec(u, L, 0); }
Expr polyExpr(Expr &&u, Expr &L) { return collectRec(u, L, 0); }
Expr polyExpr(Expr &u, Expr &&L) { return collectRec(u, L, 0); }

void includeVariable(std::vector<Expr> &vars, Expr u) {
  bool included = false;

  for (Expr k : vars) {
    if (k == (u)) {
      included = true;
      break;
    }
  }

  if (!included) {
    vars.push_back(u);
  }
}

bool isGeneralMonomial(Expr u, Expr v) {
  Expr S;
  if (v.kind() != Kind::Set) {
    S = set({v});
  } else {
    S = v;
  }

  if (exists(S, u)) {

    return true;
  } else if (u.kind() == Kind::Power) {
    Expr b = u[0];
    Expr e = u[1];

    if (exists(S, b) && e.kind() == Kind::Integer && e.value() > 1) {

      return true;
    }
  } else if (u.kind() == Kind::Multiplication) {
    for (unsigned int i = 0; i < u.size(); i++) {
      if (isGeneralMonomial(u[i], S) == false) {

        return false;
      }
    }

    return true;
  }

  return u.freeOfElementsInSet(S);
}

bool isGerenalPolynomial(Expr u, Expr v) {
  Expr S;

  if (v.kind() != Kind::Set) {
    S = set({v});
  } else {
    S = v;
  }

  if (u.kind() != Kind::Addition && u.kind() != Kind::Subtraction) {
    bool r = isGeneralMonomial(u, S);

    return r;
  }

  if (exists(S, u)) {

    return true;
  }

  for (unsigned int i = 0; i < u.size(); i++) {
    if (isGeneralMonomial(u[i], S) == false) {

      return false;
    }
  }

  return true;
}

Expr coeffVarMonomial(Expr u, Expr S) {
  if (!isGeneralMonomial(u, S))
    return undefined();

  if (isConstant(u))
    return list({u, integer(1)});

  if (exists(S, u))
    return list({integer(1), u});

  if (u.kind() == Kind::Power && exists(S, u[0]))
    return list({integer(1), u});

  if (u.kind() == Kind::Multiplication) {
    Expr C = list({});
    Expr V = list({});

    for (unsigned int i = 0; i < u.size(); i++) {
      Expr L = coeffVarMonomial(u[i], S);

      Expr CL = list({L[0]});
      Expr VL = list({L[1]});

      Expr C_ = join(C, CL);
      Expr V_ = join(V, VL);

      C = C_;
      V = V_;
    }

    Expr coefs = mul({});
    Expr vars = mul({});

    for (unsigned int i = 0; i < C.size(); i++) {
      if (C[i].kind() == Kind::Integer && C[i].value() == 1)
        continue;

      coefs.insert(C[i]);
    }

    for (unsigned int i = 0; i < V.size(); i++) {
      if (V[i].kind() == Kind::Integer && V[i].value() == 1)
        continue;
      vars.insert(V[i]);
    }

    if (coefs.size() == 0) {

      coefs = integer(1);
    } else if (coefs.size() == 1) {
      Expr coefs_ = coefs[0];

      coefs = coefs_;
    }

    if (vars.size() == 0) {

      vars = integer(1);
    } else if (vars.size() == 1) {
      Expr vars_ = vars[0];

      vars = vars_;
    }

    return list({coefs, vars});
  }

  return list({u, integer(1)});
}

Expr collectTerms(Expr u, Expr S) {
  if (u.kind() != Kind::Addition) {
    Expr L = coeffVarMonomial(u, S);
    if (L.kind() == Kind::Undefined) {

      return undefined();
    }

    return u;
  }

  if (exists(S, u)) {
    return u;
  }

  int N = 0;

  Expr T = list({});

  for (unsigned int i = 0; i < u.size(); i++) {
    Expr f = coeffVarMonomial(u[i], S);

    if (f.kind() == Kind::Undefined) {

      return undefined();
    }

    int j = 1;
    bool combined = false;

    while (!combined && j <= N) {
      int j_ = j - 1;

      if (f[1] == (T[j_][1])) {

        Expr Tj = list({add({T[j_][0], f[0]}), f[1]});

        Expr Tj_ = T[j_];
        T.remove(j_);

        T.insert(Tj, j_);

        combined = true;
      }

      j = j + 1;
    }

    if (!combined) {
      T.insert(f, N);
      N = N + 1;
    }
  }

  Expr v = add({});

  for (int j = 0; j < N; j++) {
    if (T[j][1].kind() == Kind::Integer && T[j][1].value() == 1) {
      v.insert(T[j][0]);
    } else {
      v.insert(mul({
          T[j][0],
          T[j][1],
      }));
    }
  }

  if (v.size() == 0) {

    return integer(0);
  }

  if (v.size() == 1) {
    Expr v_ = v[0];

    v = v_;
  }

  return v;
}

Expr degreeGME(Expr u, Expr v) {
  if (u.kind() == Kind::Integer && u.value() == 0)
    return Expr(Kind::MinusInfinity);

  if (isConstant(u))
    return integer(0);

  Expr S;
  if (v.kind() != Kind::Set) {
    S = set({v});
  } else {
    S = v;
  }

  if (exists(S, u)) {
    return 1;
  } else if (u.kind() == Kind::Power) {
    Expr b = u[0];
    Expr e = u[1];

    if (exists(S, b) && isConstant(e)) {
      return e;
    }

  } else if (u.kind() == Kind::Multiplication) {
    Expr deg = integer(0);
    for (unsigned int i = 0; i < u.size(); i++) {
      Expr deg_ = degreeGME(u[i], S);
      if (deg_.value() > deg.value()) {

        deg = deg_;
      } else {
      }
    }

    return deg;
  }

  return 0;
}

Expr degree(Expr u, Expr v) {
  Expr S;

  if (u.kind() == Kind::Integer && u.value() == 0) {
    return Expr(Kind::MinusInfinity);
  }

  if (v.kind() != Kind::Set) {
    S = set({v});
  } else {
    S = v;
  }

  if (u.kind() != Kind::Addition && u.kind() != Kind::Subtraction) {
    Expr r = degreeGME(u, S);

    return r;
  }

  if (exists(S, u)) {

    return integer(1);
  }

  Expr deg = integer(0);

  for (unsigned int i = 0; i < u.size(); i++) {
    Expr deg_ = degreeGME(u[i], S);

    if (deg_.value() > deg.value()) {
      deg = deg_;
    }
  }

  return deg;
}

Expr variablesRec(Expr u) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Fraction)
    return set({});

  if (u.kind() == Kind::Power) {
    Expr b = u[0];
    Expr e = u[1];

    if (e.kind() == Kind::Integer && e.value() > 1)
      return set({b});

    return set({u});
  }

  if (u.kind() == Kind::Addition || u.kind() == Kind::Subtraction ||
      u.kind() == Kind::Multiplication) {
    Expr S = set({});

    for (unsigned int i = 0; i < u.size(); i++) {

      Expr S_ = variablesRec(u[i]);
      Expr S__ = unification(S, S_);

      S = S__;
    }

    return S;
  }

  return set({u});
}

Expr variables(Expr u) {
  Expr v = variablesRec(u);
  Expr t = list({});

  for (Int i = 0; i < v.size(); i++) {
    t.insert(list({v[i], degree(u, v[i])}));
  }

  Expr a = list({});

  // TODO: optimize sorting
  for (Int i = 0; i < t.size(); i++) {
    for (Int j = i + 1; j < t.size(); j++) {
      if (t[i][1].value() < t[j][1].value()) {
        Expr tmp = t[i];
        t[i] = t[j];
        t[j] = tmp;
      } else if (t[i][1] == t[j][1]) {
        if (t[i][0].kind() == Kind::FunctionCall &&
            t[j][0].kind() != Kind::FunctionCall) {
          Expr tmp = t[i];
          t[i] = t[j];
          t[j] = tmp;
        }
      }
    }
  }

  for (Int i = 0; i < t.size(); i++) {
    a.insert(t[i][0]);
  }

  return a;
}
Expr coefficientGME(Expr u, Expr x) {
  if (u == (x)) {
    return list({1, 1});
  }

  if (u.kind() == Kind::Power) {
    Expr b = u[0];
    Expr e = u[1];

    if (b == (x) && e.kind() == Kind::Integer && e.value() > 0) {
      return list({1, e});
    }

  } else if (u.kind() == Kind::Multiplication) {
    Expr m = 0;
    Expr c = u;

    for (unsigned int i = 0; i < u.size(); i++) {
      Expr f = coefficientGME(u[i], x);

      if (f.kind() == Kind::Undefined) {

        return undefined();
      }
      if (f[1].kind() != Kind::Integer || f[1].value() != 0) {

        m = f[1];
        Expr c_ = div(u, power(x, m));
        c = algebraicExpand(c_);
      }
    }

    return list({c, m});
  }

  if (u.freeOf(x)) {
    return list({u, 0});
  }

  return undefined();
}

Expr coeff(Expr u, Expr x, Expr j) {

  if (u.kind() != Kind::Addition && u.kind() != Kind::Subtraction) {
    Expr f = coefficientGME(u, x);

    if (f.kind() == Kind::Undefined)
      return f;

    if (j == (f[1])) {
      Expr k = f[0];

      return k;
    }

    return integer(0);
  }

  if (x == (u)) {
    if (j.kind() == Kind::Integer && j.value() == 1) {
      return integer(1);
    }

    return integer(0);
  }

  Expr c = integer(0);

  for (unsigned int i = 0; i < u.size(); i++) {
    Expr f = coefficientGME(u[i], x);

    if (f.kind() == Kind::Undefined)
      return f;

    if (j == (f[1])) {
      Expr k = f[0];

      if (c.kind() == Kind::Integer && c.value() == 0) {

        c = Expr(u.kind());
        c.insert(k);
      } else {
        c.insert(k);
      }
    }
  }

  if (c.kind() != Kind::Integer && c.size() == 1) {
    Expr l = c[0];

    return l;
  }

  return c;
}

Expr leadCoeff(Expr u, Expr x) {
  Expr d = degree(u, x);
  Expr lc = coeff(u, x, d);

  return lc;
}

Expr divideGPE(Expr u, Expr v, Expr x) {
  Expr t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

  Expr q = integer(0);
  Expr r = u;

  Expr m = degree(r, x);
  Expr n = degree(v, x);

  Expr lcv = leadCoeff(v, x);

  while (m.kind() != Kind::MinusInfinity &&
         (m.kind() == Kind::Integer && n.kind() == Kind::Integer &&
          m.value() >= n.value())) {
    Expr lcr = leadCoeff(r, x);

    Expr s = div(lcr, lcv);

    t1 = power(x, m - n);
    t2 = mulPoly(s, t1);
    t3 = addPoly(q, t2);

    q = reduceAST(t3);

    t1 = power(x, m);
    t2 = mulPoly(lcr, t1);
    t3 = subPoly(r, t2);

    t4 = power(x, n);
    t5 = mulPoly(lcv, t4);
    t6 = subPoly(v, t5);

    t7 = mulPoly(t6, s);
    t8 = power(x, sub({m, n}));
    t9 = mulPoly(t7, t8);
    t10 = subPoly(t3, t9);

    r = reduceAST(t10);

    m = degree(r, x);
  }

  Expr res = list({algebraicExpand(q), algebraicExpand(r)});

  return res;
}

Expr quotientGPE(Expr u, Expr v, Expr x) { return divideGPE(u, v, x)[0]; }

Expr remainderGPE(Expr u, Expr v, Expr x) { return divideGPE(u, v, x)[1]; }

Expr expandGPE(Expr u, Expr v, Expr x, Expr t) {
  if (u == 0)
    return 0;

  Expr d = divideGPE(u, v, x);

  Expr q = d[0];
  Expr r = d[1];

  Expr expoent = add({mul({t, expandGPE(q, v, x, t)}), r});

  return algebraicExpand(expoent);
}

Expr gcdGPE(Expr u, Expr v, Expr x) {
  if (u.kind() == Kind::Integer && u.value() == 0 &&
      v.kind() == Kind::Integer && v.value() == 0) {
    return integer(0);
  }

  Expr U = u;
  Expr V = v;

  while (V.kind() != Kind::Integer ||
         (V.kind() == Kind::Integer && V.value() != 0)) {
    Expr R = remainderGPE(U, V, x);

    U = V;
    V = R;
  }

  Expr e = mul({div(integer(1), leadCoeff(U, x)), U});
  Expr res = algebraicExpand(e);

  return res;
}

Expr extendedEuclideanAlgGPE(Expr u, Expr v, Expr x) {
  if (u.kind() == Kind::Integer && u.value() == 0 &&
      v.kind() == Kind::Integer && v.value() == 0) {
    return list({integer(0), integer(0), integer(0)});
  }

  Expr U = u;
  Expr V = v;

  Expr App = 1, Ap = 0, Bpp = 0, Bp = 1;

  while (V != 0) {
    Expr d = divideGPE(U, V, x);

    Expr q = d[0];
    Expr r = d[1];

    Expr A_ = sub({App, mul({q, Ap})});

    Expr B_ = sub({Bpp, mul({q, Bp})});

    Expr A = algebraicExpand(A_);
    Expr B = algebraicExpand(B_);

    App = Ap;

    Ap = A;

    Bpp = Bp;

    Bp = B;

    U = V;

    V = r;
  }

  Expr c = leadCoeff(U, x);

  Expr App_ = quotientGPE(App, c, x);

  App = App_;

  Expr Bpp_ = quotientGPE(Bpp, c, x);

  Bpp = Bpp_;

  Expr U_ = quotientGPE(U, c, x);

  U = U_;

  return list({U, App, Bpp});
}

Expr mulPolyRec(Expr p1, Expr p2) {
  if (p1.kind() == Kind::Addition) {
    Expr res = add({});

    for (unsigned int i = 0; i < p1.size(); i++) {
      res.insert(mulPoly(p1[i], p2));
    }

    return res;
  }

  if (p2.kind() == Kind::Addition) {
    Expr res = add({});

    for (unsigned int i = 0; i < p2.size(); i++) {
      res.insert(mulPoly(p2[i], p1));
    }

    return res;
  }

  return mul({p1, p2});
}

Expr mulPoly(Expr p1, Expr p2) {
  Expr t1 = mulPolyRec(p1, p2);
  Expr t2 = reduceAST(t1);

  return t2;
}

Expr subPoly(Expr p1, Expr p2) {
  Expr s = sub({p1, p2});
  Expr p = reduceAST(s);

  return p;
}

Expr addPoly(Expr p1, Expr p2) {
  Expr s = add({p1, p2});
  Expr p = reduceAST(s);

  return p;
}

Expr recPolyDiv(Expr u, Expr v, Expr L, Expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q",
         "Field needs to be Z or Q");

  if (L.size() == 0) {
    Expr d = algebraicExpand(u / v);

    if (K.identifier() == "Z") {
      if (d.kind() == Kind::Integer) {
        return list({d, 0});
      }

      return list({0, u});
    }

    return list({d, 0});
  }

  Expr x = L[0];
  Expr r = u;

  Expr m = degree(r, x);
  Expr n = degree(v, x);

  Expr q = 0;
  Expr lcv = leadCoeff(v, x);
  Expr R = rest(L);

	while (m.kind() != Kind::MinusInfinity && m.value() >= n.value()) {
    Expr lcr = leadCoeff(r, x);

    Expr d = recPolyDiv(lcr, lcv, R, K);

    if (d[1] != 0) {
      return list({algebraicExpand(q), r});
    }

    Expr j = power(x, m - n);

    q = q + d[0] * j;

    Expr t1 = mulPoly(v, d[0]);
    Expr t2 = mulPoly(t1, j);
    Expr t3 = subPoly(r, t2);

    r = reduceAST(t3);

    m = degree(r, x);
  }

  return list({algebraicExpand(q), r});
}

Expr recQuotient(Expr u, Expr v, Expr L, Expr K) {
  return recPolyDiv(u, v, L, K)[0];
}

Expr recRemainder(Expr u, Expr v, Expr L, Expr K) {
  return recPolyDiv(u, v, L, K)[1];
}

Expr pdiv(Expr f, Expr g, Expr x) {
  assert(g != 0, "Division by zero!");

  Expr lg, k, q, r, t, m, n, j;
  Expr t1, t2, t3, t4, t5, t6;

  m = degree(f, x);
  n = degree(g, x);

  if (m.value() < n.value()) {
    return list({0, f});
  }

  if (g == 1) {
    return list({f, 0});
  }

  q = 0;
  r = f;
  t = m;

  k = m - n + 1;

  lg = leadCoeff(g, x);

  while (true) {
    t1 = leadCoeff(r, x);
    j = sub({t, n});
    k = sub({k, integer(1)});
    t3 = power(x, j);

    t2 = mulPoly(q, lg); // mul({ q, lg });
    t4 = mulPoly(t1, t3);
    q = addPoly(t2, t4);

    t4 = mulPoly(r, lg); // mul({ r, lg });
    t5 = mulPoly(g, t1); // mul({ g, t1, t3 });
    t6 = mulPoly(t5, t3);
    r = subPoly(t4, t6);

    t = degree(r, x);

    if (t.kind() == Kind::MinusInfinity || t.value() < n.value()) {
      break;
    }
  }

  q = mul({q, power(lg, k)});
  r = mul({r, power(lg, k)});

  t1 = algebraicExpand(q);
  t2 = algebraicExpand(r);

  return list({t1, t2});
}

Expr pseudoDivision(Expr u, Expr v, Expr x) {
  Expr p = 0;
  Expr s = u;

  Expr m = degree(s, x);
  Expr n = degree(v, x);

  Expr delta = max(m.value() - n.value() + 1, 0);

  Expr lcv = leadCoeff(v, x);

  Int tal = 0;

  while (m.kind() != Kind::MinusInfinity && m.value() >= n.value()) {
    Expr lcs = leadCoeff(s, x);

    Expr j = power(x, sub({m, n}));

    Expr t1 = mulPoly(lcv, p);
    Expr t2 = mulPoly(lcs, j);

    Expr t3 = addPoly(t1, t2);

    p = reduceAST(t3);

    Expr t4 = mulPoly(lcv, s);
    Expr t5 = mulPoly(lcs, v);
    Expr t6 = mulPoly(t5, j);

    s = subPoly(t4, t6);

    tal = tal + 1;

    m = degree(s, x);
  }

  Expr k = power(lcv, delta.value() - tal);

  Expr A = mulPoly(k, p);
  Expr B = mulPoly(k, s);

  Expr Q = reduceAST(A);
  Expr R = reduceAST(B);

  return list({Q, R});
}

Expr pseudoQuotient(Expr u, Expr v, Expr x) {
  Expr r = pseudoDivision(u, v, x);
  Expr q = r[0];

  return q;
}

Expr pseudoRemainder(Expr u, Expr v, Expr x) {
  Expr r = pseudoDivision(u, v, x);
  Expr q = r[1];

  return q;
}

Expr getNormalizationFactor(Expr u, Expr L, Expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q",
         "field must be Z or Q");

  if (u == (0)) {
    return integer(0);
  }

  if (isConstant(u)) {
    if (u > 0) {
      if (K.identifier() == "Z") {
        return integer(1);
      }

      return power(u, integer(-1));
    } else {
      if (K.identifier() == "Z") {
        return -1;
      }

      return -1 * power(u, -1);
    }
  }

  if (L.size() == 0) {
    return undefined();
  }

  Expr lc = leadCoeff(u, L[0]);

  Expr rL = rest(L);

  Expr cf = getNormalizationFactor(lc, rL, K);

  return cf;
}

Expr normalizePoly(Expr u, Expr L, Expr K) {
  if (u.kind() == Kind::Integer && u.value() == 0) {
    return integer(0);
  }

  Expr u__ = mul({getNormalizationFactor(u, L, K), u});

  Expr u_ = algebraicExpand(u__);

  return u_;
}

Expr unitNormal(Expr v, Expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q",
         "field must be Z or Q");

  if (K.identifier() == "Z") {
    if (v < 0) {
      return -1;
    }
    return 1;
  }

  if (K.identifier() == "Q") {
    if (v < 0) {
      return -1 * power(v, -1);
    }

    return power(v, -1);
  }

  return 1;
}

// Finds the content of u with respect to x using
// the auxiliary variables R with coeff domain K,
// with is Z or Q
Expr polynomialContent(Expr u, Expr x, Expr R, Expr K) {
  if (u == (0)) {
    return integer(0);
  }

  Expr n = degree(u, x);

  Expr g = coeff(u, x, n);

  Expr k = sub({u, mul({g, power(x, n)})});

  Expr v = algebraicExpand(k);

  if (v == (0)) {
    Expr un = unitNormal(g, K);
    Expr t = mul({un, g});

    g = reduceAST(t);

  } else {
    while (v != 0) {
      Expr d = degree(v, x);
      Expr c = leadCoeff(v, x);

      Expr t = mvPolyGCD(g, c, R, K);

      g = t;

      k = sub({v, mul({c, power(x, d)})});
      v = algebraicExpand(k);
    }
  }

  return g;
}

// Finds the content of u with respect to x using
// the auxiliary variables R with coeff domain K,
// with is Z or Q
Expr polynomialContentSubResultant(Expr u, Expr x, Expr R, Expr K) {
  if (u == (0)) {
    return integer(0);
  }

  Expr n = degree(u, x);

  Expr g = coeff(u, x, n);

  Expr v = algebraicExpand(u - g * power(x, n));

  if (v == 0) {
    g = reduceAST(unitNormal(g, K) * g);
  } else {
    while (v != 0) {
      Expr d = degree(v, x);
      Expr c = leadCoeff(v, x);
      Expr t = mvSubResultantGCD(g, c, R, K);

      g = t;

      v = algebraicExpand(v - c * power(x, d));
    }
  }
  return g;
}

Expr subResultantGCDRec(Expr u, Expr v, Expr L, Expr K) {
  if (L.size() == 0) {
    if (K.identifier() == "Z") {
      return gcd(u.value(), v.value());
    }

    if (K.identifier() == "Q") {
      return integer(1);
    }
  }

  Expr x = first(L);

  Expr du = degree(u, x);
  Expr dv = degree(v, x);

  Expr U = undefined();
  Expr V = undefined();

  if (du.value() >= dv.value()) {
    U = u;
    V = v;
  } else {
    U = v;
    V = u;
  }

  Expr R = rest(L);

  Expr contU = polynomialContentSubResultant(U, x, R, K);
  Expr contV = polynomialContentSubResultant(V, x, R, K);

  Expr d = subResultantGCDRec(contU, contV, R, K);

  Expr tmp1 = recQuotient(U, contU, L, K);
  Expr tmp2 = recQuotient(V, contV, L, K);

  U = tmp1;

  V = tmp2;

  Expr tmp3 = leadCoeff(U, x);
  Expr tmp4 = leadCoeff(V, x);

  Expr g = subResultantGCDRec(tmp3, tmp4, R, K);

  int i = 1;

  Expr delta = undefined();
  Expr y = undefined();
  Expr b = undefined();
  Expr dp = undefined();

  while (V != 0) {
    Expr r = pseudoRemainder(U, V, x);

    if (r != 0) {
      if (i == 1) {

        Expr tmp3 =
            add({degree(U, x), mul({integer(-1), degree(V, x)}), integer(1)});

        delta = algebraicExpand(tmp3);

        y = integer(-1);

        Expr tmp4 = power(integer(-1), delta);

        b = algebraicExpand(tmp4);

      } else {
        dp = delta;

        Expr tmp3 =
            add({degree(U, x), mul({integer(-1), degree(V, x)}), integer(1)});

        delta = algebraicExpand(tmp3);

        Expr f = leadCoeff(U, x);

        Expr tmp4 = power(mul({integer(-1), f}), sub({dp, integer(1)}));
        Expr tmp5 = power(y, sub({dp, integer(2)}));

        Expr tmp6 = algebraicExpand(tmp4);
        Expr tmp7 = algebraicExpand(tmp5);

        y = recQuotient(tmp6, tmp7, R, K);

        Expr tmp8 = mul({integer(-1), f, power(y, sub({delta, integer(1)}))});

        b = algebraicExpand(tmp8);
      }

      U = V;

      V = recQuotient(r, b, L, K);

      i = i + 1;
    } else {
      U = V;

      V = r;
    }
  }

  Expr tmp5 = leadCoeff(U, x);

  Expr s = recQuotient(tmp5, g, R, K);

  Expr W = recQuotient(U, s, L, K);

  Expr contW = polynomialContentSubResultant(W, x, R, K);
  Expr ppW = recQuotient(W, contW, L, K);

  Expr tmp6 = mul({d, ppW});
  Expr res = algebraicExpand(tmp6);

  return res;
}

Expr mvSubResultantGCD(Expr u, Expr v, Expr L, Expr K) {
  if (u == (0)) {
    return normalizePoly(v, L, K);
  }

  if (v == (0)) {
    return normalizePoly(u, L, K);
  }

  Expr gcd = subResultantGCDRec(u, v, L, K);

  Expr r = normalizePoly(gcd, L, K);

  return r;
}

Expr mvPolyGCDRec(Expr u, Expr v, Expr L, Expr K) {
  if (L.size() == 0) {
    if (K.identifier() == "Z") {
      return gcd(u.value(), v.value());
    }

    if (K.identifier() == "Q") {
      return integer(1);
    }
  }
	printf("\nu(x) = %s\n", u.toString().c_str());
	printf("\nv(x) = %s\n", v.toString().c_str());
  Expr x = first(L);
  Expr R = rest(L);
  Expr cont_u = polynomialContent(u, x, R, K);
  Expr cont_v = polynomialContent(v, x, R, K);

  Expr d = mvPolyGCDRec(cont_u, cont_v, R, K);

  Expr pp_u = recQuotient(u, cont_u, L, K);
  Expr pp_v = recQuotient(v, cont_v, L, K);

  while (pp_v != 0) {
    Expr r = pseudoRemainder(pp_u, pp_v, x);

    Expr pp_r = undefined();

    if (r == (0)) {
      pp_r = integer(0);
    } else {
      Expr cont_r = polynomialContent(r, x, R, K);
      pp_r = recQuotient(r, cont_r, L, K);
    }

    pp_u = pp_v;
    pp_v = pp_r;
  }

  Expr k = mul({d, pp_u});
  Expr result = algebraicExpand(k);

  return result;
}

Expr mvPolyGCD(Expr u, Expr v, Expr L, Expr K) {
	printf("\n======> GCD u(x) = %s\n", u.toString().c_str());
	printf("\n======> GCD v(x) = %s\n", v.toString().c_str());

	if (u == (0)) {
    return normalizePoly(v, L, K);
  }

  if (v == (0)) {
    return normalizePoly(u, L, K);
  }

  Expr gcd = mvPolyGCDRec(u, v, L, K);
  Expr r = normalizePoly(gcd, L, K);

  return r;
}

Expr leadMonomial(Expr u, Expr L) {
  if (L.size() == 0) {
    return u;
  }

  Expr x = first(L);
  Expr m = degree(u, x);

  Expr c = coeff(u, x, m);

  Expr restL = rest(L);

  Expr r_ = mul({power(x, m), leadMonomial(c, restL)});

  Expr r = algebraicExpand(r_);

  return r;
}

bool wasSimplified(Expr u) {
  if (u.kind() == Kind::Symbol)
    return true;

  if (u.kind() == Kind::Integer)
    return true;

  if (u.kind() == Kind::Division)
    return false;

  if (u.kind() == Kind::Multiplication) {

    for (unsigned int i = 0; i < u.size(); i++) {
      if (u[i].kind() == Kind::Fraction) {
        return false;
      }

      if (u[i].kind() == Kind::Power && u[i][1].kind() == Kind::Integer &&
          u[i][1].value() < 0) {
        return false;
      }
    }

    return true;
  }

  return false;
}

/**
 * Return summation(u[i]/v) if v divides u[i]
 */
Expr G(Expr u, Expr v) {
  if (u.kind() == Kind::Addition || u.kind() == Kind::Subtraction) {
    Expr k = Expr(u.kind());

    for (unsigned int i = 0; i < u.size(); i++) {
      Expr z_ = div(u[i], v);
      Expr z = algebraicExpand(z_);

      if (wasSimplified(z)) {
        k.insert(z);
      } else {
      }
    }

    if (k.size() == 0) {

      return integer(0);
    }

    return k;
  }

  Expr z_ = div(u, v);

  Expr z = algebraicExpand(z_);

  if (wasSimplified(z)) {
    return z;
  }

  return integer(0);
}

Expr monomialPolyDiv(Expr u, Expr v, Expr L) {
  Expr q = integer(0);
  Expr r = u;
  Expr vt = leadMonomial(v, L);

  Expr f = G(r, vt);

  while (f.kind() != Kind::Integer || f.value() != 0) {

    q = add({q, f});

    Expr r_ = sub({r, mul({f, v})});

    r = algebraicExpand(r_);

    f = G(r, vt);
  }

  Expr l = list({reduceAST(q), reduceAST(r)});

  return l;
}

// TODO

// monomialBasedPolyExpansion(a^2*b + 2*a*b^2 + b^3 + 2*a + 2*b + 3, a+b, [a,
// b], t) -> b*t^2 + 2*t + 3
Expr monomialBasedPolyExpansion(Expr u, Expr v, Expr L, Expr t) {
  if (u.kind() == Kind::Integer && u.value() == 0) {
    return integer(0);
  }

  Expr d = monomialPolyDiv(u, v, L);
  Expr q = d[0];
  Expr r = d[1];

  Expr k = add({mul({
                    t,
                    monomialBasedPolyExpansion(q, v, L, t),
                }),
                r});

  Expr x = algebraicExpand(k);

  return x;
}

// monomialPolyRem can be used for simplification, for example
// monomialPolyRem(a*i^3 + b*i^2 + c*i + d, i^2 + 1, [i]) -> -a*i - b + c*i + d
// simplification when i^2 + 1 = 0
// also
// monomialPolyRem(sin^4(x)+sin^3(x)+2*sin^2(x)cos^2(x)+cos^4(x),
// sin^2(x)+cos^2(x)-1, [cos(x), sin(x)]) -> 1 + sin^3(x)
Expr monomialPolyRem(Expr u, Expr v, Expr L) {
  Expr d = monomialPolyDiv(u, v, L);
  Expr r = d[1];

  return r;
}

Expr monomialPolyQuo(Expr u, Expr v, Expr L) {
  Expr d = monomialPolyDiv(u, v, L);
  Expr r = d[0];

  return r;
}

Expr algebraicExpandRec(Expr u);

Expr expandProduct(Expr r, Expr s) {
  // printf("IN = %s * %s\n", r.toString().c_str(), s.toString().c_str());
  if (r == 0 || s == 0)
    return 0;

  if (r.kind() == Kind::Addition && r.size() == 0)
    return 0;
  if (s.kind() == Kind::Addition && s.size() == 0)
    return 0;
  // if (r.kind() == Kind::Multiplication) r = algebraicExpand(r);

  if (r.kind() == Kind::Addition) {
    Expr f = r[0];

    Expr k = r;

    k.remove(0);

    Expr a = expandProduct(f, s);
    Expr b = expandProduct(k, s);

    bool c0 = a == 0;
    bool c1 = b == 0;

    Expr z = 0;

    if (!c0 && !c1)
      z = a + b;
    else if (!c0 && c1)
      z = a;
    else if (c0 && !c1)
      z = b;
    else
      z = 0;

    // printf("1. OUT = %s\n", z.toString().c_str());
    return z;
  }

  if (s.kind() == Kind::Addition) {
    return expandProduct(s, r);
  }

  // printf("2.OUT = %s\n", reduceAST(r * s).toString().c_str());
  return reduceAST(r * s);
}

Expr expandPower(Expr u, Expr n) {
  if (u == 1)
    return 1;
  if (n == 0)
    return u == 0 ? undefined() : 1;

  if (u.kind() == Kind::Addition) {
    Expr f = u[0];

    Expr o = reduceAST(u - f);
    Expr s = 0;

    Int N = n.value();

    for (Int k = 0; k <= N; k++) {
      Int d = fact(N) / (fact(k) * fact(N - k));

      Expr z = d * power(f, N - k);
      Expr t = expandPower(o, k);

      s = s + expandProduct(z, t);
    }

    return s;
  }

  return reduceAST(power(u, n));
}

Expr expandProductRoot(Expr r, Expr s) {
  if (r.kind() == Kind::Addition) {
    Expr f = r[0];
    Expr k = r;

    k.remove(0);

    return f * s + k * s;
  }

  if (s.kind() == Kind::Addition) {
    return expandProductRoot(s, r);
  }

  return r * s;
}

Expr expandPowerRoot(Expr u, Expr n) {
  if (u.kind() == Kind::Addition) {
    Expr f = u[0];

    Expr r = reduceAST(u - f);

    Expr s = 0;

    Int N = n.value();

    for (Int k = 0; k <= n.value(); k++) {
      Expr c = fact(N) / (fact(k) * fact(N - k));
      Expr z = reduceAST(c * power(f, N - k));

      Expr t = expandPowerRoot(r, k);

      s = s + expandProductRoot(z, t);
    }

    return s;
  }

  Expr v = power(u, n);

  return reduceAST(v);
}

Expr algebraicExpandRoot(Expr u) {
  if (u.isTerminal())
    return reduceAST(u);

  Expr u_ = reduceAST(u);

  if (u_.kind() == Kind::Addition) {
    Expr v = u_[0];
    Expr k = reduceAST(u_ - v);
    u_ = algebraicExpandRoot(v) + algebraicExpandRoot(k);
  }

  if (u_.kind() == Kind::Multiplication) {
    Expr v = u_[0];
    Expr t = reduceAST(u_ / v);
    u_ = expandProductRoot(t, v);
  }

  if (u_.kind() == Kind::Power) {

    Expr b = u_[0];
    Expr e = u_[1];

    if (e.kind() == Kind::Integer && e.value() >= 2) {
      Expr t = expandPowerRoot(b, e);
      u_ = reduceAST(t);
    }

    if (e.kind() == Kind::Integer && e.value() <= -2) {
      Expr p = reduceAST(power(u_, -1));
      Expr t = expandPowerRoot(p[0], p[1]);
      u_ = power(t, -1);
    }
  }

  Expr k = reduceAST(u_);

  return k;
}

Expr algebraicExpandRec(Expr u) {
  if (u.isTerminal())
    return u;

  if (u.kind() == Kind::Subtraction || u.kind() == Kind::Division ||
      u.kind() == Kind::Factorial) {
    u = reduceAST(u);
  }

  if (u.kind() == Kind::Power) {
    Expr t = algebraicExpand(u[0]);
    Expr z = algebraicExpand(u[1]);

    if (z.kind() == Kind::Integer) {
      u = expandPower(t, z);
    } else {
      u = power(t, z);
    }
  }

  if (u.kind() == Kind::Multiplication) {
    Expr t = algebraicExpand(u[0]);

    for (Int i = 1; i < u.size(); i++) {
      t = expandProduct(t, algebraicExpand(u[i]));
    }

    u = t;
  }

  if (u.kind() == Kind::Addition) {
    Expr t = algebraicExpand(u[0]);

    for (Int i = 1; i < u.size(); i++) {
      t = t + algebraicExpand(u[i]);
    }

    u = t;
  }

  // printf("out -> %s\n", u.toString().c_str());
  return u;
}

Expr algebraicExpand(Expr u) {
  Expr t = algebraicExpandRec(u);
  // printf("out -> %s\n", t.toString().c_str());
  Expr k = reduceAST(t);
  // printf("out -> %s\n", k.toString().c_str());
  return k;
}

Expr cont(Expr u, Expr x) {
  Expr n, c, c1, c2, tmp;

  u = algebraicExpand(u);

  if (u == 0) {
    return 0;
  }

  if (u.size() >= 2) {
    n = degree(u, x);

    c1 = coeff(u, x, n);

    u = algebraicExpand(u - c1 * power(x, n));

    n = degree(u, x);

    c2 = coeff(u, x, n);

    u = algebraicExpand(u - c2 * power(x, n));

    c = gcd(c1.value(), c2.value());

    while (u != 0) {
      n = degree(u, x);

      c1 = coeff(u, x, n);

      tmp = u - c1 * power(x, n);

      u = algebraicExpand(tmp);

      c2 = gcd(c.value(), c1.value());

      c = c2;
    }

    return c;
  }

  if (u.size() == 1) {
    return coeff(u, x, degree(u, x));
  }

  return 0;
}

Expr cont(Expr u, Expr L, Expr K) {
  return polynomialContentSubResultant(u, L[0], rest(L), K);
}

Expr pp(Expr u, Expr L, Expr K) {
  Expr c = cont(u, L, K);

  Expr p = recQuotient(u, c, L, K);

  return p;
}

Expr pp(Expr u, Expr c, Expr L, Expr K) {
  Expr R = rest(L);

  Expr p = recQuotient(u, c, L, K);

  return p;
}

Expr groundLeadCoeff(Expr u, Expr L) {
	if(L.size() == 0) return u;
	return groundLeadCoeff(leadCoeff(u, L[0]), rest(L));
}

Expr monic(Expr u, Expr L, Expr K) {
	Expr lc = groundLeadCoeff(u, L);
	return recQuotient(u, lc, L, K);
}

Expr gcdPoly(Expr u, Expr v, Expr L, Expr K) {
	if(L.size() == 0) {
		assert(u.kind() == Kind::Integer, "polynomial with zero variables should be a integer");
		assert(v.kind() == Kind::Integer, "polynomial with zero variables should be a integer");

		return gcd(u.value(), v.value());
	}

	Expr u_cnt = cont(u, L, K);
	Expr u_ppr = pp(u, u_cnt, L, K);

	Expr v_cnt = cont(v, L, K);
	Expr v_ppr = pp(v, v_cnt, L, K);

	Expr R = rest(L);

	Expr h = polyRemSeq(u_ppr, v_ppr, L, K)[0];

	printf("%s\n", h.toString().c_str());

	Expr c = gcdPoly(u_cnt, v_cnt, R, K);

	h = pp(h, L, K);
	h = mulPoly(h, c);

	return monic(h, L, K);
}


bool isZeroPolyExpr(Expr &u) {
  if (u.isTerminal()) {
    return u.kind() == Kind::Integer && u.value() == 0;
  }

  if (u.kind() == Kind::Addition && u.size() > 1) {
    return 0;
  }

  return isZeroPolyExpr(u[0]);
}

bool isConstantPolyExpr(Expr &u, Expr v) {
	if(v == 0) return isZeroPolyExpr(u);

	if (u.isTerminal()) {
    return (u.kind() == Kind::Integer || u.kind() == Kind::Fraction) && u == v;
  }

	if(u.kind() == Kind::Multiplication) {
		if(u.size() > 2) return false;
		if(u.size() == 2) {
			assert(u[1].kind() == Kind::Power, "not a poly expression");
			assert(u[1][1].kind() == Kind::Integer, "not a poly expression");

			if(u[1][1].value() != 0) {
				return false;
			}
		}
	}

	if (u.kind() == Kind::Addition && u.size() > 1) {
    return 0;
  }

  return isConstantPolyExpr(u[0], v);
}
bool isConstantPolyExpr(Expr &u) {
  if (u.isTerminal()) {
    return u.kind() == Kind::Integer || u.kind() == Kind::Fraction;
  }

  if (u.kind() == Kind::Addition && u.size() > 1) {
    return 0;
  }

  return isConstantPolyExpr(u[0]);
}

Expr mulPolyExpr(Expr &&p1, Expr &&p2) {
  if (p1.isTerminal() && p2.isTerminal()) {
    return reduceAST(p1 * p2);
  }

  if (p2.isTerminal()) {
    return mulPolyExpr(p2, p1);
  }

  if (p1.isTerminal()) {
    if (p2.kind() == Kind::Multiplication) {
      return mulPolyExpr(p1, p2[0]) * p2[1];
    }
    if (p2.kind() == Kind::Addition) {
      Expr g = Expr(Kind::Addition);

      for (size_t i = 0; i < p2.size(); i++) {
        g.insert(mulPolyExpr(p1, p2[i][0]) * p2[i][1]);
      }

      return g;
    }
  }

  // TODO: should this function perform division in the case that the other
  // operand is a division? In the case of powPolyExpr, a division is returned
  // if the expoent is < 0, in that case it will be usefful to perform the
  // division here

  // if(p2.kind() == Kind::Division) {
  // 	//assert(p2[1].kind() == Kind::Addition, "denominator of division needs
  // to be a polynomial expression"); 	Expr a = 1; 	Expr b = 1;

  // 	if(p1.kind() == Kind::Division) {
  // 		a = p1[0];
  // 		b = p1[1];
  // 	} else {
  // 		a = p1;
  // 		b = 1;
  // 	}

  // 	Expr c = p2[0];
  // 	Expr d = p2[1];

  // 	Expr e = mulPolyExpr(a, c);
  // 	Expr f = mulPolyExpr(b, d);

  // }

  Expr x = p1[0][1][0];

  std::map<Int, Expr> coeffs;

  for (size_t i = 0; i < p1.size(); ++i) {
    assert(p1[i][1][0] == x, "p1 was collected incorrectly!!!");

    Expr u = p1[i];

    for (size_t j = 0; j < p2.size(); j++) {
      assert(p2[j][1][0] == x, "p2 was collected incorrectly!!!");

      Expr v = p2[j];

      Int e = u[1][1].value() + v[1][1].value();

      Expr c = mulPolyExpr(u[0], v[0]);

      if (coeffs.count(e) == 0) {
        coeffs[e] = c;
      } else {
        coeffs[e] = addPolyExpr(coeffs[e], c);
      }
    }
  }

  Expr g = Expr(Kind::Addition);

  for (std::map<Int, Expr>::iterator it = coeffs.begin(); it != coeffs.end();
       it++) {
    g.insert(it->second * power(x, it->first));
  }

  return g;
}

Expr mulPolyExpr(Expr &p1, Expr &p2) {
  return mulPolyExpr(std::forward<Expr>(p1), std::forward<Expr>(p2));
}

Expr addColPolyRec(Expr &u, Expr &v, unsigned int i = 0, unsigned int j = 0) {
  if (i == u.size() && j == v.size())
    return 0;

  if (i == u.size()) {
    std::vector<Expr> ops;

    for (size_t t = j; t < v.size(); t++) {
      ops.push_back(v[t]);
    }

    return Expr(Kind::Addition, ops);
  }

  if (j == v.size()) {
    std::vector<Expr> ops;

    for (size_t t = i; t < u.size(); t++) {
      ops.push_back(u[t]);
    }

    return Expr(Kind::Addition, ops);
  }

  if (u[i].kind() == Kind::Multiplication && u[i].size() == 2 &&
      v[j].kind() == Kind::Multiplication && v[j].size() == 2) {

    Expr &ucoeff = u[i][0];
    Expr &upower = u[i][1];

    Expr &vcoeff = v[j][0];
    Expr &vpower = v[j][1];

    if (upower[0] == vpower[0]) {

      bool uzero = isZeroPolyExpr(ucoeff);
      bool vzero = isZeroPolyExpr(vcoeff);

      if (uzero && vzero) {
        Expr a = addColPolyRec(u, v, i + 1, j + 1);

        if (isZeroPolyExpr(a)) {
          return Expr(Kind::Addition, {0 * power(upower[0], 1)});
        }

        return a;
      }

      if (uzero) {
        return addColPolyRec(u, v, i + 1, j);
      }

      if (vzero) {
        return addColPolyRec(u, v, i, j + 1);
      }

      if (upower[1] == vpower[1]) {
        Expr a = addPolyExpr(ucoeff, vcoeff);

        if (!isZeroPolyExpr(a)) {

          if (a.kind() == Kind::Multiplication) {
            a = Expr(Kind::Addition, {a});
          }

          a = a * upower;
        }

        Expr b = addColPolyRec(u, v, i + 1, j + 1);

        if (isZeroPolyExpr(b)) {
          if (isZeroPolyExpr(a)) {
            return Expr(Kind::Addition, {0 * power(upower[0], 1)});
          }

          return Expr(Kind::Addition, {a});
        }

        if (isZeroPolyExpr(a))
          return b;

        if (b.kind() == Kind::Addition) {
          b.insert(a, 0);
          return b;
        }

        return a + b;
      }

      if (upower[1].value() < vpower[1].value()) {
        Expr a = addColPolyRec(u, v, i + 1, j);

        if (isZeroPolyExpr(a))
          return Expr(Kind::Addition, {u[i]});

        if (a.kind() == Kind::Addition) {
          a.insert(u[i], 0);
          return a;
        }

        return u[i] + a;
      }

      Expr a = addColPolyRec(u, v, i, j + 1);
      Expr c = Expr(Kind::Addition, {v[j]});

      if (isZeroPolyExpr(a)) {
        return c;
      }

      for (Int i = 0; i < a.size(); i++) {
        c.insert(a[i]);
      }

      return c;
    }
  }

  Expr a = addColPolyRec(u, v, i + 1, j + 1);

  Expr b;

  if ((u[i].kind() == Kind::Integer || u[i].kind() == Kind::Fraction) &&
      (v[j].kind() == Kind::Integer || v[j].kind() == Kind::Fraction)) {
    b = reduceAST(u[i] + v[j]);
  } else {
    b = u[i] + v[j];
  }

  if (a == 0)
    return b;

  if (a.kind() == Kind::Addition) {
    a.insert(b, 0);
    return a;
  }

  return b + a;
}

Expr addPolyExpr(Expr &u, Expr &v) {
  if ((u.kind() == Kind::Integer || u.kind() == Kind::Fraction) &&
      (v.kind() == Kind::Integer || v.kind() == Kind::Fraction)) {
    return reduceAST(u + v);
  }

  return addColPolyRec(u, v, 0, 0);
}

Expr addPolyExpr(Expr &&u, Expr &&v) {
  if ((u.kind() == Kind::Integer || u.kind() == Kind::Fraction) &&
      (v.kind() == Kind::Integer || v.kind() == Kind::Fraction)) {
    return reduceAST(u + v);
  }

  return addColPolyRec(u, v, 0, 0);
}

Expr subColPolyRec(Expr &u, Expr &v, unsigned int i = 0, unsigned int j = 0) {
  Expr k = -1;

  if (i == u.size() && j == v.size())
    return 0;

  if (i == u.size()) {
    Expr g = Expr(Kind::Addition);

    for (size_t t = j; t < v.size(); t++) {
      g.insert(v[t]);
    }

    return mulPolyExpr(k, g);
  }

  if (j == v.size()) {
    std::vector<Expr> ops;

    for (size_t t = i; t < u.size(); t++) {
      ops.push_back(u[t]);
    }

    return Expr(Kind::Addition, ops);
  }

  if (u[i].kind() == Kind::Multiplication && u[i].size() == 2 &&
      v[j].kind() == Kind::Multiplication && v[j].size() == 2) {

    Expr &ucoeff = u[i][0];
    Expr &upower = u[i][1];

    Expr &vcoeff = v[j][0];
    Expr &vpower = v[j][1];

    if (upower[0] == vpower[0]) {
      bool uzero = isZeroPolyExpr(ucoeff);
      bool vzero = isZeroPolyExpr(vcoeff);

      if (uzero && vzero) {
        Expr a = subColPolyRec(u, v, i + 1, j + 1);

        if (isZeroPolyExpr(a)) {
          return Expr(Kind::Addition, {0 * power(upower[0], 1)});
        }

        return a;
      }

      if (uzero) {
        return subColPolyRec(u, v, i + 1, j);
      }

      if (vzero) {
        return subColPolyRec(u, v, i, j + 1);
      }

      if (upower[1] == vpower[1]) {
        Expr a = subPolyExpr(ucoeff, vcoeff);

        if (!isZeroPolyExpr(a)) {

          if (a.kind() == Kind::Multiplication) {
            a = Expr(Kind::Addition, {a});
          }

          a = a * upower;
        }

        Expr b = subColPolyRec(u, v, i + 1, j + 1);

        if (isZeroPolyExpr(b)) {

          if (isZeroPolyExpr(a)) {
            return Expr(Kind::Addition, {0 * power(upower[0], 1)});
          }

          return Expr(Kind::Addition, {a});
        }

        if (isZeroPolyExpr(a))
          return b;

        if (b.kind() == Kind::Addition) {
          b.insert(a, 0);
          return b;
        }

        return a + b;
      }

      if (upower[1].value() < vpower[1].value()) {
        Expr a = subColPolyRec(u, v, i + 1, j);

        if (isZeroPolyExpr(a)) {
          return Expr(Kind::Addition, {u[i]});
        }

        if (a.kind() == Kind::Addition) {
          a.insert(u[i], 0);
          return a;
        }

        return u[i] + a;
      }

      Expr a = subColPolyRec(u, v, i, j + 1);

      Expr c = mulPolyExpr(-1, Expr(Kind::Addition, {v[j]}));

      if (isZeroPolyExpr(a)) {
        return c;
      }
      for (Int i = 0; i < a.size(); i++) {
        c.insert(a[i]);
      }

      return c;
    }
  }

  Expr a = subColPolyRec(u, v, i + 1, j + 1);

  Expr b;

  if ((u[i].kind() == Kind::Integer || u[i].kind() == Kind::Fraction) &&
      (v[j].kind() == Kind::Integer || v[j].kind() == Kind::Fraction)) {
    b = reduceAST(u[i] + -v[j]);
  } else {
    b = u[i] + mulPolyExpr(k, v[j]);
  }

  if (isZeroPolyExpr(a)) {
    return b;
  }

  if (a.kind() == Kind::Addition) {
    a.insert(b, 0);
    return a;
  }

  return b + mulPolyExpr(k, a);
}

Expr subPolyExpr(Expr &u, Expr &v) {
  if ((u.kind() == Kind::Integer || u.kind() == Kind::Fraction) &&
      (v.kind() == Kind::Integer || v.kind() == Kind::Fraction)) {
    return reduceAST(u - v);
  }

  return subColPolyRec(u, v, 0, 0);
}

Expr subPolyExpr(Expr &&u, Expr &&v) {
  if ((u.kind() == Kind::Integer || u.kind() == Kind::Fraction) &&
      (v.kind() == Kind::Integer || v.kind() == Kind::Fraction)) {
    return reduceAST(u - v);
  }

  return subColPolyRec(u, v, 0, 0);
}

Expr powPolyExpr(Expr &u, Int v, Expr &L) {
  if (v < 0) {
    return 1 / powPolyExpr(u, -v, L);
  }

  Expr g = polyExpr(1, L);

  Expr x = u;

  while (v) {
    if (v % 2 == 1)
      g = mulPolyExpr(g, x);

    v = v / 2;

    x = mulPolyExpr(x, x);
  }

  return g;
}

Expr powPolyExpr(Expr &&x, Int v, Expr &L) {
  if (v < 0) {
    return 1 / powPolyExpr(x, -v, L);
  }

  Expr g = polyExpr(1, L);

  while (v) {
    if (v % 2 == 1)
      g = mulPolyExpr(g, x);

    v = v / 2;

    x = mulPolyExpr(x, x);
  }

  return g;
}

Expr leadCoeffPolyExpr(Expr &u) {
  if (u.isTerminal())
    return u;

  assert(u.kind() == Kind::Addition, "u wasn't collected correctly!!!");

  return u[u.size() - 1][0];
}

Expr degreePolyExpr(Expr &u) {
  if (isZeroPolyExpr(u))
    return -inf();

  if (u.kind() == Kind::Integer)
    return 0;
  if (u.kind() == Kind::Fraction)
    return 0;

  assert(u.kind() == Kind::Addition, "u wasn't collected correctly!!!");

  return u[u.size() - 1][1][1];
}

Expr raiseToPolyExpr(Expr &x, Int e, Expr &L) {
  return polyExpr(power(x, e), L);
}

Expr raiseToPolyExpr(Expr &&x, Int e, Expr &&L) {
  return polyExpr(power(x, e), L);
}

Expr raisePolyExpr(Expr &u, Int exp, Expr &x) {
  return add({u * power(x, exp)});
}

Expr raisePolyExpr(Expr &&u, Int exp, Expr &x) {
  return add({u * power(x, exp)});
}

Expr divPolyExpr(Expr &&u, Expr &&v, Expr &L, Expr &K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q",
         "Field needs to be Z or Q");

  if (L.size() == 0) {
    Expr d = reduceAST(u / v);

    if (K.identifier() == "Z") {
      if (d.kind() == Kind::Integer) {
        return list({d, 0});
      }

      return list({0, u});
    }

    return list({d, 0});
  }

  Expr r = u;

  Expr m = degreePolyExpr(r);
  Expr n = degreePolyExpr(v);

  Expr q = polyExpr(0, L);

  Expr lcv = leadCoeffPolyExpr(v);

  Expr R = rest(L);

  while (m != -inf() && m.value() >= n.value()) {
    Expr lcr = leadCoeffPolyExpr(r);

    Expr d = divPolyExpr(lcr, lcv, R, K);

    if (!isZeroPolyExpr(d[1])) {
      return list({q, r});
    }

    Expr k = raisePolyExpr(d[0], m.value() - n.value(), L[0]);

    q = addColPolyRec(q, k);

    Expr g = raisePolyExpr(d[0], 0, L[0]);

    Expr j = raiseToPolyExpr(L[0], m.value() - n.value(), L);

    Expr t1 = mulPolyExpr(v, g);
    Expr t2 = mulPolyExpr(t1, j);

    r = subPolyExpr(r, t2);

    m = degreePolyExpr(r);
  }

  return list({q, r});
}

Expr divPolyExpr(Expr &u, Expr &v, Expr &L, Expr &K) {
  return divPolyExpr(std::forward<Expr>(u), std::forward<Expr>(v), L, K);
}

Expr quoPolyExpr(Expr &u, Expr &v, Expr &L, Expr &K) {
  return divPolyExpr(std::forward<Expr>(u), std::forward<Expr>(v), L, K)[0];
}

Expr remPolyExpr(Expr &u, Expr &v, Expr &L, Expr &K) {
  return divPolyExpr(std::forward<Expr>(u), std::forward<Expr>(v), L, K)[1];
}

Expr pseudoDivPolyExpr(Expr &&u, Expr &&v, Expr &L) {
  Expr p = polyExpr(0, L);

  Expr s = u;

  Expr m = degreePolyExpr(s);
  Expr n = degreePolyExpr(v);

  Int t = 0;

  Int d = max(m.value() - n.value() + 1, 0);

  Expr lv = raisePolyExpr(leadCoeffPolyExpr(v), 0, L[0]);

  while (m.kind() != Kind::MinusInfinity && m.value() >= n.value()) {
    Expr ls = raisePolyExpr(leadCoeffPolyExpr(s), 0, L[0]);
    Expr jx = raiseToPolyExpr(L[0], m.value() - n.value(), L);

    Expr t1 = mulPolyExpr(lv, p);
    Expr t2 = mulPolyExpr(ls, jx);

    p = addPolyExpr(t1, t2);

    Expr t4 = mulPolyExpr(lv, s);
    Expr t5 = mulPolyExpr(ls, v);
    Expr t6 = mulPolyExpr(t5, jx);

    s = subPolyExpr(t4, t6);

    t = t + 1;

    m = degreePolyExpr(s);
  }

  Expr k = powPolyExpr(lv, d - t, L);

  Expr Q = mulPolyExpr(k, p);
  Expr R = mulPolyExpr(k, s);

  return list({Q, R});
}

Expr pseudoDivPolyExpr(Expr &u, Expr &v, Expr &L) {
  return pseudoDivPolyExpr(std::forward<Expr>(u), std::forward<Expr>(v), L);
}

Expr pseudoQuoPolyExpr(Expr &u, Expr &v, Expr &L) {
  return pseudoDivPolyExpr(std::forward<Expr>(u), std::forward<Expr>(v), L)[0];
}

Expr pseudoRemPolyExpr(Expr &u, Expr &v, Expr &L) {
  return pseudoDivPolyExpr(std::forward<Expr>(u), std::forward<Expr>(v), L)[1];
}

Expr colPolyGCDRec(Expr &u, Expr &v, Expr &L, Expr &K) {
  if (L.size() == 0) {
    if (K.identifier() == "Z") {
      return polyExpr(gcd(u.value(), v.value()), L);
    }

    if (K.identifier() == "Q") {
      return polyExpr(1, L);
    }
  }

  Expr R = rest(L);

  Expr ct_u = contPolyExpr(u, L, K);
  Expr ct_v = contPolyExpr(v, L, K);

  Expr d = colPolyGCDRec(ct_u, ct_v, R, K);

  ct_v = raisePolyExpr(ct_v, 0, L[0]);
  ct_u = raisePolyExpr(ct_u, 0, L[0]);

  Expr pp_u = quoPolyExpr(u, ct_u, L, K);
  Expr pp_v = quoPolyExpr(v, ct_v, L, K);

  Expr pp_r = 0;
  Expr ct_r = 0;

  Expr r = 0;

  while (!isZeroPolyExpr(pp_v)) {
    r = pseudoRemPolyExpr(pp_u, pp_v, L);

    if (isZeroPolyExpr(r)) {
      pp_r = polyExpr(0, L);
    } else {
      ct_r = contPolyExpr(r, L, K);
      ct_r = raisePolyExpr(ct_r, 0, L[0]);

      pp_r = quoPolyExpr(r, ct_r, L, K);
    }

    pp_u = pp_v;
    pp_v = pp_r;
  }

  d = raisePolyExpr(d, 0, L[0]);

  return mulPolyExpr(d, pp_u);
}

// bool isConstantPolyExpr(Expr& u) {
// 	if(u.kind() == Kind::Integer) {
// 		return true;
// 	}

// 	if(u.kind() == Kind::Addition && u.size() > 1) {
// 		return false;
// 	}

// 	if(u.kind() == Kind::Multiplication) {
// 		assert(u.size() == 2, "u should have only two operands");
// 		return u[1][1] == 0 && isConstantPolyExpr(u[0]);
// 	}

// 	return isConstantPolyExpr(u[0]);
// }

Expr getColPolyNormFactor(Expr &u, Expr &L, Expr &K) {
  if (isZeroPolyExpr(u)) {
    return Expr(u);
  }

  if (u.kind() == Kind::Integer || u.kind() == Kind::Fraction) {
    if (u.value() > 0) {
      if (K.identifier() == "Z") {
        return 1;
      }

      return u;
    } else {
      if (K.identifier() == "Z") {
        return -1;
      }

      return u;
    }
  }

  if (L.size() == 0) {
    return 1;
  }

  Expr lc = leadCoeffPolyExpr(u);

  Expr rL = rest(L);

  return Expr(Kind::Addition,
              {getColPolyNormFactor(lc, rL, K) * power(L[0], 0)});
}

Expr normalizePolyExpr(Expr &u, Expr &L, Expr &K) {
  if (isZeroPolyExpr(u))
    return Expr(u);
  Expr k = getColPolyNormFactor(u, L, K);
  return quoPolyExpr(u, k, L, K);
}

Expr gcdPolyExpr(Expr &u, Expr &v, Expr &L, Expr &K) {
  if (isZeroPolyExpr(u)) {
    return normalizePolyExpr(v, L, K);
  }

  if (isZeroPolyExpr(v)) {
    return normalizePolyExpr(u, L, K);
  }

  Expr g = colPolyGCDRec(u, v, L, K);

  return normalizePolyExpr(g, L, K);
}

Expr unitNormalColPoly(Expr v, Expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q",
         "field must be Z or Q");

  if (v.kind() == Kind::Addition) {
    return 1;
  }

  short sign = 1;

  if (v.kind() == Kind::Integer) {
    sign = v.value() < 0 ? -1 : 1;
  }
  if (v.kind() == Kind::Fraction) {
    sign = (v[0].value() < 0) ^ (v[1].value() < 0) ? -1 : 1;
  }

  if (K.identifier() == "Z") {
    if (sign < 0) {
      return -1;
    }

    return 1;
  }

  if (K.identifier() == "Q") {
    if (sign < 0) {
      if (v.kind() == Kind::Fraction) {
        return fraction(-1 * v[1].value(), v[0].value());
      }

      return fraction(-1, v);
    }

    if (v.kind() == Kind::Fraction) {
      return fraction(v[1].value(), v[0].value());
    }

    return fraction(1, v);
  }

  return 1;
}

Expr contPolyExpr(Expr &&u, Expr &L, Expr &K) {
  if (isZeroPolyExpr(u)) {
    return Expr(u);
  }

  Expr R = rest(L);

  long i = u.size() - 1;

  Expr g = u[i][0];

  if (i == 0) {
    Expr t = unitNormalColPoly(g, K);
    return mulPolyExpr(t, g);
  } else {
    i = i - 1;

    while (i >= 0) {
      Expr ui = u[i][0];
      g = gcdPolyExpr(g, ui, R, K);
      i = i - 1;
    }
  }

  return g;
}

Expr contPolyExpr(Expr &u, Expr &L, Expr &K) {
  return contPolyExpr(std::forward<Expr>(u), L, K);
}

Expr ppPolyExpr(Expr &&u, Expr &L, Expr &K) {
  Expr c = raisePolyExpr(contPolyExpr(u, L, K), 0, L[0]);
  return quoPolyExpr(u, c, L, K);
}

Expr ppPolyExpr(Expr &u, Expr &L, Expr &K) {
  return ppPolyExpr(std::forward<Expr>(u), L, K);
}

Expr contAndPpPolyExpr(Expr &&u, Expr &L, Expr &K) {
  Expr c = contPolyExpr(u, L, K);
  Expr C = raisePolyExpr(c, 0, L[0]);
  return list({c, quoPolyExpr(u, C, L, K)});
}

Expr contAndPpPolyExpr(Expr &u, Expr &L, Expr &K) {
  return contAndPpPolyExpr(std::forward<Expr>(u), L, K);
}

Expr expandPolyExpr(Expr &&u) {
	if(u.isTerminal() || u.kind() == Kind::FunctionCall) {
		return Expr(u);
	}

	if(u.kind() == Kind::Power) {
		if(u[1] == 0) return 1;
		if(u[1] == 1) return u[0];
		return Expr(u);
	}

  if(u.kind() == Kind::Multiplication) {
		Expr g = Expr(Kind::Multiplication);

		for(Int i = 0; i < u.size(); i++) {
			Expr c = expandPolyExpr(u[i]);

			if(c == 0) return 0;

			if(c == 1) continue;

			if(c.kind() == Kind::Multiplication) {
				for(Int j = 0; j < c.size(); j++) {
					g.insert(c[j]);
				}
			} else {
				g.insert(c);
			}
		}

		Expr t = expandPolyExpr(u[0]);

		if(t == 0) return 0;

		Expr p = u[1];

		if(p[1] == 0) return t;

		if(t.kind() == Kind::Addition) {
			Expr g = Expr(Kind::Addition);

			for(Int i = 0; i < u.size(); i++) {
				if(p[1] == 1) {
					g.insert(t[i]*p[0]);
				} else {
					g.insert(t[i]*p);
				}
			}

			return g;
		}

		if(p[1] == 1) return t*p[0];

		return t*p;
	}

	assert(u.kind() == Kind::Addition, "not a poly expr");

	Expr g = Expr(Kind::Addition);

  for (Int i = 0; i < u.size(); i++) {
    Expr c = expandPolyExpr(u[i]);

		if (c == 0)
      continue;

		if (c.kind() == Kind::Addition) {
			for (Int j = 0; j < c.size(); j++) {
				g.insert(c[j]);
			}
		} else {
			g.insert(c);
		}

  }

  if(g.size() == 0) return 0;
	if(g.size() == 1) return g[0];

	return g;
}

Expr expandPolyExpr(Expr& u) {
	return expandPolyExpr(std::forward<Expr>(u));
}

Expr diffPolyExprRec(Expr &u, Expr& x, bool* was_diff = nullptr) {
	if(u.kind() == Kind::Integer || u.kind() == Kind::Fraction) {
		return u;
	}

	if(u.kind() == Kind::Multiplication) {
		assert(u.size() == 2, "not a poly expr");

		if(u[1][0] == x) {
			Expr d = u[1][1];

			assert(d.kind() == Kind::Integer, "not a poyl expr");

			if (d == 0) return raisePolyExpr(0, 0, x);

			Expr c = mulPolyExpr(d, u[0]);

			if(was_diff) *was_diff = true;

			return c * power(x, d.value() - 1);
		}

		bool cond = false;

		Expr udx = diffPolyExprRec(u[0], x, &cond);

		if(cond == false) return raisePolyExpr(0, 0, x);

		return udx*u[1];
	}

  Expr g = Expr(Kind::Addition);

	for (Int i = 0; i < u.size(); i++) {
		Expr c = diffPolyExprRec(u[i], x, was_diff);

		if(isZeroPolyExpr(c)) continue;

		g.insert(c);
  }

	if(g.size() == 0) {
		return raisePolyExpr(0, 0, x);
	}

	return g;
}

Expr diffPolyExpr(Expr &u, Expr& x) {
	return diffPolyExprRec(u, x);
}

} // namespace polynomial
