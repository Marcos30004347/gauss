#include "Multiplication.hpp"
#include "Addition.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Rational/Rational.hpp"
#include "Power.hpp"
#include "Rationals.hpp"

#include <tuple>
#include <vector>

using namespace ast;
// using namespace expand;
using namespace algebra;
using namespace rational;

namespace simplification {

Expr simplifyProductRec(Expr L);

Expr mergeProducts(Expr p, Expr q) {
  if (p.size() == 0)
    return q;

  if (q.size() == 0)
    return p;

  Expr L = list({p[0], q[0]});

  Expr H = simplifyProductRec(L);

  if (H.size() == 0) {
    return mergeProducts(rest(p), rest(q));
  }

  if (H.size() == 1) {
    Expr R = mergeProducts(rest(p), rest(q));
    return adjoin(H[0], R, simplifyProductRec);
  }

  if (H[0] == p[0]) {
    Expr mer = mergeProducts(rest(p), q);
    Expr res = adjoin(p[0], mer, simplifyProductRec);

    return res;
  }

  Expr mer = mergeProducts(p, rest(q));
  Expr res = adjoin(q[0], mer, simplifyProductRec);

  return res;
}

Expr simplifyProductRec(Expr L) {

  if (L.size() == 2 && L[0].kind() != Kind::Multiplication &&
      L[1].kind() != Kind::Multiplication) {
    Expr u1 = L[0];
    Expr u2 = L[1];

    if (isConstant(u1) && isConstant(u2)) {
      Expr P_ = mul({u1, u2});
      Expr P = reduceRNEAST(P_);

      if (P.kind() == Kind::Integer && P.value() == 1) {

        return list({});
      }

      return list({P});
    }

    if (u1.kind() == Kind::Infinity) {
      if (u2.kind() == Kind::Integer && u2.value() == 0) {
        return list({undefined()});
      } else if (u2.kind() == Kind::Integer && u2.value() == -1) {
        return list({Expr(Kind::MinusInfinity)});
      } else {
        return list({Expr(Kind::Infinity)});
      }
    }

    if (u1.kind() == Kind::MinusInfinity) {
      if (u2.kind() == Kind::Integer && u2.value() == 0) {
        return list({undefined()});
      } else if (u2.kind() == Kind::Integer && u2.value() == -1) {
        return list({Expr(Kind::Infinity)});
      } else {
        return list({Expr(Kind::MinusInfinity)});
      }
    }

    if (u2.kind() == Kind::Infinity) {
      if (u1.kind() == Kind::Integer && u1.value() == 0) {
        return list({undefined()});
      } else if (u1.kind() == Kind::Integer && u1.value() == -1) {
        return list({Expr(Kind::MinusInfinity)});
      } else {
        return list({Expr(Kind::Infinity)});
      }
    }

    if (u2.kind() == Kind::MinusInfinity) {
      if (u1.kind() == Kind::Integer && u1.value() == 0) {
        return list({undefined()});
      } else if (u1.kind() == Kind::Integer && u1.value() == -1) {
        return list({Expr(Kind::Infinity)});
      } else {
        return list({Expr(Kind::MinusInfinity)});
      }
    }

    if (u1.kind() == Kind::Integer && u1.value() == 1) {
      return list({u2});
    }

    if (u1.kind() == Kind::Integer && u1.value() == 0) {
      return list({integer(0)});
    }

    if (u2.kind() == Kind::Integer && u2.value() == 1) {
      return list({u1});
    }

    if (u2.kind() == Kind::Integer && u2.value() == 0) {
      return list({integer(0)});
    }

    Expr base_u1 = base(u1);
    Expr base_u2 = base(u2);

    if (base_u1 == base_u2) {
      Expr O = expoent(u1) + expoent(u2);
      Expr P = reducePowerExpr(power(base(u1), reduceAdditionAST(O)));

      if (P == 1) {
        return list({});
      }

      return list({P});
    }

    if (orderRelation(u2, u1)) {
      return list({u2, u1});
    }

    return list({u1, u2});
  }

  if (L.size() == 2 && (L[0].kind() == Kind::Multiplication ||
                        L[1].kind() == Kind::Multiplication)) {
    Expr u1 = L[0];
    Expr u2 = L[1];

    if (u1.kind() == Kind::Multiplication &&
        u2.kind() == Kind::Multiplication) {

      Expr U1 = list({});
      Expr U2 = list({});

      for (unsigned int i = 0; i < u1.size(); i++)
        U1.insert(u1[i]);

      for (unsigned int i = 0; i < u2.size(); i++)
        U2.insert(u2[i]);

      return mergeProducts(U1, U2);
    }

    if (u1.kind() == Kind::Multiplication) {
      Expr U1 = list({});
      Expr U2 = list({});

      for (unsigned int i = 0; i < u1.size(); i++)
        U1.insert(u1[i]);

      U2.insert(u2);

      return mergeProducts(U1, U2);
    }

    if (u2.kind() == Kind::Multiplication) {
      Expr U1 = list({});
      Expr U2 = list({});

      for (unsigned int i = 0; i < u2.size(); i++)
        U2.insert(u2[i]);

      U1.insert(u1);

      return mergeProducts(U1, U2);
    }
  }

  Expr u1 = L[0];

  Expr w = simplifyProductRec(rest(L));

  if (u1.kind() == Kind::Multiplication) {
    Expr U1 = list({});

    for (unsigned int i = 0; i < u1.size(); i++)
      U1.insert(u1[i]);

    return mergeProducts(U1, w);
  }

  return mergeProducts(list({u1}), w);
}

Expr reduceMultiplicationAST(Expr u) {

  if (u == undefined())
    return undefined();
  for (unsigned int i = 0; i < u.size(); i++) {
    if (u[i] == 0)
      return u[i];
  }

  if (u.size() == 1)
    return u[0];

  Expr L = list({});

  for (unsigned int i = 0; i < u.size(); i++)
    L.insert(u[i]);

  Expr R = simplifyProductRec(L);

  if (R.size() == 1) {
    Expr r = R[0];

    return r;
  }

  if (R.size() == 0) {

    return integer(1);
  }

  Expr res = Expr(Kind::Multiplication);

  for (unsigned int i = 0; i < R.size(); i++) {
    res.insert(R[i]);
  }

  return res;
}

Expr constMultiplication(Expr &u, Expr &v) {
  if (u == 1)
    return v;
  if (v == 1)
    return u;

  if (u == 0)
    return 0;
  if (v == 0)
    return 0;

  if (u.kind() == Kind::Integer && v.kind() == Kind::Integer)
    return u.value() * v.value();

  Expr a = numerator(u);
  Expr b = numerator(v);
  Expr c = denominator(u);
  Expr d = denominator(v);

  Expr num = constMultiplication(a, b);
  Expr den = constMultiplication(c, d);

  Int g = gcd(num.value(), den.value());

  return fraction(num.value() / g, den.value() / g);
}

Expr constSummation(Expr &u, Expr &v) {
  if (u == 0)
    return v;
  if (v == 0)
    return u;

  if (u.kind() == Kind::Integer && v.kind() == Kind::Integer)
    return u.value() + v.value();

  Expr a = numerator(u);
  Expr b = numerator(v);
  Expr c = denominator(u);
  Expr d = denominator(v);

  Expr e = constMultiplication(a, d);
  Expr f = constMultiplication(b, c);

  Expr num = constSummation(e, f);
  Expr den = constMultiplication(c, d);

  Int g = gcd(num.value(), den.value());

  return fraction(num.value() / g, den.value() / g);
}

bool mulConstants(std::vector<Expr> &L, long l, long r) {
  if (L[l] == 0 && L[r] == 0)
    return false;

  Expr P = constMultiplication(L[l], L[r]);

  L[l] = P;
  L[r] = 1;

  return true;
}

Expr split(Expr &u) {
  if (isConstant(u))
    return list({u, 1});

  if (u.kind() == Kind::Addition)
    return list({1, u});

  if (u.kind() == Kind::Multiplication) {

    Expr c = 1;
    Expr n = 1;

    for (long i = 0; i < u.size(); i++) {
      if (isConstant(u[i]))
        c = constMultiplication(c, u[i]);
      else
        n = n == 1 ? u[i] : n * u[i];
    }

    return list({c, n});
  }

  return list({1, u});
}

bool mulNonConstans(std::vector<Expr> &L, long l, long r) {

  Expr t0 = split(L[l]);
  Expr t1 = split(L[r]);

  if (base(t0[1]) == base(t1[1])) {
    Expr a = expoent(t0[1]);
    Expr b = expoent(t1[1]);

    Expr c = constMultiplication(t0[0], t1[0]);
    Expr d = constSummation(a, b);
    Expr e = base(t0[1]);

    if (d == 0)
      L[l] = 1;
    if (d == 1)
      L[l] = e;
    else if (c == 0)
      L[l] = 0;
    else if (c == 1)
      L[l] = power(e, d);
    else
      L[l] = c * power(e, d);

    L[r] = 1;

    return true;
  }

  return false;
}

void flatMultiplications(Expr &u, std::vector<Expr> &L) {
  if (u.kind() != Kind::Multiplication) {
    L.push_back(u);
    return;
  }

  for (long i = 0; i < u.size(); i++) {
    if (u[i].kind() == Kind::Multiplication) {
      flatMultiplications(u[i], L);
    } else
      L.push_back(u[i]);
  }
}

Expr reduceMultiplicationExpr(Expr &&u) {
  if (u.kind() != Kind::Multiplication)
    return Expr(u);

  if (u.size() == 0)
    return 0;
  if (u.size() == 1)
    return u[0];

  std::vector<Expr> L;

  flatMultiplications(u, L);

  sort(L, true);

  size_t left = 0;
  size_t righ = 1;

  size_t size = L.size();

  bool have_pos_inf = false;
  bool have_neg_inf = false;
  bool have_zero = false;
  while (righ < size) {
    bool merged = true;

    while (merged && righ < size) {
      merged = false;

      if (L[left] == inf() || L[righ] == inf()) {
        have_pos_inf = true;

        break;
      }

      if (L[left] == -inf() || L[righ] == -inf()) {
        have_neg_inf = true;

        break;
      }

      if (L[left] == 0 || L[righ] == 0) {
        have_zero = true;

        break;
      }

      bool c0 = isConstant(L[left]);
      bool c1 = isConstant(L[righ]);

      if (c0 && c1) {
        merged = mulConstants(L, left, righ);
      }

      if (!c0 && !c1) {
        merged = mulNonConstans(L, left, righ);
      }

      righ += merged;
    }

    left = righ;
    righ = righ + 1;
  }

  if (have_pos_inf && have_neg_inf)
    return undefined();
  if (have_pos_inf)
    return +inf();
  if (have_neg_inf)
    return -inf();
  if (have_zero)
    return 0;

  Expr S = Expr(Kind::Multiplication);

  for (size_t i = 0; i < L.size(); i++) {
    if (L[i] != 1)
      S.insert(L[i]);
  }

  if (S.size() == 0)
    return 1;
  if (S.size() == 1)
    return S[0];

  return S;
}

Expr reduceMultiplicationExpr(Expr &u) {
  if (u.kind() != Kind::Multiplication)
    return Expr(u);
  return reduceMultiplicationExpr(std::forward<Expr>(u));
}

} // namespace simplification
