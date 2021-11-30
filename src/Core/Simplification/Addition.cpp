#include "Addition.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Rational/Rational.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"
#include "Rationals.hpp"
#include <cstddef>
#include <cstdio>
#include <ios>
#include <list>
#include <utility>
#include <vector>

using namespace ast;
using namespace algebra;
using namespace rational;

namespace simplification {

Expr simplifyAdditionRec(Expr L);

Expr mergeAdditions(Expr p, Expr q) {
  if (p.size() == 0) {
    return q;
  }

  if (q.size() == 0) {
    return p;
  }

  Expr L = list({p[0], q[0]});
  Expr H = simplifyAdditionRec(L);

  if (H.size() == 0) {
    return mergeAdditions(rest(p), rest(q));
  }

  if (H.size() == 1) {
    Expr R = mergeAdditions(rest(p), rest(q));
    return adjoin(H[0], R, simplifyAdditionRec);
  }

  if (H[0] == p[0]) {
    Expr mer = mergeAdditions(rest(p), q);
    Expr res = adjoin(p[0], mer, simplifyAdditionRec);

    return res;
  }

  Expr mer = mergeAdditions(p, rest(q));
  Expr res = adjoin(q[0], mer, simplifyAdditionRec);

  return res;
}

Expr nonConstantCoefficient(Expr a) {
  if (a.kind() == Kind::FunctionCall) {
    bool non_constant = false;

    for (long i = 0; i < a.size(); i++) {
      if (!isConstant(a[i])) {
        non_constant = true;
        break;
      }
    }

    if (non_constant) {
      return a;
    }

    return undefined();
  }

  if (a.kind() == Kind::Power) {
    if (!isConstant(a[0]) || !isConstant(a[1]))
      return a;
  }

  Expr res = Expr(Kind::Multiplication);

  for (long i = 0; i < a.size(); i++) {
    if (a[i].kind() == Kind::FunctionCall) {
      Expr k = nonConstantCoefficient(a[i]);
      if (k != undefined()) {
        res.insert(k);
      }
    } else if (!isConstant(a[i])) {
      res.insert(a[i]);
    }
  }

  if (res.size() == 0) {

    return undefined();
  }

  if (res.size() == 1) {
    Expr r = res[0];

    return r;
  }

  return res;
}

Expr constantCoefficient(Expr a) {
  if (a.kind() == Kind::FunctionCall) {
    bool non_constant = false;

    for (long i = 0; i < a.size(); i++) {
      if (!isConstant(a[i])) {
        non_constant = true;
        break;
      }
    }

    if (non_constant) {
      return integer(1);
    }

    return a;
  }

  if (a.kind() == Kind::Power) {
    if (!isConstant(a[0]) || !isConstant(a[1]))
      return integer(1);

    return a;
  }

  Expr res = Expr(Kind::Multiplication);

  for (long i = 0; i < a.size(); i++) {
    if (a[i].kind() == Kind::FunctionCall) {
      Expr k = constantCoefficient(a[i]);
      if (k.kind() == Kind::Integer && k.value() == 1) {

      } else {
        res.insert(k);
      }
    } else if (isConstant(a[i])) {
      res.insert(a[i]);
    }
  }

  if (res.size() == 0) {

    return integer(1);
  }

  if (res.size() == 1) {
    Expr r = res[0];

    return r;
  }

  if (res.size() > 1) {
    Expr old = res;
    res = reduceRNEAST(res);
  }

  return res;
}

Expr simplifyAdditionRec(Expr L) {
  if (L.size() == 2 && L[0].kind() != Kind::Addition &&
      L[1].kind() != Kind::Addition) {
    Expr u1 = L[0];
    Expr u2 = L[1];

    if (isConstant(u1) && isConstant(u2)) {
      Expr K = u1 + u2;

      Expr P = reduceRNEAST(K);

      if (P == 0) {
        return list({});
      }

      return list({P});
    }

    if (u2 == inf()) {

      if (u1 == -inf()) {
        return {undefined()};
      }

      return list({inf()});
    }

    if (u2 == -inf()) {
      if (u1 == inf()) {
        return list({undefined()});
      }

      return list({-inf()});
    }

    if (u1 == 0) {
      return list({u2});
    }

    if (u2 == 0) {
      return list({u1});
    }

    Expr nc_u1 = nonConstantCoefficient(u1);
    Expr nc_u2 = nonConstantCoefficient(u2);

    if (nc_u1 == nc_u2) {
      Expr A = constantCoefficient(u1) + constantCoefficient(u2);
      Expr B = reduceAdditionAST(A) * nonConstantCoefficient(u1);

      Expr P = reduceMultiplicationAST(B);

      if (P == 0) {
        return list({});
      }

      return list({P});
    }

    if (orderRelation(u2, u1)) {
      return list({u2, u1});
    }

    return list({u1, u2});
  }

  if (L.size() == 2 &&
      (L[0].kind() == Kind::Addition || L[1].kind() == Kind::Addition)) {

    Expr u1 = L[0];
    Expr u2 = L[1];

    if (u1.kind() == Kind::Addition && u2.kind() == Kind::Addition) {
      Expr U1 = list({});
      Expr U2 = list({});

      for (long i = 0; i < u1.size(); i++)
        U1.insert(u1[i]);

      for (long i = 0; i < u2.size(); i++)
        U2.insert(u2[i]);

      Expr L_ = mergeAdditions(U1, U2);

      return L_;
    }

    if (u1.kind() == Kind::Addition) {
      Expr U1 = list({});
      Expr U2 = list({});

      for (long i = 0; i < u1.size(); i++)
        U1.insert(u1[i]);

      U2.insert(u2);

      return mergeAdditions(U1, U2);
    }

    if (u2.kind() == Kind::Addition) {
      Expr U1 = list({});
      Expr U2 = list({});

      for (long i = 0; i < u2.size(); i++)
        U2.insert(u2[i]);

      U1.insert(u1);

      return mergeAdditions(U1, U2);
    }
  }

  Expr u1 = L[0];

  Expr restL = rest(L);

  Expr w = simplifyAdditionRec(restL);

  if (u1.kind() == Kind::Addition) {
    Expr U1 = list({});

    for (long i = 0; i < u1.size(); i++)
      U1.insert(u1[i]);

    return mergeAdditions(U1, w);
  }

  return mergeAdditions(list({u1}), w);
}

Expr reduceAdditionAST(Expr u) {
  if (u.kind() == Kind::Undefined) {
    return undefined();
  }

  if (u.size() == 1) {
    return u[0];
  }

  Expr L = list({});

  for (long i = 0; i < u.size(); i++) {
    L.insert(u[i]);
  }

  Expr R = simplifyAdditionRec(L);

  if (R.size() == 0) {

    return integer(0);
  }

  if (R.size() == 1) {
    return R[0];
  }

  Expr res = Expr(Kind::Addition);

  for (Int i = 0; i < R.size(); i++) {
    res.insert(R[i]);
  }

  return res;
}

Expr constProduct(Expr &u, Expr &v) {
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

  Expr num = constProduct(a, b);
  Expr den = constProduct(c, d);

  Int g = gcd(num.value(), den.value());

  return fraction(num.value() / g, den.value() / g);
}

Expr constAddition(Expr &u, Expr &v) {
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

  Expr e = constProduct(a, d);
  Expr f = constProduct(b, c);

  Expr num = constAddition(e, f);
  Expr den = constProduct(c, d);

  Int g = gcd(num.value(), den.value());

  return fraction(num.value() / g, den.value() / g);
}

Expr splitTerm(Expr &u) {
  if (isConstant(u))
    return list({u, 1});
  if (u.kind() == Kind::Addition)
    return list({1, u});

  if (u.kind() == Kind::Multiplication) {
    Expr c = 1;
    Expr n = 1;

    for (long i = 0; i < u.size(); i++) {
      if (isConstant(u[i]))
        c = constProduct(c, u[i]);
      else
        n = n == 1 ? u[i] : n * u[i];
    }

    return list({c, n});
  }

  return list({1, u});
}

bool addConstants(std::vector<Expr> &L, long l, long r) {
  if (L[l] == 0 && L[r] == 0)
    return false;

  Expr P = constAddition(L[l], L[r]);

  L[l] = P;
  L[r] = 0;

  return true;
}

bool addNonConstans(std::vector<Expr> &L, long l, long r) {
  Expr t0 = splitTerm(L[l]);
  Expr t1 = splitTerm(L[r]);

  if (t0[1] == t1[1]) {
    Expr c = constAddition(t0[0], t1[0]);

    if (c == 0)
      L[l] = 0;
    else if (c == 1)
      L[l] = t0[1];
    else
      L[l] = c * t0[1];

    L[r] = 0;

    return true;
  }

  return false;
}

void print_terms(std::vector<Expr> &L) {
  for (size_t i = 0; i < L.size(); i++) {
    printf("* %s ", L[i].toString().c_str());
  }
  printf("\n");
}
/*
void mergeAdditionExpr(std::vector<Expr> &L, std::vector<Expr> &temp, long l,
                       long m, long &r) {

  // printf("\n******\n***** merging %li %li %li\n******\n", l ,m , r);
  size_t left_pos = l;
  size_t left_end = m;

  size_t temp_pos = l;

  size_t righ_end = r;
  size_t righ_pos = m + 1;

  while (left_pos <= left_end && righ_pos <= righ_end) {

    if (orderRelation(L[left_pos], L[righ_pos])) {
      temp[temp_pos++] = std::move(L[left_pos++]);
    } else {
      temp[temp_pos++] = std::move(L[righ_pos++]);
    }
  }

  while (left_pos <= left_end) {
    temp[temp_pos++] = std::move(L[left_pos++]);
  }

  while (righ_pos <= righ_end) {
    temp[temp_pos++] = std::move(L[righ_pos++]);
  }

  size_t num = r - l + 1;

  for (size_t i = 0; i < num; i++) {
    L[righ_end] = std::move(temp[righ_end]);
    righ_end--;
  }
}

void sortAdditionExpr(std::vector<Expr> &L, std::vector<Expr> &tmp, long l,
                      long r) {
  if (l < r) {
    long m = l + (r - l) / 2;

    sortAdditionExpr(L, tmp, l, m);
    sortAdditionExpr(L, tmp, m + 1, r);

    mergeAdditionExpr(L, tmp, l, m, r);
  }
}
*/
void flatAddition(Expr &u, std::vector<Expr> &L) {
  if (u.kind() != Kind::Addition) {
    L.push_back(u);
    return;
  }

  for (long i = 0; i < u.size(); i++) {
    if (u[i].kind() == Kind::Addition) {
      flatAddition(u[i], L);
    } else
      L.push_back(u[i]);
  }
}

Expr reduceAdditionExpr(Expr &&u) {
  if (u.kind() != Kind::Addition)
    return Expr(u);

  if (u.size() == 0)
    return 0;

  if (u.size() == 1)
    return u[0];

  std::vector<Expr> L;

  flatAddition(u, L);

  sort(L);

  size_t left = 0;
  size_t righ = 1;

  size_t size = L.size();

  bool have_pos_inf = false;
  bool have_neg_inf = false;

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

      bool c0 = isConstant(L[left]);
      bool c1 = isConstant(L[righ]);

      if (c0 && c1) {
        merged = addConstants(L, left, righ);
      }

      if (!c0 && !c1) {
        merged = addNonConstans(L, left, righ);
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

  Expr S = Expr(Kind::Addition);

  for (size_t i = 0; i < L.size(); i++) {
    if (L[i] != 0)
      S.insert(L[i]);
  }

	if(S.size() == 0) {
		return 0;
	}

	if (S.size() == 1)
    return S[0];

  return S;
}

Expr reduceAdditionExpr(Expr &u) {
  if (u.kind() != Kind::Addition)
    return Expr(u);
  return reduceAdditionExpr(std::forward<Expr>(u));
}

} // namespace simplification
