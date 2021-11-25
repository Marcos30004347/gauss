#include "Addition.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/List.hpp"
//#include "Core/Expand/Expand.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"
#include "Rationals.hpp"

#include <vector>

using namespace ast;
// using namespace expand;
using namespace algebra;

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

} // namespace simplification
