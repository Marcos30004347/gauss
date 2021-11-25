#include "Multiplication.hpp"
#include "Addition.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/List.hpp"
//#include "Core/Expand/Expand.hpp"
#include "Power.hpp"
#include "Rationals.hpp"

#include <vector>

using namespace ast;
// using namespace expand;
using namespace algebra;

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
      Expr P = reducePowerAST(power(base(u1), reduceAdditionAST(O)));

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

} // namespace simplification
