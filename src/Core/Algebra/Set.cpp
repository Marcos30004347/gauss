#include "Set.hpp"
#include "Core/Debug/Assert.hpp"

using namespace ast;

namespace algebra {

Expr set(std::vector<Expr> e) {
  Expr S = Expr(Kind::Set);

  for (long unsigned int i = 0; i < e.size(); i++) {
    S.insert(e[i]);
  }

  return S;
}

void combUtil(Expr* ans, Expr tmp, Expr n, Int left, Int k) {
  if (k == 0) {
    ans->insert(tmp);
    return;
  }

  for (Int i = left; i < n.size(); ++i) {
    tmp.insert(n[i]);
    combUtil(ans, tmp, n, i + Int(1), k - Int(1));
    tmp.remove(tmp.size() - 1);
  }
}

Expr combination(Expr n, Expr k) {
  assert(k.kind() == Kind::Integer, "k needs to be an integer!");

  Expr ans = set({});
  Expr tmp = set({});

  combUtil(&ans, tmp, n, 0, k.value());

  return ans;
}

Expr difference(Expr L, Expr M) {
  assert(L.kind() == Kind::Set, "L is not a Set!\n");
  assert(M.kind() == Kind::Set, "M is not a Set!\n");

  Expr S = Expr(Kind::Set);

  for (unsigned int i = 0; i < L.size(); i++) {
    bool inc = false;

    for (unsigned int j = 0; j < M.size(); j++) {
      if (L[i] == M[j]) {
        inc = true;
        break;
      }
    }

    if (inc)
      continue;
    S.insert(L[i]);
  }

  return S;
}

Expr unification(Expr L, Expr M) {
  Expr D = difference(M, L);
  Expr S = Expr(Kind::Set);

  for (unsigned int i = 0; i < L.size(); i++) {
    S.insert(L[i]);
  }

  for (unsigned int i = 0; i < D.size(); i++) {
    S.insert(D[i]);
  }

  return S;
}

Expr intersection(Expr L, Expr M) {
  return difference(L, difference(L, M));
}

bool exists(Expr L, Expr e) {
  for (unsigned int i = 0; i < L.size(); i++) {
    if (L[i] == e)
      return true;
  }

  return false;
}

Expr cleanUp(Expr C, Expr t) {

  assert(C.kind() == Kind::Set, "C is not a set!");
  assert(t.kind() == Kind::Set, "t is not a set!");

  Expr S = Expr(Kind::Set);

  for (unsigned int i = 0; i < C.size(); i++) {
    bool inc = false;

    Expr s = C[i];

    for (unsigned int j = 0; j < s.size(); j++) {
      for (unsigned int k = 0; k < t.size(); k++) {
        if (t[k] == s[j]) {
          inc = true;
          break;
        }
      }

      if (inc) {
        break;
      }
    }

    if (!inc) {
      S.insert(s);
    }
  }

  return S;
}

} // namespace algebra
