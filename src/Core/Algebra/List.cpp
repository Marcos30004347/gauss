#include "List.hpp"
#include "Core/Debug/Assert.hpp"

using namespace ast;

namespace algebra {

Expr list(std::vector<Expr> e) { return Expr(Kind::List, e); }

Expr remove(Expr L, Expr M) {
  assert(L.kind() == Kind::List, "L is not a List!\n");
  assert(M.kind() == Kind::List, "M is not a List!\n");

  Expr l = Expr(Kind::List);

  int p = 0;

  for (unsigned int i = 0; i < L.size(); i++) {
    bool inc = false;

    for (unsigned int j = p; j < M.size(); j++) {
      if (L[i] == (M[j])) {
        p = j + 1;
        inc = true;
        break;
      }
    }

    if (inc)
      continue;

    l.insert(L[i]);
  }

  return l;
}

Expr join(Expr L, Expr M) {
  assert(L.kind() == Kind::List, "L is not a list!\n");
  assert(M.kind() == Kind::List, "M is not a list!\n");

  Expr l = Expr(Kind::List);

  for (unsigned int i = 0; i < L.size(); i++) {
    l.insert(L[i]);
  }

  for (unsigned int j = 0; j < M.size(); j++) {
    bool inc = false;

    for (unsigned int i = 0; i < L.size(); i++) {
      if (L[i] == (M[j])) {
        inc = true;
        break;
      }
    }

    if (inc)
      continue;

    l.insert(M[j]);
  }

  return l;
}

Expr append(Expr L, Expr M) {
  assert(L.kind() == Kind::List, "L is not a list!\n");
  assert(M.kind() == Kind::List, "M is not a list!\n");

  Expr l = Expr(Kind::List);

  for (unsigned int i = 0; i < L.size(); i++) {
    l.insert(L[i]);
  }

  for (unsigned int j = 0; j < M.size(); j++) {
    l.insert(M[j]);
  }

  return l;
}

Expr adjoin(Expr x, Expr L, Expr (*f)(Expr const)) {
  assert(L.kind() == Kind::List, "L is not a list!\n");
  if (f) {
    if (L.size() > 0) {
      Expr K = list({x, L[0]});
      Expr r = f(K);
      Expr l = Expr(Kind::List);

      if (r.kind() == Kind::List) {
        for (unsigned int i = 0; i < r.size(); i++) {
          l.insert(r[i]);
        }
        for (unsigned int i = 1; i < L.size(); i++) {
          l.insert(L[i]);
        }

        return l;
      }

      l.insert(r);

      for (unsigned int i = 1; i < L.size(); i++) {
        l.insert(L[i]);
      }

      return l;
    }

    Expr l = Expr(Kind::List);
    l.insert(x);

    return l;
  }

  Expr l = Expr(Kind::List);

  l.insert(x);

  for (unsigned int i = 0; i < L.size(); i++) {
    l.insert(L[i]);
  }

  return l;
}

Expr first(Expr L) { return L[0]; }

Expr rest(Expr L, int i) {
  std::vector<Expr> ops = L.operands();
	int j = 0;

	while (j < i) {
    ops.erase(ops.begin());
		j++;
	}

  return list(ops);
}

} // namespace algebra
