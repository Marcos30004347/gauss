#include <assert.h>

#include "Multiplication.hpp"

using namespace ast;
using namespace algebra;

namespace expand {

Expr expandProductOfSums(Expr a, Expr b) {
  assert(a.kind() == Kind::Addition || b.kind() == Kind::Addition);

  if (a.kind() == Kind::Addition && b.kind() == Kind::Addition) {
    Expr u = Expr(Kind::Addition);

    for (unsigned int i = 0; i < a.size(); i++) {
      for (unsigned int j = 0; j < b.size(); j++) {
        u.insert(mul({a[i], b[j]}));
      }
    }

    return {u};
  }

  if (a.kind() != Kind::Addition && b.kind() == Kind::Addition) {
    Expr u = Expr(Kind::Addition);

    for (unsigned int j = 0; j < b.size(); j++) {
      u.insert(mul({a, b[j]}));
    }

    return u;
  }

  Expr u = Expr(Kind::Addition);

  for (unsigned int j = 0; j < a.size(); j++) {
    u.insert(mul({a[j], b}));
  }

  return u;
}

Expr expandMultiplicationAST(Expr u) {

  if (u.size() == 1) {
    Expr res = u[0];
    return res;
  }

  Expr expanded = u;

  for (unsigned int i = 0; i < expanded.size(); i++) {
    if (expanded[i].kind() == Kind::Addition) {
      signed long no = expanded.size();

      unsigned int k = (i + 1) % no;

      if (i == k)
        continue;

      Expr a = expanded[i];
      Expr b = expanded[k];

      if (k > i) {
        expanded.remove(k);
        expanded.remove(i);
      } else {
        expanded.remove(i);
        expanded.remove(k);
      }

      i = -1;

      expanded.insert(expandProductOfSums(a, b), 0);
    }
  }

  return expanded;
}

} // namespace expand
