#include "Polynomial.hpp"
#include <assert.h>

#include "MathSystem/Simplification/Simplification.hpp"
#include "MathSystem/Simplification/Subtraction.hpp"

using namespace ast;
using namespace algebra;
using namespace simplification;

namespace expand {

void getBinomialsParameters(Int n, Int sum, Int *out, Int index,
                            std::vector<std::vector<Int>> &res) {
  if (index > n || sum < 0)
    return;

  if (index == n) {
    if (sum == 0) {

      res.push_back(std::vector<Int>());

      for (long long i = 0; i < n; i++) {
        res[res.size() - 1].push_back(out[i]);
      }
    }
    return;
  }

  for (Int i = 0; i <= 9; i++) {
    out[index.longValue()] = i;
    getBinomialsParameters(n, sum - i, out, index + Int(1), res);
  }
}

std::vector<std::vector<Int>> findNDigitNums(Int n, Int sum) {
  std::vector<std::vector<Int>> res;
  Int out[n.longValue() + 1];

  for (Int i = 0; i <= 9; i++) {
    out[0] = i;
    getBinomialsParameters(n, sum - i, out, 1, res);
  }

  return res;
}

Expr expandMultinomial(Expr b, Expr e) {
  assert(b.kind() == Kind::Addition || b.kind() == Kind::Subtraction);
  assert(e.kind() == Kind::Integer && e.value() > 0);

  Expr k;

  if (b.kind() == Kind::Subtraction) {
    k = reduceSubtractionAST(b);
  } else {
    k = b;
  }

  long long m = k.size();
  Int n = e.value();

  Expr r = Expr(Kind::Addition);

  std::vector<Int> s = std::vector<Int>(m + 1, 0);

  s[m] = n;

  std::vector<std::vector<Int>> binomials = findNDigitNums(m, n);
  for (std::vector<Int> bi : binomials) {
    Expr p = Expr(Kind::Multiplication);

    p.insert(binomial(n, bi));

    for (int i = 0; i < m; i++) {
      p.insert(power(k[i], integer(bi[i])));
    }
    r.insert(p);
  }

  return reduceAST(r);
}

Expr expandMultinomialAST(Expr u) {
  Expr b = base(u);
  Expr e = expoent(u);

  if ((b.kind() == Kind::Addition || b.kind() == Kind::Subtraction) &&
      b.size() > 1 && e.kind() == Kind::Integer && e.value() > 0) {
    Expr res = expandMultinomial(b, e);

    return res;
  }

  return u;
}

} // namespace expand
