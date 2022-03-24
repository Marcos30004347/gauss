#include "Polynomial.hpp"

#include "gauss/Algebra/Sorting.hpp"
#include "gauss/Algebra/Reduction.hpp"
#include "gauss/Algebra/Expression.hpp"
#include "gauss/Error/error.hpp"
#include "gauss/Factorization/Wang.hpp"
#include "gauss/GaloisField/GaloisField.hpp"
#include "Resultant.hpp"

#include <climits>
#include <cstddef>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

using namespace alg;
using namespace galoisField;

namespace poly {

 expr degreeGME(expr u, set S) {
  if (u.kind() == kind::INT && u.value() == 0)
    return -inf();

  if (is(&u, kind::CONST))
    return 0;

  if (exists(S, u)) {
    return 1;
  } else if (u.kind() == kind::POW) {
    expr b = u[0];
    expr e = u[1];

    if (exists(S, b) && is(&e, kind::CONST)) {
      return e;
    }

  } else if (u.kind() == kind::MUL) {
    expr deg = 0;
    for (unsigned int i = 0; i < u.size(); i++) {
      expr deg_ = degreeGME(u[i], S);
      if (deg_.value() > deg.value()) {

        deg = deg_;
      } else {
      }
    }

    return deg;
  }

  return 0;
}

expr degree(expr u, expr v) {
  if (u.kind() == kind::INT && u.value() == 0) {
    return -inf();
  }

	set S = set({ v });

  if (u.kind() != kind::ADD && u.kind() != kind::SUB) {
    expr r = degreeGME(u, S);

    return r;
  }

  if (exists(S, u)) {
    return 1;
  }

  expr deg = 0;

  for (unsigned int i = 0; i < u.size(); i++) {
    expr deg_ = degreeGME(u[i], S);

    if (deg_.value() > deg.value()) {
      deg = deg_;
    }
  }

  return deg;
}

list coefficientGME(expr u, expr x) {
	if (u == x) {
		return list({ 1, 1 });
  }

  if (u.kind() == kind::POW) {
    expr b = u[0];
    expr e = u[1];

    if (b == (x) && e.kind() == kind::INT && e.value() > 0) {
      return list({ 1, e });
    }
  }

	if (u.kind() == kind::MUL) {
    expr m = 0;
    expr c = u;

		for (size_t i = 0; i < u.size(); i++) {
      list f = coefficientGME(u[i], x);

			if (f.size() == 0) {
        return f;
      }

      if (f[1] != 0) {
        m = f[1];
				c = expand(u / pow(x, m));
      }
    }

    return list({c, m});
  }

  if (u.freeOf(x)) {
		return list({ u, 0 });
  }

  return list({});
}

expr coeff(expr u, expr x, expr d) {
	if (u.kind() != kind::ADD && u.kind() != kind::SUB) {
		list f = coefficientGME(u, x);

		if (f.size() == 0) {
			return undefined();
		}

    if (d == f[1]) {
      return f[0];
    }

    return 0;
  }

  if (x == u) {
    if (d == 1) {
      return 1;
    }

    return 0;
  }

  expr c = undefined();

	for (size_t i = 0; i < u.size(); i++) {
		list f = coefficientGME(u[i], x);

    if (f.size() == 0) continue;

    if (d == f[1]) {
      if (c == undefined()) {
        c = expr(u.kind());
      }

			c.insert(f[0]);
    }
  }

  return reduce(c);
}


long int sortSplit(list *a, list* d, long l, long r) {
  long int i = l - 1;

  expr *p = &d->members[r];

  for (long int j = l; j < r; j++) {
    if (compare(&d->members[j], p, kind::ADD) < 0) {
			i = i + 1;

			std::swap(d->members[i], d->members[j]);
			std::swap(a->members[i], a->members[j]);
    }
  }

	i = i + 1;

	std::swap(a->members[i], a->members[r]);
	std::swap(d->members[i], d->members[r]);

  return i;
}

void sortFreeVarsByDegrees(list *a, list* d, long int l, long int r) {
  if (l < r) {
    long int m = sortSplit(a, d, l, r);

    sortFreeVarsByDegrees(a, d, l, m - 1);
    sortFreeVarsByDegrees(a, d, m + 1, r);
  }
}

list getVariableListForPolyExpr(expr a) {
	list f = freeVariables(a);
	list d = {};

	for(size_t i = 0; i < f.size(); i++) {
		expr dg = degree(a, f[i]);

		assert(is(&dg, kind::FRAC | kind::INT));
		d.insert(dg);
	}

	sortFreeVarsByDegrees(&f, &d, 0, f.size() - 1);

	return f;
}

Int collectDegree(expr &u, expr &x) {
  if (u.kind() == kind::INT || u.kind() == kind::FRAC) {
    return 0;
  }

  if (u.kind() == kind::SYM) {
    if (u.identifier() == x.identifier()) {
      return 1;
    }

    return 0;
  }

  if (u.kind() == kind::POW) {
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

expr collectCoeff(expr &u, expr &x, Int d) {
  if (u == x && d == 1)
    return 1;
  if (u.kind() == x.kind()) {
    if (d == 0 && u != x) {
      return u;
    }

    return 0;
  }

  if (u.kind() == kind::POW) {
    if (d == 0) {
      if (u[0] == x)
        return u[1] == 0 ? 1 : 0;
      return u;
    }
    return u[0] == x && u[1] == d ? 1 : 0;
  }

  if (u.kind() == kind::MUL) {
    expr c = expr(kind::MUL);
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

expr collectRec(expr &u, expr &L, Int i) {
  if (! is(&L, kind::LIST)) {
		raise(error(ErrorCode::ARG_IS_NOT_LIST_EXPR, 0));
  }

  if (is(&u, kind::DIV)) {
		raise(error(ErrorCode::ARG_IS_NOT_POLY_EXPR, 0));
  }

  if (i == L.size()) {
    if (i == 0 && u.kind() == kind::POW) {
      return expr(kind::ADD, {1 * u});
    }

    if (i == 0 && u.kind() == kind::MUL) {
      return expr(kind::ADD, {u});
    }

    if (i == 0 && u.kind() == kind::SYM) {
      return expr(kind::ADD, {pow(u, 0)});
    }

    if (i == 0 && u.kind() == kind::FUNC) {
      return expr(kind::ADD, {pow(u, 0)});
    }

    return u;
  }

  if (u.kind() == kind::MUL && u.size() == 2 && u[1].kind() == kind::POW &&
      u[1][0] == L[i]) {
    return expr(kind::ADD, {u});
  }

  expr c = 1;

  if (u.kind() == kind::POW && u[0] == L[i]) {
    return expr(kind::ADD, {collectRec(c, L, i + 1) * u});
  }

  if (u.kind() == kind::SYM || u.kind() == kind::FUNC) {
    if (u == L[i]) {
      return create(kind::ADD, {collectRec(c, L, i + 1) * pow(L[i], 1)});
    }

    return create(kind::ADD, {collectRec(u, L, i + 1) * pow(L[i], 0)});
  }

  Int d = collectDegree(u, L[i]);

  if (u.kind() == kind::MUL) {
    Int k = collectDegree(u, L[i]);
    expr c = collectCoeff(u, L[i], k);

    return expr(kind::ADD, {collectRec(c, L, i + 1) * pow(L[i], k)});
  }

  if (u.kind() == kind::ADD) {
    std::map<Int, expr> coeffs;

    for (size_t j = 0; j < u.size(); j++) {
      Int k = collectDegree(u[j], L[i]);
      expr c = collectCoeff(u[j], L[i], k);
      if (c == 0)
        continue;

      if (coeffs.count(k) == 0)
        coeffs.insert(std::pair<Int, expr>({k, c}));
      else
        coeffs[k] = coeffs[k] + c;
    }

    expr g = expr(kind::ADD);

    for (std::map<Int, expr>::iterator it = coeffs.begin(); it != coeffs.end();
         it++) {
      if (it->second != 0) {
        g.insert(collectRec(it->second, L, i + 1) * pow(L[i], it->first));
      }
    }

    return g;
  }

  return expr(kind::ADD, {collectRec(u, L, i + 1) * pow(L[i], d)});
}

expr polyExpr(expr &&u, expr &&L) { return collectRec(u, L, 0); }

expr polyExpr(expr &u, expr &L) { return collectRec(u, L, 0); }

expr polyExpr(expr &&u, expr &L) { return collectRec(u, L, 0); }

expr polyExpr(expr &u, expr &&L) { return collectRec(u, L, 0); }

expr groundLeadCoeffPolyExpr(expr u) {
  if (u.kind() == kind::INT || u.kind() == kind::FRAC) {
    return u;
  }

  return groundLeadCoeffPolyExpr(leadCoeffPolyExpr(u));
}

expr monicPolyExpr(expr u, expr L, expr K) {
  expr lc = groundLeadCoeffPolyExpr(u);
  return quoPolyExpr(u, lc, L, K);
}

expr constDenLcmPolyExprRec(expr u, expr L, unsigned int j) {
  if (L.size() == j) {
    assert(u.kind() == kind::INT || u.kind() == kind::FRAC);

    if (u.kind() == kind::FRAC) {
      return u[1];
    }

    return 1;
  }

  expr g = 1;

  for (Int i = 0; i < u.size(); i++) {
    expr t = constDenLcmPolyExprRec(u[i][0], L, j + 1);
    g = lcm(g, t);
  }

  return g;
}

expr constDenLcmPolyExpr(expr u, expr L) {
  return constDenLcmPolyExprRec(u, L, 0);
}

expr removeDenominatorsPolyExpr(expr u, expr L, expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q");
  expr f = constDenLcmPolyExpr(u, L);
  expr t = polyExpr(f, L);
  expr k = mulPolyExpr(u, t);
  return list({t, k});
}

expr igcdPolyExpr(expr u, expr v, expr L, expr K) {
  if (L.size() == 0) {
    assert(u.kind() == kind::INT);
    assert(v.kind() == kind::INT);

    return abs(gcd(u.value(), v.value()));
  }

  expr U = contAndPpPolyExpr(u, L, K);

  expr u_cnt = U[0];
  expr u_ppr = U[1];

  expr V = contAndPpPolyExpr(v, L, K);

  expr v_cnt = V[0];
  expr v_ppr = V[1];

  expr R = rest(L);

  expr h = remSeqPolyExpr(u_ppr, v_ppr, L, K)[0];

  h = ppPolyExpr(h, L, K);

  expr c = raisePolyExpr(gcdPolyExpr(u_cnt, v_cnt, R, K), 0, L[0]);
  h = mulPolyExpr(h, c);

  expr M = monicPolyExpr(h, L, K);

  return M;
}

expr gcdPolyExpr(expr u, expr v, expr L, expr K) {
  expr a = 1;
  expr b = 1;

  if (K.identifier() == "Q") {
    u = removeDenominatorsPolyExpr(u, L, K)[1];
    v = removeDenominatorsPolyExpr(v, L, K)[1];
  }

  // if(isZeroPolyExpr(u)) return v;
  // if(isZeroPolyExpr(v)) return u;

  expr Z = expr("Z");
  expr H = heuristicGcdPolyExpr(u, v, L, Z);

  if (H != fail()) {
    if (K == Z) {
      return H[0];
    }

    return monicPolyExpr(H[0], L, K);
  }
  expr g = igcdPolyExpr(u, v, L, Z);

  if (K == Z) {
    return g;
  }

  return monicPolyExpr(g, L, K);
}

bool isZeroPolyExpr(expr &u) {
  if (is(&u, kind::TERMINAL)) {
    return u.kind() == kind::INT && u.value() == 0;
  }

  if (u.kind() == kind::ADD && u.size() > 1) {
    return 0;
  }

  return isZeroPolyExpr(u[0]);
}

bool isConstantPolyExpr(expr &u, expr v) {
  if (v == 0)
    return isZeroPolyExpr(u);

  if (is(&u, kind::TERMINAL)) {
    return (u.kind() == kind::INT || u.kind() == kind::FRAC) && u == v;
  }

  if (u.kind() == kind::MUL) {
    if (u.size() > 2)
      return false;
    if (u.size() == 2) {
      assert(u[1].kind() == kind::POW);
      assert(u[1][1].kind() == kind::INT);

      if (u[1][1].value() != 0) {
        return false;
      }
    }
  }

  if (u.kind() == kind::ADD && u.size() > 1) {
    return 0;
  }

  return isConstantPolyExpr(u[0], v);
}

bool isConstantPolyExpr(expr &u) {
  if (is(&u, kind::TERMINAL)) {
    return u.kind() == kind::INT || u.kind() == kind::FRAC;
  }

  if (u.kind() == kind::ADD && u.size() > 1) {
    return false;
  }

  if (u.kind() == kind::MUL) {
    assert(u.size() == 2);

    expr p = degree(u[1]);

    if (p != 0 && p != -inf())
      return false;
  }

  return isConstantPolyExpr(u[0]);
}

expr mulPolyExpr(expr &&p1, expr &&p2) {
  if (is(&p1, kind::CONST) && is(&p2, kind::CONST)) {
    return reduce(p1 * p2);
  }

  if (is(&p2, kind::CONST)) {
    return mulPolyExpr(p2, p1);
  }

  if (is(&p1, kind::CONST)) {
    if (p2.kind() == kind::MUL) {
      return mulPolyExpr(p1, p2[0]) * p2[1];
    }

    if (p2.kind() == kind::ADD) {
      expr g = expr(kind::ADD);

      for (size_t i = 0; i < p2.size(); i++) {
        g.insert(mulPolyExpr(p1, p2[i][0]) * p2[i][1]);
      }

      return g;
    }
  }

  expr x = p1[0][1][0];
  std::map<Int, expr> coeffs;

  for (size_t i = 0; i < p1.size(); ++i) {
    assert(p1[i][1][0] == x);

    expr u = p1[i];

    for (size_t j = 0; j < p2.size(); j++) {
      assert(p2[j][1][0] == x);

      expr v = p2[j];

      Int e = u[1][1].value() + v[1][1].value();

      expr c = mulPolyExpr(u[0], v[0]);

      if (coeffs.count(e) == 0) {
        coeffs[e] = c;
      } else {
        coeffs[e] = addPolyExpr(coeffs[e], c);
      }
    }
  }

  expr g = expr(kind::ADD);

  for (std::map<Int, expr>::iterator it = coeffs.begin(); it != coeffs.end();
       it++) {
    if (!isZeroPolyExpr(it->second) || g.size() == 0) {
      g.insert(it->second * pow(x, it->first));
    }
  }

  return g;
}

expr mulPolyExpr(expr &p1, expr &p2) {
  return mulPolyExpr(std::forward<expr>(p1), std::forward<expr>(p2));
}

expr addColPolyRec(expr &u, expr &v, unsigned int i = 0, unsigned int j = 0) {
  if (is(&u, kind::CONST) && is(&v, kind::CONST)) {
    return reduce(u + v);
  }

  if ((is(&u, kind::ADD) && i == u.size()) &&
      (is(&v, kind::ADD) && j == v.size())) {
    return 0;
  }

  if (is(&u, kind::ADD) && u.size() == 0) {
    return v;
  }

  if (is(&v, kind::ADD) && v.size() == 0) {
    return u;
  }

  if (is(&u, kind::ADD) && i == u.size()) {
    expr r = create(kind::ADD);

    for (size_t t = j; t < v.size(); t++) {
      r.insert(v[t]);
    }

    return r;
  }

  if (is(&v, kind::ADD) && j == v.size()) {
    expr r = create(kind::ADD);

    for (size_t t = i; t < u.size(); t++) {
      r.insert(u[t]);
    }

    return r;
  }

  if (u[i].kind() == kind::MUL && u[i].size() == 2 &&
      v[j].kind() == kind::MUL && v[j].size() == 2) {

    expr &ucoeff = u[i][0];
    expr &upower = u[i][1];

    expr &vcoeff = v[j][0];
    expr &vpower = v[j][1];

    if (upower[0] == vpower[0]) {

      bool uzero = isZeroPolyExpr(ucoeff);
      bool vzero = isZeroPolyExpr(vcoeff);

      if (uzero && vzero) {
        expr a = addColPolyRec(u, v, i + 1, j + 1);

        if (isZeroPolyExpr(a)) {
          return expr(kind::ADD, {0 * pow(upower[0], 0)});
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
        expr a = addPolyExpr(ucoeff, vcoeff);

        if (!isZeroPolyExpr(a)) {

          if (a.kind() == kind::MUL) {
            a = expr(kind::ADD, {a});
          }

          a = a * upower;
        }

        expr b = addColPolyRec(u, v, i + 1, j + 1);

        if (isZeroPolyExpr(b)) {
          if (isZeroPolyExpr(a)) {
            return expr(kind::ADD, {0 * pow(upower[0], 1)});
          }

          return expr(kind::ADD, {a});
        }

        if (isZeroPolyExpr(a))
          return b;

        if (b.kind() == kind::ADD) {
          b.insert(a, 0);
          return b;
        }

        return a + b;
      }

      if (upower[1].value() < vpower[1].value()) {
        expr a = addColPolyRec(u, v, i + 1, j);

        if (isZeroPolyExpr(a))
          return expr(kind::ADD, {u[i]});

        if (a.kind() == kind::ADD) {
          a.insert(u[i], 0);
          return a;
        }

        return u[i] + a;
      }

      expr a = addColPolyRec(u, v, i, j + 1);
      expr c = expr(kind::ADD, {v[j]});

      if (isZeroPolyExpr(a)) {
        return c;
      }

      for (Int i = 0; i < a.size(); i++) {
        c.insert(a[i]);
      }

      return c;
    }
  }

  expr a = addColPolyRec(u, v, i + 1, j + 1);

  expr b;

  if ((u[i].kind() == kind::INT || u[i].kind() == kind::FRAC) &&
      (v[j].kind() == kind::INT || v[j].kind() == kind::FRAC)) {
    b = reduce(u[i] + v[j]);
  } else {
    b = u[i] + v[j];
  }

  if (isZeroPolyExpr(a))
    return b;

  if (a.kind() == kind::ADD) {
    a.insert(b, 0);
    return a;
  }

  return b + a;
}

expr addPolyExpr(expr &u, expr &v) { return addColPolyRec(u, v, 0, 0); }

expr addPolyExpr(expr &&u, expr &&v) {
  if (is(&u, kind::CONST) && is(&v, kind::CONST)) {
    return reduce(u + v);
  }

  return addColPolyRec(u, v, 0, 0);
}

expr subColPolyRec(expr &u, expr &v, unsigned int i = 0, unsigned int j = 0) {
  expr k = -1;

  if (i == u.size() && j == v.size())
    return 0;

  if (i == u.size()) {
    expr g = expr(kind::ADD);

    for (size_t t = j; t < v.size(); t++) {
      g.insert(v[t]);
    }

    return mulPolyExpr(k, g);
  }

  if (j == v.size()) {
    expr q = create(kind::ADD);

    for (size_t t = i; t < u.size(); t++) {
      q.insert(u[t]);
    }

    return q;
  }

  if (u[i].kind() == kind::MUL && u[i].size() == 2 &&
      v[j].kind() == kind::MUL && v[j].size() == 2) {

    expr &ucoeff = u[i][0];
    expr &upower = u[i][1];

    expr &vcoeff = v[j][0];
    expr &vpower = v[j][1];

    if (upower[0] == vpower[0]) {
      bool uzero = isZeroPolyExpr(ucoeff);
      bool vzero = isZeroPolyExpr(vcoeff);

      if (uzero && vzero) {
        expr a = subColPolyRec(u, v, i + 1, j + 1);

        if (isZeroPolyExpr(a)) {
          return expr(kind::ADD, {0 * pow(upower[0], 1)});
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
        expr a = subPolyExpr(ucoeff, vcoeff);

        if (!isZeroPolyExpr(a)) {

          if (a.kind() == kind::MUL) {
            a = expr(kind::ADD, {a});
          }

          a = a * upower;
        }

        expr b = subColPolyRec(u, v, i + 1, j + 1);

        if (isZeroPolyExpr(b)) {

          if (isZeroPolyExpr(a)) {
            return expr(kind::ADD, {0 * pow(upower[0], 1)});
          }

          return expr(kind::ADD, {a});
        }

        if (isZeroPolyExpr(a))
          return b;

        if (b.kind() == kind::ADD) {
          b.insert(a, 0);
          return b;
        }

        return a + b;
      }

      if (upower[1].value() < vpower[1].value()) {
        expr a = subColPolyRec(u, v, i + 1, j);

        if (isZeroPolyExpr(a)) {
          return expr(kind::ADD, {u[i]});
        }

        if (a.kind() == kind::ADD) {
          a.insert(u[i], 0);
          return a;
        }

        return u[i] + a;
      }

      expr a = subColPolyRec(u, v, i, j + 1);

      expr c = mulPolyExpr(-1, expr(kind::ADD, {v[j]}));

      if (isZeroPolyExpr(a)) {
        return c;
      }
      for (Int i = 0; i < a.size(); i++) {
        c.insert(a[i]);
      }

      return c;
    }
  }

  expr a = subColPolyRec(u, v, i + 1, j + 1);

  expr b;

  if (is(&u[i], kind::CONST) && is(&v[j], kind::CONST)) {
    b = reduce(u[i] + -v[j]);
  } else {
    b = u[i] + mulPolyExpr(k, v[j]);
  }

  if (isZeroPolyExpr(a)) {
    return b;
  }

  if (a.kind() == kind::ADD) {
    a.insert(b, 0);
    return a;
  }

  return b + mulPolyExpr(k, a);
}

expr subPolyExpr(expr &u, expr &v) {

  if (is(&u, kind::CONST) && is(&v, kind::CONST)) {
    return reduce(u - v);
  }

  return subColPolyRec(u, v, 0, 0);
}

expr subPolyExpr(expr &&u, expr &&v) {

  if (is(&u, kind::CONST) && is(&v, kind::CONST)) {
    return reduce(u - v);
  }

  return subColPolyRec(u, v, 0, 0);
}

expr raiseToExpression(expr c, expr u) {
  if (is(&u, kind::CONST)) {
    return c;
  }

  if (u.kind() != kind::ADD || u.size() == 0) {
    return u;
  }

  expr p = u[0][0];
  expr k = u[0][1];

  return expr(kind::ADD, {raiseToExpression(c, p) * pow(base(k), 0)});
}

expr normalizeToPolyExprs(expr u, expr v) {
	expr t = reduce(expand(u) + expand(v));

	expr L = getVariableListForPolyExpr(t);

	return list({ L, polyExpr(u, L), polyExpr(v, L) });
}


expr powPolyExpr(expr &u, Int v) {
  if (v < 0) {
    return 1 / powPolyExpr(u, -v);
  }

  expr g = raiseToExpression(1, u);

  expr x = u;

  while (v) {
    if (v % 2 == 1) {
      g = mulPolyExpr(g, x);
    }

    v = v / 2;

    x = mulPolyExpr(x, x);
  }

  return g;
}

expr powPolyExpr(expr &&x, Int v) {
  if (v < 0) {
    return 1 / powPolyExpr(x, -v);
  }

  expr g = 1;

  while (v) {
    if (v % 2 == 1)
      g = mulPolyExpr(g, x);

    v = v / 2;

    x = mulPolyExpr(x, x);
  }

  return g;
}

expr leadCoeffPolyExpr(expr &u) {
  if (is(&u, kind::TERMINAL))
    return u;

  assert(u.kind() == kind::ADD);

  return u[u.size() - 1][0];
}

expr degreePolyExpr(expr &u) {
  if (isZeroPolyExpr(u))
    return -inf();

  if (u.kind() == kind::INT)
    return 0;
  if (u.kind() == kind::FRAC)
    return 0;

  assert(u.kind() == kind::ADD);

  return u[u.size() - 1][1][1];
}

expr raiseToPolyExpr(expr &x, Int e, expr &L) { return polyExpr(pow(x, e), L); }

expr raiseToPolyExpr(expr &&x, Int e, expr &&L) {
  return polyExpr(pow(x, e), L);
}

expr raisePolyExpr(expr &u, Int exp, expr &x) {
  return create(kind::ADD, {u * pow(x, exp)});
}

expr raisePolyExpr(expr &&u, Int exp, expr &x) {
  return create(kind::ADD, {u * pow(x, exp)});
}

expr divPolyExpr(expr &&u, expr &&v, expr &L, expr &K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q");

  if (L.size() == 0) {
    expr d = reduce(u / v);
    if (K.identifier() == "Z") {
      if (d.kind() == kind::INT) {
        return list({d, 0});
      }

      return list({0, u});
    }

    return list({d, 0});
  }

  expr r = u;

  expr m = degreePolyExpr(r);
  expr n = degreePolyExpr(v);

  expr q = polyExpr(0, L);

  expr lcv = leadCoeffPolyExpr(v);

  expr R = rest(L);

  while (m != -inf() && m.value() >= n.value()) {
    expr lcr = leadCoeffPolyExpr(r);

    expr d = divPolyExpr(lcr, lcv, R, K);

    if (!isZeroPolyExpr(d[1])) {
      return list({q, r});
    }

    expr k = raisePolyExpr(d[0], m.value() - n.value(), L[0]);

    q = addColPolyRec(q, k);

    expr g = raisePolyExpr(d[0], 0, L[0]);

    expr j = raiseToPolyExpr(L[0], m.value() - n.value(), L);

    expr t1 = mulPolyExpr(v, g);
    expr t2 = mulPolyExpr(t1, j);

    r = subPolyExpr(r, t2);

    m = degreePolyExpr(r);
  }

  return list({q, r});
}

expr divPolyExpr(expr &u, expr &v, expr &L, expr &K) {
  return divPolyExpr(std::forward<expr>(u), std::forward<expr>(v), L, K);
}

expr quoPolyExpr(expr &u, expr &v, expr &L, expr &K) {
  return divPolyExpr(std::forward<expr>(u), std::forward<expr>(v), L, K)[0];
}

expr remPolyExpr(expr &u, expr &v, expr &L, expr &K) {
  return divPolyExpr(std::forward<expr>(u), std::forward<expr>(v), L, K)[1];
}

expr pseudoDivPolyExpr(expr &&u, expr &&v, expr &L) {
  expr p = polyExpr(0, L);

  expr s = u;

  expr m = degreePolyExpr(s);
  expr n = degreePolyExpr(v);

  Int t = 0;

  Int d = max(m.value() - n.value() + 1, 0);

  expr lv = raisePolyExpr(leadCoeffPolyExpr(v), 0, L[0]);

  while (m != -inf() && m.value() >= n.value()) {
    expr ls = raisePolyExpr(leadCoeffPolyExpr(s), 0, L[0]);
    expr jx = raiseToPolyExpr(L[0], m.value() - n.value(), L);

    expr t1 = mulPolyExpr(lv, p);
    expr t2 = mulPolyExpr(ls, jx);

    p = addPolyExpr(t1, t2);

    expr t4 = mulPolyExpr(lv, s);
    expr t5 = mulPolyExpr(ls, v);
    expr t6 = mulPolyExpr(t5, jx);

    s = subPolyExpr(t4, t6);

    t = t + 1;

    m = degreePolyExpr(s);
  }

  expr k = powPolyExpr(lv, d - t);

  expr Q = mulPolyExpr(k, p);
  expr R = mulPolyExpr(k, s);

  return list({Q, R});
}

expr pseudoDivPolyExpr(expr &u, expr &v, expr &L) {
  return pseudoDivPolyExpr(std::forward<expr>(u), std::forward<expr>(v), L);
}

expr pseudoQuoPolyExpr(expr &u, expr &v, expr &L) {
  return pseudoDivPolyExpr(std::forward<expr>(u), std::forward<expr>(v), L)[0];
}

expr pseudoRemPolyExpr(expr &u, expr &v, expr &L) {
  return pseudoDivPolyExpr(std::forward<expr>(u), std::forward<expr>(v), L)[1];
}

expr getColPolyNormFactor(expr &u, expr &L, expr &K) {
  if (isZeroPolyExpr(u)) {
    return expr(u);
  }

  if (u.kind() == kind::INT || u.kind() == kind::FRAC) {
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

  expr lc = leadCoeffPolyExpr(u);

  expr rL = rest(L);

  return expr(kind::ADD, {getColPolyNormFactor(lc, rL, K) * pow(L[0], 0)});
}

expr normalizePolyExpr(expr &u, expr &L, expr &K) {
  if (isZeroPolyExpr(u))
    return expr(u);
  expr k = getColPolyNormFactor(u, L, K);
  return quoPolyExpr(u, k, L, K);
}

expr unitNormalColPoly(expr v, expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q");

  if (v.kind() == kind::ADD) {
    return 1;
  }

  short sign = 1;

  if (v.kind() == kind::INT) {
    sign = v.value() < 0 ? -1 : 1;
  }
  if (v.kind() == kind::FRAC) {
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
      if (v.kind() == kind::FRAC) {
        return fraction(-1 * v[1].value(), v[0].value());
      }

      return reduce(-1 / v);
    }

    if (v.kind() == kind::FRAC) {
      return fraction(v[1].value(), v[0].value());
    }

    return reduce(1 / v);
  }

  return 1;
}

expr contPolyExpr(expr &&u, expr &L, expr &K) {
  if (isZeroPolyExpr(u)) {
    return polyExpr(0, rest(L));
  }

  expr R = rest(L);

  long i = u.size() - 1;

  expr g = u[i][0];

  if (i == 0) {
    expr t = unitNormalColPoly(g, K);
    return mulPolyExpr(t, g);
  } else {
    i = i - 1;

    while (i >= 0) {
      expr ui = u[i][0];
      g = gcdPolyExpr(g, ui, R, K);
      i = i - 1;
    }
  }

  return g;
}

expr contPolyExpr(expr &u, expr &L, expr &K) {
  return contPolyExpr(std::forward<expr>(u), L, K);
}

expr ppPolyExpr(expr &&u, expr &L, expr &K) {
  if (L.size() == 0)
    return 1;

  expr c = raisePolyExpr(contPolyExpr(u, L, K), 0, L[0]);
  return quoPolyExpr(u, c, L, K);
}

expr ppPolyExpr(expr &u, expr &L, expr &K) {
  return ppPolyExpr(std::forward<expr>(u), L, K);
}

expr contAndPpPolyExpr(expr &&u, expr &L, expr &K) {
  expr c = contPolyExpr(u, L, K);
  expr C = raisePolyExpr(c, 0, L[0]);
  return list({c, quoPolyExpr(u, C, L, K)});
}

expr contAndPpPolyExpr(expr &u, expr &L, expr &K) {
  return contAndPpPolyExpr(std::forward<expr>(u), L, K);
}

expr expandPolyExpr(expr &&u) {
  if (is(&u, kind::TERMINAL) || u.kind() == kind::FUNC) {
    return expr(u);
  }

  if (u.kind() == kind::POW) {
    if (u[1] == 0)
      return 1;
    if (u[1] == 1)
      return u[0];
    return expr(u);
  }

  if (u.kind() == kind::MUL) {
    expr g = 1;

    for (Int i = 0; i < u.size(); i++) {
      expr c = expandPolyExpr(u[i]);

      if (c == 0) {
        return 0;
      }

      if (c == 1) {
        continue;
      }

      if (g == 1) {
        g = c;
        continue;
      }

      if (c.kind() == kind::POW || c.kind() == kind::INT ||
          c.kind() == kind::SYM || c.kind() == kind::FRAC ||
          c.kind() == kind::FUNC) {
        if (g.kind() == kind::ADD) {
          expr t = expr(kind::ADD);

          for (Int j = 0; j < g.size(); j++) {
            t = t + g[j] * c;
          }

          g = t;
        } else {
          g = g * c;
        }

        continue;
      }

      // Multiplication * Multiplication
      if (c.kind() == kind::MUL && g.kind() == kind::MUL) {
        for (Int j = 0; j < c.size(); j++) {
          g.insert(c[j]);
        }

        continue;
      }

      // Addition * Addition
      if (c.kind() == kind::ADD && g.kind() == kind::ADD) {
        expr t = expr(kind::ADD);

        for (Int j = 0; j < g.size(); j++) {
          for (Int k = 0; k < c.size(); k++) {
            t = t + g[j] * c[k];
          }
        }

        g = t;

        continue;
      }

      // Addition * Multiplication
      expr a = c.kind() == kind::ADD ? c : g;
      expr b = c.kind() == kind::ADD ? g : c;

      assert(a.kind() == kind::ADD);
      assert(b.kind() == kind::MUL);

      expr t = expr(kind::ADD);

      for (Int j = 0; j < a.size(); j++) {
        t = t + b * a[j];
      }

      g = t;
    }

    if (g.kind() == kind::ADD && g.size() == 1)
      return g[0];
    if (g.kind() == kind::MUL && g.size() == 1)
      return g[0];
    return g;
  }

  assert(u.kind() == kind::ADD);

  expr g = expr(kind::ADD);

  for (Int i = 0; i < u.size(); i++) {
    expr c = expandPolyExpr(u[i]);

    if (c != 0) {
      if (c.kind() == kind::ADD) {
        for (Int j = 0; j < c.size(); j++) {
          g = g + c[j];
        }
      } else {
        g = g + c;
      }
    }
  }

  if (g.size() == 0)
    return 0;
  if (g.size() == 1)
    return g[0];

  return g;
}

expr expandPolyExpr(expr &u) { return expandPolyExpr(std::forward<expr>(u)); }

expr diffPolyExprRec(expr &u, expr &x, bool *was_diff = nullptr) {
  if (u.kind() == kind::INT || u.kind() == kind::FRAC) {
    return u;
  }

  if (u.kind() == kind::MUL) {
    assert(u.size() == 2);

    if (u[1][0] == x) {
      expr d = u[1][1];

      assert(d.kind() == kind::INT);

      if (d == 0)
        return raisePolyExpr(0, 0, x);

      expr c = mulPolyExpr(d, u[0]);

      *was_diff = true;

      return c * pow(x, d.value() - 1);
    }

    bool cond = false;

    expr udx = diffPolyExprRec(u[0], x, &cond);

    if (cond == false)
      return raisePolyExpr(0, 0, x);

    *was_diff = true;

    return udx * u[1];
  }

  expr g = expr(kind::ADD);

  for (Int i = 0; i < u.size(); i++) {
    expr c = diffPolyExprRec(u[i], x, was_diff);

    if (isZeroPolyExpr(c))
      continue;

    g.insert(c);
  }

  if (g.size() == 0) {
    return raiseToExpression(0, u);
  }

  return g;
}

expr diffPolyExpr(expr &u, expr &x) {
  bool was_diff = false;

  return diffPolyExprRec(u, x, &was_diff);
}

Int normPolyExpr(expr u, expr L, expr K, size_t i = 0) {
  if (i == L.size()) {
    assert(u.kind() == kind::INT);

    return u.value();
  }

  Int k = 0;

  for (Int j = 0; j < u.size(); j++) {
    k = max(abs(normPolyExpr(u[j][0], L, K, i + 1)), k);
  }

  return k;
}

expr replaceAndReduce(expr u, expr x, expr c) {
  expr g = replace(u, x, c);
  return reduce(g);
}

expr evalPolyExpr(expr u, expr x, Int c) {
  if (is(&u, kind::CONST)) {
    return u;
  }

  expr g = expr(kind::ADD);

  for (Int i = 0; i < u.size(); i++) {
    if (base(u[i][1]) == x) {
      expr e = degree(u[i][1]);

      expr k = pow(c, e.value());

      expr t = mulPolyExpr(u[i][0], k);

      g = addPolyExpr(g, t);
    } else {
      g = g + (evalPolyExpr(u[i][0], x, c)) * u[i][1];
    }
  }

  if (!is(&g, kind::TERMINAL) && g.size() == 0) {
    return 0;
  }

  if (g.size() == 1 && is(&g[0], kind::TERMINAL | kind::FRAC)) {
    return g[0];
  }

  return g;
}

expr evalPolyExpr(expr u, expr x, expr c) {
  if (is(&u, kind::CONST)) {
    return u;
  }

  expr g = expr(kind::ADD);

  for (Int i = 0; i < u.size(); i++) {
    if (base(u[i][1]) == x) {
      expr e = degree(u[i][1]);
      expr k = powPolyExpr(c, e.value());
      expr t = mulPolyExpr(u[i][0], k);

      g = addPolyExpr(g, t);
    } else {
      g = g + evalPolyExpr(u[i][0], x, c) * u[i][1];
    }
  }

  if (is(&g, kind::ADD) && g.size() == 0) {
    return 0;
  }

  if (is(&g, kind::ADD) && g.size() == 1 && is(&g[0], kind::TERMINAL)) {
    return g[0];
  }

  return g;
}

expr evalTailPolyExpr(expr u, expr L, expr A, Int from) {
  assert(L.kind() == kind::LIST);
  assert(A.kind() == kind::LIST);

  expr t = u;

  for (Int i = from; i < L.size(); i++) {
    t = evalPolyExpr(t, L[i], A[i - from]);
  }

  return t;
}

expr invertPolyExpr(expr f) {
  expr i = -1;
  return mulPolyExpr(f, i);
}

expr interpolatePolyExpr(expr h, Int p, expr x, expr R, expr K) {

  expr f = expr(kind::ADD);

  expr y = polyExpr(p, R);

  Int i = 0;

  while (!isZeroPolyExpr(h)) {
    expr g = gfPolyExpr(h, p, true);

    if (!isZeroPolyExpr(g))
      f.insert(g * pow(x, i));

    h = subPolyExpr(h, g);
    h = quoPolyExpr(h, y, R, K);
    i = i + 1;
  }

  if (f.size() == 0) {
    return raisePolyExpr(polyExpr(0, R), 0, x);
  }

  expr lc = groundLeadCoeffPolyExpr(f);

  if (lc.value() < 0) {
    return invertPolyExpr(f);
  }

  return f;
}

expr groundContPolyExprRec(expr f) {
  if (f.kind() == kind::INT) {
    return f.value();
  }

  if (f.kind() == kind::MUL) {
    assert(f.size() == 2);
    return groundContPolyExprRec(f[0]);
  }

  expr g = 0;

  expr r, u, e, t, x, R;

  for (Int i = 0; i < f.size(); i++) {
    expr l = leadCoeffPolyExpr(f[i][0]);
    expr t = groundContPolyExprRec(l);

    g = gcd(g, t);
  }

  return g;
}

expr groundContPolyExpr(expr f) { return groundContPolyExprRec(f); }

expr groundDivPolyExpr(expr u, expr v) {
  assert(v.kind() == kind::INT);

  if (isZeroPolyExpr(u))
    return u;

  if (u.kind() == kind::INT) {
    return u.value() / v.value();
  }

  expr g = expr(kind::ADD);

  for (Int i = 0; i < u.size(); i++) {
    expr k = groundDivPolyExpr(u[i][0], v);

    if (!isZeroPolyExpr(k)) {
      g.insert(k * u[i][1]);
    }
  }

  if (g.size() == 0) {
    return raiseToExpression(0, u);
  }

  return g;
}

expr groundMulPolyExpr(expr u, Int v) {
  if (isZeroPolyExpr(u))
    return u;

  if (u.kind() == kind::INT) {
    return u.value() * v;
  }

  expr g = expr(kind::ADD);

  for (Int i = 0; i < u.size(); i++) {
    expr k = groundMulPolyExpr(u[i][0], v);

    if (!isZeroPolyExpr(k)) {
      g.insert(k * u[i][1]);
    }
  }

  if (g.size() == 0) {
    return raiseToExpression(0, u);
  }

  return g;
}

expr groundPPPolyExpr(expr u) {
  expr c = groundContPolyExpr(u);
  return groundDivPolyExpr(u, c);
}

expr heuristicGcdPolyExpr(expr u, expr v, expr L, expr K) {
  // References:
  // [1] Liao, H.-C., & Fateman, R. J.(1995). Evaluation of the heuristic
  // polynomial GCD

  assert(K == expr("Z"));

  if ((isConstantPolyExpr(u) && isConstantPolyExpr(v)) || L.size() == 0) {
    Int a = groundLeadCoeffPolyExpr(u).value();
    Int b = groundLeadCoeffPolyExpr(v).value();

    Int c = abs(gcd(a, b));

    if (L.size()) {
      return list({polyExpr(c, L), polyExpr(a / c, L), polyExpr(b / c, L)});
    }

    return list({c, a / c, b / c});
  }

  expr ucont = groundContPolyExpr(u);
  expr vcont = groundContPolyExpr(v);

  expr g = polyExpr(gcd(ucont, vcont), L);

  u = quoPolyExpr(u, g, L, K);
  v = quoPolyExpr(v, g, L, K);

  Int un = normPolyExpr(u, L, K);
  Int vn = normPolyExpr(v, L, K);

  Int b = 2 * min(un, vn) + 29;

  Int uc = groundLeadCoeffPolyExpr(u).value();
  Int vc = groundLeadCoeffPolyExpr(v).value();

  Int x = max(min(b, 99 * isqrt(b)), 2 * min(un / uc, vn / vc) + 2);

  for (short i = 0; i < 6; i++) {
    expr ux = evalPolyExpr(u, L[0], x);
    expr vx = evalPolyExpr(v, L[0], x);

    expr R = rest(L);

    if (!isZeroPolyExpr(ux) && !isZeroPolyExpr(vx)) {

      expr HGCD = heuristicGcdPolyExpr(ux, vx, R, K);
      expr h = HGCD[0];
      expr a = HGCD[1];
      expr b = HGCD[2];
      expr cu, cv, ru, rv;

      h = interpolatePolyExpr(h, x, L[0], R, K);

      h = groundPPPolyExpr(h);

      expr U = divPolyExpr(u, h, L, K);

      cu = U[0];
      ru = U[1];

      if (isZeroPolyExpr(ru)) {
        expr V = divPolyExpr(v, h, L, K);

        cv = V[0];
        rv = V[1];

        if (isZeroPolyExpr(rv)) {
          return list({h, cu, cv});
        }
      }

      a = interpolatePolyExpr(a, x, L[0], R, K);
      // a = raisePolyExpr(a, 0, L[0]);

      U = divPolyExpr(u, a, L, K);
      h = U[0];
      ru = U[1];

      if (isZeroPolyExpr(ru)) {
        expr V = divPolyExpr(v, h, L, K);

        cv = V[0];
        rv = V[1];

        if (isZeroPolyExpr(rv)) {
          h = mulPolyExpr(h, g);

          return list({h, a, cv});
        }
      }

      b = interpolatePolyExpr(b, x, L[0], R, K);
      // b = raisePolyExpr(b, 0, L[0]);

      expr V = divPolyExpr(v, b, L, K);

      h = V[0];
      rv = V[1];

      if (isZeroPolyExpr(rv)) {
        U = divPolyExpr(u, h, L, K);
        expr c = U[0];
        ru = U[1];

        if (isZeroPolyExpr(ru)) {
          h = mulPolyExpr(h, g);
          return list({h, c, b});
        }
      }
    }
    x = 73794 * x * isqrt(isqrt(x)) / 27011;
  }

  return fail();
}

expr groundInvertPolyExpr(expr p) {
  expr k = -1;
  return mulPolyExpr(k, p);
}

expr insertSymbolPolyExpr(expr u, expr x, Int d, Int level, Int i) {
  if (is(&u, kind::TERMINAL)) {
    return expr(kind::ADD, {u * pow(x, d)});
  }

  assert(u.kind() == kind::ADD);

  if (i == level) {

    expr g = expr(kind::ADD);

    for (Int j = 0; j < u.size(); j++) {
      g.insert(expr(kind::ADD, {u[j][0] * pow(x, d)}) * u[j][1]);
    }

    return g;
  }

  expr g = expr(kind::ADD);

  for (Int j = 0; j < u.size(); j++) {
    g.insert(insertSymbolPolyExpr(u[j][0], x, d, level, i + 1) * u[j][1]);
  }

  return g;
}

expr insertSymbolPolyExpr(expr u, expr x, Int d, Int level) {
  return insertSymbolPolyExpr(u, x, d, level, 0);
}

expr insertSymbolsPolyExpr(expr u, expr L, Int d, Int level) {
  expr g = u;

  for (Int i = 0; i < L.size(); i++) {
    g = insertSymbolPolyExpr(g, L[i], d, level + i);
  }

  return g;
}

expr degreePolyExpr(expr u, expr x) {
  if (is(&u, kind::TERMINAL))
    return -inf();
  if (u.size() == 0)
    return -inf();

  if (base(u[u.size() - 1][1]) == x) {
    return degree(u[u.size() - 1][1]);
  }

  expr m = -inf();

  for (Int i = 0; i < u.size(); i++) {
    expr d = degreePolyExpr(u[i][0], x);

    assert(d.kind() == kind::INT);

    if (m == -inf()) {
      m = d;
    } else {
      m = max(m.value(), d.value());
    }
  }

  return m;
}

bool isPolynomial(expr &a) {
	//TODO: return an error message detailing the reason, for instance, "Infinity is not a polynomial"
	// or "expressions with a symbolic degree are not a valid polynomial";

  if (is(&a, kind::DIV | kind::FACT | kind::UNDEF | kind::FAIL | kind::FUNC |
                 kind::INF | kind::LIST | kind::SET | kind::ROOT)) {
    return false;
  }

	if(is(&a, kind::FRAC | kind::INT | kind::SYM)) {
		return true;
	}

	assert(is(&a, kind::ADD | kind::SUB | kind::MUL | kind::POW));

	if(is(&a, kind::POW)) {
		assert(is(&a[1], kind::INT | kind::FRAC));
	}

	for(size_t i = 0; i < a.size(); i++) {
		if(!isPolynomial(a[i])) {
			return false;
		}
	}

	return true;
}


expr factorPolyExprAndExpand(expr u, expr L, expr K) {
 	assert(K.identifier() == "Z" || K.identifier() == "Q");

	expr c = 1;
	expr f = u;

	if (K.identifier() == "Q") {

		expr T = removeDenominatorsPolyExpr(u, L, K);

		c = T[0];
		f = T[1];
	}

	expr Z = expr("Z");

	expr F = factorization::factorsPolyExpr(f, L, Z);

	c = c * F[0];

	expr v = 1;

	for(size_t i = 0; i < F[1].size(); i++) {
		expr t = contAndPpPolyExpr(F[1][i][0], L, Z);

		c = c * t[0];
		v = v * pow(t[1], F[1][i][1]);
	}

	c = reduce(c);
	v = reduce(v);
	return c == 1 ? v : c * v;
}

expr lcmPolyExpr(expr u, expr v, expr L, expr K) {
	expr a = gcdPolyExpr(u, v, L, K);
	expr b = quoPolyExpr(u, a, L, K);

	return mulPolyExpr(v, b);
}

} // namespace polynomial
