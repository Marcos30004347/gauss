#include "Polynomial.hpp"
#include "MathSystem/Algebra/Expression.hpp"
#include "Resultant.hpp"
#include "MathSystem/GaloisField/GaloisField.hpp"

#include <climits>
#include <cstddef>
#include <map>
#include <numeric>
#include <utility>
#include <vector>

using namespace alg;
using namespace galoisField;

namespace polynomial {

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

		if (u.kind() == kind::MUL && u.size() == 2 &&
      u[1].kind() == kind::POW && u[1][0] == L[i]) {
    return expr(kind::ADD, {u});
  }

  expr c = 1;

  if (u.kind() == kind::POW && u[0] == L[i]) {
    return expr(kind::ADD, {collectRec(c, L, i + 1) * u});
  }

  if (u.kind() == kind::SYM || u.kind() == kind::FUNC) {
    if (u == L[i]) {
      return create(kind::ADD, { collectRec(c, L, i + 1) * pow(L[i], 1) });
    }

    return create(kind::ADD, {collectRec(u, L, i + 1) * pow(L[i], 0)});
  }

  Int d = collectDegree(u, L[i]);

  if (u.kind() == kind::MUL) {
    Int k = collectDegree(u, L[i]);
    expr c = collectCoeff(u, L[i], k);

    return expr(kind::ADD,
                {collectRec(c, L, i + 1) * pow(L[i], k)});
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

void includeVariable(std::vector<expr> &vars, expr u) {
  bool included = false;

  for (expr k : vars) {
    if (k == (u)) {
      included = true;
      break;
    }
  }

  if (!included) {
    vars.push_back(u);
  }
}

// bool isGeneralMonomial(expr u, expr v) {
//   expr S = set({});

// 	if (v.kind() != kind::SET) {
//     S = set({v});
//   } else {
//     S = v;
//   }

//   if (exists(S, u)) {

//     return true;
//   } else if (u.kind() == kind::POW) {
//     expr b = u[0];
//     expr e = u[1];

//     if (exists(S, b) && e.kind() == kind::INT && e.value() > 1) {

//       return true;
//     }
//   } else if (u.kind() == kind::MUL) {
//     for (unsigned int i = 0; i < u.size(); i++) {
//       if (isGeneralMonomial(u[i], S) == false) {

//         return false;
//       }
//     }

//     return true;
//   }

// 	return u.freeOfElementsInSet(S);
// }

// bool isGerenalPolynomial(expr u, expr v) {
//   expr S;

//   if (v.kind() != kind::SET) {
//     S = set({v});
//   } else {
//     S = v;
//   }

//   if (u.kind() != kind::ADD && u.kind() != kind::SUB) {
//     bool r = isGeneralMonomial(u, S);

//     return r;
//   }

//   if (exists(S, u)) {

//     return true;
//   }

//   for (unsigned int i = 0; i < u.size(); i++) {
//     if (isGeneralMonomial(u[i], S) == false) {

//       return false;
//     }
//   }

//   return true;
// }

list coeffVarMonomial(expr u, set S) {
  // if (!isGeneralMonomial(u, S))
  //   return undefined();

  if (is(&u, kind::CONST))
    return list({u, 1});

  if (exists(S, u))
    return list({1, u});

  if (u.kind() == kind::POW && exists(S, u[0]))
    return list({1, u});

  if (u.kind() == kind::MUL) {
    expr C = list({});
    expr V = list({});

    for (unsigned int i = 0; i < u.size(); i++) {
      list L = coeffVarMonomial(u[i], S);

		  C = join(C, list({L[0]}));
			V = join(V, list({L[1]}));
    }

    expr coef = create(kind::MUL);
    expr vars = create(kind::MUL);

    for (unsigned int i = 0; i < C.size(); i++) {
      if (C[i].kind() == kind::INT && C[i].value() == 1)
        continue;

      coef.insert(C[i]);
    }

    for (unsigned int i = 0; i < V.size(); i++) {
      if (V[i].kind() == kind::INT && V[i].value() == 1)
        continue;
      vars.insert(V[i]);
    }

    if (coef.size() == 0) {
      coef = 1;
    } else if (coef.size() == 1) {
      coef = coef[0];
    }

    if (vars.size() == 0) {
      vars = 1;
    } else if (vars.size() == 1) {
      vars = vars[0];
    }

    return list({coef, vars});
  }

  return list({u, 1});
}

expr collectTerms(expr u, set S) {
  if (u.kind() != kind::ADD) {

    list L = coeffVarMonomial(u, S);
    if (L.size() == 0) {
      return undefined();
    }

    return u;
  }

  if (exists(S, u)) {
    return u;
  }

  int N = 0;

  list T({});

  for (size_t i = 0; i < u.size(); i++) {
    list f = coeffVarMonomial(u[i], S);

    if (f.size() == 0) {
      return undefined();
    }

    int j = 1;

		bool combined = false;

    while (!combined && j <= N) {
      int k = j - 1;

      if (f[1] == (T[k][1])) {
        T[k] = list({create(kind::ADD, { T[k][0], f[0] }), f[1]});

				combined = true;
      }

      j = j + 1;
    }


    if (!combined) {
      T.insert(f, N++);
    }
  }

  expr v = create(kind::ADD, {});

  for (int j = 0; j < N; j++) {
    if (T[j][1].kind() == kind::INT && T[j][1].value() == 1) {
      v.insert(T[j][0]);
    } else {
      v.insert(create(kind::MUL, {
          T[j][0],
          T[j][1],
      }));
    }
  }

  if (v.size() == 0) {

    return 0;
  }

  if (v.size() == 1) {
    expr v_ = v[0];

    v = v_;
  }

  return v;
}

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

set variablesRec(expr u) {
  if (u.kind() == kind::INT || u.kind() == kind::FRAC) {
		return set({});
	}

  if (u.kind() == kind::POW) {
    expr b = u[0];
    expr e = u[1];

    if (e.kind() == kind::INT && e.value() > 1) {
      return set({ b });
		}

    return set({ u });
  }

  if (u.kind() == kind::ADD || u.kind() == kind::SUB ||
      u.kind() == kind::MUL) {
    set S = set({});

    for (unsigned int i = 0; i < u.size(); i++) {
			set T = variablesRec(u[i]);

			S = unification(S, T);
    }

    return S;
  }
	return set({u});
}

set variables(expr u) {
  return variablesRec(u);

}


list coefficientGME(expr u, expr x) {
	if (u == x) {
		// printf("aa\n");
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
			// printf("----> %s  %s\n",  to_string(u).c_str(), to_string(x).c_str());
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

  expr c = 0;
  for (size_t i = 0; i < u.size(); i++) {

		// printf("coeffGME %s %s\n", to_string(u[i]).c_str(), to_string(x).c_str());
		list f = coefficientGME(u[i], x);
		// printf("f = %s\n",  to_string(f).c_str());
    if (f.size() == 0) {
			// printf("----> %s  %s\n",  to_string(u[i]).c_str(), to_string(x).c_str());
			return undefined();
		}

    if (d == f[1]) {
      expr k = f[0];

      if (c == 0) {

        c = expr(u.kind());
        c.insert(k);
      } else {
        c.insert(k);
      }
    }
  }

  if (c.kind() != kind::INT && c.size() == 1) {
    return c[0];
  }

  return c;
}

expr leadCoeff(expr u, expr x) {
  return coeff(u, x, degree(u, x));
}

expr divideGPE(expr u, expr v, expr x) {
  expr t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

  expr q = 0;
  expr r = u;

  expr m = degree(r, x);
  expr n = degree(v, x);

  // printf("-> lead %s\n", to_string(v).c_str());
  expr lcv = leadCoeff(v, x);
  // printf("<- lead %s\n", to_string(lcv).c_str());

  while (m != -inf() && (m.kind() == kind::INT && n.kind() == kind::INT &&
                         m.value() >= n.value())) {
    // printf("-> lead %s\n", to_string(r).c_str());
    expr lcr = leadCoeff(r, x);
    // printf("<- lead %s\n", to_string(lcr).c_str());

    expr s = lcr / lcv;
    // printf("a -------\n");
    t1 = pow(x, m - n);

    t2 = mulPoly(s, t1);
    q = addPoly(q, t2);

    // printf("b -------\n");
    t1 = pow(x, m);
    t2 = mulPoly(lcr, t1);
    // printf("(%s) * (%s) = %s\n", to_string(lcr).c_str(), to_string(t1).c_str(),
    //        to_string(t2).c_str());

    t3 = subPoly(r, t2);

    // printf("(%s) - (%s) = %s\n", to_string(r).c_str(), to_string(t2).c_str(),
    //        to_string(t3).c_str());

    // printf("c -------\n");

    t4 = pow(x, n);
    t5 = mulPoly(lcv, t4);

    // printf("(%s) * (%s) = %s\n", to_string(lcv).c_str(), to_string(t4).c_str(),
    //        to_string(t5).c_str());

    t6 = subPoly(v, t5);

    // printf("(%s) - (%s) = %s\n", to_string(v).c_str(), to_string(t5).c_str(),
    //        to_string(t6).c_str());

    // printf("d -------\n");
    t7 = mulPoly(t6, s);

    // printf("(%s) * (%s) = %s\n", to_string(t6).c_str(), to_string(s).c_str(),
    //        to_string(t7).c_str());

    t8 = pow(x, m - n);
    t9 = mulPoly(t7, t8);

    // printf("(%s) * (%s) = %s\n", to_string(t7).c_str(), to_string(t8).c_str(),
    //        to_string(t9).c_str());

    // printf("********** e -------\n");
    // printf("********** e -------\n");
    // printf("********** e -------\n");
    // printf("********** e -------\n");
    r = subPoly(t3, t9);
    // printf("(%s) - (%s) = %s\n", to_string(t3).c_str(), to_string(t9).c_str(),
    //        to_string(r).c_str());

    // printf("f -------\n");
    m = degree(r, x);
    // printf("deg(%s) = %s\n", to_string(r).c_str(), to_string(m).c_str());

    // printf("g -------\n");
  }

  // printf("aa\n");
  return list({(q), (r)});
}

expr quotientGPE(expr u, expr v, expr x) { return divideGPE(u, v, x)[0]; }

expr remainderGPE(expr u, expr v, expr x) {
	return divideGPE(u, v, x)[1];
}

expr expandGPE(expr u, expr v, expr x, expr t) {
  if (u == 0)
    return 0;

  expr d = divideGPE(u, v, x);

  expr q = d[0];
  expr r = d[1];

  expr expoent = ((t * expandGPE(q, v, x, t)) + r);

  return expand(expoent);
}

expr gcdGPE(expr u, expr v, expr x) {
  if (u == 0 && v == 0) {
    return 0;
  }

  expr U = u;
  expr V = v;

	while (V != 0) {
    expr R = remainderGPE(U, V, x);

    U = V;
    V = R;
  }

  return expand((1 / leadCoeff(U, x)) * U);
}

expr extendedEuclideanAlgGPE(expr u, expr v, expr x) {
  if (u == 0 && v == 0) {
    return list({0, 0, 0});
  }

  expr U = u;
  expr V = v;

  expr App = 1, Ap = 0, Bpp = 0, Bp = 1;

  while (V != 0) {
    expr d = divideGPE(U, V, x);

    expr q = d[0];
    expr r = d[1];

    expr A_ = (App - (q * Ap));
    expr B_ = (Bpp - (q * Bp));

    expr A = expand(A_);
    expr B = expand(B_);

    App = Ap;

    Ap = A;

    Bpp = Bp;

    Bp = B;

    U = V;

    V = r;
  }

  expr c = leadCoeff(U, x);

  expr App_ = quotientGPE(App, c, x);

  App = App_;

  expr Bpp_ = quotientGPE(Bpp, c, x);

  Bpp = Bpp_;

  expr U_ = quotientGPE(U, c, x);

  U = U_;

  return list({U, App, Bpp});
}

expr mulPolyRec(expr p1, expr p2) {
  if (p1.kind() == kind::ADD) {
    expr res = create(kind::ADD, {});

    for (unsigned int i = 0; i < p1.size(); i++) {
      res.insert(mulPoly(p1[i], p2));
    }

    return res;
  }

  if (p2.kind() == kind::ADD) {
    expr res = create(kind::ADD);

    for (unsigned int i = 0; i < p2.size(); i++) {
      res.insert(mulPoly(p2[i], p1));
    }

    return res;
  }

  return (p1 * p2);
}

expr mulPoly(expr p1, expr p2) {
  expr t1 = mulPolyRec(p1, p2);
  expr t2 = reduce(t1);

  return t2;
}

expr subPoly(expr p1, expr p2) {
  expr s = (p1 - p2);
  expr p = reduce(s);

  return p;
}

expr addPoly(expr p1, expr p2) {
  expr s = (p1 + p2);
  expr p = reduce(s);

  return p;
}

expr recPolyDiv(expr u, expr v, expr L, expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q");

  if (L.size() == 0) {

		expr d = expand(u / v);

		//printf("(%s) / (%s) =  %s\n", to_string(u).c_str(), to_string(v).c_str(), to_string(d).c_str());
		if (K.identifier() == "Z") {
      if (d.kind() == kind::INT) {
        return list({d, 0});
      }

      return list({0, u});
    }

    return list({d, 0});
  }

  expr x = L[0];
  expr r = u;

  expr m = degree(r, x);
  expr n = degree(v, x);

  expr q = 0;

  expr lcv = leadCoeff(v, x);
  expr R = rest(L);

	while (m != -inf() && m.value() >= n.value()) {
    expr lcr = leadCoeff(r, x);
		//printf("lcr %s\n", to_string(lcr).c_str());
    expr d = recPolyDiv(lcr, lcv, R, K);
		//printf("d2 %s\n", to_string(d).c_str());

    if (d[1] != 0) {
      return list({expand(q), r});
    }

    expr j = pow(x, m - n);
		//printf("j %s\n", to_string(j).c_str());

    q = q + d[0] * j;

		//printf("q %s\n", to_string(q).c_str());
    expr t1 = mulPoly(v, d[0]);
    expr t2 = mulPoly(t1, j);
    expr t3 = subPoly(r, t2);

    r = reduce(t3);

		//printf("r %s\n", to_string(r).c_str());
    m = degree(r, x);
  }

  return list({expand(q), r});
}

expr recQuotient(expr u, expr v, expr L, expr K) {
  return recPolyDiv(u, v, L, K)[0];
}

expr recRemainder(expr u, expr v, expr L, expr K) {
  return recPolyDiv(u, v, L, K)[1];
}

expr pdiv(expr f, expr g, expr x) {
  assert(g != 0);

  expr lg, k, q, r, t, m, n, j;
  expr t1, t2, t3, t4, t5, t6;

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

		j = (t - n);
    k = (k - 1);

		t3 = pow(x, j);

    t2 = mulPoly(q, lg); // mul({ q, lg });
    t4 = mulPoly(t1, t3);
    q = addPoly(t2, t4);

    t4 = mulPoly(r, lg); // mul({ r, lg });
    t5 = mulPoly(g, t1); // mul({ g, t1, t3 });
    t6 = mulPoly(t5, t3);
    r = subPoly(t4, t6);

    t = degree(r, x);

    if (t == -inf() || t.value() < n.value()) {
      break;
    }
  }

  q = q * pow(lg, k);
  r = r * pow(lg, k);

  t1 = expand(q);
  t2 = expand(r);

  return list({t1, t2});
}

expr pseudoDivision(expr u, expr v, expr x) {
  expr p = 0;
  expr s = u;

  expr m = degree(s, x);
  expr n = degree(v, x);

  expr delta = max(m.value() - n.value() + 1, 0);

  expr lcv = leadCoeff(v, x);

  Int tal = 0;

  while (m != -inf() && m.value() >= n.value()) {
    expr lcs = leadCoeff(s, x);

    expr j = pow(x, (m - n));

    expr t1 = mulPoly(lcv, p);
    expr t2 = mulPoly(lcs, j);

    expr t3 = addPoly(t1, t2);

    p = reduce(t3);

    expr t4 = mulPoly(lcv, s);
    expr t5 = mulPoly(lcs, v);
    expr t6 = mulPoly(t5, j);

    s = subPoly(t4, t6);

    tal = tal + 1;

    m = degree(s, x);
  }

  expr k = pow(lcv, delta.value() - tal);

  expr A = mulPoly(k, p);
  expr B = mulPoly(k, s);

  expr Q = reduce(A);
  expr R = reduce(B);

  return list({Q, R});
}

expr pseudoQuotient(expr u, expr v, expr x) {
  expr r = pseudoDivision(u, v, x);
  expr q = r[0];

  return q;
}

expr pseudoRemainder(expr u, expr v, expr x) {
  expr r = pseudoDivision(u, v, x);
  expr q = r[1];

  return q;
}

expr getNormalizationFactor(expr u, expr L, expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q");

  if (u == 0) {
    return 0;
  }

  if (is(&u, kind::CONST)) {
		if(is(&u, kind::INT)) {
			if (u.value() > 0) {
				if (K.identifier() == "Z") {
					return 1;
				}
			}
		}

		if(is(&u, kind::FRAC)) {
			if(u[0].value() > 0) {
				if (K.identifier() == "Z") {
					return 1;
				}
			}
		}

		return pow(u, -1);

	}

	if (L.size() == 0) {
    return undefined();
  }

	expr lc = leadCoeff(u, L[0]);

  return getNormalizationFactor(lc, rest(L), K);
}

expr normalizePoly(expr u, expr L, expr K) {
	if (u == 0) {
    return 0;
  }

  return expand(getNormalizationFactor(u, L, K) * u);
}

expr unitNormal(expr v, expr K) {
  assert(K.identifier() == "Z" || K.identifier() == "Q");

  if (K.identifier() == "Z") {
		if(is(&v, kind::INT)) {
			if (v.value() < 0) {
				return -1;
			}
		}

		if(is(&v, kind::FRAC)) {
			if(v[0].value() < 0) {
				if (K.identifier() == "Z") {
					return -1;
				}
			}
		}

    return 1;
  }

  if (K.identifier() == "Q") {

		if(is(&v, kind::INT)) {
			if (v.value() < 0) {
				return reduce(-1 * pow(v, -1));
			}
		}

		if(is(&v, kind::FRAC)) {
			if(v[0].value() < 0) {
				if (K.identifier() == "Z") {
					return reduce(-1 * pow(v, -1));
				}
			}
		}

    return pow(v, -1);
  }

  return 1;
}

// Finds the content of u with respect to x using
// the auxiliary variables R with coeff domain K,
// with is Z or Q
expr polynomialContent(expr u, expr x, expr R, expr K) {
  if (u == (0)) {
    return 0;
  }

  expr n = degree(u, x);

  expr g = coeff(u, x, n);

  expr k = (u - (g * pow(x, n)));

  expr v = expand(k);

  if (v == (0)) {
    expr un = unitNormal(g, K);

		expr t = (un * g);

    g = reduce(t);

  } else {
    while (v != 0) {
      expr d = degree(v, x);
      expr c = leadCoeff(v, x);

      expr t = gcdPoly(g, c, R, K);

      g = t;

      k = (v - (c * pow(x, d)));
      v = expand(k);
    }
  }

  return g;
}

// Finds the content of u with respect to x using
// the auxiliary variables R with coeff domain K,
// with is Z or Q
expr polynomialContentSubResultant(expr u, expr x, expr R, expr K) {
  if (u == (0)) {
    return 0;
  }
	//printf("u =  %s\n", to_string(u).c_str());
  expr n = degree(u, x);

  expr g = coeff(u, x, n);

  expr v = expand(u - g * pow(x, n));

  if (v == 0) {
    g = reduce(unitNormal(g, K) * g);
  } else {
		while (v != 0) {
			//printf("v = %s\n", to_string(v).c_str());
      expr d = degree(v, x);
			//printf("d3 = %s\n", to_string(d).c_str());
      expr c = coeff(v, x, d);
			//printf("c = %s\n", to_string(c).c_str());

			g = gcdPoly(g, c, R, K);

			v = expand(v - c * pow(x, d));
			//printf("v = %s\n", to_string(v).c_str());
    }
  }
  //printf("cont = %s\n", to_string(g).c_str());
  return g;
}

expr subResultantGCDRec(expr u, expr v, expr L, expr K) {
  if (L.size() == 0) {
    if (K.identifier() == "Z") {
      return gcd(u.value(), v.value());
    }

    if (K.identifier() == "Q") {
      return 1;
    }
  }

  expr x = first(L);

  expr du = degree(u, x);
  expr dv = degree(v, x);

  expr U = undefined();
  expr V = undefined();

  if (du.value() >= dv.value()) {
    U = u;
    V = v;
  } else {
    U = v;
    V = u;
  }

  expr R = rest(L);

  expr contU = polynomialContentSubResultant(U, x, R, K);
  expr contV = polynomialContentSubResultant(V, x, R, K);

  expr d = subResultantGCDRec(contU, contV, R, K);

  expr tmp1 = recQuotient(U, contU, L, K);
  expr tmp2 = recQuotient(V, contV, L, K);

  U = tmp1;

  V = tmp2;

  expr tmp3 = leadCoeff(U, x);
  expr tmp4 = leadCoeff(V, x);

  expr g = subResultantGCDRec(tmp3, tmp4, R, K);

  int i = 1;

  expr delta = undefined();
  expr y = undefined();
  expr b = undefined();
  expr dp = undefined();

  while (V != 0) {
    expr r = pseudoRemainder(U, V, x);

    if (r != 0) {
      if (i == 1) {

        expr tmp3 = (degree(U, x) + -degree(V, x) + 1);

        delta = expand(tmp3);

        y = -1;

        expr tmp4 = pow(-1, delta);

        b = expand(tmp4);

      } else {
        dp = delta;

        expr tmp3 = (degree(U, x) + -degree(V, x) + 1);

        delta = expand(tmp3);

        expr f = leadCoeff(U, x);

        expr tmp4 = pow(-f, (dp - 1));
        expr tmp5 = pow(y, (dp - 2));

        expr tmp6 = expand(tmp4);
        expr tmp7 = expand(tmp5);

        y = recQuotient(tmp6, tmp7, R, K);

        expr tmp8 = -f * pow(y, delta - 1);

        b = expand(tmp8);
      }

      U = V;

      V = recQuotient(r, b, L, K);

      i = i + 1;
    } else {
      U = V;

      V = r;
    }
  }

  expr tmp5 = leadCoeff(U, x);

  expr s = recQuotient(tmp5, g, R, K);

  expr W = recQuotient(U, s, L, K);

  expr contW = polynomialContentSubResultant(W, x, R, K);
  expr ppW = recQuotient(W, contW, L, K);

  expr tmp6 = (d * ppW);
  expr res = expand(tmp6);

  return res;
}

expr mvSubResultantGCD(expr u, expr v, expr L, expr K) {
  if (u == (0)) {
    return normalizePoly(v, L, K);
  }

  if (v == (0)) {
    return normalizePoly(u, L, K);
  }

  expr gcd = subResultantGCDRec(u, v, L, K);

  expr r = normalizePoly(gcd, L, K);

  return r;
}

expr mvPolyGCDRec(expr u, expr v, expr L, expr K) {
  if (L.size() == 0) {
    if (K.identifier() == "Z") {
      return gcd(u.value(), v.value());
    }

    if (K.identifier() == "Q") {
      return 1;
    }
  }

	expr x = first(L);
  expr R = rest(L);
  expr cont_u = polynomialContent(u, x, R, K);
  expr cont_v = polynomialContent(v, x, R, K);

  expr d = mvPolyGCDRec(cont_u, cont_v, R, K);

  expr pp_u = recQuotient(u, cont_u, L, K);
  expr pp_v = recQuotient(v, cont_v, L, K);

  while (pp_v != 0) {
    expr r = pseudoRemainder(pp_u, pp_v, x);

    expr pp_r = undefined();

    if (r == (0)) {
      pp_r = 0;
    } else {
      expr cont_r = polynomialContent(r, x, R, K);
      pp_r = recQuotient(r, cont_r, L, K);
    }

    pp_u = pp_v;
    pp_v = pp_r;
  }

  expr k = (d * pp_u);
  expr result = expand(k);

  return result;
}

expr mvPolyGCD(expr u, expr v, expr L, expr K) {

	if (u == (0)) {
    return normalizePoly(v, L, K);
  }

  if (v == (0)) {
    return normalizePoly(u, L, K);
  }

  expr gcd = mvPolyGCDRec(u, v, L, K);
  expr r = normalizePoly(gcd, L, K);

  return r;
}

expr leadMonomial(expr u, expr L) {
	if (L.size() == 0) {
		return u;
  }

  expr x = first(L);
  expr m = degree(u, x);
  expr c = coeff(u, x, m);

  return expand(pow(x, m) * leadMonomial(c, rest(L)));
}

bool wasSimplified(expr u) {
  if (u.kind() == kind::SYM)
    return true;

  if (u.kind() == kind::INT)
    return true;

  if (u.kind() == kind::DIV)
    return false;

  if (u.kind() == kind::MUL) {

    for (unsigned int i = 0; i < u.size(); i++) {
      if (u[i].kind() == kind::FRAC) {
        return false;
      }

      if (u[i].kind() == kind::POW && u[i][1].kind() == kind::INT &&
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
expr G(expr u, expr v) {
  // //printf("A\n");
  if (u.kind() == kind::ADD || u.kind() == kind::SUB) {

		expr k = create(u.kind());

    for (unsigned int i = 0; i < u.size(); i++) {
      // //printf("B %s\n", to_string(u[i]).c_str());
      // //printf("B %s\n", to_string(v).c_str());

      expr z = expand(u[i] / v);

      // //printf("C\n");
      if (wasSimplified(z)) {
        // //printf("D\n");
        k.insert(z);
      }
    }

    if (k.size() == 0) {
      return 0;
    }

    // //printf("D\n");
    return k;
  }

  expr z = expand(u / v);

  if (wasSimplified(z)) {
    return z;
  }

  return 0;
}

expr monomialPolyDiv(expr u, expr v, expr L) {
  expr q = 0;

	expr t = leadMonomial(v, L);
  expr f = G(u, t);

	while (f != 0) {
    q = reduce(q + f);
    u = expand(u - (f * v));
    f = G(u, t);
  }

  return list({q, u});
}

// TODO

// monomialBasedPolyExpansion(a^2*b + 2*a*b^2 + b^3 + 2*a + 2*b + 3, a+b, [a,
// b], t) -> b*t^2 + 2*t + 3
expr monomialBasedPolyExpansion(expr u, expr v, expr L, expr t) {
  if (u.kind() == kind::INT && u.value() == 0) {
    return 0;
  }

  expr d = monomialPolyDiv(u, v, L);
  expr q = d[0];
  expr r = d[1];

  expr k = ((
                    t *
                    monomialBasedPolyExpansion(q, v, L, t)
                ) +
                r);

  expr x = expand(k);

  return x;
}

// monomialPolyRem can be used for simplification, for example
// monomialPolyRem(a*i^3 + b*i^2 + c*i + d, i^2 + 1, [i]) -> -a*i - b + c*i + d
// simplification when i^2 + 1 = 0
// also
// monomialPolyRem(sin^4(x)+sin^3(x)+2*sin^2(x)cos^2(x)+cos^4(x),
// sin^2(x)+cos^2(x)-1, [cos(x), sin(x)]) -> 1 + sin^3(x)
expr monomialPolyRem(expr u, expr v, expr L) {
  expr d = monomialPolyDiv(u, v, L);
  expr r = d[1];

  return r;
}

expr monomialPolyQuo(expr u, expr v, expr L) {
  expr d = monomialPolyDiv(u, v, L);
  expr r = d[0];

  return r;
}

expr algebraicExpandRec(expr u);

expr expandProduct(expr r, expr s) {
  if (r == 0 || s == 0)
    return 0;

  if (r.kind() == kind::ADD && r.size() == 0)
    return 0;
  if (s.kind() == kind::ADD && s.size() == 0)
    return 0;
  // if (r.kind() == kind::MUL) r = expand(r);

  if (r.kind() == kind::ADD) {
    expr f = r[0];

    expr k = r;

    k.remove(0);

    expr a = expandProduct(f, s);
    expr b = expandProduct(k, s);

    bool c0 = a == 0;
    bool c1 = b == 0;

    expr z = 0;

    if (!c0 && !c1)
      z = a + b;
    else if (!c0 && c1)
      z = a;
    else if (c0 && !c1)
      z = b;
    else
      z = 0;

    return z;
  }

  if (s.kind() == kind::ADD) {
    return expandProduct(s, r);
  }

  return reduce(r * s);
}

expr expandPow(expr u, expr n) {
  if (u == 1)
    return 1;
  if (n == 0)
    return u == 0 ? undefined() : 1;

  if (u.kind() == kind::ADD) {
    expr f = u[0];

    expr o = reduce(u - f);
    expr s = 0;

    Int N = n.value();

    for (Int k = 0; k <= N; k++) {
      Int d = fact(N) / (fact(k) * fact(N - k));

      expr z = d * pow(f, N - k);
      expr t = expandPow(o, k);

      s = s + expandProduct(z, t);
    }

    return s;
  }

  return reduce(pow(u, n));
}

expr expandProductRoot(expr r, expr s) {
  if (r.kind() == kind::ADD) {
    expr f = r[0];
    expr k = r;

    k.remove(0);

    return f * s + k * s;
  }

  if (s.kind() == kind::ADD) {
    return expandProductRoot(s, r);
  }

  return r * s;
}

expr expandPowerRoot(expr u, expr n) {
  if (u.kind() == kind::ADD) {
    expr f = u[0];

    expr r = reduce(u - f);

    expr s = 0;

    Int N = n.value();

    for (Int k = 0; k <= n.value(); k++) {
      expr c = fact(N) / (fact(k) * fact(N - k));
      expr z = reduce(c * pow(f, N - k));

      expr t = expandPowerRoot(r, k);

      s = s + expandProductRoot(z, t);
    }

    return s;
  }

  expr v = pow(u, n);

  return reduce(v);
}

expr algebraicExpandRoot(expr u) {
  if (is(&u, kind::TERMINAL))
    return reduce(u);

  expr u_ = reduce(u);

  if (u_.kind() == kind::ADD) {
    expr v = u_[0];
    expr k = reduce(u_ - v);
    u_ = algebraicExpandRoot(v) + algebraicExpandRoot(k);
  }

  if (u_.kind() == kind::MUL) {
    expr v = u_[0];
    expr t = reduce(u_ / v);
    u_ = expandProductRoot(t, v);
  }

  if (u_.kind() == kind::POW) {

    expr b = u_[0];
    expr e = u_[1];

    if (e.kind() == kind::INT && e.value() >= 2) {
      expr t = expandPowerRoot(b, e);
      u_ = reduce(t);
    }

    if (e.kind() == kind::INT && e.value() <= -2) {
      expr p = reduce(pow(u_, -1));
      expr t = expandPowerRoot(p[0], p[1]);
      u_ = pow(t, -1);
    }
  }

  expr k = reduce(u_);

  return k;
}


expr cont(expr u, expr x) {
  expr n, c, c1, c2, tmp;

  u = expand(u);

  if (u == 0) {
    return 0;
  }

  if (u.size() >= 2) {
    n = degree(u, x);

    c1 = coeff(u, x, n);

    u = expand(u - c1 * pow(x, n));

    n = degree(u, x);

    c2 = coeff(u, x, n);

    u = expand(u - c2 * pow(x, n));

    c = gcd(c1.value(), c2.value());

    while (u != 0) {
      n = degree(u, x);

      c1 = coeff(u, x, n);

      tmp = u - c1 * pow(x, n);

      u = expand(tmp);

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

expr cont(expr u, expr L, expr K) {
  return polynomialContentSubResultant(u, L[0], rest(L), K);
}

expr pp(expr u, expr L, expr K) {
  expr c = cont(u, L, K);

  expr p = recQuotient(u, c, L, K);

  return p;
}

expr pp(expr u, expr c, expr L, expr K) {
  expr R = rest(L);

  expr p = recQuotient(u, c, L, K);

  return p;
}

expr groundLeadCoeffPoly(expr u, expr L) {
	if(L.size() == 0) {
		assert(u.kind() == kind::INT || u.kind() == kind::FRAC);
		return u;
	}

	return groundLeadCoeffPoly(leadCoeff(u, L[0]), rest(L));
}

expr groundLeadCoeffPolyExpr(expr u) {
	if(u.kind() == kind::INT || u.kind() == kind::FRAC) {
		return u;
	}

	return groundLeadCoeffPolyExpr(leadCoeffPolyExpr(u));
}

expr monic(expr u, expr L, expr K) {
	expr lc = groundLeadCoeffPoly(u, L);
	return recQuotient(u, lc, L, K);
}

expr monicPolyExpr(expr u, expr L, expr K) {
	expr lc = groundLeadCoeffPolyExpr(u);
	return quoPolyExpr(u, lc, L, K);
}

expr igcdPoly(expr u, expr v, expr L, expr K) {
	if(L.size() == 0) {
		assert(u.kind() == kind::INT);
		assert(v.kind() == kind::INT);

		return abs(gcd(u.value(), v.value()));
	}
	expr heu = heuristicGcdPoly(u, v, L, K);

	if(heu != fail()) {
		return heu[0];
	}

	expr u_cnt = cont(u, L, K);
	expr u_ppr = pp(u, u_cnt, L, K);

	expr v_cnt = cont(v, L, K);
	expr v_ppr = pp(v, v_cnt, L, K);

	expr R = rest(L);

	expr h = polyRemSeq(u_ppr, v_ppr, L, K)[0];

	expr tmp = cont(h, L, K);

	h = pp(h, L, K);

	expr c = gcdPoly(u_cnt, v_cnt, R, K);

	h = mulPoly(h, c);

	return monic(h, L, K);
}

expr constDenLcmPolyExprRec(expr u, expr L, unsigned int j) {
	if(L.size() == j) {
		assert(u.kind() == kind::INT || u.kind() == kind::FRAC);

		if(u.kind() == kind::FRAC) {
			return u[1];
		}

		return 1;
	}

	expr g = 1;

	for(Int i = 0; i < u.size(); i++) {
		expr t = constDenLcmPolyExprRec(u[i][0], L, j + 1);
		g = lcm(g, t);
	}

	return g;
}

expr constDenLcmPolyExpr(expr u, expr L) {
	return constDenLcmPolyExprRec(u, L, 0);
}

expr constDenLcmPoly(expr u, expr L) {
	if(L.size() == 0) {
		if(u.kind() == kind::FRAC) {
			return u[1];
		}

		return 1;
	}

	expr g = 1;

	expr R = rest(L);

	expr d = degree(u, L[0]);

	for(Int i = 0; i <= d.value(); i++) {
		expr t = constDenLcmPoly(coeff(u, L[0], i), R);
		g = lcm(g, t);
	}

	return g;
}

expr removeDenominatorsPoly(expr u, expr L, expr K) {
	assert(K.identifier() == "Z" || K.identifier() == "Q");

	expr f = constDenLcmPoly(u, L);
	return list({ f , mulPoly(f, u) });
}

expr removeDenominatorsPolyExpr(expr u, expr L, expr K) {
	assert(K.identifier() == "Z" || K.identifier() == "Q");
	expr f = constDenLcmPolyExpr(u, L);
	expr t = polyExpr(f, L);
	expr k =  mulPolyExpr(u, t);
	return list({ t , k });
}

expr igcdPolyExpr(expr u, expr v, expr L, expr K) {
	if(L.size() == 0) {
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

	if(K.identifier() == "Q") {
		u = removeDenominatorsPolyExpr(u, L, K)[1];
		v = removeDenominatorsPolyExpr(v, L, K)[1];
	}

	// if(isZeroPolyExpr(u)) return v;
	// if(isZeroPolyExpr(v)) return u;

	expr Z = expr("Z");
	expr H = heuristicGcdPolyExpr(u, v, L, Z);

	if(H != fail()) {
		if(K == Z) {
			return H[0];
		}

		return monicPolyExpr(H[0], L, K);
	}
	expr g = igcdPolyExpr(u, v, L, Z);

	if(K == Z) {
		return g;
	}

	return monicPolyExpr(g, L, K);
}

expr gcdPoly(expr u, expr v, expr L, expr K) {
	expr a = 1;
	expr b = 1;

	if(K.identifier() == "Q") {
		u = removeDenominatorsPoly(u, L, K)[1];
		v = removeDenominatorsPoly(v, L, K)[1];
	}

	expr Z = expr("Z");
	expr H = heuristicGcdPoly(u, v, L, Z);

	if(H != fail()) {
		if(K == Z) {
			return H[0];
		}

		return monic(H[0], L, K);
	}

	expr g = igcdPoly(u, v, L, Z);

	// //printf("--> %s\n",to_string(g).c_str() );

	if(K == Z) {
		return g;
	}
	return monic(g, L, K);
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
	if(v == 0) return isZeroPolyExpr(u);

	if (is(&u, kind::TERMINAL)) {
    return (u.kind() == kind::INT || u.kind() == kind::FRAC) && u == v;
  }

	if(u.kind() == kind::MUL) {
		if(u.size() > 2) return false;
		if(u.size() == 2) {
			assert(u[1].kind() == kind::POW);
			assert(u[1][1].kind() == kind::INT);

			if(u[1][1].value() != 0) {
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

  if(u.kind() == kind::MUL) {
		assert(u.size() == 2);

		expr p = degree(u[1]);

		if(p != 0 && p != -inf()) return false;
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

  for(std::map<Int, expr>::iterator it = coeffs.begin(); it != coeffs.end(); it++) {
		if(!isZeroPolyExpr(it->second) || g.size() == 0) {
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

	if(is(&u, kind::ADD) && u.size() == 0) {
		return v;
	}

	if(is(&v, kind::ADD) && v.size() == 0) {
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

expr addPolyExpr(expr &u, expr &v) {


  return addColPolyRec(u, v, 0, 0);
}

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
	if(is(&u, kind::CONST)) {
		return c;
	}

	if(u.kind() != kind::ADD || u.size() == 0) {
		return u;
	}

	expr p = u[0][0];
	expr k = u[0][1];

	return expr(kind::ADD, { raiseToExpression(c, p)*pow(base(k), 0) });
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

expr raiseToPolyExpr(expr &x, Int e, expr &L) {
  return polyExpr(pow(x, e), L);
}

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

// expr colPolyGCDRec(expr &u, expr &v, expr &L, expr &K) {
//   if (L.size() == 0) {
//     if (K.identifier() == "Z") {
//       return polyExpr(gcd(u.value(), v.value()), L);
//     }

//     if (K.identifier() == "Q") {
//       return polyExpr(1, L);
//     }
//   }

//   expr R = rest(L);

//   expr ct_u = contPolyExpr(u, L, K);
//   expr ct_v = contPolyExpr(v, L, K);

//   expr d = colPolyGCDRec(ct_u, ct_v, R, K);

//   ct_v = raisePolyExpr(ct_v, 0, L[0]);
//   ct_u = raisePolyExpr(ct_u, 0, L[0]);

//   expr pp_u = quoPolyExpr(u, ct_u, L, K);
//   expr pp_v = quoPolyExpr(v, ct_v, L, K);

//   expr pp_r = 0;
//   expr ct_r = 0;

//   expr r = 0;

//   while (!isZeroPolyExpr(pp_v)) {
//     r = pseudoRemPolyExpr(pp_u, pp_v, L);

//     if (isZeroPolyExpr(r)) {
//       pp_r = polyExpr(0, L);
//     } else {
//       ct_r = contPolyExpr(r, L, K);
//       ct_r = raisePolyExpr(ct_r, 0, L[0]);

//       pp_r = quoPolyExpr(r, ct_r, L, K);
//     }

//     pp_u = pp_v;
//     pp_v = pp_r;
//   }

//   d = raisePolyExpr(d, 0, L[0]);

//   return mulPolyExpr(d, pp_u);
// }

// bool isConstantPolyExpr(expr& u) {
// 	if(u.kind() == kind::INT) {
// 		return true;
// 	}

// 	if(u.kind() == kind::ADD && u.size() > 1) {
// 		return false;
// 	}

// 	if(u.kind() == kind::MUL) {
// 		assert(u.size() == 2, "u should have only two operands");
// 		return u[1][1] == 0 && isConstantPolyExpr(u[0]);
// 	}

// 	return isConstantPolyExpr(u[0]);
// }

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

  return expr(kind::ADD,
              {getColPolyNormFactor(lc, rL, K) * pow(L[0], 0)});
}

expr normalizePolyExpr(expr &u, expr &L, expr &K) {
  if (isZeroPolyExpr(u))
    return expr(u);
  expr k = getColPolyNormFactor(u, L, K);
  return quoPolyExpr(u, k, L, K);
}

// expr gcdPolyExpr(expr &u, expr &v, expr &L, expr &K) {
//   if (isZeroPolyExpr(u)) {
//     return normalizePolyExpr(v, L, K);
//   }

//   if (isZeroPolyExpr(v)) {
//     return normalizePolyExpr(u, L, K);
//   }

//   expr g = colPolyGCDRec(u, v, L, K);

//   return normalizePolyExpr(g, L, K);
// }

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
	if(L.size() == 0) return 1;

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
	if(is(&u, kind::TERMINAL) || u.kind() == kind::FUNC) {
		return expr(u);
	}

	if(u.kind() == kind::POW) {
		if(u[1] == 0) return 1;
		if(u[1] == 1) return u[0];
		return expr(u);
	}

  if(u.kind() == kind::MUL) {
		expr g = 1;

		for(Int i = 0; i < u.size(); i++) {
			expr c = expandPolyExpr(u[i]);

			if(c == 0) {
				return 0;
			}

			if(c == 1) {
				continue;
			}

			if(g == 1) {
				g = c;
				continue;
			}

			if(
				 c.kind() == kind::POW ||
				 c.kind() == kind::INT ||
				 c.kind() == kind::SYM ||
				 c.kind() == kind::FRAC ||
				 c.kind() == kind::FUNC
			 ) {
				if(g.kind() == kind::ADD) {
					expr t = expr(kind::ADD);

					for(Int j = 0; j < g.size(); j++) {
						t = t + g[j]*c;
					}

					g = t;
				} else {
					g = g*c;
				}

				continue;
			}

			// Multiplication * Multiplication
			if(c.kind() == kind::MUL && g.kind() == kind::MUL) {
				for(Int j = 0; j < c.size(); j++) {
					g.insert(c[j]);
				}

				continue;
			}

			// Addition * Addition
			if(c.kind() == kind::ADD && g.kind() == kind::ADD) {
				expr t = expr(kind::ADD);

				for(Int j = 0; j < g.size(); j++) {
					for(Int k = 0; k < c.size(); k++) {
						t = t + g[j]*c[k];
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

			for(Int j = 0; j < a.size(); j++) {
				t = t + b*a[j];
			}

			g = t;
		}

		if(g.kind() == kind::ADD && g.size() == 1) return g[0];
		if(g.kind() == kind::MUL && g.size() == 1) return g[0];
		return g;
	}

	assert(u.kind() == kind::ADD);

	expr g = expr(kind::ADD);

  for (Int i = 0; i < u.size(); i++) {
		expr c = expandPolyExpr(u[i]);

		if (c != 0) {
			if(c.kind() == kind::ADD) {
				for(Int j = 0; j < c.size(); j++) {
					g = g + c[j];
				}
			} else {
				g = g + c;
			}
		}
  }


	if(g.size() == 0) return 0;
	if(g.size() == 1) return g[0];

	return g;
}

expr expandPolyExpr(expr& u) {
	return expandPolyExpr(std::forward<expr>(u));
}

expr diffPolyExprRec(expr &u, expr& x, bool* was_diff = nullptr) {
	if(u.kind() == kind::INT || u.kind() == kind::FRAC) {
		return u;
	}

	if(u.kind() == kind::MUL) {
		assert(u.size() == 2);

		if(u[1][0] == x) {
			expr d = u[1][1];

			assert(d.kind() == kind::INT);

			if (d == 0) return raisePolyExpr(0, 0, x);

			expr c = mulPolyExpr(d, u[0]);

			*was_diff = true;

			return c * pow(x, d.value() - 1);
		}

		bool cond = false;

		expr udx = diffPolyExprRec(u[0], x, &cond);

		if(cond == false) return raisePolyExpr(0, 0, x);

		*was_diff = true;

		return udx*u[1];
	}

  expr g = expr(kind::ADD);

	for (Int i = 0; i < u.size(); i++) {
		expr c = diffPolyExprRec(u[i], x, was_diff);

		if(isZeroPolyExpr(c)) continue;

		g.insert(c);
  }

	if(g.size() == 0) {
		return raiseToExpression(0, u);
	}

	return g;
}

expr diffPolyExpr(expr &u, expr& x) {
	bool was_diff = false;

	return diffPolyExprRec(u, x, &was_diff);
}

Int norm(expr u, expr L, expr K, size_t i = 0) {
	if(u == 0 || i == L.size()) {
		assert(u.kind() == kind::INT);
		return u.value();
	}

	Int k = 0;

	expr q, p, t, c, n;

	n = degree(u, L[i]);

	if(n == -inf()) {
		return 0;
	}

	p = expand(u);

	for(Int e = n.value(); e >= 0; e--)
	{
		c = coeff(u, L[i], e);

		k = max(abs(norm(c, L, K, i + 1)), k);

		t = c * pow(L[i], e);

		q = subPoly(p, t);

		p = expand(q);
	}

	return k;
}

Int normPolyExpr(expr u, expr L, expr K, size_t i = 0)
{
	if(i == L.size())
	{
		assert(u.kind() == kind::INT);

		return u.value();
	}

	Int k = 0;

	for(Int j =0; j < u.size(); j++) {
		k = max(abs(normPolyExpr(u[j][0], L, K, i+1)), k);
	}

	return k;
}

expr replaceAndReduce(expr u, expr x, expr c) {
	expr g = replace(u, x, c);
	return  reduce(g);
}

expr evalPolyExpr(expr u, expr x, Int c) {
	if(is(&u, kind::CONST)) {
		return u;
	}

	expr g = expr(kind::ADD);

	for(Int i = 0; i < u.size(); i++) {
		if(base(u[i][1]) == x) {
			expr e = degree(u[i][1]);

			expr k = pow(c, e.value());

			expr t = mulPolyExpr(u[i][0], k);

			g = addPolyExpr(g, t);
		} else {
			g = g + (evalPolyExpr(u[i][0], x, c))*u[i][1];
		}
	}

	if(!is(&g, kind::TERMINAL) && g.size() == 0) {
		return 0;
	}

	if(g.size() == 1 && is(&g[0], kind::TERMINAL | kind::FRAC)) {
		return g[0];
	}

	return g;
}

expr evalPolyExpr(expr u, expr x, expr c) {
	if(is(&u, kind::CONST)) {
		return u;
	}

	expr g = expr(kind::ADD);

	for(Int i = 0; i < u.size(); i++) {
		if(base(u[i][1]) == x) {
			expr e = degree(u[i][1]);
			expr k = powPolyExpr(c, e.value());
			expr t = mulPolyExpr(u[i][0], k);

			g = addPolyExpr(g, t);
		} else {
			g = g + evalPolyExpr(u[i][0], x, c)*u[i][1];
		}
	}

	if(is(&g, kind::ADD) && g.size() == 0) {
		return 0;
	}

	if(is(&g, kind::ADD) && g.size() == 1 && is(&g[0], kind::TERMINAL)) {
		return g[0];
	}


	return g;
}





expr evalTailPolyExpr(expr u, expr L, expr A, Int from) {
	assert(L.kind() == kind::LIST);
	assert(A.kind() == kind::LIST);

	// if(L.size() == 0) return u;
	// if(A.size() == 0) return u;

	// if(A.size() <= from) return u;
	// if(L.size() <= from) return u;

	// if(u.kind() == kind::INT || u.kind() == kind::FRAC) {
	// 	return u;
	// }

	expr t = u;

	for(Int i = from; i < L.size(); i++) {
		t = evalPolyExpr(t, L[i], A[i - from]);
	}

	return t;
	// expr g = expr(kind::ADD);

	// for(Int i = 0; i < u.size(); i++) {
	// 	if(u[i][1][0] == L[from]) {
	// 		expr e = expoent(u[i][1]);
	// 		expr k = pow(A[from].value(), e.value());
	// 		expr t = evalTailPolyExpr(u[i][0], L, A, from + 1);
	// 		expr f = mulPolyExpr(t, k);

	// 		g = addPolyExpr(g, f);
	// 	} else {
	// 		g = g + evalTailPolyExpr(u[i][0], L, A, from)*u[i][1];
	// 	}
	// }

	// if(g.size() == 0) return 0;
	// if(g.size() == 1 && g[0].isTerminal()) return g[0];

	// return g;
}



expr interpolate(expr h, Int p, expr x, expr R, expr K) {
  expr f = 0;

  Int i = 0;

  while (h != 0) {
    expr g = gf(h, p, true);

    f = g * pow(x, i) + f;

    h = subPoly(h, g);
    h = recQuotient(h, p, R, K);

    i = i + 1;
  }

	f = reduce(f);

  if (f == 0) {
    return 0;
	}

	expr lc = groundLeadCoeffPoly(leadCoeff(f, x), R);

  if (lc.value() < 0) {
    f = mulPoly(f, -1);
  }

  return f;
}

expr invertPolyExpr(expr f) {
	expr i = -1;
	return mulPolyExpr(f, i);
}

expr interpolatePolyExpr(expr h, Int p, expr x, expr R, expr K) {

	expr f = expr(kind::ADD);

	expr y = polyExpr(p, R);

	Int i = 0;

	while(!isZeroPolyExpr(h)) {
		expr g = gfPolyExpr(h, p, true);

		if(!isZeroPolyExpr(g)) f.insert(g*pow(x, i));

		h = subPolyExpr(h, g);
		h = quoPolyExpr(h, y, R, K);
		i = i + 1;
	}

	if(f.size() == 0) {
		return raisePolyExpr(polyExpr(0, R), 0, x);
	}

	expr lc = groundLeadCoeffPolyExpr(f);

	if(lc.value() < 0) {
		return invertPolyExpr(f);
	}

	return f;
}

expr groundContRec(expr f, expr L, expr K) {
  if (is(&f, kind::CONST)) {
    return f;
  }

	expr p = f;

  expr g = 0;

  expr r, u, e, t, x, R, z;

  x = L[0];

  R = rest(L);

  expr d = degree(f, x);

  while (p != 0) {
	  // printf("p = %s\n", to_string(p).c_str());
		t = leadCoeff(p, x);
	  // printf("t = %s\n", to_string(t).c_str());

		z = groundContRec(t, R, K);
	  // printf("z = %s\n", to_string(z).c_str());

		g = gcd(g, z);

		// printf("gcd(%s, %s) = %s\n", to_string(g).c_str(), to_string(z).c_str(), to_string(g).c_str());

    e = pow(x, degree(p, x));

		// printf("e = %s\n", to_string(e).c_str());

		u = mulPoly(t, e);

		// printf("---> (%s) * (%s) = %s\n", to_string(t).c_str(), to_string(e).c_str(),  to_string(u).c_str());

		// printf("**************************\n");
		// printf("**************************\n");
		// printf("**************************\n");
		// printf("**************************\n");
		t = subPoly(p, u);

		// printf("**************************\n");
		// printf("**************************\n");
		// printf("**************************\n");
		// printf("**************************\n");
		// printf("\n---> (%s) - (%s) = %s\n", to_string(p).c_str(), to_string(u).c_str(),  to_string(t).c_str());

    p = t;
  }

  return g;
}


expr groundContPoly(expr f, expr L, expr K) {
	return groundContRec(f, L, K);
}

expr groundContPolyExprRec(expr f) {
  if (f.kind() == kind::INT) {
    return f.value();
  }

	if(f.kind() == kind::MUL) {
		assert(f.size() == 2);
		return groundContPolyExprRec(f[0]);
	}

  expr g = 0;

  expr r, u, e, t, x, R;

	for(Int i = 0; i < f.size(); i++) {
		expr l = leadCoeffPolyExpr(f[i][0]);
		expr t = groundContPolyExprRec(l);

		g = gcd(g, t);
	}

  return g;
}

expr groundContPolyExpr(expr f) {
  return groundContPolyExprRec(f);
}

expr groundPPPoly(expr u, expr L, expr K) {
	expr c = groundContPoly(u, L, K);
	return recQuotient(u, c, L, K);
}

expr heuristicGcdPoly(expr u, expr v, expr L, expr K) {
	// References:
	// [1] Liao, H.-C., & Fateman, R. J.(1995). Evaluation of the heuristic polynomial GCD
	assert(K == expr("Z"));

	if ((u.kind() == kind::INT && v.kind() == kind::INT) || L.size() == 0) {
    assert(u.kind() == kind::INT);
    assert(v.kind() == kind::INT);

		Int g = abs(gcd(u.value(), v.value()));

		return list({ g, u.value() / g, v.value() / g});
  }
	// printf("aaa  %s\n", to_string(u).c_str());
  expr&& ucont = groundContPoly(u, L, K);
	// printf("aaa  %s\n", to_string(v).c_str());
  expr&& vcont = groundContPoly(v, L, K);

	// printf("gc = %s\n", to_string(ucont).c_str());
	// printf("vc = %s\n", to_string(vcont).c_str());

  expr&& g = gcd(ucont, vcont);

	// printf("g = %s\n", to_string(g).c_str());

  u = recQuotient(u, g, L, K);
  v = recQuotient(v, g, L, K);

	// printf("u = %s\n", to_string(u).c_str());
	// printf("v = %s\n", to_string(v).c_str());

	Int un = norm(u, L, K);
  Int vn = norm(v, L, K);

  Int b = 2 * min(un, vn) + 29;
  Int uc = groundLeadCoeffPoly(u, L).value();
  Int vc = groundLeadCoeffPoly(v, L).value();

  Int x = max(min(b, 99 * isqrt(b)), 2 * min(un / uc, vn / vc) + 2);

	// printf("uc = %s\n", to_string(uc).c_str());
	// printf("vc = %s\n", to_string(vc).c_str());


  for (short i = 0; i < 6; i++) {
		////printf("----> u = %s\n", to_string(u).c_str());
		// printf("x = %s\n", to_string(x).c_str());
		expr&& ux = replaceAndReduce(u, L[0], x);
    expr&& vx = replaceAndReduce(v, L[0], x);
		// printf("ux = %s\n", to_string(ux).c_str());
		// printf("vx = %s\n", to_string(vx).c_str());

		expr&& R = rest(L);

    if (ux != 0 && vx != 0) {
      expr&& HGCD = heuristicGcdPoly(ux, vx, R, K);

      expr& h = HGCD[0];

			expr& a = HGCD[1];
      expr& b = HGCD[2];

      //expr cv, rv;

      h = interpolate(h, x, L[0], R, K);
      h = groundPPPoly(h, L, K);

      expr U = recPolyDiv(u, h, L, K);

			expr& cu = U[0];
      expr& ru = U[1];

      if (ru == 0) {
        expr&& V = recPolyDiv(v, h, L, K);

        expr& cv = V[0];
        expr& rv = V[1];

        if (rv == 0) {

          return list({h, cu, cv});
        }
      }

      a = interpolate(a, x, L[0], R, K);

      U = recPolyDiv(u, a, L, K);

			h = U[0];
      ru = U[1];

			if (ru == 0) {
        expr&& V = recPolyDiv(v, h, L, K);

        expr& cv = V[0];
        expr& rv = V[1];

        if (rv == 0) {
          h = mulPoly(h, g);

          return list({h, a, cv});
        }
      }

      b = interpolate(b, x, L[0], R, K);

      expr&& V = recPolyDiv(v, b, L, K);

      h = V[0];
      expr& rv = V[1];

      if (rv == 0) {
        U = recPolyDiv(u, h, L, K);
        expr& c = U[0];
        ru = U[1];

        if (ru == 0) {
          h = mulPoly(h, g);
          return list({h, c, b});
        }
      }
    }

		x = 73794 * x * isqrt(isqrt(x)) / 27011;
	}

	return fail();
}

expr groundDivPolyExpr(expr u, expr v) {
	assert(v.kind() == kind::INT);

	if(isZeroPolyExpr(u)) return u;

	if(u.kind() == kind::INT) {
		return u.value() / v.value();
	}

	expr g = expr(kind::ADD);

	for(Int i = 0; i < u.size(); i++) {
		expr k = groundDivPolyExpr(u[i][0], v);

		if(!isZeroPolyExpr(k)) {
			g.insert(k*u[i][1]);
		}
	}

	if(g.size() == 0) {
		return raiseToExpression(0, u);
	}

	return g;
}


expr groundMulPolyExpr(expr u, Int v) {
	if(isZeroPolyExpr(u)) return u;

	if(u.kind() == kind::INT) {
		return u.value() * v;
	}

	expr g = expr(kind::ADD);

	for(Int i = 0; i < u.size(); i++) {
		expr k = groundMulPolyExpr(u[i][0], v);

		if(!isZeroPolyExpr(k)) {
			g.insert(k*u[i][1]);
		}
	}

	if(g.size() == 0) {
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
	// [1] Liao, H.-C., & Fateman, R. J.(1995). Evaluation of the heuristic polynomial GCD

	assert(K == expr("Z"));

	if ((isConstantPolyExpr(u) && isConstantPolyExpr(v)) || L.size() == 0) {
		Int a = groundLeadCoeffPolyExpr(u).value();
		Int b = groundLeadCoeffPolyExpr(v).value();

		Int c = abs(gcd(a, b));

		if(L.size()) {
			return list({ polyExpr(c, L), polyExpr(a/c, L), polyExpr(b/c, L) });
		}

		return list({ c, a/c, b/c });
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
			//a = raisePolyExpr(a, 0, L[0]);

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
			//b = raisePolyExpr(b, 0, L[0]);

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
expr groundInvert(expr p) {
  return mulPoly(p, -1);
}

expr groundInvertPolyExpr(expr p) {
	expr k = -1;
	return mulPolyExpr(k, p);
}


expr insertSymbolPolyExpr(expr u, expr x, Int d, Int level, Int i) {
	if(is(&u, kind::TERMINAL)) {
		return expr(kind::ADD, { u*pow(x, d) });
	}

	assert(u.kind() == kind::ADD);

	if(i == level) {

		expr g = expr(kind::ADD);

		for(Int j = 0; j < u.size(); j++) {
			g.insert(expr(kind::ADD, { u[j][0]*pow(x, d) })*u[j][1]);
		}

		return g;
	}

	expr g = expr(kind::ADD);

	for(Int j = 0; j < u.size(); j++) {
		g.insert(insertSymbolPolyExpr(u[j][0], x, d, level, i + 1)*u[j][1]);
	}

	return g;
}

expr insertSymbolPolyExpr(expr u, expr x, Int d, Int level) {
	return insertSymbolPolyExpr(u, x, d, level, 0);
}

expr insertSymbolsPolyExpr(expr u, expr L, Int d, Int level) {
	expr g = u;

	for(Int i = 0; i < L.size(); i++) {
		g = insertSymbolPolyExpr(g, L[i], d, level + i);
	}

	return g;
}


expr degreePolyExpr(expr u, expr x) {
	if(is(&u, kind::TERMINAL)) return -inf();
	if(u.size() == 0) return -inf();

	if(base(u[u.size() - 1][1]) == x) {
		return degree(u[u.size() - 1][1]);
	}

	expr m = -inf();

	for(Int i = 0; i < u.size(); i++) {
		expr d = degreePolyExpr(u[i][0], x);

		assert(d.kind() == kind::INT);

		if(m == -inf()) {
			m = d;
		} else {
			m = max(m.value(), d.value());
		}
	}

	return m;
}

} // namespace polynomial
