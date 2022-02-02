/**
 * @file GaloisField.cpp
 * @author Marcos Vincius Moreira Santos (marcos30004347@gmail.com)
 * @brief This file implement some of the Finite field methods used in the
 * library
 * @version 0.1
 * @date 2021-10-11
 *
 * @ref Michael Monagan - In-place Arithmetic for Polynomials over Zn
 *
 * @copyright Copyright (c) 2021
 */

#include "GaloisField.hpp"

#include "Core/AST/Integer.hpp"
#include "Core/Debug/Assert.hpp"

#include "Core/Polynomial/Polynomial.hpp"

#include <climits>
#include <cstddef>
#include <limits>
#include <random>

using namespace alg;
using namespace polynomial;

namespace galoisField {

Int mod(Int a, Int b, bool symmetric) {
  Int n = (b + (a % b)) % b;

  if (symmetric) {
    if (0 <= n && n <= b / 2) {
      return n;
    }

    return n - b;
  }

  return n;
}

Int randomGf(Int p, bool symmetric) {
  std::random_device dev;

  std::mt19937 rng(dev());

  std::uniform_int_distribution<std::mt19937::result_type> dist(
      std::numeric_limits<long long>::min(),
      std::numeric_limits<long long>::max());

  return mod((long long)dist(rng), p, symmetric);
}
// Function for extended Euclidean Algorithm
Int gcdExtended(Int a, Int b, Int* x, Int* y, bool symmetric)
{

    // Base Case
    if (a == 0)
    {
        *x = 0, *y = 1;
        return b;
    }

    Int x1, y1;

    Int gcd = gcdExtended(mod(b, a, symmetric), a, &x1, &y1, symmetric);

    // Update x and y using results of recursive
    // call
    *x = y1 - (b / a) * x1;
    *y = x1;

    return gcd;
}
// Function to find modulo inverse of a
Int _inverseGf(Int a, Int m, bool symmetric)
{
    Int x, y;
    Int g = gcdExtended(a, m, &x, &y, symmetric);
    if (g != 1) {

			printf("%s have no inverse mod %s\n", a.to_string().c_str(),
						 m.to_string().c_str());
			exit(1);
			//TODO: error
    }

		return mod(x, m, symmetric);// (x % m + m) % m;
}



Int inverseGf(Int a, Int b, bool symmetric) {
  Int t, nt, r, nr, q, tmp;

	if(b < 0) b = -b;

	// all the computations are made on the
	// non-symetric representation and
	// converted back at the end to its
	// right representation
	a = mod(a, b, false);

	t = 0;
  nt = 1;
  r = b;

  nr = mod(a, b, false);

  while (nr != 0) {
    q = r / nr;
    tmp = nt;
    nt = t - q * nt;
    t = tmp;
    tmp = nr;
    nr = r - q * nr;
    r = tmp;
  }

  if (r > 1) {
		// TODO: better error handling
		printf("%s have no inverse mod %s\n", a.to_string().c_str(),
           b.to_string().c_str());

    exit(1);
  }

  return mod(t, b, symmetric);
}

Int quoGf(Int s, Int t, Int p, bool symmetric) {
  return mod((s * inverseGf(t, p, symmetric)), p, symmetric);
}

// expr gf(expr u, expr x, Int s, bool symmetric) {
//   expr k = expand(u);

//   if (k.kind() == Kind::FAIL || k.kind() == Kind::UNDEF) {
//     return k;
//   }

//   if (k.kind() == Kind::MinusINF || k.kind() == Kind::INF) {
//     return undefined();
//   }

//   if (k.kind() == Kind::Integer) {
//     Int p = k.value();

//     return integer(mod(p, s, symmetric));
//   }

//   if (k.kind() == Kind::SYM) {
//     return k;
//   }

//   if (k.kind() == Kind::FRAC) {
//     assert(k[0].kind() == Kind::Integer,
//            "numerator of a fraction needs to be a integer");

//     assert(k[1].kind() == Kind::Integer,
//            "denominator of a fraction needs to be a integer");

//     Int n = k[0].value();
//     Int d = k[1].value();

//     return mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s,
//     symmetric);
//   }

//   if (k.kind() == Kind::Derivative) {
//     expr p = expr(Kind::Derivative, {gf(k[0], x, s, symmetric), k[1]});

//     return p;
//   }

//   if (k.kind() == Kind::Integral) {
//     expr p = expr(Kind::Integral, {gf(k[0], x, s, symmetric), k[1]});

//     return p;
//   }

//   if (k.kind() == Kind::FACT) {
//     if (k[0].kind() == Kind::Integer) {
//       expr f = reduce(k);
//       return gf(f, x, s, symmetric);
//     }

//     return k;
//   }

//   if (k.kind() == Kind::Division) {
//     expr p = div(gf(k[0], x, s, symmetric), gf(k[1], x, s, symmetric));

//     return gf(reduce(p), x, s, symmetric);
//   }

//   if (k.kind() == Kind::POW) {
//     expr p = pow(gf(k[0], x, s, symmetric), k[1]);

//     return p;
//   }

//   if (k.kind() == Kind::FunctionCall) {
//     return k;
//   }

//   if (k.kind() == Kind::MUL) {
//     expr p = expr(Kind::MUL);

//     for (long i = 0; i < k.size(); i++) {
//       p.insert(gf(k[i], x, s, symmetric));
//     }

//     return p;
//   }

//   if (k.kind() == Kind::ADD || k.kind() == Kind::SUB) {
//     expr p = expr(k.kind());

//     expr d = degree(k, x);

//     for (Int n = 0; n <= d.value(); n++) {
//       expr c = coeff(k, x, n);
//       p.insert(gf(c, x, s, symmetric) * pow(x, n));
//     }

//     return reduce(p);
//   }

//   return k;
// }

expr groundGf(expr u, Int s, bool symmetric) {
  expr k = u;

  if (k.kind() == kind::FAIL || k.kind() == kind::UNDEF) {
    return k;
  }

  // if (k.kind() == kind::MinusINF || k.kind() == kind::INF) {
  //   return undefined();
  // }

  if (k.kind() == kind::INT) {
    return integer(mod(k.value(), s, symmetric));
  }

  if (k.kind() == kind::SYM) {
    return k;
  }

  if (k.kind() == kind::FRAC) {
    assert(k[0].kind() == kind::INT,
           "numerator of a fraction needs to be a integer");

    assert(k[1].kind() == kind::INT,
           "denominator of a fraction needs to be a integer");

    Int n = k[0].value();
    Int d = k[1].value();

    return mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric);
  }

  // if (k.kind() == kind::Derivative) {
  //   expr p = expr(kind::Derivative, {groundGf(k[0], s, symmetric), k[1]});

  //   return p;
  // }

  // if (k.kind() == kind::Integral) {
  //   return expr(kind::Integral, {groundGf(k[0], s, symmetric), k[1]});
  // }

  if (k.kind() == kind::FACT) {
    if (k[0].kind() == kind::INT) {
      return groundGf(fact(k.value()), s, symmetric);
    }

    return k;
  }

  if (k.kind() == kind::DIV) {
    expr p = groundGf(k[0], s, symmetric) / groundGf(k[1], s, symmetric);

		expr t = reduce(p);

    return groundGf(t, s, symmetric);
  }

  if (k.kind() == kind::POW) {
    return pow(groundGf(k[0], s, symmetric), k[1]);
  }

  if (k.kind() == kind::FUNC) {
    return k;
  }

  if (k.kind() == kind::MUL) {
    expr p = expr(kind::MUL);

    for (size_t i = 0; i < k.size(); i++) {
      p.insert(groundGf(k[i], s, symmetric));
    }

    return reduce(p);
  }

  if (k.kind() == kind::ADD || k.kind() == kind::SUB) {
    expr p = expr(k.kind());

    for (size_t i = 0; i < k.size(); i++) {
      p.insert(groundGf(k[i], s, symmetric));
    }

    expr r = reduce(p);

    return r;
  }

  return k;
}

expr gf(expr u, Int s, bool symmetric) {
  expr k = expand(u);

  if (k.kind() == kind::FAIL || k.kind() == kind::UNDEF) {
    return k;
  }

  // if (k.kind() == kind::MinusINF || k.kind() == kind::INF) {
  //   return undefined();
  // }

  if (k.kind() == kind::INT) {
    Int p = k.value();

    return mod(p, s, symmetric);
  }

  if (k.kind() == kind::SYM) {
    return k;
  }

  if (k.kind() == kind::FRAC) {
    assert(k[0].kind() == kind::INT,
           "numerator of a fraction needs to be a integer");

    assert(k[1].kind() == kind::INT,
           "denominator of a fraction needs to be a integer");

    Int n = k[0].value();
    Int d = k[1].value();

    return integer(
        mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric));
  }

  // if (k.kind() == kind::Derivative) {
  //   expr p = expr(kind::Derivative, {gf(k[0], s, symmetric), k[1]});

  //   return p;
  // }

  // if (k.kind() == kind::Integral) {
  //   expr p = expr(kind::Integral, {gf(k[0], s, symmetric), k[1]});

  //   return p;
  // }

  if (k.kind() == kind::FACT) {
    if (k[0].kind() == kind::INT) {
      expr f = reduce(k);

      expr p = gf(f, s, symmetric);

      return p;
    }

    return k;
  }

  if (k.kind() == kind::DIV) {
    expr p = gf(k[0], s, symmetric) / gf(k[1], s, symmetric);
    expr t = reduce(p);
    expr r = gf(t, s, symmetric);

    return r;
  }

  if (k.kind() == kind::POW) {
    expr p = pow(gf(k[0], s, symmetric), k[1]);

    return p;
  }

  if (k.kind() == kind::FUNC) {
    return k;
  }

  if (k.kind() == kind::MUL) {
    expr p = expr(kind::MUL);
    for (size_t i = 0; i < k.size(); i++) {
      p.insert(gf(k[i], s, symmetric));
    }

    return reduce(p);
  }

  if (k.kind() == kind::ADD || k.kind() == kind::SUB) {
    expr p = expr(k.kind());
    for (size_t i = 0; i < k.size(); i++) {
      p.insert(gf(k[i], s, symmetric));
    }

    return reduce(p);
  }

  return k;
}

expr divPolyGf(expr a, expr b, expr x, Int p, bool symmetric) {

	expr da = degree(a, x);
  expr db = degree(b, x);

	if(db == -inf()) {
		// TODO: throw division by zero
		return undefined();
	}

  if (da == -inf() || da.value() < db.value()) {
    return list({0, a});
  }

  long long k, j, d;
  Int s, e, lb;

  expr dq, dr, q, r;
  expr t1, t2, t3, ex;

  std::vector<expr> A = std::vector<expr>(da.value().longValue() + 1);
  std::vector<expr> B = std::vector<expr>(db.value().longValue() + 1);

  for (k = da.value().longValue(); k >= 0; k--) {
    A[k] = coeff(a, x, k);
  }

  for (k = db.value().longValue(); k >= 0; k--) {
    B[k] = coeff(b, x, k);
  }

  dq = da.value() - db.value();

  dr = db.value() - 1;

  t1 = leadCoeff(b, x);
	//printf("lc = %s\n", t1.toString().c_str());

	//printf("lc = %s\n", t1.value().to_string().c_str());
	// TODO: remove the false from here
	lb = inverseGf(t1.value(), p, false);

	//printf("lb = %s\n", lb.to_string().c_str());

	//	printf("lb = %s\n", mod(t1.value()*lb, p, symmetric).to_string().c_str());
  for (k = da.value().longValue(); k >= 0; k--) {
    t1 = A[k];

    s = max(0, k - dq.value());
    e = min(dr.value(), k);
		//printf("h = ");
		//for(int i = 0; i < A.size(); i++) {
		//	printf("%s ", A[i].toString().c_str());
		//}
		//printf("\n");
    for (j = s.longValue(); j <= e; j++) {
			//printf("%s %s\n", A[k - j + db.value().longValue()].toString().c_str(), B[j].toString().c_str());
			t2 = mulPoly(B[j], A[k - j + db.value().longValue()]);
      t1 = subPoly(t1, t2);
    }

    t1 = reduce(t1);
		// printf("%s\n", p.to_string().c_str());
		// printf("%s\n", t1.toString().c_str());
    t1 = mod(t1.value(), p, symmetric);

		//printf("%s\n", t1.toString().c_str());
    // if (t1.value() < 0) {
    //   t1 = mod(t1.value() + p, p, symmetric);
    // }

    if (da.value() - k <= dq.value()) {
      t1 = reduce(mulPoly(t1, lb));
			//printf("%s\n", t1.toString().c_str());
			//printf("%s\n", lb.to_string().c_str());
			// TODO: remove the false from here
			t1 = mod(t1.value(), p, false);
    }

    A[k] = t1;
  }

  q = create(kind::ADD, {0});
  r = create(kind::ADD, {0});

  d = 0;

  for (k = da.value().longValue() - dq.value().longValue(); k <= da.value();
       k++) {
    q.insert(A[k] * pow(x, integer(d)));
    d = d + 1;
  }

  d = 0;

  for (k = 0; k <= dr.value(); k++) {
    r.insert(A[k] * pow(x, integer(d)));

    d = d + 1;
  }
  t1 = gf(q, p, symmetric);
  t2 = gf(r, p, symmetric);

  return list({t1, t2});
}

expr remPolyGf(expr a, expr b, expr x, Int p, bool symmetric) {
  return divPolyGf(a, b, x, p, symmetric)[1];
}

expr quoPolyGf(expr a, expr b, expr x, Int p, bool symmetric) {
  return divPolyGf(a, b, x, p, symmetric)[0];
}

expr monicPolyGf(expr f, expr x, Int p, bool symmetric) {
  if (f == 0) {
    return 0;
  }

  expr lc = leadCoeff(f, x);
  expr F = quoPolyGf(f, lc, x, p, symmetric);

  return list({lc, F});
}

expr gcdPolyGf(expr a, expr b, expr x, Int p, bool symmetric) {
  expr da = degree(a, x);
  expr db = degree(b, x);

  if (da == -inf() || db.value() > da.value()) {
    return gcdPolyGf(b, a, x, p, symmetric);
  }

  expr t;
  while (b != 0 && db != -inf() && db.value() >= 0) {
		t = a;
    a = b;

		b = remPolyGf(t, b, x, p, symmetric);
		db = degree(b, x);
  }

  b = monicPolyGf(a, x, p, symmetric);
  return b[1];
}

expr addPolyGf(expr f, expr g, expr x, Int p, bool symmetric) {
  expr u = addPoly(f, g);

  return gf(u, p, symmetric);
}

expr subPolyGf(expr f, expr g, expr x, Int p, bool symmetric) {
  expr t, u;

  u = subPoly(f, g);

  t = gf(u, p, symmetric);

  return t;
}

expr mulPolyGf(expr f, expr g, expr x, Int p, bool symmetric) {
  expr t, u;

  u = mulPoly(f, g);

  t = gf(u, p, symmetric);

  return t;
}

expr powModPolyGf(expr f, expr g, expr x, Int n, Int p, bool symmetric) {
  if (n == 0)
    return 1;

  expr a = f;
  expr b = 1;
  expr t = 0;

  while (n > 1) {
    if (n % 2 == 0) {
      t = mulPolyGf(a, a, x, p, symmetric);
			a = remPolyGf(t, g, x, p, symmetric);
      n = n / 2;
    } else {
      t = mulPolyGf(a, b, x, p, symmetric);
      b = remPolyGf(t, g, x, p, symmetric);
      t = mulPolyGf(a, a, x, p, symmetric);
      a = remPolyGf(t, g, x, p, symmetric);
      n = (n - 1) / 2;
    }
  }

  t = mulPolyGf(a, b, x, p, symmetric);
  return remPolyGf(t, g, x, p, symmetric);
}

expr randPolyGf(Int d, expr x, Int p, bool symmetric) {
  Int k = 0;

  if (d == 0) {
    return randomGf(p, symmetric);
  }

  if (d == 1) {
    k = randomGf(p, symmetric);

    if (k == 0)
      return x;

    return x + k;
  }

  expr r = expr(kind::ADD, {pow(x, d)});

  for (Int i = d - 1; i >= 2; i--) {
    k = randomGf(p, symmetric);

    if (k != 0) {
      r = r + k * pow(x, i);
    }
  }

  k = randomGf(p, symmetric);

  if (k != 0) {
    r = r + k * x;
  }

  k = randomGf(p, symmetric);

  if (k != 0) {
    r = r + k;
  }

  return r;
}

expr extendedEuclidGf(expr f, expr g, expr x, Int p, bool sym) {
  if (f == 0 || g == 0) {
    return list({1, 0, 0});
  }

  expr t, s, i, lc, k1, t0, t3, s0, s1, Q, R, T;
  expr t1 = monicPolyGf(f, x, p, sym);
  expr t2 = monicPolyGf(g, x, p, sym);
  expr p0 = t1[0];
  expr r0 = t1[1];

  expr p1 = t2[0];
  expr r1 = t2[1];

  if (f == 0) {
    return list({0, inverseGf(p1.value(), p, sym), r1});
  }

  if (g == 0) {
    return list({inverseGf(p0.value(), p, sym), 0, r0});
  }


  s0 = inverseGf(p0.value(), p, sym);
  s1 = 0;
  t0 = 0;
  t1 = inverseGf(p1.value(), p, sym);

	while (true) {
    T = divPolyGf(r0, r1, x, p, sym);

    Q = T[0];
    R = T[1];

    if (R == 0) {
      break;
    }

    T = monicPolyGf(R, x, p, sym);

    r0 = r1;

    lc = T[0];
    r1 = T[1];

    i = inverseGf(lc.value(), p, sym);

    k1 = mulPolyGf(s1, Q, x, p, sym);
    s = subPolyGf(s0, k1, x, p, sym);

    k1 = mulPolyGf(t1, Q, x, p, sym);
    t = subPolyGf(t0, k1, x, p, sym);

    s0 = s1;
    t0 = t1;

    s1 = mulPolyGf(s, i, x, p, sym);
    t1 = mulPolyGf(t, i, x, p, sym);
  }

  return list({r1, s1, t1});
}

expr gfPolyExpr(expr u, Int p, bool symmetric) {
  if (u.kind() == kind::INT) {
    return mod(u.value(), p, symmetric);
  }

  if (u.kind() == kind::MUL) {
    assert(u.size() == 2, "not a polynomial expr");
    return gfPolyExpr(u[0], p, symmetric) * u[1];
  }

  assert(u.kind() == kind::ADD, "not a polynomial expr");

  expr g = expr(kind::ADD);
	expr x = 0;

  for (Int i = 0; i < u.size(); i++) {
    assert(u[i].kind() == kind::MUL && u[i].size() == 2,
           "not a polynomial expr");

    expr c = gfPolyExpr(u[i][0], p, symmetric);

		x = u[i][1][0];

		if (!isZeroPolyExpr(c)) {
      g.insert(c * u[i][1]);
    }
	}

  if (g.size() == 0) {
		return raiseToExpression(0, u);
	}

  return g;
}

expr addPolyExprGf(expr f, expr g, Int p, bool sym) {
  return gfPolyExpr(addPolyExpr(f, g), p, sym);
}

expr subPolyExprGf(expr f, expr g, Int p, bool sym) {
  return gfPolyExpr(subPolyExpr(f, g), p, sym);
}

expr mulPolyExprGf(expr f, expr g, Int p, bool sym) {
  return gfPolyExpr(mulPolyExpr(f, g), p, sym);
}

expr divPolyExprGf(expr a, expr b, expr L, Int p, bool symmetric) {
  assert(L.kind() == kind::LIST && L.size() == 1, "not a univariate poly expr");

	assert(a.kind() == kind::ADD, "not a poly expr");
  assert(b.kind() == kind::ADD, "not a poly expr");

	if(isZeroPolyExpr(a)) {
		return list({polyExpr(0, L), polyExpr(0, L)});
	}

	if(isZeroPolyExpr(b)) {
		// TODO: handle error division by zero
		printf("division by zero\n");
		abort();
	}


	expr x = L[0];

	expr da = degreePolyExpr(a);
  expr db = degreePolyExpr(b);

  assert(da.kind() == kind::INT,
         "degree of polynomial should be an integer\n");
  assert(db.kind() == kind::INT,
         "degree of polynomial should be an integer\n");

  if (da.value() < db.value()) {
    return list({polyExpr(0, L), a});
  }

  //long long j;
  Int s, e;

  Int dq, dr;

  Int t1, t2, t3, lb, d;

  std::vector<Int> A = std::vector<Int>(da.value().longValue() + 1, 0);
  std::vector<Int> B = std::vector<Int>(db.value().longValue() + 1, 0);

  for (size_t k = 0; k < a.size(); k++) {
    assert(a[k].kind() == kind::MUL && a[k].size() == 2,
           "not a poly expr");

    assert(a[k][0].kind() == kind::INT, "not a univariate poly expr");

    expr d = a[k][1][1];

    assert(d.kind() == kind::INT,
           "poly expr should have only integers as degrees");

    assert(a[k][0].kind() == kind::INT,
           "poly expr should have only integers as degrees");

    A[d.value().longValue()] = a[k][0].value();
  }

  for (size_t k = 0; k < b.size(); k++) {
    assert(b[k].kind() == kind::MUL && b[k].size() == 2,
           "not a poly expr");

    assert(b[k][0].kind() == kind::INT, "not a univariate poly expr");

    expr d = b[k][1][1];

    assert(d.kind() == kind::INT,
           "poly expr should have only integers as degrees");
    assert(b[k][0].kind() == kind::INT,
           "poly expr should have only integers as degrees");
    B[d.value().longValue()] = b[k][0].value();
  }

  dq = da.value() - db.value();
  dr = db.value() - 1;

  t1 = leadCoeffPolyExpr(b).value();

	// TODO: remove the false from here
	//printf("t1 = %s\n", t1.to_string().c_str());
  lb = inverseGf(t1, p, false);
	//printf("lb = %s\n", lb.to_string().c_str());
	//printf("mb = %s\n", mod(lb*t1, p, symmetric).to_string().c_str());

	assert(mod(lb*t1, p, symmetric) == 1, "inverse mod p is wrong");

	for (long long k = da.value().longValue(); k >= 0; k--) {

    t1 = A[k];
    s = max(0, k - dq);
    e = min(dr, k);

    for (long long j = s.longValue(); j <= e; j++) {
      t2 = B[j] * A[k - j + db.value().longValue()];
      t1 = t1 - t2;
    }

    t1 = mod(t1, p, symmetric);

    // if (t1 < 0) {
    //   t1 = mod(t1 + p, p, symmetric);
    // }

    if (da.value() - k <= dq) {
			// TODO: remove false from here
      t1 = mod(t1 * lb, p, false);
    }

    A[k] = t1;
  }

  expr q = create(kind::ADD, {});
  expr r = create(kind::ADD, {});

  d = 0;

  for (long long k = da.value().longValue() - dq.longValue(); k <= da.value(); k++) {
    if (A[k] != 0) {
      q.insert(mod(A[k], p, symmetric) * pow(x, d));
    }

    d = d + 1;
  }

  d = 0;

  for (long long  k = 0; k <= dr; k++) {
    if (A[k] != 0) {
      r.insert(mod(A[k], p, symmetric) * pow(x, d));
    }

    d = d + 1;
  }

	if(q.size() == 0) {
		q.insert(0*pow(x, 0));
	}

	if(r.size() == 0) {
		r.insert(0*pow(x, 0));
	}

  return list({q, r});
}


expr quoPolyExprGf(expr a, expr b, expr L, Int p, bool symmetric) {
  return divPolyExprGf(a, b, L, p, symmetric)[0];
}

expr remPolyExprGf(expr a, expr b, expr L, Int p, bool symmetric) {
  return divPolyExprGf(a, b, L, p, symmetric)[1];
}

expr monicPolyExprGf(expr f, expr L, Int p, bool symmetric) {
	expr x = L[0];

  if (isZeroPolyExpr(f)) {
    return expr(kind::ADD, {0 * pow(x, 0)});
  }

  expr lc = expr(kind::ADD, {leadCoeffPolyExpr(f) * pow(x, 0)});

  expr F = quoPolyExprGf(f, lc, L, p, symmetric);

  return list({lc, F});
}


// expr randPolyExprGf(Int d, expr L, Int p, bool symmetric) {
// 	expr x = L[0];

// 	Int k = 0;

//   if (d == 0) {
//     return randomGf(p, symmetric);
//   }

//   if (d == 1) {
//     k = randomGf(p, symmetric);

//     if (k == 0)
//       return x;

//     return polyExpr(x + k, L);
//   }

//   expr r = expr(kind::ADD, { pow(x, d) });

//   for (Int i = d - 1; i >= 2; i--) {
//     k = randomGf(p, symmetric);

//     if (k != 0) {
//       r = r + k * pow(x, i);
//     }
//   }

//   k = randomGf(p, symmetric);

//   if (k != 0) {
//     r = r + k * x;
//   }

//   k = randomGf(p, symmetric);

//   if (k != 0) {
//     r = r + k;
//   }

//   return r;
// }


expr randPolyExprGf(Int d, expr L, Int p, bool symmetric) {
	assert(L.kind() == kind::LIST && L.size() <= 1, "L should be a list with one element");
	expr r = expr(kind::ADD);

  expr x = L[0];

  for (Int i = 0; i < d; i++) {
    Int k = randomGf(p, symmetric);

    if (k != 0) {
			r.insert(k*pow(x, i));
    }
  }

	r.insert(1*pow(x, d));

  return r;
}

expr powModPolyExprGf(expr f, expr g, expr L, Int n, Int p, bool symmetric) {
	assert(L.kind() == kind::LIST, "L should be a list");
	expr b = expr(kind::ADD, {1 * pow(L[0], 0)});

  if (n == 0)
    return b;

  expr a = f;
  expr t = 0;

  while (n > 1) {
    if (n % 2 == 0) {
      t = mulPolyExprGf(a, a, p, symmetric);
      a = remPolyExprGf(t, g, L, p, symmetric);
      n = n / 2;

    } else {
      t = mulPolyExprGf(a, b, p, symmetric);
      b = remPolyExprGf(t, g, L, p, symmetric);

      t = mulPolyExprGf(a, a, p, symmetric);
      a = remPolyExprGf(t, g, L, p, symmetric);
      n = (n - 1) / 2;
    }
  }

  t = mulPolyExprGf(a, b, p, symmetric);
  return remPolyExprGf(t, g, L, p, symmetric);
}

expr gcdPolyExprGf(expr a, expr b, expr L, Int p, bool symmetric) {
  expr da = degreePolyExpr(a);
  expr db = degreePolyExpr(b);

  if (da == -inf() || db.value() > da.value()) {
    return gcdPolyExprGf(b, a, L, p, symmetric);
  }

  expr t;

  while (!isZeroPolyExpr(b) && db != -inf() &&
         db.value() >= 0) {
    t = a;
    a = b;

    b = remPolyExprGf(t, b, L, p, symmetric);
    db = degreePolyExpr(b);
  }

  b = monicPolyExprGf(a, L, p, symmetric);

  return b[1];
}

expr extendedEuclidPolyExprGf(expr f, expr g, expr L, Int p, bool sym) {

  if (f == 0 || g == 0) {
    return list({
        polyExpr(1, L),
        polyExpr(0, L),
        polyExpr(0, L),
    });
  }

  expr t, s, i, lc, k1, t0, t3, s0, s1, Q, R, T;
  expr t1 = monicPolyExprGf(f, L, p, sym);
  expr t2 = monicPolyExprGf(g, L, p, sym);

	expr p0 = leadCoeffPolyExpr(t1[0]);
  expr r0 = t1[1];

  expr p1 = leadCoeffPolyExpr(t2[0]);
  expr r1 = t2[1];

  if (isZeroPolyExpr(f)) {
    return list({
				polyExpr(0, L),
				polyExpr(inverseGf(p1.value(), p, sym), L),
				r1,
			});
  }

  if (isZeroPolyExpr(g)) {
    return list({
				polyExpr(inverseGf(p0.value(), p, sym), L),
				polyExpr(0, L),
				r0
			});
  }

  s0 = polyExpr(inverseGf(p0.value(), p, sym), L);
	s1 = polyExpr(0, L);

	t0 = polyExpr(0, L);
	t1 = polyExpr(inverseGf(p1.value(), p, sym), L);
  while (true) {
		T = divPolyExprGf(r0, r1, L, p, sym);

    Q = T[0];
    R = T[1];

    if (isZeroPolyExpr(R)) {
      break;
    }

    T = monicPolyExprGf(R, L, p, sym);

    r0 = r1;

    lc = leadCoeffPolyExpr(T[0]);
    r1 = T[1];

    assert(lc.kind() == kind::INT, "lc of univariate should be a integer");

    i = polyExpr(inverseGf(lc.value(), p, sym), L);

    k1 = mulPolyExprGf(s1, Q, p, sym);
		s = subPolyExprGf(s0, k1, p, sym);

    k1 = mulPolyExprGf(t1, Q, p, sym);
    t = subPolyExprGf(t0, k1, p, sym);

    s0 = s1;
    t0 = t1;

		s1 = mulPolyExprGf(s, i, p, sym);
    t1 = mulPolyExprGf(t, i, p, sym);
	}

	return list({ r1, s1, t1 });
}

expr groundQuoPolyExprGf(expr u, Int v, Int p, bool sym) {
	Int k = inverseGf(v, p, sym);

	return gf(groundMulPolyExpr(u, k), p, sym);
}

} // namespace galoisField
