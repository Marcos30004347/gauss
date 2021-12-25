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

#include "Core/AST/AST.hpp"
#include "Core/AST/Integer.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Debug/Assert.hpp"

#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"

#include <climits>
#include <limits>
#include <random>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

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

	// if(symmetric && a < 0) {
	// 	a += b/2;
	//	}

	// if (b < 0) b = -b;
	// if (a < 0) a = b - (-a % b);

  t = 0;
  nt = 1;
  r = b;

  nr = mod(a, b, symmetric);
	// nr = a % b;
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
    printf("%s have no inverse mod %s\n", a.to_string().c_str(),
           b.to_string().c_str());

    exit(1);
  }

  //if (t < 0)
  //   t += b;

  return mod(t, b, symmetric);
}

Int quoGf(Int s, Int t, Int p, bool symmetric) {
  return mod((s * inverseGf(t, p, symmetric)), p, symmetric);
}

// Expr gf(Expr u, Expr x, Int s, bool symmetric) {
//   Expr k = algebraicExpand(u);

//   if (k.kind() == Kind::Fail || k.kind() == Kind::Undefined) {
//     return k;
//   }

//   if (k.kind() == Kind::MinusInfinity || k.kind() == Kind::Infinity) {
//     return undefined();
//   }

//   if (k.kind() == Kind::Integer) {
//     Int p = k.value();

//     return integer(mod(p, s, symmetric));
//   }

//   if (k.kind() == Kind::Symbol) {
//     return k;
//   }

//   if (k.kind() == Kind::Fraction) {
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
//     Expr p = Expr(Kind::Derivative, {gf(k[0], x, s, symmetric), k[1]});

//     return p;
//   }

//   if (k.kind() == Kind::Integral) {
//     Expr p = Expr(Kind::Integral, {gf(k[0], x, s, symmetric), k[1]});

//     return p;
//   }

//   if (k.kind() == Kind::Factorial) {
//     if (k[0].kind() == Kind::Integer) {
//       Expr f = reduceAST(k);
//       return gf(f, x, s, symmetric);
//     }

//     return k;
//   }

//   if (k.kind() == Kind::Division) {
//     Expr p = div(gf(k[0], x, s, symmetric), gf(k[1], x, s, symmetric));

//     return gf(reduceAST(p), x, s, symmetric);
//   }

//   if (k.kind() == Kind::Power) {
//     Expr p = power(gf(k[0], x, s, symmetric), k[1]);

//     return p;
//   }

//   if (k.kind() == Kind::FunctionCall) {
//     return k;
//   }

//   if (k.kind() == Kind::Multiplication) {
//     Expr p = Expr(Kind::Multiplication);

//     for (long i = 0; i < k.size(); i++) {
//       p.insert(gf(k[i], x, s, symmetric));
//     }

//     return p;
//   }

//   if (k.kind() == Kind::Addition || k.kind() == Kind::Subtraction) {
//     Expr p = Expr(k.kind());

//     Expr d = degree(k, x);

//     for (Int n = 0; n <= d.value(); n++) {
//       Expr c = coeff(k, x, n);
//       p.insert(gf(c, x, s, symmetric) * power(x, n));
//     }

//     return reduceAST(p);
//   }

//   return k;
// }

Expr groundGf(Expr u, Int s, bool symmetric) {
  Expr k = u;

  if (k.kind() == Kind::Fail || k.kind() == Kind::Undefined) {
    return k;
  }

  if (k.kind() == Kind::MinusInfinity || k.kind() == Kind::Infinity) {
    return undefined();
  }

  if (k.kind() == Kind::Integer) {
    return integer(mod(k.value(), s, symmetric));
  }

  if (k.kind() == Kind::Symbol) {
    return k;
  }

  if (k.kind() == Kind::Fraction) {
    assert(k[0].kind() == Kind::Integer,
           "numerator of a fraction needs to be a integer");

    assert(k[1].kind() == Kind::Integer,
           "denominator of a fraction needs to be a integer");

    Int n = k[0].value();
    Int d = k[1].value();

    return mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric);
  }

  if (k.kind() == Kind::Derivative) {
    Expr p = Expr(Kind::Derivative, {groundGf(k[0], s, symmetric), k[1]});

    return p;
  }

  if (k.kind() == Kind::Integral) {
    return Expr(Kind::Integral, {groundGf(k[0], s, symmetric), k[1]});
  }

  if (k.kind() == Kind::Factorial) {
    if (k[0].kind() == Kind::Integer) {
      return groundGf(fact(k.value()), s, symmetric);
    }

    return k;
  }

  if (k.kind() == Kind::Division) {
    Expr p = div(groundGf(k[0], s, symmetric), groundGf(k[1], s, symmetric));
    Expr t = reduceAST(p);

    return groundGf(t, s, symmetric);
  }

  if (k.kind() == Kind::Power) {
    return power(groundGf(k[0], s, symmetric), k[1]);
  }

  if (k.kind() == Kind::FunctionCall) {
    return k;
  }

  if (k.kind() == Kind::Multiplication) {
    Expr p = Expr(Kind::Multiplication);

    for (long i = 0; i < k.size(); i++) {
      p.insert(groundGf(k[i], s, symmetric));
    }

    return reduceAST(p);
  }

  if (k.kind() == Kind::Addition || k.kind() == Kind::Subtraction) {
    Expr p = Expr(k.kind());

    for (long i = 0; i < k.size(); i++) {
      p.insert(groundGf(k[i], s, symmetric));
    }

    Expr r = reduceAST(p);

    return r;
  }

  return k;
}

Expr gf(Expr u, Int s, bool symmetric) {
  Expr k = algebraicExpand(u);

  if (k.kind() == Kind::Fail || k.kind() == Kind::Undefined) {
    return k;
  }

  if (k.kind() == Kind::MinusInfinity || k.kind() == Kind::Infinity) {
    return undefined();
  }

  if (k.kind() == Kind::Integer) {
    Int p = k.value();

    return integer(mod(p, s, symmetric));
  }

  if (k.kind() == Kind::Symbol) {
    return k;
  }

  if (k.kind() == Kind::Fraction) {
    assert(k[0].kind() == Kind::Integer,
           "numerator of a fraction needs to be a integer");

    assert(k[1].kind() == Kind::Integer,
           "denominator of a fraction needs to be a integer");

    Int n = k[0].value();
    Int d = k[1].value();

    return integer(
        mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric));
  }

  if (k.kind() == Kind::Derivative) {
    Expr p = Expr(Kind::Derivative, {gf(k[0], s, symmetric), k[1]});

    return p;
  }

  if (k.kind() == Kind::Integral) {
    Expr p = Expr(Kind::Integral, {gf(k[0], s, symmetric), k[1]});

    return p;
  }

  if (k.kind() == Kind::Factorial) {
    if (k[0].kind() == Kind::Integer) {
      Expr f = reduceAST(k);

      Expr p = gf(f, s, symmetric);

      return p;
    }

    return k;
  }

  if (k.kind() == Kind::Division) {
    Expr p = div(gf(k[0], s, symmetric), gf(k[1], s, symmetric));
    Expr t = reduceAST(p);
    Expr r = gf(t, s, symmetric);

    return r;
  }

  if (k.kind() == Kind::Power) {
    Expr p = power(gf(k[0], s, symmetric), k[1]);

    return p;
  }

  if (k.kind() == Kind::FunctionCall) {
    return k;
  }

  if (k.kind() == Kind::Multiplication) {
    Expr p = Expr(Kind::Multiplication);

    for (long i = 0; i < k.size(); i++) {
      p.insert(gf(k[i], s, symmetric));
    }

    return p;
  }

  if (k.kind() == Kind::Addition || k.kind() == Kind::Subtraction) {
    Expr p = Expr(k.kind());

    for (long i = 0; i < k.size(); i++) {
      p.insert(gf(k[i], s, symmetric));
    }

    Expr r = reduceAST(p);

    return r;
  }

  return k;
}

Expr divPolyGf(Expr a, Expr b, Expr x, Int p, bool symmetric) {

	Expr da = degree(a, x);
  Expr db = degree(b, x);

	if(db == -inf()) {
		// TODO: throw division by zero
		return undefined();
	}

  if (da == -inf() || da.value() < db.value()) {
    return list({0, a});
  }

  long long k, j, d;
  Int s, e, lb;

  Expr dq, dr, q, r;
  Expr t1, t2, t3, ex;

  std::vector<Expr> A = std::vector<Expr>(da.value().longValue() + 1);
  std::vector<Expr> B = std::vector<Expr>(db.value().longValue() + 1);

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

    t1 = reduceAST(t1);
		// printf("%s\n", p.to_string().c_str());
		// printf("%s\n", t1.toString().c_str());
    t1 = mod(t1.value(), p, symmetric);

		//printf("%s\n", t1.toString().c_str());
    // if (t1.value() < 0) {
    //   t1 = mod(t1.value() + p, p, symmetric);
    // }

    if (da.value() - k <= dq.value()) {
      t1 = reduceAST(mulPoly(t1, lb));
			//printf("%s\n", t1.toString().c_str());
			//printf("%s\n", lb.to_string().c_str());
			t1 = mod(t1.value(), p, false);
    }

    A[k] = t1;
  }

  q = add({0});
  r = add({0});

  d = 0;

  for (k = da.value().longValue() - dq.value().longValue(); k <= da.value();
       k++) {
    q.insert(A[k] * power(x, integer(d)));
    d = d + 1;
  }

  d = 0;

  for (k = 0; k <= dr.value(); k++) {
    r.insert(A[k] * power(x, integer(d)));

    d = d + 1;
  }
  t1 = gf(q, p, symmetric);
  t2 = gf(r, p, symmetric);

  return list({t1, t2});
}

Expr remPolyGf(Expr a, Expr b, Expr x, Int p, bool symmetric) {
  return divPolyGf(a, b, x, p, symmetric)[1];
}

Expr quoPolyGf(Expr a, Expr b, Expr x, Int p, bool symmetric) {
  return divPolyGf(a, b, x, p, symmetric)[0];
}

Expr monicPolyGf(Expr f, Expr x, Int p, bool symmetric) {
  if (f == 0) {
    return 0;
  }

  Expr lc = leadCoeff(f, x);

  Expr F = quoPolyGf(f, lc, x, p, symmetric);

  return list({lc, F});
}

Expr gcdPolyGf(Expr a, Expr b, Expr x, Int p, bool symmetric) {
  Expr da = degree(a, x);
  Expr db = degree(b, x);

  if (da.kind() == Kind::MinusInfinity || db.value() > da.value()) {
    return gcdPolyGf(b, a, x, p, symmetric);
  }

  Expr t;
  while (b != 0 && db.kind() != Kind::MinusInfinity && db.value() >= 0) {
		printf("a = %s\n", a.toString().c_str());
		printf("b = %s\n", b.toString().c_str());

		t = a;
    a = b;

		Expr t1 = divPolyGf(t, b, x, p, symmetric);
		printf("q = %s\n", t1.toString().c_str());
		Expr t2 = mulPolyGf(b, t1[0], x, p, symmetric);
		Expr t3 = addPolyGf(t2, t1[1], x, p, symmetric);
		printf("t = %s\n", t3.toString().c_str());
		//assert(t == t3, "");
		b = remPolyGf(t, b, x, p, symmetric);


		db = degree(b, x);
  }

	printf("A = %s\n", a.toString().c_str());
  b = monicPolyGf(a, x, p, symmetric);
	printf("b = %s\n", b.toString().c_str());
  return b[1];
}

Expr addPolyGf(Expr f, Expr g, Expr x, Int p, bool symmetric) {
  Expr u = addPoly(f, g);

  return gf(u, p, symmetric);
}

Expr subPolyGf(Expr f, Expr g, Expr x, Int p, bool symmetric) {
  Expr t, u;

  u = subPoly(f, g);

  t = gf(u, p, symmetric);

  return t;
}

Expr mulPolyGf(Expr f, Expr g, Expr x, Int p, bool symmetric) {
  Expr t, u;

  u = mulPoly(f, g);

  t = gf(u, p, symmetric);

  return t;
}

Expr powModPolyGf(Expr f, Expr g, Expr x, Int n, Int p, bool symmetric) {
  if (n == 0)
    return 1;

  Expr a = f;
  Expr b = 1;
  Expr t = 0;

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

Expr randPolyGf(Int d, Expr x, Int p, bool symmetric) {
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

  Expr r = Expr(Kind::Addition, {power(x, d)});

  for (Int i = d - 1; i >= 2; i--) {
    k = randomGf(p, symmetric);

    if (k != 0) {
      r = r + k * power(x, i);
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

Expr extendedEuclidGf(Expr f, Expr g, Expr x, Int p, bool sym) {
  if (f == 0 || g == 0) {
    return list({1, 0, 0});
  }

  Expr t, s, i, lc, k1, t0, t3, s0, s1, Q, R, T;

  Expr t1 = monicPolyGf(f, x, p, sym);
  Expr t2 = monicPolyGf(g, x, p, sym);

  Expr p0 = t1[0];
  Expr r0 = t1[1];

  Expr p1 = t2[0];
  Expr r1 = t2[1];

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

Expr gfPolyExpr(Expr u, Int p, bool symmetric) {
  if (u.kind() == Kind::Integer) {
    return mod(u.value(), p, symmetric);
  }

  if (u.kind() == Kind::Multiplication) {
    assert(u.size() == 2, "not a polynomial expr");
    return gfPolyExpr(u[0], p, symmetric) * u[1];
  }

  assert(u.kind() == Kind::Addition, "not a polynomial expr");

  Expr g = Expr(Kind::Addition);
	Expr x = 0;

  for (Int i = 0; i < u.size(); i++) {
    assert(u[i].kind() == Kind::Multiplication && u[i].size() == 2,
           "not a polynomial expr");

    Expr c = gfPolyExpr(u[i][0], p, symmetric);

		x = u[i][1][0];

		if (!isZeroPolyExpr(c)) {
      g.insert(c * u[i][1]);
    }
	}

  if (g.size() == 0) {
		if(x == 0) return 0;
		g.insert(0*power(x, 0));
	}

  return g;
}

Expr addPolyExprGf(Expr f, Expr g, Int p, bool sym) {
  return gfPolyExpr(addPolyExpr(f, g), p, sym);
}

Expr subPolyExprGf(Expr f, Expr g, Int p, bool sym) {
  return gfPolyExpr(subPolyExpr(f, g), p, sym);
}

Expr mulPolyExprGf(Expr f, Expr g, Int p, bool sym) {
  return gfPolyExpr(mulPolyExpr(f, g), p, sym);
}

Expr divPolyExprGf(Expr a, Expr b, Expr L, Int p, bool symmetric) {
  assert(a.kind() == Kind::Addition, "not a poly expr");
  assert(b.kind() == Kind::Addition, "not a poly expr");
  assert(L.kind() == Kind::List && L.size() == 1, "not a univariate poly expr");

  Expr x = L[0];

  Expr da = degreePolyExpr(a);
  Expr db = degreePolyExpr(b);

  assert(da.kind() == Kind::Integer,
         "degree of polynomial should be an integer\n");
  assert(db.kind() == Kind::Integer,
         "degree of polynomial should be an integer\n");

  if (da.value() < db.value()) {
    return list({0, a});
  }

  long long k, j;
  Int s, e;

  Int dq, dr;

  Int t1, t2, t3, lb, d;

  std::vector<Int> A = std::vector<Int>(da.value().longValue() + 1, 0);
  std::vector<Int> B = std::vector<Int>(db.value().longValue() + 1, 0);

  for (k = 0; k < a.size(); k++) {
    assert(a[k].kind() == Kind::Multiplication && a[k].size() == 2,
           "not a poly expr");

    assert(a[k][0].kind() == Kind::Integer, "not a univariate poly expr");

    Expr d = a[k][1][1];

    assert(d.kind() == Kind::Integer,
           "poly expr should have only integers as degrees");

    assert(a[k][0].kind() == Kind::Integer,
           "poly expr should have only integers as degrees");

    A[d.value().longValue()] = a[k][0].value();
  }

  for (k = 0; k < b.size(); k++) {
    assert(b[k].kind() == Kind::Multiplication && b[k].size() == 2,
           "not a poly expr");

    assert(b[k][0].kind() == Kind::Integer, "not a univariate poly expr");

    Expr d = b[k][1][1];

    assert(d.kind() == Kind::Integer,
           "poly expr should have only integers as degrees");
    assert(b[k][0].kind() == Kind::Integer,
           "poly expr should have only integers as degrees");
    B[d.value().longValue()] = b[k][0].value();
  }

  dq = da.value() - db.value();
  dr = db.value() - 1;

  t1 = leadCoeffPolyExpr(b).value();

  lb = inverseGf(t1, p, symmetric);

  for (k = da.value().longValue(); k >= 0; k--) {

    t1 = A[k];
    s = max(0, k - dq);
    e = min(dr, k);

    for (j = s.longValue(); j <= e; j++) {
      t2 = B[j] * A[k - j + db.value().longValue()];
      t3 = t1 - t2;
      t1 = t3;
    }

    t1 = mod(t1, p, symmetric);

    if (t1 < 0) {
      t1 = mod(t1 + p, p, symmetric);
    }

    if (da.value() - k <= dq) {
      t1 = mod(t1 * lb, p, symmetric);
    }

    A[k] = t1;
  }

  Expr q = add({});
  Expr r = add({});

  d = 0;

  for (k = da.value().longValue() - dq.longValue(); k <= da.value(); k++) {
    if (A[k] != 0) {
      q.insert(A[k] * power(x, d));
    }

    d = d + 1;
  }

  d = 0;

  for (k = 0; k <= dr; k++) {
    if (A[k] != 0) {
      r.insert(A[k] * power(x, d));
    }

    d = d + 1;
  }

	if(q.size() == 0) {
		q.insert(0*power(x, 0));
	}

	if(r.size() == 0) {
		r.insert(0*power(x, 0));
	}

  return list({q, r});
}


Expr quoPolyExprGf(Expr a, Expr b, Expr L, Int p, bool symmetric) {
  return divPolyExprGf(a, b, L, p, symmetric)[0];
}

Expr remPolyExprGf(Expr a, Expr b, Expr L, Int p, bool symmetric) {
  return divPolyExprGf(a, b, L, p, symmetric)[1];
}

Expr monicPolyExprGf(Expr f, Expr L, Int p, bool symmetric) {
	Expr x = L[0];

  if (isZeroPolyExpr(f)) {
    return Expr(Kind::Addition, {0 * power(x, 0)});
  }

  Expr lc = Expr(Kind::Addition, {leadCoeffPolyExpr(f) * power(x, 0)});

  Expr F = quoPolyExprGf(f, lc, L, p, symmetric);

  return list({lc, F});
}

Expr randPolyExprGf(Int d, Expr L, Int p, bool symmetric) {
  Expr r = Expr(Kind::Addition, {});
	Expr x = L[0];

  for (Int i = 0; i < d; i++) {
    Int k = randomGf(p, symmetric);

    if (k != 0) {
      r = r + k * power(x, i);
    }
  }

  r = r + 1 * power(x, d);

  return r;
}

Expr powModPolyExprGf(Expr f, Expr g, Expr L, Int n, Int p, bool symmetric) {
  Expr b = Expr(Kind::Addition, {1 * power(L[0], 0)});

  if (n == 0)
    return b;

  Expr a = f;
  Expr t = 0;

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

Expr gcdPolyExprGf(Expr a, Expr b, Expr L, Int p, bool symmetric) {
  Expr da = degreePolyExpr(a);
  Expr db = degreePolyExpr(b);

  if (da.kind() == Kind::MinusInfinity || db.value() > da.value()) {
    return gcdPolyExprGf(b, a, L, p, symmetric);
  }

  Expr t;

  while (!isZeroPolyExpr(b) && db.kind() != Kind::MinusInfinity &&
         db.value() >= 0) {
    t = a;
    a = b;

    b = remPolyExprGf(t, b, L, p, symmetric);
    db = degreePolyExpr(b);
  }

  b = monicPolyExprGf(a, L, p, symmetric);

  return b[1];
}

Expr extendedEuclidPolyExprGf(Expr f, Expr g, Expr L, Int p, bool sym) {
	Expr x = L[0];

	if (f == 0 || g == 0) {
    return list({
				raisePolyExpr(1, 0, x),
				raisePolyExpr(0, 0, x),
				raisePolyExpr(0, 0, x),
			});
  }

  Expr t, s, i, lc, k1, t0, t3, s0, s1, Q, R, T;

  Expr t1 = monicPolyExprGf(f, L, p, sym);
  Expr t2 = monicPolyExprGf(g, L, p, sym);
  Expr p0 = t1[0][0][0];
  Expr r0 = t1[1];

  Expr p1 = t2[0][0][0];
  Expr r1 = t2[1];

  if (f == 0) {
    return list({
				raisePolyExpr(0, 0, x),
				raisePolyExpr(inverseGf(p1.value(), p, sym), 0, x),
				r1,
			});
  }

  if (g == 0) {
    return list({
				raisePolyExpr(inverseGf(p0.value(), p, sym), 0, x),
				raisePolyExpr(0, 0, x),
				r0
			});
  }

  s0 = raisePolyExpr(inverseGf(p0.value(), p, sym), 0, x);

  s1 = raisePolyExpr(0, 0, x);
	t0 = raisePolyExpr(0, 0, x);

	t1 = raisePolyExpr(inverseGf(p1.value(), p, sym), 0, x);

  while (true) {
		T = divPolyExprGf(r0, r1, L, p, sym);

    Q = T[0];
    R = T[1];

    if (isZeroPolyExpr(R)) {
      break;
    }
    T = monicPolyExprGf(R, L, p, sym);

    r0 = r1;

    lc = T[0][0][0];
    r1 = T[1];

    assert(lc.kind() == Kind::Integer, "lc of univariate should be a integer");

    i = raisePolyExpr(inverseGf(lc.value(), p, sym), 0, x);

    k1 = mulPolyExprGf(s1, Q, p, sym);
		s = subPolyExprGf(s0, k1, p, sym);

    k1 = mulPolyExprGf(t1, Q, p, sym);
    t = subPolyExprGf(t0, k1, p, sym);

    s0 = s1;
    t0 = t1;
		s1 = mulPolyExprGf(s, i, p, sym);
    t1 = mulPolyExprGf(t, i, p, sym);

	}

	return list({r1, s1, t1});
}

} // namespace galoisField
