// reference: A Formalization of the Berlekamp–Zassenhaus Factorization
// Algorithm

#include "Berlekamp.hpp"

#include "MathSystem/Polynomial/Polynomial.hpp"
#include "MathSystem/GaloisField/GaloisField.hpp"
#include <cstddef>

using namespace alg;
using namespace polynomial;
using namespace galoisField;

namespace factorization {

void swapRows(expr& M, Int j, Int t) {
  for (Int i = 0; i < M.size(); i++) {
    expr Mji = M[j][i];
    expr Mti = M[t][i];

		M[j][i] = Mti;
		M[t][i] = Mji;
  }
}

void addFreeVariableToBase(expr& v, Int n, long var_idx) {
  v.insert(list({}));

  for (long k = 0; k < n; k++) {
    if (k == var_idx) {
      v[v.size() - 1].insert(integer(1));
    } else {
      v[v.size() - 1].insert(integer(0));
    }
  }
}

expr buildBerlekampBasis(expr A, expr w, bool symmetric) {
  expr M;
  Int lead, row_count, col_count, n, x, q;

  q = w.value();

  n = A.size();

  M = list({});
  // M = (A - I)'
  for (size_t i = 0; i < n; i++) {
    M.insert(list({}));

    for (size_t j = 0; j < n; j++) {
      if (i == j) {
        M[i].insert(integer(mod(A[j][i].value() - 1, q, symmetric)), j);
      } else {
        M[i].insert(integer(A[j][i].value()), j);
      }
    }
  }

  lead = 0;

  row_count = n;
  col_count = n;

  for (size_t r = 0; r < row_count; r++) {
    if (col_count <= lead) {
      break;
    }

    size_t i = r;

    while (M[i][lead] == 0) {
      i = i + 1;

      if (row_count == i) {
        i = r;
        lead = lead + 1;

        if (col_count == lead) {
          break;
        }
      }
    }
    if (col_count == lead) {
      break;
    }

    swapRows(M, i, r);

    if (M[r][lead] != (0)) {
      Int x = M[r][lead].value();

      for (size_t j = 0; j < n; j++) {
        Int v = M[r][j].value();
        M[r][j] = quoGf(v, x, q, symmetric);
      }
    }

    for (i = 0; i < row_count; i++) {
      if (i != r) {
        Int x = M[i][lead].value();

        for (size_t j = 0; j < n; j++) {
          Int v = M[r][j].value();
					Int t = M[i][j].value();

          M[i][j] = mod(t - mod(x * v, q, symmetric), q, symmetric);
        }
      }
    }
  }

  expr v = list({});

  size_t k = 0;

  for (size_t i = 0; i < n; i++) {
    addFreeVariableToBase(v, n, -1);
  }

  k = 0;

  for (size_t i = 0; i < n; i++) {
    while (k < n && M[i][k] == 0) {
			v[k][k] = 1;
			// v[k].remove(k);
      // v[k].insert(integer(1), k);
      k++;
    }

    if (k < n) {
      for (size_t j = 0; j < n; j++) {
        if (j != k) {
          v[j][k] = mod(-1 * M[i][j].value(), q, symmetric);
          // v[j].remove(k);
          // v[j].insert(integer(x), k);
        }
      }
    }

    k = k + 1;
  }

  for (size_t i = 0; i < v.size(); i++) {
    bool rem = true;

    for (size_t j = 0; j < n && rem; j++) {
      if (v[i][j] != 0) {
        rem = false;
      }
    }

    if (rem) {
      v.remove(i);
      i = i - 1;
    }
  }

  return v;
}

expr initBerkelampBasisMatrix(expr n) {
  expr Q = list({});
  Q.insert(list({}));

  Q[0].insert(1, 0);

  for (long i = 1; i < n.value(); i++) {
    Q[0].insert(0);
  }

  return Q;
}

void addVecToBasisMatrix(expr& Q, expr r, expr x, long i, Int n) {
  expr ex, ri;

  Q.insert(list({}), i);

	for (long k = 0; k < n; k++) {
    Q[i].insert(coeff(r, x, k));
  }

}

expr buildBerkelampMatrix(expr ax, expr x, expr p, bool symmetric) {

	expr n = degree(ax, x);
  expr Q = initBerkelampBasisMatrix(n);

  expr r0 = powModPolyGf(x, ax, x, p.value(), p.value(), symmetric);

  addVecToBasisMatrix(Q, r0, x, 1, n.value());

  expr rx = r0;

  expr zx = 0;

  for (long m = 2; m < n.value(); m++) {
		zx = mulPolyGf(r0, rx, x, p.value(), symmetric);
    rx = remPolyGf(zx, ax, x, p.value(), symmetric);

		addVecToBasisMatrix(Q, rx, x, m, n.value());
  }

  return Q;
}

expr buildBerlekampBasisPolynomials(expr B, expr x, expr n) {
  expr basis = list({});

  // Build polynomial basis ignoring the '1' basis
  for (size_t i = 1; i < B.size(); i++) {
    expr bx = 0;

    for (size_t j = 0; j < B[i].size(); j++) {
      if (B[i][j] != 0) {
        bx = bx + (B[i][j] * pow(x, integer(n.value() - j - 1)));
      }
    }

    basis.insert(reduce(bx));
  }

  return basis;
}

expr berlekampFactors(expr sfx, expr x, expr p, bool symmetric) {
  expr Q, B, lc, fx, n, H, F, h, v, f;

  lc = leadCoeff(sfx, x);

  fx = quoPolyGf(sfx, lc, x, p.value(), symmetric);

  n = degree(fx, x);

  Q = buildBerkelampMatrix(fx, x, p, symmetric);
  B = buildBerlekampBasis(Q, p, symmetric);

  if (B.size() == 1) {
    return sfx;
  }

  size_t r = B.size();

  H = buildBerlekampBasisPolynomials(B, x, n);

  size_t i = 0;

  F = list({fx});

  f = F[0];

  for (size_t k = 0; k < r; k++) {
    if (F.size() == r || f == 1) {
      break;
    }

    f = F[i];
    h = H[k];

    for (size_t s = 0; s < p.value(); s++) {
      expr g = integer(s);

      v = subPolyGf(h, g, x, p.value(), symmetric);

      g = gcdPolyGf(f, v, x, p.value(), symmetric);

      if (g != 1 && g != f) {
        v = quoPolyGf(f, g, x, p.value(), symmetric);

        F.remove(i);

        F.insert(g);
        F.insert(v);

        i = F.size() - 1;

        f = F[i];
      }
    }
  }

  F.insert(lc, 0);

  return F;
}

} // namespace factorization
