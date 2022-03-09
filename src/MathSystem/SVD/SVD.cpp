#include "SVD.hpp"
#include "MathSystem/Algebra/Matrix.hpp"
#include <cmath>
#include <tuple>

using namespace alg;

double norm(double *v, long long n) {
  double norm = 0;
  for (long long i = 0; i < n; i++)
    norm += v[i] * v[i];
  return sqrt(norm);
}

// Given a real n-vector x, the Householder matrix
// P = I - u*u'/tal can be constructed so that P is
// orthogonal and P*x = += |x| * [1, 0 ... 0]
// The function outputs z = Py given x
void householderOne(double *x, double *y, double *z, long long n) {
  long long i;
  double sxx = 0, suy = 0, s = 0, tal = 0, d = 0;

  double *u = (double *)malloc(sizeof(double) * n);

  for (i = 0; i < n; i++)
    sxx += x[i] * x[i];

  if (sxx == 0) {
    for (i = 0; i < n; i++)
      z[i] = y[i];
  } else {
    s = sqrt(sxx);
    u[0] = x[0] + s;

    for (i = 1; i < n; i++)
      u[i] = x[i];

    tal = sxx + s * x[0];

    for (i = 0; i < n; i++)
      suy = suy + u[i] * y[i];

    d = suy / tal;

    for (i = 0; i < n; i++)
      z[i] = y[i] - u[i] * d;
  }

  delete u;
}

// Given a real n-vector x, the Householder matrix
// P = I - u*u'/tal can be constructed so that P is
// orthogonal and P*x = += |x| * [1, 0 ... 0]
// The function outputs z = Py given x
void householderTwo(double *x, double *y, double *z, long long n) {
  long long i;

  double sxx = 0, sxy = 0, s = 0, d = 0;

  double *u = (double *)malloc(sizeof(double) * n);

  for (i = 0; i < n; i++)
    sxx += x[i] * x[i];

  if (sxx == 0) {
    for (i = 0; i < n; i++)
      z[i] = y[i];
  } else {
    s = sqrt(sxx);

    u[0] = x[0] + s;

    for (i = 0; i < n; i++)
      sxy = sxy + x[i] * y[i];

    z[0] = -sxy / s;

    d = (y[0] - z[0]) / u[0];

    for (i = 1; i < n; i++)
      z[i] = y[i] - x[i] * d;
  }

  delete u;
}

// Given x[0:n] this function computes a = alpha
// and the vector u[0:n] such that (I - u*u')x = alpha [1,0...0]
void housegen(double *x, double *u, double *a, long long n) {
  long long i;

  double p, s, v;

  for (i = 0; i < n; i++)
    u[i] = x[i];

  v = norm(x, n);

  if (v == 0) {
    u[0] = sqrt(2);
    *a = v;
  } else {
    if (u[0] >= 0.0) {
      p = u[0] / fabs(u[0]);
    } else {
      p = 1.0;
    }

    for (i = 0; i < n; i++) {
      u[i] = (p / v) * u[i];
    }

    u[0] = 1 + u[0];

    s = sqrt(u[0]);

    for (i = 0; i < n; i++) {
      u[i] = u[i] / s;
    }

    *a = -p * v;
  }
}

// Calculate X[:, k:n]*H given that H*z = -alhpa*e1
double applyHouseholderRight(matrix &X, double *z, long long k, long long n,
                             long long m, double *Vk) {
  long long i, j;

  double alpha;

  double *u = (double *)malloc(sizeof(double) * n - k);
  double *b = (double *)malloc(sizeof(double) * m);

  housegen(z, u, &alpha, n - k);

  // A*H = A(I - u*u') = A - (A*u)*u

  // We now have (I - u*u') = V

  // Vk = (I - u*u')'
  for (i = 0; i < n - k; i++) {
    for (j = 0; j < n - k; j++) {
      if (i == j)
        Vk[i * (n - k) + j] = 1.;
      else
        Vk[i * (n - k) + j] = 0.;
    }
  }

  for (i = 0; i < n - k; i++)
    for (j = 0; j < n - k; j++)
      Vk[i * (n - k) + j] = Vk[i * (n - k) + j] - u[i] * u[j];

  // we need X[:][k:n]*V'

  // given that V is hermitian,
  // X[:][k:n]*V'
  //	= X[:][k:n] * (I - u*u')
  //	= X[:][k:n] - (X[:][k:n]*u)*u'

  // X[:][k:n]*u
  for (i = 0; i < m; i++) {
    b[i] = 0;
    for (j = k; j < n; j++) {
      b[i] += X[i][j] * u[j - k];
    }
  }

  // X[:][k:n] - b*u'
  for (i = 0; i < m; i++) {
    for (j = k; j < n; j++) {
      X[i][j] = X[i][j] - b[i] * u[j - k];

      // Inverted results here
      // X[i][j] = -X[i][j];
    }
  }

  delete[] u;
  delete[] b;

  return alpha;
}

// arrays e and q are gonna hold the values of the diagonals
// e is the upper diagonal and q the lower.
void golubReinschHouseholderBidiagonalization(matrix &a, double *q, double *e,
                                              matrix &u, matrix &v,
                                              double tol = 2.22e-16,
                                              double *x_ = nullptr) {
  long long i, j, k, l;
  double f, g, h, s, x, y;

  // m >= n assumed
  long long m = a.lines();
  long long n = a.columns();

  // Copy matrix a into u
  // u = matrix(a);
  u = matrix(m, m);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      u[i][j] = a[i][j];

  v = matrix(n, n);

  g = 0;
  x = 0;

  // Householder reduction t o bidiagonal form
  for (i = 0; i < n; i++) {

    e[i] = g;

    s = 0;
    l = i + 1;

    for (j = i; j < m; j++)
      s = s + (u[j][i] * u[j][i]);

    if (s < tol)
      g = 0;
    else {
      f = u[i][i];
      g = f < 0 ? sqrt(s) : -sqrt(s);

      h = f * g - s;

      u[i][i] = f - g;

      for (j = l; j < n; j++) {
        s = 0;

        for (k = i; k < m; k++)
          s = s + (u[k][i] * u[k][j]);

        f = s / h;

        for (k = i; k < m; k++)
          u[k][j] = u[k][j] + (f * u[k][i]);
      }
    }

    q[i] = g;
    s = 0;

    for (j = l; j < n; j++)
      s = s + (u[i][j] * u[i][j]);

    if (s < tol)
      g = 0;
    else {

      f = u[i][i + 1];
      g = f < 0 ? sqrt(s) : -sqrt(s);

      h = f * g - s;
      u[i][i + 1] = f - g;

      for (j = l; j < n; j++)
        e[j] = u[i][j] / h;

      for (j = l; j < m; j++) {
        s = 0;

        for (k = l; k < n; k++)
          s = s + (u[j][k] * u[i][k]);

        for (k = l; k < n; k++)
          u[j][k] = u[j][k] + (s * e[k]);
      }
    }
    y = fabs(q[i]) + fabs(e[i]);

    if (y > x)
      x = y;
  }

  // accumulation of right-hand transformations.
  for (i = n - 1; i >= 0; i--) {
    if (g != 0) {
      h = u[i][i + 1] * g;

      for (j = l; j < n; j++)
        v[j][i] = u[i][j] / h;

      for (j = l; j < n; j++) {
        s = 0;

        for (k = l; k < n; k++)
          s = s + (u[i][k] * v[k][j]);

        for (k = l; k < n; k++)
          v[k][j] = v[k][j] + (s * v[k][i]);
      }
    }

    for (j = l; j < n; j++) {
      v[i][j] = 0;
      v[j][i] = 0;
    }

    v[i][i] = 1;

    g = e[i];

    l = i;
  }
  // accumulation of the left-hand transformations.
  for (i = n; i < m; i++) {
    for (j = n; j > m; j++)
      u[i][j] = 0;

    u[i][i] = 1;
  }

  for (i = n - 1; i >= 0; i--) {
    l = i + 1;
    g = q[i];

    for (j = l; j < m; j++)
      u[i][j] = 0;

    if (g != 0) {
      h = u[i][i] * g;
      for (j = l; j < m; j++) {
        s = 0;
        for (k = l; k < m; k++)
          s = s + (u[k][i] * u[k][j]);

        f = s / h;

        for (k = i; k < m; k++)
          u[k][j] = u[k][j] + (f * u[k][i]);
      }

      for (j = i; j < m; j++)
        u[j][i] = u[j][i] / g;
    } else {
      for (j = i; j < m; j++)
        u[j][i] = 0;
    }

    u[i][i] = u[i][i] + 1;
  }

  if (x_)
    *x_ = x;
}

double house(double *x, double *v, long long n, double tol) {
  long long i;
  double mu, beta, sigma, v0;

  // std::cout << "x:\n";
  // for(i=0; i<n; i++)
  // 	std::cout << x[i] << " ";
  // std::cout << "\n";

  sigma = 0.0;

  for (i = 1; i < n; i++) {
    sigma = sigma + (x[i] * x[i]);
  }

  v[0] = 1;

  for (i = 1; i < n; i++) {
    v[i] = x[i];
  }

  if (fabs(sigma) <= tol) {
    beta = 0;
  } else {
    mu = sqrt((x[0] * x[0]) + sigma);

    if (x[0] <= 0) {
      v[0] = x[0] - mu;
    } else {
      v[0] = -sigma / (x[0] + mu);
    }

    beta = (2 * (v[0] * v[0])) / (sigma + (v[0] * v[0]));

    v0 = v[0];
    // std::cout << "house\n";
    // for(i = 0; i < n; i++)
    // 	std::cout << v[i] << " ";
    // std::cout << "\n";

    for (i = 0; i < n; i++) {
      v[i] = v[i] / v0;
    }
  }

  return beta;
}

// returns Px where P = (I - u*u')
void applyHouseholderToVector(double *u, double beta, double *x, double *Px,
                              long long n) {
  // P*x = (I * beta*u*u')*x = x - beta*u(u'*x)
  long long i;
  double ux = 0;

  // u'*x
  for (i = 0; i < n; i++) {
    ux = ux + u[i] * x[i];
  }

  // beta * u * (u'*x)
  for (i = 0; i < n; i++) {
    Px[i] = beta * u[i] * ux;
  }

  // P*x = x - u(u'*x)
  for (i = 0; i < n; i++) {
    Px[i] = x[i] - Px[i];
  }
}

// Computes [v, beta] = house(A[j:m, p]) such that (I - v*v')A[j:m, p] = beta*e1
double houseCol(matrix &A, long long m, long long j, long long p, double *v,
                double tol) {
  long long i;
  double beta;

  double *x = new double[m - j]; // (double*)malloc(sizeof(double)*(m-j));

  for (i = 0; i < m - j; i++) {
    x[i] = A[j + i][p];
  }

  beta = house(x, v, m - j, tol);

  delete[] x;

  return beta;
}

// Computes [v, beta] = house(A[p, j:n]') such that A[j, j:n]*(I - v*v')=
// beta*e1
double houseRow(matrix &A, long long n, long long p, long long j, double *v,
                double tol) {
  long long i;

  double beta;

  double *x = new double[n - j];

  for (i = 0; i < n - j; i++)
    x[i] = 0;

  for (i = 0; i < n - j; i++) {
    x[i] = A[p][i + j];
  }

  beta = house(x, v, n - j, tol);

  delete[] x;

  return beta;
}

/**
 * @brief Compute P*A where P = (I - beta*u[0:m - p]*u'[0:m - p])
 * @obs: The essential part of u is inserted where the zeros where introduced in
 * column k.
 *
 * This function calculates
 * 				= (I - beta * u[0:m - p] * u[0:m - p]') *
 * A[p:m][k:n] = A[p:m][k:n] - beta * u[1:m - p] * (u[1:m - p]' * A[p:m][k:n])
 */
double preHouseholdermatrix(double *u, double beta, matrix &A, unsigned m,
                            unsigned n, unsigned p, unsigned k) {
  double *uA;
  long long i, j;

  uA = new double[n - k]; // (double*)malloc(sizeof(double)*(n - k));

  // uA[1:n - k] = u[1:m - p]' * A[p:m][k:n]
  for (i = 0; i < n - k; i++) {
    uA[i] = 0;
    for (j = 0; j < m - p; j++) {
      uA[i] = uA[i] + (u[j] * A[j + p][i + k]);
    }
  }

  // A[p:m][k:n] - u[1:m - p] * (beta * u[1:m - p]' * A[p:m][k:n])
  //		= A[p:m][k:n] - beta * u[1:m - p] * uA[1:n - k]
  for (i = 0; i < m - p; i++) {
    for (j = 0; j < n - k; j++) {
      A[p + i][k + j] = A[p + i][k + j] - (beta * u[i] * uA[j]);
    }
  }

  // Store esential part of u where new zeros where introduced.
  for (i = 0; i < m - (p + 1); i++) {
    A[i + (p + 1)][k] = u[i + 1];
  }

  delete[] uA;

  return A[p][k];
}

/**
 * @brief Compute A*P where P = (I - beta*u[0:n - k]*u'[0:n - k])
 * @obs: The essential part of u is inserted where the zeros where introduced in
 * column k.
 *
 * This function calculates
 * 				= A[k:m][p:n] * (I - u[1:m - k]*(u[1:n - p]')
 * 				= A[k:m][p:n] - (A[k:m][p:n] * u[1:n - p]) *
 * (beta
 * * u[1:n - p]')
 */
double posHouseholdermatrix(double *u, double beta, matrix &A, unsigned m,
                            unsigned n, unsigned k, unsigned p) {
  double *Au;
  long long i, j;

  Au = new double[m - k]; // (double*)malloc(sizeof(double)*(m - k));

  // Au = A[k:m][p:n] * u[1:n - p]
  for (i = 0; i < m - k; i++) {
    Au[i] = 0;

    for (j = 0; j < n - p; j++) {
      Au[i] = Au[i] + (A[i + k][j + p] * u[j]);
    }
  }

  // A[k:m][p:n] - Au[1:m - k] * (beta*u[1:n-p]).T
  for (i = 0; i < m - k; i++) {
    for (j = 0; j < n - p; j++) {
      A[k + i][p + j] = A[k + i][p + j] - (Au[i] * (beta * u[j]));
    }
  }

  // Store esential part of u where new zeros where introduced.
  for (i = 0; i < n - (p + 1); i++) {
    A[k][i + p + 1] = u[i + 1];
  }

  delete[] Au;

  return A[k][p];
}

void householderBidiagonalization(matrix &A, double *diag, double *sdiag,
                                  matrix &u, matrix &v, double tol) {
  long long m, n, i, j, k, q;

  double *w, *wQ, *beta, *gamma;

  m = A.lines();
  n = A.columns();

  w = new double[m];
  wQ = new double[m];

  beta = new double[n];
  gamma = new double[n];

  gamma[0] = 0;

  u = identity(m, m);
  v = identity(n, n);

  for (j = 0; j < n; j++) {
    beta[j] = houseCol(A, m, j, j, w, tol);

    diag[j] = preHouseholdermatrix(w, beta[j], A, m, n, j, j);

    if (j < n - 1) {
      gamma[j + 1] = houseRow(A, n, j, j + 1, w, tol);
      sdiag[j + 1] = posHouseholdermatrix(w, gamma[j + 1], A, m, n, j, j + 1);
    }
  }

  // Accumulate left transformations
  k = m;

  for (j = n - 1; j >= 0; j--) {

    w[j] = 1;
    for (i = 1; i < m - j; i++) {
      w[j + i] = A[j + i][j];
    }

    // w[j:m]' * u[j:m][j:k]
    for (q = 0; q < k - j; q++) {
      wQ[j + q] = 0;
      for (i = 0; i < m - j; i++) {
        wQ[j + q] = wQ[j + q] + (w[j + i] * u[j + i][j + q]);
      }
    }

    for (q = 0; q < k - j; q++) {
      for (i = 0; i < m - j; i++) {
        u[j + q][j + i] = u[j + q][j + i] - (beta[j] * w[j + q] * wQ[j + i]);
      }
    }
  }

  // Accumulate right transformations
  k = n;
  for (j = n - 2; j >= 0; j--) {
    w[j] = 1;
    for (i = 1; i < n - (j + 1); i++) {
      w[j + i] = A[j][j + i + 1];
    }

    for (q = 0; q < k - (j + 1); q++) {
      wQ[j + q] = 0;
      for (i = 0; i < n - (j + 1); i++) {
        wQ[j + q] = wQ[j + q] + (w[j + i] * v[j + i + 1][j + q + 1]);
      }
    }

    for (q = 0; q < k - (j + 1); q++) {
      for (i = 0; i < n - (j + 1); i++) {
        v[j + q + 1][j + i + 1] =
            v[j + q + 1][j + i + 1] - (gamma[j + 1] * w[j + q] * wQ[j + i]);
      }
    }
  }

  delete[] beta;
  delete[] gamma;
  delete[] wQ;
  delete[] w;
}

long long sign(double v) {
  if (v > 0)
    return 1;
  if (v == 0)
    return 0;
  return -1;
}

void sortUV(double *s, unsigned n, matrix &u, matrix &v) {
  long long i, k, rows, j, i_last;
  double s_last, tmp;

  for (i = 0; i < n; ++i) {
    s_last = s[i];
    i_last = i;

    for (j = i + 1; j < n; ++j) {
      if (fabs(s[j]) > fabs(s_last)) {
        s_last = s[j];
        i_last = j;
      }
    }

    if (i_last != i) {
      tmp = s[i];

      s[i] = s[i_last];
      s[i_last] = tmp;

      rows = u.lines();

      for (k = 0; k < rows; ++k) {
        tmp = u[k][i];
        u[k][i] = u[k][i_last];
        u[k][i_last] = tmp;
      }

      rows = v.lines();

      for (k = 0; k < rows; ++k) {
        tmp = v[i][k];
        v[i][k] = v[i_last][k];
        v[i_last][k] = tmp;
      }
    }
  }
}

// GvL pg. 216 : algo 5.1.3
void givens(double a, double b, double *c, double *s) {
  // Computes scalars c and s such that
  //   [c, s; -s, c].T * [a, b] = [r, 0]

  if (b == 0) {
    *c = 1.0;
    *s = 0.0;
  } else {
    double r = hypot(a, b);
    *c = a / r;
    *s = -b / r;
  }
}

// GvL pg 216 section 5.1.9
// Computes [c, s; -s,c].T * A[:, [i, k]]
void leftGivens(matrix &A, double c, double s, long long i, long long k) {
  // Apply [c, s; -s, c].T to A[[i,k], :]
  long long j, n;
  double t1, t2;

  n = A.columns();

  for (j = 0; j < n; j++) {
    t1 = A[i][j];
    t2 = A[k][j];

    A[i][j] = c * t1 - s * t2;
    A[k][j] = s * t1 + c * t2;

    // if(fabs(A[i][j]) <= tol)
    // 	A[i][j] = 0.0;

    // if(fabs(A[k][j]) <= tol)
    // 	A[k][j] = 0.0;
  }
}

// GvL pg 216 section 5.1.9
// Computes A[:, [i, k]] * [c, s; -s,c]
void rightGivens(matrix &A, double c, double s, long long i, long long k) {
  // Apply [c, s; -s, c].T to A[[i,k], :]
  long long j, m;
  double t1, t2;

  m = A.lines();

  for (j = 0; j < m; j++) {
    t1 = A[j][i];
    t2 = A[j][k];

    A[j][i] = c * t1 - s * t2;
    A[j][k] = s * t1 + c * t2;
  }
}

// delete extra diagonals in diags
// pointer diags[2] and diags[3]
// stays allocated
void freeDiagmatrix(double **diags) {
  diags[1] = nullptr;
  diags[2] = nullptr;

  delete[] diags[0];
  delete[] diags[3];
  delete[] diags[4];
  delete[] diags;
}

double trailing2x2Eigenvalue(matrix &B, long long q) {
  double a, b, c, d, mu, t11, t21, t22;
  // 1. Find y and z
  //	 Let u be the eigenvalue of the trailing 2x2 submatrix
  // 	 of T = B'*B that is closer to t_nn

  //	T[0:2, 0:2] = B[:,-2:].T*B[:,-2:]
  //  B[:, -2] = {{0,..., 0, a, c}, { 0,..., 0, b, d}}.T

  if (q >= 3) {
    a = B[q - 3][q - 2];
    b = B[q - 2][q - 1];
    c = B[q - 2][q - 2];
    d = B[q - 1][q - 1];
  } else {
    a = B[q - 2][q - 2];
    b = B[q - 2][q - 1];
    c = B[q - 1][q - 2];
    d = B[q - 1][q - 1];
  }

  t11 = a * a + c * c;
  t21 = b * c;
  t22 = b * b + d * d;

  // compute wilkinson shift value "µ" ("mu")
  d = (t11 - t22) / 2.0;
  mu = t22 - (sign(d) * (t21 * t21)) / (fabs(d) + sqrt((d * d) + (t21 * t21)));

  return mu;
}

void golubKahanStep(matrix &B, long long p, long long q, matrix &uT,
                    matrix &v) {
  long long k;
  double y, z, c, s, mu, t11, t12;

  t11 = B[p][p] * B[p][p];
  t12 = B[p][p] * B[p][p + 1];
  mu = trailing2x2Eigenvalue(B, q);
  y = t11 - mu;
  z = t12;

  givens(y, z, &c, &s);
  rightGivens(B, c, s, p, p + 1);
  leftGivens(v, c, s, p, p + 1);

  for (k = p; k < q - 2; k++) {

    y = B[k][k];
    z = B[k + 1][k];

    givens(y, z, &c, &s);
    leftGivens(B, c, s, k, k + 1);
    rightGivens(uT, c, s, k, k + 1);

    y = B[k][k + 1];
    z = B[k][k + 2];

    givens(y, z, &c, &s);
    rightGivens(B, c, s, k + 1, k + 2);
    leftGivens(v, c, s, k + 1, k + 2);
  }

  y = B[q - 2][q - 2];
  z = B[q - 1][q - 2];
  givens(y, z, &c, &s);
  leftGivens(B, c, s, q - 2, q - 1);
  rightGivens(uT, c, s, q - 2, q - 1);
}

void clean(matrix &B, double tol) {
  long long i, n;
  n = B.columns();
  for (i = 0; i < n; i++) {
    if (fabs(B[i][i]) < tol) {
      B[i][i] = 0.0;
    }
  }
  for (i = 0; i < n - 1; i++) {
    if (fabs(B[i][i + 1]) < tol * (fabs(B[i][i]) + fabs(B[i + 1][i + 1]))) {
      B[i][i + 1] = 0.0;
    }
  }
}

void findLimits(matrix &B, long long *p, long long *q) {
  long long i, n;

  n = B.columns();

  *q = -1;
  *p = -1;

  for (i = n - 1; i >= 1; i--) {
    // Pick largest q such that B33 is diagonal

    if (*q == -1 && B[i - 1][i]) {
      *p = *q = i - 1;
    } else if (*q != -1) {
      if (B[i - 1][i]) {
        *p = i - 1;
      } else {
        break;
      }
    }
  }

  if (*q == -1) {
    *p = 0;
    *q = 0;
  } else {
    *p = *p;
    *q = *q + 1;
  }
}

double getDiagMat(double **diags, long long l, long long c) {
  if ((l + 2 - c >= 0) && (l + 2 - c < 5)) {
    return diags[l + 2 - c][c];
  }

  return 0;
}

void walkBlemishOutRight(matrix &B, long long q, long long row,
                         long long start_col, matrix &uT) {
  matrix B_ = matrix(B);

  long long col, n;
  double c, s, old;

  n = B.columns();

  for (col = start_col; col < q; col++) {
    givens(B[col][col], B[row][col], &c, &s);
    old = B[row][col];

    B[row][col] = 0;
    B[col][col] = B[col][col] * c - old * s;

    if (col < n - 1) {
      B[row][col + 1] = s * B[col][col + 1];
      B[col][col + 1] = c * B[col][col + 1];
    }

    rightGivens(uT, c, s, col, row);
  }
}

// start row is kk-1 and col kk sup diag od k
void walkBlemishOutUp(matrix &B, long long p, long long start_row,
                      long long col, matrix &v) {
  long long row;
  double c, s, old;

  for (row = start_row; row >= p; row--) {
    givens(B[row][row], B[row][col], &c, &s);

    old = B[row][col];

    B[row][col] = 0;

    B[row][row] = B[row][row] * c - old * s;

    if (row > 0) {
      B[row - 1][col] = s * B[row - 1][row];

      B[row - 1][row] = c * B[row - 1][row];
    }

    leftGivens(v, c, s, row, col);
  }
}

// GvL pg 454
bool doZeroDiag(matrix &B, long long p, long long q, matrix &uT, matrix &v,
                double tol) {
  long long i = 0;

  bool zeroed = false;

  // Get zeros idx in diag
  for (i = p; i < q; i++) {

    if (fabs(B[i][i]) <= tol) {
      zeroed = true;

      if (i < q - 1) {
        walkBlemishOutRight(B, q, i, i + 1, uT);
      } else if (i == q - 1) {
        walkBlemishOutUp(B, p, i - 1, i, v);
      }
    }
  }

  return zeroed;
}

void golubKahanSVD(double *b_diag, double *b_sdiag, long long m, long long n,
                   matrix &uT, matrix &v, double tol) {
  long long i, q = 0, p = 0;

  matrix B = bidiag(b_diag, b_sdiag, m, n);

  v = transpose(v);

  clean(B, tol);

  findLimits(B, &p, &q);

  while (q - p) {
    bool zeroed = doZeroDiag(B, p, q + 1, uT, v, tol);
    if (!zeroed) {
      golubKahanStep(B, p, q + 1, uT, v);
    }

    clean(B, tol);
    findLimits(B, &p, &q);
  }

  for (i = 0; i < n; i++) {
    b_diag[i] = B[i][i];
  }

  sortUV(b_diag, n, uT, v);
}

std::tuple<matrix, matrix, matrix> alg::svd(matrix a, double tol) {
  matrix A;

  if (a.columns() > a.lines()) {
    A = transpose(a);
  } else {
    A = matrix(a);
  }

  unsigned m = A.lines();
  unsigned n = A.columns();

  matrix u = identity(m, m);
  matrix v = identity(n, n);

  double *phi = new double[n];

  phi[0] = 0;

  double *s = new double[n];

  householderBidiagonalization(A, s, phi, u, v, tol);

  golubKahanSVD(s, phi, m, n, u, v, tol);

  delete[] phi;

  matrix d = identity(m, n);

  for (unsigned i = 0; i < n; i++) {
    d[i][i] = s[i];
  }

  delete[] s;

  if (a.columns() != a.lines() && A.lines() == a.columns() &&
      A.columns() == a.lines()) {
    matrix tmp = u;
    u = transpose(v);
    v = transpose(tmp);
    d = transpose(d);
  }

  return std::make_tuple(u, d, v);
}
