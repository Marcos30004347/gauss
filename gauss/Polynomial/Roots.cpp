
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stack>
#include <vector>

#include "Polynomial.hpp"
#include "Roots.hpp"
#include "gauss/Algebra/Reduction.hpp"

using namespace alg;
using namespace poly;

// Internal Complex number representation.
class complex {
public:
  double real;
  double imag;

  complex(double r, double i);
  complex(double r);
  complex();

  double abs();

  complex conj() const;

  complex operator*(double other);
  complex operator+(double other);
  complex operator-(double other);
  complex operator/(double other);

  complex operator*(complex &other);
  complex operator+(complex &other);
  complex operator-(complex &other);
  complex operator/(complex &other);

  complex operator*(const complex &other);
  complex operator+(const complex &other);
  complex operator-(const complex &other);
  complex operator/(const complex &other);
};

// Internal RealPolynomialRepresentation
class RealPoly {
public:
  double *p_coefficients;

  long p_power;

  RealPoly();
  RealPoly(const RealPoly &other);
  RealPoly(unsigned long long power);
  RealPoly(unsigned long long power, double *coefs);
  RealPoly(unsigned long long power, std::initializer_list<double> coefs);
  ~RealPoly();

  long power() const;

  RealPoly &operator=(const RealPoly &other);
  RealPoly operator/(const RealPoly &other);
  RealPoly operator/(const double other);
  RealPoly operator+(const RealPoly &other);
  RealPoly operator*(const RealPoly &other);
  RealPoly operator*(const double other);
  RealPoly operator-(const RealPoly &other);

  bool operator==(const RealPoly &other);

  RealPoly cauchy();
  RealPoly normalized();
  RealPoly dx();
  RealPoly abs();
  RealPoly reverse();
  RealPoly noLeadingTerm();
  RealPoly noConstantTerm();

  double &operator[](unsigned long long i) const;
  double eval(double x);
  complex eval(complex x);
  RealPoly root(std::vector<complex> &roots, double tol = 0.00001);
  // void roots(std::vector<complex>& roots, double tol = 0.00001);
};

double newtonRaphson(RealPoly &p, double tol = 0.0001);

RealPoly addPoly(const RealPoly &A, const RealPoly &B);
RealPoly subPoly(const RealPoly &A, const RealPoly &B);
RealPoly mulPoly(const RealPoly &A, const double B);
RealPoly mulPoly(const RealPoly &A, const RealPoly &B);
std::vector<complex> polyRoots(RealPoly p, double tol = std::numeric_limits<double>::epsilon());

void divPoly(const RealPoly p, const RealPoly &d, RealPoly &q, RealPoly &r);

bool quadraShift(std::vector<complex> &xs, RealPoly &p, RealPoly &K, complex &x,
                 unsigned long long &M, unsigned long long &L, bool &qflag,
                 bool &lflag, double tol);
bool linearShift(std::vector<complex> &xs, RealPoly &p, RealPoly &K, complex &x,
                 unsigned long long M, unsigned long long L, bool &qflag,
                 bool &lflag, double tol);

complex::complex(double r, double i) {
  this->real = r;
  this->imag = i;
}

complex::complex(double r) {
  this->real = r;
  this->imag = 0;
}

complex::complex() {
  this->real = 0;
  this->imag = 0;
}

complex complex::conj() const { return complex(this->real, -1 * this->imag); }

complex complex::operator*(complex &other) {
  // 1 = i^2
  return complex((this->real * other.real) + (-1 * this->imag * other.imag),
                 (this->real * other.imag) + (this->imag * other.real));
}

complex complex::operator+(complex &other) {
  return complex(this->real + other.real, this->imag + other.imag);
}

complex complex::operator-(complex &other) {
  return complex(this->real - other.real, this->imag - other.imag);
}

complex complex::operator/(complex &other) {
  complex a = *this;

  complex b = other;

  complex n = a * other.conj();
  complex d = b * other.conj();

  return complex(n.real / d.real, n.imag / d.real);
}

complex complex::operator*(const complex &other) {
  // 1 = i^2
  return complex((this->real * other.real) + (-1 * this->imag * other.imag),
                 (this->real * other.imag) + (this->imag * other.real));
}

complex complex::operator+(const complex &other) {
  return complex(this->real + other.real, this->imag + other.imag);
}

complex complex::operator-(const complex &other) {
  return complex(this->real - other.real, this->imag - other.imag);
}

complex complex::operator/(const complex &other) {
  complex a = *this;

  complex b = other;

  complex n = a * other.conj();
  complex d = b * other.conj();

  return complex(n.real / d.real, n.imag / d.real);
}

double complex::abs() { return std::hypot(this->real, this->imag); }

complex complex::operator*(double other) {
  return complex(this->real * other, this->imag * other);
}

complex complex::operator+(double other) {
  return complex(this->real + other, this->imag);
}

complex complex::operator-(double other) {
  return complex(this->real - other, this->imag);
}

complex complex::operator/(double other) {
  return complex(this->real / other, this->imag / other);
}

void divPoly(const RealPoly p, const RealPoly &d, RealPoly &q, RealPoly &r) {
  long i, j;
  long power = p.power() - d.power();
  double ratio;

  if (power < 0)
    return;

  q = RealPoly(power, {});

  r = RealPoly(p);

  for (i = p.power(); i >= d.power(); i--) {
    q[i - d.power()] = ratio = r[i] / d[d.p_power];

    // std::cout << "quotient: " << r[i]<< std::endl;
    // std::cout << "quotient: " << d[d.p_power] << std::endl;

    r[i] = 0;

    for (j = 0; j < d.power(); j++) {
      r[i - d.power() + j] -= d[j] * ratio;
    }
  }
  while (r.p_power >= 0 && r[--r.p_power] == 0)
    ;
}

void divPoly(const RealPoly p, const double d, RealPoly &q) {
  q = RealPoly(p);
  for (int j = 0; j <= p.power(); j++)
    q[j] = p[j] / d;
}

RealPoly addPoly(const RealPoly &A, const RealPoly &B) {
  RealPoly X = A.power() > B.power() ? RealPoly(A) : RealPoly(B);
  RealPoly Y = A.power() > B.power() ? RealPoly(B) : RealPoly(A);

  for (long i = 0; i <= Y.power(); i++) {
    X[i] += Y[i];
  }

  while (X[X.p_power] == 0) {
    X.p_power--;
  }
  return X;
}

RealPoly subPoly(const RealPoly &A, const RealPoly &B) {
  RealPoly B_ = mulPoly(B, -1.0);

  RealPoly X = A.power() > B.power() ? RealPoly(A) : RealPoly(B_);
  RealPoly Y = A.power() > B.power() ? RealPoly(B_) : RealPoly(A);

  for (long i = 0; i <= Y.power(); i++) {
    X[i] += Y[i];
  }

  while (X[X.p_power] == 0) {
    X.p_power--;
  }

  return X;
}

RealPoly mulPoly(const RealPoly &A, const RealPoly &B) {
  unsigned long long m = A.power() + 1;
  unsigned long long n = B.power() + 1;

  RealPoly AB(m + n - 2, {});

  for (unsigned long long i = 0; i < m; i++) {
    for (unsigned long long j = 0; j < n; j++) {
      AB[i + j] += A[i] * B[j];
    }
  }
  return AB;
}

RealPoly mulPoly(const RealPoly &A, double B) {
  RealPoly AB(A);
  for (int j = 0; j <= AB.power(); j++)
    AB[j] = A[j] * B;
  return AB;
}

RealPoly::RealPoly() {
  this->p_power = 0;
  this->p_coefficients = new double[this->p_power + 1];
  this->p_coefficients[0] = 0.0;
}

RealPoly::RealPoly(const RealPoly &other) {
  this->p_power = other.p_power;
  this->p_coefficients = new double[this->p_power + 1];

  std::copy(other.p_coefficients, other.p_coefficients + other.p_power + 1,
            this->p_coefficients);
}

RealPoly::RealPoly(unsigned long long power,
                   std::initializer_list<double> coefs) {
  this->p_power = power;
  this->p_coefficients = new double[++power];

  if (coefs.size() == 0) {
    std::fill(&this->p_coefficients[0], &this->p_coefficients[power - 1] + 1,
              0);
  } else {
    std::copy(coefs.begin(), coefs.end(), this->p_coefficients);
  }
}

RealPoly::RealPoly(unsigned long long power) {
  this->p_power = power;
  this->p_coefficients = new double[++power];
  std::fill(&this->p_coefficients[0], &this->p_coefficients[power - 1] + 1, 0);
}

RealPoly::RealPoly(unsigned long long power, double *coefs) {
  this->p_power = power;
  this->p_coefficients = new double[++power];

  if (coefs == nullptr) {
    std::fill(&this->p_coefficients[0], &this->p_coefficients[power - 1] + 1,
              0);
  } else {
    std::copy(&coefs[0], &coefs[power - 1] + 1, this->p_coefficients);
  }
}

RealPoly::~RealPoly() { delete this->p_coefficients; }

long RealPoly::power() const { return this->p_power; }

RealPoly RealPoly::operator/(const RealPoly &other) {
  RealPoly q, r;
  divPoly(*this, other, q, r);
  return q;
}

RealPoly RealPoly::operator/(const double other) {
  RealPoly q;
  divPoly(*this, other, q);
  return q;
}

RealPoly RealPoly::operator*(const double other) {
  return mulPoly(*this, other);
}

RealPoly RealPoly::operator*(const RealPoly &other) {
  return mulPoly(*this, other);
}

RealPoly RealPoly::operator+(const RealPoly &other) {
  return addPoly(*this, other);
}

RealPoly RealPoly::operator-(const RealPoly &other) {
  return subPoly(*this, other);
}

RealPoly &RealPoly::operator=(const RealPoly &other) {
  this->p_power = other.p_power;
  this->p_coefficients = new double[this->p_power + 1];
  std::copy(other.p_coefficients, other.p_coefficients + other.p_power + 1,
            this->p_coefficients);
  return *this;
}

bool RealPoly::operator==(const RealPoly &other) {
  if (other.p_power != this->p_power)
    return false;

  for (long i = 0; i <= this->p_power; i++) {
    if (this->p_coefficients[i] != other.p_coefficients[i])
      return false;
  }

  return true;
}

double &RealPoly::operator[](unsigned long long i) const {
  return this->p_coefficients[i];
}

// horners method
double RealPoly::eval(double x) {
  double px = this->p_coefficients[this->power()];

  for (long i = this->power() - 1; i >= 0; i--) {
    px = px * x + this->p_coefficients[i];
  }

  return px;
}

double syntheticDivision(RealPoly &p, double x, RealPoly &q) {
  q = RealPoly(p.power() - 1, {});

  q[p.power() - 1] = p[p.power()];

  for (int i = p.power() - 2; i >= 0; i--) {
    q[i] = p[i + 1] + q[i + 1] * x;
  }

  return p[0] + q[0] * x;
}

complex pow(complex c, unsigned i) {
  complex p = c;

  for (unsigned j = 0; j < i; j++) {
    p = p * c;
  }

  return p;
}

// horners method
complex RealPoly::eval(complex x) {
  complex px(0, 0);

  for (long i = this->power(); i >= 0; i--) {
    px = px * x + this->p_coefficients[i];
  }

  return px;
}

RealPoly RealPoly::dx() {
  RealPoly d(this->power() - 1, {});

  for (long i = 1; i <= this->power(); i++) {
    d[i - 1] = this->p_coefficients[i] * i;
  }

  return d;
}

RealPoly RealPoly::abs() {
  RealPoly d(this->power(), {});

  for (long i = 0; i <= this->power(); i++) {
    d[i] = this->p_coefficients[i] >= 0 ? this->p_coefficients[i]
                                        : -this->p_coefficients[i];
  }

  return d;
}

RealPoly RealPoly::reverse() {
  RealPoly d(this->power(), {});

  for (long i = 0; i <= this->power(); i++) {
    d[this->power() - i] = this->p_coefficients[i];
  }

  return d;
}

RealPoly RealPoly::noLeadingTerm() {
  RealPoly d(this->power() - 1, {});

  for (long i = 0; i < this->power(); i++) {
    d[i] = this->p_coefficients[i];
  }

  return d;
}

RealPoly RealPoly::noConstantTerm() {
  RealPoly d(this->power(), {});

  for (long i = 0; i <= this->power(); i++) {
    d[i] = this->p_coefficients[i];
  }

  d[0] = 0;

  return d;
}

RealPoly RealPoly::normalized() {

  RealPoly n(this->power(), {});

  if (n[n.power()] == 1)
    return n;

  for (long i = 0; i <= this->power(); i++) {
    n[i] = this->p_coefficients[i] / this->p_coefficients[this->p_power];
  }

  return n;
}

RealPoly RealPoly::cauchy() {
  RealPoly p = this->normalized();

  for (long i = 0; i <= p.power(); i++) {
    p[i] = fabs(p[i]);
  }

  p[p.power()] *= -1;

  return p;
}

// double randomBetweenZeroOne()
// {
// 	std::mt19937_64 rng;
// 	uint64_t timeSeed =
// std::chrono::high_resolution_clock::now().time_since_epoch().count();
// 	std::seed_seq ss{uint32_t(timeSeed & 0xffffffff),
// uint32_t(timeSeed>>32)}; 	rng.seed(ss); 	std::uniform_real_distribution<double>
// unif(0, 1); 	return unif(rng);
// }

double newtonRaphson(RealPoly &p, double x0, unsigned max_it = 300) {
  RealPoly pdx = p.dx();

  double x = x0;

  double prev = 0;

  for (unsigned i = 0; i < max_it; i++) {
    prev = x;

    x = x - p.eval(x) / pdx.eval(x);

    if (fabs(prev - x) < std::numeric_limits<double>::epsilon()) {
      break;
    }
  }

  return x;
}

double rootRadius(RealPoly &p) {
  RealPoly p_ = p.abs();

  p_[0] *= -1;

  // printPoly(p_);

  return newtonRaphson(p_, 1.0, 100);
}

void jenkinsTraubPhaseOne(RealPoly &p, RealPoly &K, unsigned long long M) {
  RealPoly z = RealPoly(1, {0, 1});

  K = p.dx() / (p.power() + 1);

  for (unsigned long long i = 1; i < M; i++) {
    K = (K.noConstantTerm() + p.noConstantTerm() * (-K[0] / p[0])) / z;
  }
}

RealPoly nextSigma(RealPoly &p, RealPoly &sigma, RealPoly &K, double &a,
                   double &b, double &c, double &d) {
  double u = sigma[1];
  double v = sigma[0];

  double b1 = -K[0] / p[0];
  double b2 = -(K[1] + b1 * p[1]) / p[0];

  double a1 = b * c - a * d;
  double a2 = a * c + u * a * d + v * b * d;
  double c2 = b1 * a2;
  double c3 = b1 * b1 * (a * a + u * a * b + v * b * b);
  double c4 = v * b2 * a1 - c2 - c3;
  double c1 =
      c * c + u * c * d + v * d * d + b1 * (a * c + u * b * c + v * b * d) - c4;
  double du = -(u * (c2 + c3) + v * (b1 * a1 + b2 * a2)) / c1;
  double dv = v * c4 / c1;

  return RealPoly(2, {v + dv, u + du, 1});
}

RealPoly nextKQuaraticShift(RealPoly K, RealPoly &sigma, RealPoly &p_q,
                            RealPoly &k_q, double &a, double &b, double &c,
                            double &d) {
  double t = (a * a + sigma[1] * a * b + sigma[0] * b * b) / (b * c - a * d);

  // printPoly(sigma);

  RealPoly lin(1, {});

  lin[0] = -(a * c + sigma[1] * a * d + sigma[0] * b * d) / (b * c - a * d);
  lin[1] = 1;

  K = (k_q * t) + (lin * p_q);

  K[0] += b;

  return K;
}

int fixedShiftKPoly(RealPoly &p, RealPoly &sigma, RealPoly &K, complex &x,
                    double &a, double &b, double &c, double &d,
                    unsigned max_it) {
  RealPoly p_r, p_q, k_q, k_r;

  sigma[0] = x.real * x.real + x.imag * x.imag;
  sigma[1] = -2.0 * x.real;
  sigma[2] = 1;

  divPoly(p, sigma, p_q, p_r);

  b = p_r[p_r.power()];
  a = p_r[1] - b * sigma[1];

  complex p_x = x.conj() * -b + a;

  // std::cout << p_x.real << " + " << p_x.imag << "i\n";

  complex t_l[3] = {complex(0, 0), complex(0, 0), complex(0, 0)};
  double s_l[3] = {0, 0, 0};

  for (unsigned j = 0; j < max_it; j++) {
    K = K.normalized();

    divPoly(K, sigma, k_q, k_r);

    d = k_r[k_r.power()];
    c = k_r[k_r.power() - 1] - d * sigma[sigma.power() - 1];

    RealPoly s = nextSigma(p, sigma, K, a, b, c, d);

    complex k_x = x.conj() * -d + c;

    t_l[0] = t_l[1];
    t_l[1] = t_l[2];

    s_l[0] = s_l[1];
    s_l[1] = s_l[2];

    t_l[2] = x - p_x / k_x;
    s_l[2] = s[0];

    if (std::fabs(s_l[1] - s_l[0]) < std::fabs(s_l[0]) / 2.0 &&
        std::fabs(s_l[2] - s_l[1]) < std::fabs(s_l[1]) / 2.0)
      return 2;

    if ((t_l[1] - t_l[0]).abs() < (t_l[0]).abs() / 2.0 &&
        (t_l[2] - t_l[1]).abs() < (t_l[1]).abs() / 2.0)
      return 1;

    K = nextKQuaraticShift(K, sigma, p_q, k_q, a, b, c, d);
  }

  return 0;
}

int jenkinsTraubPhaseTwo(RealPoly &p, RealPoly &K, complex &x,
                         unsigned long long L) {
  unsigned i;

  double a, b, d, c;

  double deg2rad = M_PI / 180.0;
  double phi = 49.0 * deg2rad;

  double radius = rootRadius(p);

  RealPoly sigma(2, {});

  for (i = 0; i < L; i++) {
    x = complex(radius, 0) * complex(cos(phi), sin(phi));

    int j = fixedShiftKPoly(p, sigma, K, x, a, b, c, d, L * (i + 1));

    if (j != 0) {
      return j;
    }

    phi += 94.0 * deg2rad;
  }

  return 0;
}

bool hasConverged(std::vector<complex> &roots, double tol) {
  if (roots.size() != 3)
    return false;

  double e[2];

  e[1] = (roots[2] - roots[1]).abs();
  e[0] = (roots[1] - roots[0]).abs();

  double r = roots[1].abs();

  if (e[1] <= e[0]) {
    if (r <= tol) {
      return e[1] < tol;
    } else {
      return e[1] / r <= tol;
    }
  }

  return false;
}

bool hasConverged(std::vector<double> &roots, double tol) {
  if (roots.size() != 3)
    return false;

  double e[2];

  e[1] = fabs(roots[2] - roots[1]);
  e[0] = fabs(roots[1] - roots[0]);

  double r = fabs(roots[1]);

  if (e[1] <= e[0]) {
    if (r <= tol) {
      return e[1] < tol;
    } else {
      return e[1] / r <= tol;
    }
  }

  return false;
}
void findQuadraticRoots(RealPoly &p, complex &x1, complex &x2) {
  double a = p[2], b = p[1], c = p[0];

  double d = b * b - 4 * a * c;

  if (d >= 0) {
    if (b >= 0) {
      x1 = complex((-b - sqrt(fabs(d))) / (2.0 * a), 0);
      x2 = complex((2.0 * c) / (-b - sqrt(fabs(d))), 0);
    } else {
      x1 = complex((2.0 * c) / (-b + sqrt(fabs(d))), 0);
      x2 = complex((-b + sqrt(fabs(d))) / (2.0 * a), 0);
    }
  } else {
    x1 = complex(-b / (2.0 * a), sqrt(fabs(d)) / (2.0 * a));
    x2 = complex(-b / (2.0 * a), -sqrt(fabs(d)) / (2.0 * a));
  }
}

bool linearShift(std::vector<complex> &zeros, RealPoly &p, RealPoly &K,
                 complex &x, unsigned long long M, unsigned long long L,
                 bool &qflag, bool &lflag, double tol) {
  if (lflag) {
    return false;
  }

  std::vector<double> roots;

  RealPoly defl_p, defl_k;

  double p_x = 0, px = 0, kx = 0, d = 0;

  // s^(L) = Re(s1 - P(s1)/K^(L)(s1))
  double s = (x - p.eval(x) / K.eval(x)).real;

  roots.push_back(s);

  for (unsigned long long i = 0; i < L; i++) {
    if (hasConverged(roots, tol)) {
      zeros.push_back(complex(roots[1], 0));
      p = defl_p;
      return true;
    }

    p_x = px;

    px = syntheticDivision(p, s, defl_p);

    if (fabs(px) <= tol) {
      zeros.push_back(complex(s, 0));
      p = defl_p;

      return true;
    }

    kx = syntheticDivision(K, s, defl_k);

    K = defl_k + defl_p * -kx / px;

    K = K.normalized();

    kx = K.eval(s);

    d = px / kx;

    s = s - px / kx;

    roots.push_back(s);

    if (roots.size() > 3) {
      roots.erase(roots.begin());
    }

    if (i >= 2 && fabs(d) < 0.001 * fabs(s) && fabs(p_x) < fabs(p_x)) {
      complex root(s, 0);
      return quadraShift(zeros, p, K, root, M, L, qflag, lflag, tol);
    }
  }

  lflag = true;
  return quadraShift(zeros, p, K, x, M, L, qflag, lflag, tol);
}

bool quadraShift(std::vector<complex> &zeros, RealPoly &p, RealPoly &K,
                 complex &x, unsigned long long &M, unsigned long long &L,
                 bool &qflag, bool &lflag, double tol) {
  if (qflag) {
    return false;
  }

  complex x0, x1;

  double step = 0.01;
  double a, b, c, d;
  double px = 0, p_x = 0, tmp = 0;

  std::vector<complex> roots[2];

  RealPoly p_q, p_r, k_q, k_r, sigma(2, {});

  sigma[0] = x.real * x.real + x.imag * x.imag;
  sigma[1] = -2.0 * x.real;
  sigma[2] = 1;

  roots[0].push_back(x);
  roots[1].push_back(x.conj());

  bool fixed_shift = false;

  for (int i = 0; i < 20; i++) {
    if (hasConverged(roots[0], tol) && hasConverged(roots[1], tol)) {
      zeros.push_back(roots[0][1]);
      zeros.push_back(roots[1][1]);

      p = p_q;

      return true;
    }

    divPoly(p, sigma, p_q, p_r);

    b = p_r[p_r.power()];
    a = p_r[p_r.power() - 1] - b * sigma[sigma.power() - 1];

    findQuadraticRoots(sigma, x0, x1);

    if (fabs(fabs(x0.real) - fabs(x1.real)) > 0.01 * fabs(x1.real)) {
      return linearShift(zeros, p, K, x, M, L, qflag, lflag, tol);
    }

    px = fabs(a - x0.real * b) + fabs(x0.imag * b);

    if (!fixed_shift && fabs(sigma[0] - tmp) / sigma[0] < step && p_x > px) {
      fixed_shift = true;
      fixedShiftKPoly(p, sigma, K, x0, a, b, c, d, M);
    }

    divPoly(K, sigma, k_q, k_r);

    d = k_r[k_r.power()];
    c = k_r[k_r.power() - 1] - d * sigma[sigma.power() - 1];

    tmp = sigma[0];
    sigma = nextSigma(p, sigma, K, a, b, c, d);

    K = nextKQuaraticShift(K, sigma, p_q, k_q, a, b, c, d);
    K = K.normalized();
    p_x = px;

    roots[0].push_back(x0);
    roots[1].push_back(x1);

    if (roots[0].size() > 3) {
      roots[0].erase(roots[0].begin());
      roots[1].erase(roots[1].begin());
    }
  }

  qflag = true;

  return linearShift(zeros, p, K, x, M, L, qflag, lflag, tol);
}
bool jenkinsTraubPhaseThree(std::vector<complex> &roots, RealPoly &p,
                            RealPoly &K, complex &x, unsigned long long M,
                            unsigned long long L, unsigned cnv,
                            bool &quadratic_flag, bool &linear_flag,
                            double tol) {
  /*
   * Stage 3: Find the root with variable shift iterations on the K-polynomial.
   */

  // quadratic
  if (cnv == 2) {
    return quadraShift(roots, p, K, x, M, L, quadratic_flag, linear_flag, tol);
  }

  // linear
  if (cnv == 1) {
    return linearShift(roots, p, K, x, M, L, quadratic_flag, linear_flag, tol);
  }

  return false;
}

void findLinearRoots(RealPoly &p, complex &x) { x = complex(-p[0] / p[1], 0); }

int jenkinsTraub(std::vector<complex> &roots, RealPoly &p, double tol) {
  if (p.power() == 0) {
    return 0;
  }

  if (p.power() == 1) {
    complex x;
    findLinearRoots(p, x);
    roots.push_back(x);
    return 0;
  }

  if (p.power() == 2) {
    complex x1, x2;
    findQuadraticRoots(p, x1, x2);
    roots.push_back(x1);
    roots.push_back(x2);
    return 0;
  }

  RealPoly K;

  complex x(0.0, 0.0);

  unsigned long long M = 20;
  unsigned long long L = 100;

  jenkinsTraubPhaseOne(p, K, M);

  int cnv = jenkinsTraubPhaseTwo(p, K, x, L);

  if (cnv == 0) {
    return 1;
  }

  bool qflag = false;
  bool lflag = false;

  jenkinsTraubPhaseThree(roots, p, K, x, M, L, cnv, qflag, lflag, tol);

  return 0;
}

RealPoly RealPoly::root(std::vector<complex> &roots, double tol) {
  RealPoly p = this->normalized();

  if (jenkinsTraub(roots, p, tol) != 0) {
    // TODO: error
  }

  return p;
}

RealPoly removeZeroRoots(std::vector<complex> &roots, RealPoly &p) {
  int i = 0;

  while (p[i] == 0) {
    i++;
    roots.push_back(complex(0, 0));
  }

  RealPoly p_ = RealPoly(p.power() - i, {});

  for (int j = 0; j <= p.power() - i; j++) {
    p_[p.power() - j] = p[p.power() - j];
  }

  return p_;
}

std::vector<complex> polyRoots(RealPoly p, double tol) {
  std::vector<complex> roots;

  long n = p.power();

  while ((long long)roots.size() < n) {
    p = p.normalized();
    p = removeZeroRoots(roots, p);
    p = p.root(roots, tol);
  }

  return roots;
}

expr poly::realPolyRoots(expr P) {
  expr p = expand(P);

  expr L = getVariableListForPolyExpr(p);

  if (L.size() != 1) {
    return error("Polynomial needs to be univariate to have roots!");
  }

  if (L[0] == symbol("i")) {
    return error("Polynomial can't be imaginary!");
  }

  assert(is(&p, kind::ADD | kind::MUL));

  std::vector<double> coeffs;

  expr d = degree(p, L[0]);

  if (d.kind() != kind::INT) {
    return error("Polynomial can only have integer coefficients!");
  }

  expr t = d;

  for (Int i = d.value(); i >= 0; i--) {
    if (d.kind() != kind::INT) {
      return error("Polynomial can only have integer coefficients!");
    }

    if (i != d.value()) {
      coeffs.push_back(0);
      continue;
    }

    expr c = coeff(p, L[0], Int(i));

    if (!is(&c, kind::FRAC | kind::INT)) {
      return error("Polynomial is not a valid real polynomial!");
    }

    if (c.kind() == kind::FRAC) {
      double t = numerator(c).value().doubleValue() /
                 denominator(c).value().doubleValue();

      coeffs.push_back(t);
    }

    if (c.kind() == kind::INT) {
      double t = c.value().doubleValue();

      coeffs.push_back(t);
    }

    p = reduce(p - c * pow(L[0], d));

    d = degree(p, L[0]);
  }

  std::reverse(coeffs.begin(), coeffs.end());

  RealPoly r(coeffs.size() - 1, coeffs.data());

  std::vector<complex> roots = polyRoots(r);

  expr R = list({});

  for (size_t i = 0; i < roots.size(); i++) {
    R.insert(reduce(fromDouble(roots[i].real) + fromDouble(roots[i].imag) * symbol("i")));
  }

  return R;
}
