#ifndef INTEGER_H
#define INTEGER_H

#include "Int.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

struct Int {
  Int(bint<30> *v);
  bint<30> *val;
	friend class expression;

public:

  Int();

  Int(long int);
  Int(long long);
  Int(unsigned long long);
  Int(unsigned int);
  Int(int);
	Int(double);
  Int(const Int &);
  Int(Int &&);
  ~Int();

  inline Int operator+(const Int &other) const {
    bint<30> *tmp = new bint<30>();
    bint<30>::add(this->val, other.val, tmp);
    return Int(tmp);
  }

  inline Int operator+(const Int &&other) const {
    bint<30> *tmp = new bint<30>();
    bint<30>::add(this->val, other.val, tmp);
    return Int(tmp);
  }

	inline Int operator+(const int a) const {
		bint<30>* tmp = bint<30>::from(a);
		bint<30>* res = new bint<30>();
		bint<30>::add(this->val, tmp, res);
		return Int(res);
	}

  inline Int operator-(const Int &other) const {
    bint<30> *tmp = new bint<30>();
    bint<30>::sub(this->val, other.val, tmp);
    return Int(tmp);
  }

  inline Int operator-(const Int &&other) const {
    bint<30> *tmp = new bint<30>();
    bint<30>::sub(this->val, other.val, tmp);
		return Int(tmp);
  }
	inline Int operator-(const int a) const {
		bint<30>* tmp = bint<30>::from(a);
		bint<30>* res = new bint<30>();
		bint<30>::sub(this->val, tmp, res);
		return Int(res);
	}

  inline Int operator*(const Int &other) const {
    bint<30> *tmp = new bint<30>();
    bint<30>::mul(this->val, other.val, tmp);
    return Int(tmp);
  }

  inline Int operator*(const Int &&other) const {
    bint<30> *tmp = new bint<30>();
    bint<30>::mul(this->val, other.val, tmp);
    return Int(tmp);
  }
	inline Int operator*(const int a) const {
		bint<30>* tmp = bint<30>::from(a);
		bint<30>* res = new bint<30>();
		bint<30>::mul(this->val, tmp, res);
		return Int(res);
	}

  inline Int operator/(const Int &other) const {
    bint<30> *quo = new bint<30>();
    bint<30> *rem = new bint<30>();

    bint<30>::div(this->val, other.val, quo, rem);

    delete rem;

    return Int(quo);
  }

  inline Int operator/(const Int &&other) const {
    bint<30> *quo = new bint<30>();
    bint<30> *rem = new bint<30>();

    bint<30>::div(this->val, other.val, quo, rem);

    delete rem;

    return Int(quo);
  }

	inline Int operator/(const int a) const {
		bint<30>* tmp = bint<30>::from(a);
		bint<30>* quo = new bint<30>();
		bint<30>* rem = new bint<30>();
		bint<30>::div(this->val, tmp, quo, rem);

		delete rem;

		return Int(quo);
	}

  inline Int operator%(const Int &other) const {
    bint<30> *quo = new bint<30>();
    bint<30> *rem = new bint<30>();

    bint<30>::div(this->val, other.val, quo, rem);

    delete quo;

    return Int(rem);
  }

  inline Int operator%(const Int &&other) const {
    bint<30> *quo = new bint<30>();
    bint<30> *rem = new bint<30>();

    bint<30>::div(this->val, other.val, quo, rem);

    delete quo;

    return Int(rem);
  }
	inline Int operator%(const int a) const {
		bint<30>* tmp = bint<30>::from(a);
		bint<30>* quo = new bint<30>();
		bint<30>* rem = new bint<30>();
		bint<30>::div(this->val, tmp, quo, rem);

		delete quo;

		return Int(rem);
	}

  inline Int operator+() { return *this; }
  inline Int operator-() { return this->val->sign * -1; }

  inline void operator++() { *this = *this + 1; }
  inline void operator--() { *this = *this - 1; }

  inline void operator++(int) { *this = *this + 1; }
  inline void operator--(int) { *this = *this - 1; }

  inline bool operator==(const Int &other) const {
    return bint<30>::compare(this->val, other.val) == 0;
  }

  inline bool operator==(const Int &&other) const {
    return bint<30>::compare(this->val, other.val) == 0;
  }

  inline bool operator<(const Int &other) const {
    return bint<30>::compare(this->val, other.val) < 0;
  }
  inline bool operator<(const Int &&other) const {
    return bint<30>::compare(this->val, other.val) < 0;
  }
  inline bool operator<=(const Int &other) const {
    return bint<30>::compare(this->val, other.val) <= 0;
  }
  inline bool operator<=(const Int &&other)  const {
    return bint<30>::compare(this->val, other.val) <= 0;
  }
  inline bool operator>(const Int &other) const {
    return bint<30>::compare(this->val, other.val) > 0;
  }
  inline bool operator>(const Int &&other) const {
    return bint<30>::compare(this->val, other.val) > 0;
  }
  inline bool operator>=(const Int &other) const {
    return bint<30>::compare(this->val, other.val) >= 0;
  }
  inline bool operator>=(const Int &&other) const {
    return bint<30>::compare(this->val, other.val) >= 0;
  }
  inline bool operator!=(const Int &other) const {
    return bint<30>::compare(this->val, other.val) != 0;
  };
  inline bool operator!=(const Int &&other) const {
    return bint<30>::compare(this->val, other.val) != 0;
  };

  inline Int operator=(const Int &other) {
    this->val = other.val->copy();
    return *this;
  }

  inline Int operator=(Int &&other) {
    this->val = other.val;
    other.val = nullptr;
    return *this;
  }

  std::string to_string();

  inline long long longValue() {
		long long v = 0;
		// TODO: if to_long returned -1, that means an overflow
		bint<30>::to_long(this->val, &v);
		return v;
	}

  inline long long doubleValue() {
		double v = 0.0;
		// TODO: if to_long returned -1, that means an overflow
		bint<30>::to_double(this->val, &v);
		return v;
	}

  friend Int gcd(const Int &a, const Int &b) {
    return Int(bint<30>::gcd(a.val, b.val));
  }
  friend Int gcd(const Int &&a, const Int &&b) {
    return Int(bint<30>::gcd(a.val, b.val));
  }

  friend Int lcm(const Int &a, const Int &b) {
    return Int(bint<30>::lcm(a.val, b.val));
  }
  friend Int lcm(const Int &&a, const Int &&b) {
    return Int(bint<30>::lcm(a.val, b.val));
  }

  inline void operator/=(Int v) { *this = *this / v; }
  inline void operator+=(Int v) { *this = *this + v; }
  inline void operator-=(Int v) { *this = *this - v; }
  inline void operator*=(Int v) { *this = *this * v; }

  explicit operator bool() const { return this->val->size > 0; }

  friend Int abs(const Int &&a) { return Int(bint<30>::abs(a.val)); }
  friend Int abs(const Int &a) { return Int(bint<30>::abs(a.val)); }

  friend Int fact(const Int &&a) { return Int(bint<30>::fact(a.val)); }
  friend Int fact(const Int &a) { return Int(bint<30>::fact(a.val)); }

  friend Int max(const Int &&a, Int &&b) { return Int(bint<30>::max(a.val, b.val)); }
  friend Int max(const Int &a, Int &b) { return Int(bint<30>::max(a.val, b.val)); }
  friend Int max(const Int &a, Int &&b) { return Int(bint<30>::max(a.val, b.val)); }

  friend Int min(const Int &&a, const Int &&b) { return Int(bint<30>::min(a.val, b.val)); }
  friend Int min(const Int &a, const Int &b) { return Int(bint<30>::min(a.val, b.val)); }

  friend Int pow(const Int &&a, const Int &&b) { return Int(bint<30>::pow(a.val, b.val)); }
  friend Int pow(const Int &a, const Int &b) { return Int(bint<30>::pow(a.val, b.val)); }

  friend double pow(const Int &&a, const double b) { return bint<30>::pow(a.val, b); }
  friend double pow(const Int &a, const double b) { return bint<30>::pow(a.val, b); }

  friend bool operator<(const unsigned int &a, const Int &v) {
    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) < 0;
    delete tmp;
    return res;
  }

  friend Int operator*(const int &a, const Int &v)  {
    bint<30>* tmp = bint<30>::from(a);
    bint<30>* res = new bint<30>();
		bint<30>::mul(tmp, v.val, res);
    delete tmp;
    return res;
  }


  friend bool operator>(const unsigned int &a, const Int &v)  {
    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) > 0;
    delete tmp;
    return res;
  }


  friend bool operator<=(const unsigned int &a, const Int &v)  {
    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) <= 0;
    delete tmp;
    return res;
  }


  friend bool operator>=(const unsigned int &a, const Int &v)  {
    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) >= 0;
    delete tmp;
    return res;
  }


  friend Int operator+(const long long a, const Int &v) {
    bint<30>* tmp = bint<30>::from(a);
    bint<30>* res = new bint<30>();

		bint<30>::add(tmp, v.val, res);

		delete tmp;

		return res;
	}

  friend Int operator-(const long long a, const Int &v) {
    bint<30>* tmp = bint<30>::from(a);
    bint<30>* res = new bint<30>();

		bint<30>::sub(tmp, v.val, res);

		delete tmp;

		return res;
	}
};
/*
Int abs(Int& a);
Int abs(Int &&a);

Int fact(Int i);
Int pow(Int base, Int exp);
Int max(Int a, Int b);
Int min(Int a, Int b);
*/

/*
struct Int {
  std::vector<int> a;

  static const int base = 1000000000;
  static const int base_digits = 9;

  int sign;

  int size();

  Int operator^(const Int &v);

  std::string to_string();

  int sumof();

  Int();
  Int(long long v);
  Int(const std::string &s);

  void operator++();
  void operator--();

  void operator++(int);
  void operator--(int);

  // Int& operator--();

  Int operator+(const Int &v) const;
  Int operator+(const long long v) const;
  Int operator-(const long long v) const;
  Int operator-(const Int &v) const;
  Int operator-() const;
  Int operator/(const Int &v) const;
  Int operator%(const Int &v) const;
  Int operator*(const Int &v) const;
  Int operator/(int v) const;
  Int operator*(int v) const;
  Int operator*(long long v) const;
  void operator*=(int v);
  void operator*=(long long v);
  void operator/=(int v);
  int operator%(int v) const;
  void operator+=(const Int &v);
  void operator-=(const Int &v);
  void operator*=(const Int &v);
  void operator/=(const Int &v);
  bool operator<(const Int &v) const;
  bool operator>(const Int &v) const;
  bool operator<=(const Int &v) const;
  bool operator>=(const Int &v) const;
  bool operator==(const Int &v) const;
  bool operator!=(const Int &v) const;

  friend std::pair<Int, Int> divmod(const Int &a1, const Int &b1) {
    int norm = Int::base / (b1.a.back() + 1);

    Int a = a1.abs() * norm;
    Int b = b1.abs() * norm;

    Int q, r;

    q.a.resize(a.a.size());

    for (int i = a.a.size() - 1; i >= 0; i--) {
      r *= Int::base;

      r += a.a[i];

      int s1 = r.a.size() <= b.a.size() ? 0 : r.a[b.a.size()];
      int s2 = r.a.size() <= b.a.size() - 1 ? 0 : r.a[b.a.size() - 1];

      int d = ((long long)Int::base * s1 + s2) / b.a.back();

      r -= b * d;

      while (r < 0) {
        r += b, --d;
      }

      q.a[i] = d;
    }

    q.sign = a1.sign * b1.sign;
    r.sign = a1.sign;

    q.trim();
    r.trim();

    return std::make_pair(q, r / norm);
  }

  void trim();

  bool isZero() const;

  Int abs() const;

  long long longValue() const;

  friend Int gcd(const Int &a, const Int &b) {
    // if((a > 0 && b < 0) || (a < 0 && b > 0))
    // {
    // 	return gcd(a < 0 ? -a : a, b < 0 ? -b : b);
    // }

    return b.isZero() ? a : gcd(b, a % b);
  }

  friend Int lcm(const Int &a, const Int &b) { return a / gcd(a, b) * b; }

  void read(const std::string &s);

  friend std::istream &operator>>(std::istream &stream, Int &v) {
    std::string s;
    stream >> s;
    v.read(s);
    return stream;
  }

  friend std::ostream &operator<<(std::ostream &stream, const Int &v) {
    if (v.sign == -1) {
      stream << '-';
    }

    stream << (v.a.empty() ? 0 : v.a.back());

    for (int i = (int)v.a.size() - 2; i >= 0; --i) {
      stream << std::setw(Int::base_digits) << std::setfill('0') << v.a[i];
    }

    return stream;
  }

  static std::vector<int> convert_base(const std::vector<int> &a,
                                       int old_digits, int new_digits);
  static std::vector<long long>
  karatsubaMultiply(const std::vector<long long> &a,
                    const std::vector<long long> &b);
};

Int abs(Int a);
Int fact(Int i);
Int pow(Int base, Int exp);
Int max(Int a, Int b);
Int min(Int a, Int b);

bool operator<(const unsigned int &a, const Int &v);
Int operator*(const int &a, const Int &v);
bool operator>(const unsigned int &a, const Int &v);
bool operator<=(const unsigned int &a, const Int &v);
bool operator>=(const unsigned int &a, const Int &v);

Int operator+(const long long a, const Int &v);
Int operator-(const long long a, const Int &v);
*/
#endif
