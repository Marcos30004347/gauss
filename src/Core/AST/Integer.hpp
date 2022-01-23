#ifndef INTEGER_H
#define INTEGER_H

#include "Int.hpp"
#include <climits>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

struct Int2 {
  Int2(bint<30> *v);
  bint<30> *val = nullptr;

  friend class expression;

public:
  friend struct Int;

  Int2();

  Int2(long int);
  Int2(long long);
  Int2(unsigned long long);
  Int2(unsigned int);
  Int2(int);
  Int2(double);
  Int2(const Int2 &);
  Int2(Int2 &&);
  ~Int2();


  inline Int2 operator+(const Int2 &other) const;
  inline Int2 operator+(const Int2 &&other) const;
  inline Int2 operator+(const int a) const;

  inline Int2 operator-(const Int2 &&other) const;
  inline Int2 operator-(const Int2 &other) const;
  inline Int2 operator-(const int a) const;
  inline Int2 operator*(const Int2 &other) const;
  inline Int2 operator*(const Int2 &&other) const;
  inline Int2 operator*(const int a) const;
  inline Int2 operator/(const Int2 &other) const;
  inline Int2 operator/(const Int2 &&other) const;
  inline Int2 operator/(const int a) const;
  inline Int2 operator%(const Int2 &other) const;
  inline Int2 operator%(const Int2 &&other) const;
  inline Int2 operator%(const int a) const;
  inline Int2 operator+();
  inline Int2 operator-();
  inline Int2 operator++();
  inline Int2 operator--();
  inline Int2 operator++(int);
  inline Int2 operator--(int);
  inline bool operator==(const Int2 &other) const;
  inline bool operator==(const Int2 &&other) const;
  inline bool operator<(const Int2 &other) const;
  inline bool operator<(const Int2 &&other) const;
  inline bool operator<=(const Int2 &other) const;
  inline bool operator<=(const Int2 &&other) const;
  inline bool operator>(const Int2 &other) const;
  inline bool operator>(const Int2 &&other) const;
  inline bool operator>=(const Int2 &other) const;
  inline bool operator>=(const Int2 &&other) const;
  inline bool operator!=(const Int2 &other) const;
  inline bool operator!=(const Int2 &&other) const;
  inline Int2 operator=(const Int2 &other);
  inline Int2 operator=(Int2 &&other);

  std::string to_string();
  inline Int2 ceil_log2();
  inline long long longValue();

  inline double doubleValue();

  friend Int2 gcd(const Int2 &a, const Int2 &b) {
    return Int2(bint<30>::gcd(a.val, b.val));
  }
  friend Int2 gcd(const Int2 &&a, const Int2 &&b) {
    return Int2(bint<30>::gcd(a.val, b.val));
  }

  friend Int2 lcm(const Int2 &a, const Int2 &b) {
    return Int2(bint<30>::lcm(a.val, b.val));
  }

  friend Int2 lcm(const Int2 &&a, const Int2 &&b) {
    return Int2(bint<30>::lcm(a.val, b.val));
  }

  inline void operator/=(Int2 v) { *this = *this / v; }
  inline void operator+=(Int2 v) { *this = *this + v; }
  inline void operator-=(Int2 v) { *this = *this - v; }
  inline void operator*=(Int2 v) { *this = *this * v; }

  friend Int2 abs(const Int2 &&a) { return Int2(bint<30>::abs(a.val)); }
  friend Int2 abs(const Int2 &a) { return Int2(bint<30>::abs(a.val)); }

  friend Int2 fact(const Int2 &&a) { return Int2(bint<30>::fact(a.val)); }
  friend Int2 fact(const Int2 &a) { return Int2(bint<30>::fact(a.val)); }

  friend Int2 max(const Int2 &&a, Int2 &&b) {
    return Int2(bint<30>::max(a.val, b.val));
  }

  friend Int2 max(const Int2 &a, Int2 &b) {
    return Int2(bint<30>::max(a.val, b.val));
  }

  friend Int2 max(const Int2 &a, Int2 &&b) {
    return Int2(bint<30>::max(a.val, b.val));
  }

  friend Int2 min(const Int2 &&a, const Int2 &&b) {
    return Int2(bint<30>::min(a.val, b.val));
  }

  friend Int2 min(const Int2 &a, const Int2 &b) {
    return Int2(bint<30>::min(a.val, b.val));
  }

  friend Int2 pow(const Int2 &&a, const Int2 &&b) {
    return Int2(bint<30>::pow(a.val, b.val));
  }

  friend Int2 pow(const Int2 &a, const Int2 &b) {
    return Int2(bint<30>::pow(a.val, b.val));
  }

  friend double pow(const Int2 &&a, const double b) {
    return bint<30>::pow(a.val, b);
  }

  friend double pow(const Int2 &a, const double b) {
    return bint<30>::pow(a.val, b);
  }

  friend Int2 isqrt(const Int2 &a) {
    bint<30> *res = new bint<30>();
    bint<30>::isqrt(a.val, res, nullptr);
    return res;
  }

  // Int2 gcd(const Int2 &a, const Int2 &b);
  // inline Int2 gcd(const Int2 &&a, const Int2 &&b);
  // inline Int2 lcm(const Int2 &a, const Int2 &b);
  // inline Int2 lcm(const Int2 &&a, const Int2 &&b);

  // inline void operator/=(Int2 v);
  // inline void operator+=(Int2 v);
  // inline void operator-=(Int2 v);
  // inline void operator*=(Int2 v);

  // explicit operator bool() const { return this->val->size > 0; }

  // inline Int2 abs(const Int2 &&a);
  // inline Int2 abs(const Int2 &a);

  // inline Int2 fact(const Int2 &&a);
  // inline Int2 fact(const Int2 &a);
  // inline Int2 max(const Int2 &&a, Int2 &&b);

  // inline Int2 max(const Int2 &a, Int2 &b);
  // inline Int2 max(const Int2 &a, Int2 &&b);
  // inline Int2 min(const Int2 &&a, const Int2 &&b);

  // inline Int2 min(const Int2 &a, const Int2 &b);
  // inline Int2 pow(const Int2 &&a, const Int2 &&b);
  // inline Int2 pow(const Int2 &a, const Int2 &b);

  // inline double pow(const Int2 &&a, const double b);
  // inline double pow(const Int2 &a, const double b);

  // inline Int2 isqrt(const Int2 &a);

  friend bool operator<(const unsigned int &a, const Int2 &v) {
    bint<30> *tmp = bint<30>::from(a);

    bool res = bint<30>::compare(tmp, v.val) < 0;

    delete tmp;

    return res;
  }

  friend Int2 operator*(const int &a, const Int2 &v) {
    bint<30> *tmp = bint<30>::from(a);
    bint<30> *res = new bint<30>();
    bint<30>::mul(tmp, v.val, res);
    delete tmp;
    return res;
  }

  friend bool operator>(const unsigned int &a, const Int2 &v) {
    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) > 0;
    delete tmp;
    return res;
  }

  friend bool operator<=(const unsigned int &a, const Int2 &v) {
    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) <= 0;
    delete tmp;
    return res;
  }

  friend bool operator>=(const unsigned int &a, const Int2 &v) {
    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) >= 0;
    delete tmp;
    return res;
  }

  friend Int2 operator+(const long long a, const Int2 &v) {
    bint<30> *tmp = bint<30>::from(a);
    bint<30> *res = new bint<30>();

    bint<30>::add(tmp, v.val, res);

    delete tmp;

    return res;
  }

  friend Int2 operator-(const long long a, const Int2 &v) {
    bint<30> *tmp = bint<30>::from(a);
    bint<30> *res = new bint<30>();

    bint<30>::sub(tmp, v.val, res);

    delete tmp;

    return res;
  }
};

inline long long safe_mul(long long a, long long b) {
	if(a == 0 || b == 0) return 0;

	long sign = (a < 0 ? -1 : +1) * (b < 0 ? -1 : +1);

	// printf("---> %lli * %lli\n", a, b);
	unsigned long long A = std::abs(a);
	unsigned long long B = std::abs(b);

	// printf("---> %llu * %llu\n", A, B);
	unsigned long long C = A * B;

	// printf("---> %llu * %llu = %llu   rem    %llu\n", A, B, C, C / B);
	//long long c = (unsigned long long)(sa*a) * (unsigned long long)(sb*b);



	if (A == C / B)	{
		if (C >= ((unsigned long long)LONG_LONG_MAX) - 1) {
			return LONG_LONG_MAX;
		}

		return sign * ((long long)C);
	}


	return LONG_LONG_MAX;
}

inline long long safe_add(long long a, long long b) {
  if (a > 0 && b > LONG_LONG_MAX - a) {
    return LONG_LONG_MAX;
  } else if (a < 0 && b < LONG_LONG_MIN - a) {
    return LONG_LONG_MAX;
  }

  return a + b;
}

inline long long safe_pow(long long a, long long b) {
  long long r = 1;

  while (b) {
    if (b % 2 == 1) {
      r = safe_mul(r, a);

      if (r == LONG_LONG_MAX) {
        return LONG_LONG_MAX;
      }
    }

    b = b >> 1;

    a = safe_mul(a, a);

    if (a == LONG_LONG_MAX) {
      return LONG_LONG_MAX;
    }
  }

  return r;
}

inline long long safe_sub(long long a, long long b) { return a - b; }

inline long long safe_div(long long a, long long b) { return a / b; }

inline long long safe_mod(long long a, long long b) { return (b + (a % b)) % b; }
inline long long safe_gcd(long long x, long long y) {
	if(y ==0) return x;
	long long s = x < 0 ? -1 : + 1;
	long long a = std::abs(x);
	long long b = std::abs(y);
	return safe_gcd(y, s*(a % b));
	// long long a, b, r;

  // if (x > y) {
  //   a = x;
  //   b = y;
  // } else {
  //   a = y;
  //   b = x;
  // }

  // while (b != 0) {
  //   r = a % b;
  //   a = b;
  //   b = r;
  // }

  // return a;
}

inline long long safe_lcm(long long a, long long b) {
  return (a / safe_gcd(a, b)) * b;
}

inline long long safe_fact(long long a) {
  if (a == LONG_LONG_MAX)
    return LONG_LONG_MAX;
  if (a == 0 || a == 1)
    return 1;

  if (a < 0)
    a = -a;

  return a * safe_fact(a - 1);
}

inline long long safe_isqrt(long long x) {
    long long q = 1, r = 0;

		while (q <= x) {
        q <<= 2;
    }

		long long t;

		while (q > 1) {
        q >>= 2;
        t = x - r - q;
        r >>= 1;
        if (t >= 0) {
            x = t;
            r += q;
        }
    }

    return r;
}

struct Int {
  Int(bint<30> *v);


  bool flag;

  bint<30> *val = nullptr;
  long long x;
  friend class expression;

	void to_long_if_small() {
		if(!this->flag || this->val->size > 2) return;

		long long v = 0;

		bint<30>::to_long(this->val, &v);

		delete this->val;

		this->flag = 0;
		this->x = v;
	}

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
    if (!this->flag && !other.flag) {
      long long r = safe_add(this->x, other.x);

      if (r != LONG_LONG_MAX) {
				assert(Int(r).to_string() == ( Int2(this->x) + Int2(other.x)).to_string());

				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *c = new bint<30>;

      bint<30>::add(a, b, c);

      delete a;
      delete b;

			assert(Int(c->copy()).to_string() == ( Int2(this->x) + Int2(other.x)).to_string());

      return c;
    }

    bint<30> *tmp = new bint<30>();

    if (this->flag && !other.flag) {
      bint<30> *k = bint<30>::from(other.x);

      bint<30>::add(this->val, k, tmp);

      delete k;

			assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) + Int2(other.x)).to_string());

      return tmp;
    }

    if (!this->flag && other.flag) {
      bint<30> *k = bint<30>::from(this->x);

      bint<30>::add(k, other.val, tmp);

      delete k;

			assert(Int(tmp->copy()).to_string() == ( Int2(this->x) + Int2(other.val->copy())).to_string());

      return tmp;
    }

    bint<30>::add(this->val, other.val, tmp);

		assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) + Int2(other.val->copy())).to_string());

    return tmp;
  }

  inline Int operator+(const Int &&other) const {
        if (!this->flag && !other.flag) {
      long long r = safe_add(this->x, other.x);

      if (r != LONG_LONG_MAX) {
				assert(Int(r).to_string() == ( Int2(this->x) + Int2(other.x)).to_string());

				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *c = new bint<30>();

      bint<30>::add(a, b, c);

      delete a;
      delete b;

			assert(Int(c->copy()).to_string() == ( Int2(this->x) + Int2(other.x)).to_string());

      return c;
    }

    bint<30> *tmp = new bint<30>();

    if (this->flag && !other.flag) {
      bint<30> *k = bint<30>::from(other.x);

      bint<30>::add(this->val, k, tmp);

      delete k;

			assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) + Int2(other.x)).to_string());

      return tmp;
    }

    if (!this->flag && other.flag) {
      bint<30> *k = bint<30>::from(this->x);

      bint<30>::add(k, other.val, tmp);

      delete k;

			assert(Int(tmp->copy()).to_string() == ( Int2(this->x) + Int2(other.val->copy())).to_string());

      return tmp;
    }

    bint<30>::add(this->val, other.val, tmp);

		assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) + Int2(other.val->copy())).to_string());

    return tmp;
  }

  inline Int operator+(const int z) const {
    bint<30> *res = new bint<30>();

    if (!this->flag) {
      long long r = safe_add(this->x, z);

      if (r != LONG_LONG_MAX) {
        return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from<int>(z);

      bint<30>::add(a, b, res);

      delete a;
      delete b;

			assert(Int(res->copy()).to_string() == ( Int2(this->x) + z).to_string());
      return res;
    }

    bint<30> *tmp = bint<30>::from(z);
    bint<30>::add(this->val, tmp, res);

    delete tmp;

		assert(Int(res->copy()).to_string() == ( Int2(this->val->copy()) + z).to_string());

		return res;
  }

  inline Int operator-(const Int &other) const {
    if (!this->flag && !other.flag) {
      long long r = safe_sub(this->x, other.x);

      if (r != LONG_LONG_MAX) {

				assert(Int(r).to_string() == ( Int2(this->x) - Int2(other.x)).to_string());
				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *c = new bint<30>();

      bint<30>::sub(a, b, c);

      delete a;
      delete b;

      if (c->size == 2) {
        bint<30>::to_long(c, &r);

        delete c;

				assert(Int(r).to_string() == ( Int2(this->x) - Int2(other.x)).to_string());

				return r;
      }

			assert(Int(c->copy()).to_string() == ( Int2(this->x) - Int2(other.x)).to_string());

      return c;
    }

    bint<30> *tmp = new bint<30>();

    if (this->flag && !other.flag) {
      long long r;

      bint<30> *k = bint<30>::from(other.x);

      bint<30>::sub(this->val, k, tmp);

      delete k;

      if (tmp->size == 2) {
        bint<30>::to_long(tmp, &r);

        delete tmp;

				assert(Int(r).to_string() == ( Int2(this->val->copy()) - Int2(other.x)).to_string());

				return r;
      }

			assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) - Int2(other.x)).to_string());
      return tmp;
    }

    if (!this->flag && other.flag) {

      long long r;

      bint<30> *k = bint<30>::from(this->x);

      bint<30>::sub(k, other.val, tmp);

      delete k;

      if (tmp->size == 2) {
        bint<30>::to_long(tmp, &r);

        delete tmp;

				assert(Int(r).to_string() == ( Int2(this->x) - Int2(other.val->copy())).to_string());

        return r;
      }

			assert(Int(tmp->copy()).to_string() == ( Int2(this->x) - Int2(other.val->copy())).to_string());
      return tmp;
    }

    bint<30>::sub(this->val, other.val, tmp);

		assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) - Int2(other.val->copy())).to_string());

		return tmp;
  }

  inline Int operator-(const Int &&other) const {
    if (!this->flag && !other.flag) {
      long long r = safe_sub(this->x, other.x);

      if (r != LONG_LONG_MAX) {

				assert(Int(r).to_string() == ( Int2(this->x) - Int2(other.x)).to_string());

				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *c = new bint<30>;

      bint<30>::sub(a, b, c);

      delete a;
      delete b;

      if (c->size == 2) {
        bint<30>::to_long(c, &r);

        delete c;

				assert(Int(r).to_string() == ( Int2(this->x) - Int2(other.x)).to_string());
        return r;
      }

			assert(Int(c->copy()).to_string() == ( Int2(this->x) - Int2(other.x)).to_string());
      return c;
    }

    bint<30> *tmp = new bint<30>();

    if (this->flag && !other.flag) {
      long long r;

      bint<30> *k = bint<30>::from(other.x);

      bint<30>::sub(this->val, k, tmp);

      delete k;

      if (tmp->size == 2) {
        bint<30>::to_long(tmp, &r);

        delete tmp;

				assert(Int(r).to_string() == ( Int2(this->val->copy()) - Int2(other.x)).to_string());

				return r;
      }

			assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) - Int2(other.x)).to_string());

			return tmp;
    }

    if (!this->flag && other.flag) {

      long long r;

      bint<30> *k = bint<30>::from(this->x);

      bint<30>::sub(k, other.val, tmp);

      delete k;

      if (tmp->size == 2) {
        bint<30>::to_long(tmp, &r);

        delete tmp;

				assert(Int(r).to_string() == ( Int2(this->x) - Int2(other.val->copy())).to_string());

        return r;
      }

			assert(Int(tmp->copy()).to_string() == ( Int2(this->x) - Int2(other.val->copy())).to_string());

			return tmp;
    }

    bint<30>::sub(this->val, other.val, tmp);

		assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) - Int2(other.val)).to_string());

		return tmp;
  }

  inline Int operator-(const int z) const {

    bint<30> *res = new bint<30>();

    if (!this->flag) {
      long long r = safe_sub(this->x, z);

      if (r != LONG_LONG_MAX) {

				assert(Int(r).to_string() == ( Int2(this->x) - Int2(z)).to_string());
				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(z);

      bint<30>::sub(a, b, res);

      delete a;
      delete b;

      if (res->size == 2) {
        long long r;

        bint<30>::to_long(res, &r);

        delete res;

				assert(Int(r).to_string() == ( Int2(this->x) - Int2(z)).to_string());

				return r;
      }

			assert(Int(res->copy()).to_string() == ( Int2(this->x) - Int2(z)).to_string());
      return res;
    }

    bint<30> *tmp = bint<30>::from(z);
    bint<30>::sub(this->val, tmp, res);

    delete tmp;

    if (res->size == 2) {
      long long r;

      bint<30>::to_long(res, &r);

      delete res;

			assert(Int(r).to_string() == ( Int2(this->val->copy()) - Int2(z)).to_string());

			return r;
    }

		assert(Int(res->copy()).to_string() == ( Int2(this->val->copy()) - Int2(z)).to_string());

		return res;
  }

  inline Int operator*(const Int &other) const {

    if (!this->flag && !other.flag) {
      long long r = safe_mul(this->x, other.x);

      if (r != LONG_LONG_MAX) {

				assert(Int(r).to_string() == ( Int2(this->x) * Int2(other.x)).to_string());

				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *c = new bint<30>();

      bint<30>::mul(a, b, c);

      delete a;
      delete b;

			assert(Int(c->copy()).to_string() == ( Int2(this->x) * Int2(other.x)).to_string());

			return c;
    }

    bint<30> *tmp = new bint<30>();

    if (this->flag && !other.flag) {
      bint<30> *k = bint<30>::from(other.x);
      bint<30>::mul(this->val, k, tmp);
      delete k;

			assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) * Int2(other.x)).to_string());

			return tmp;
    }

    if (!this->flag && other.flag) {
      bint<30> *k = bint<30>::from(this->x);
      bint<30>::mul(k, other.val, tmp);
      delete k;

			assert(Int(tmp->copy()).to_string() == ( Int2(this->x) * Int2(other.val->copy())).to_string());

			return tmp;
    }

    bint<30>::mul(this->val, other.val, tmp);

		assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) * Int2(other.val->copy())).to_string());

		return tmp;
  }

  inline Int operator*(const Int &&other) const {
		if (!this->flag && !other.flag) {
      long long r = safe_mul(this->x, other.x);

			// printf("%lli * %lli = %lli\n", this->x, other.x, r);

			// printf("a = %s\n", Int(r).to_string().c_str());
			// printf("b = %s\n", (Int2(this->x) * Int2(other.x)).to_string().c_str());

			if (r != LONG_LONG_MAX) {
				assert(Int(r).to_string() == ( Int2(this->x) * Int2(other.x)).to_string());

				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);
			// printf("aa\n");
			// printf("%s     %s\n", a->to_string().c_str(), b->to_string().c_str());
      bint<30> *c = new bint<30>();
			// printf("bb\n");

			// printf("cc\n");

			bint<30>::mul(a, b, c);

      delete a;
      delete b;

			// printf("c = %s\n", Int(c->copy()).to_string().c_str());

			assert(Int(c->copy()).to_string() == ( Int2(this->x) * Int2(other.x)).to_string());

			return c;
    }

    bint<30> *tmp = new bint<30>();

    if (this->flag && !other.flag) {
      bint<30> *k = bint<30>::from(other.x);
      bint<30>::mul(this->val, k, tmp);
      delete k;

			assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) * Int2(other.x)).to_string());

			return tmp;
    }

    if (!this->flag && other.flag) {
      bint<30> *k = bint<30>::from(this->x);
      bint<30>::mul(k, other.val, tmp);
      delete k;

			assert(Int(tmp->copy()).to_string() == ( Int2(this->x) * Int2(other.val->copy())).to_string());

			return tmp->copy();
    }

    bint<30>::mul(this->val, other.val, tmp);

		assert(Int(tmp->copy()).to_string() == ( Int2(this->val->copy()) * Int2(other.val->copy())).to_string());

		return tmp;
  }

  inline Int operator*(const int z) const {

    bint<30> *res = new bint<30>();

    if (!this->flag) {
      long long r = safe_mul(this->x, z);

      if (r != LONG_LONG_MAX) {

				assert(Int(r).to_string() == ( Int2(this->x) * Int2(z)).to_string());

				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(z);

      bint<30>::mul(a, b, res);

      delete a;
      delete b;

			assert(Int(res->copy()).to_string() == ( Int2(this->x) * Int2(z)).to_string());

			return res;
    }

    bint<30> *tmp = bint<30>::from(z);

    bint<30>::mul(this->val, tmp, res);

    delete tmp;

		assert(Int(res->copy()).to_string() == ( Int2(this->val->copy()) * Int2(z)).to_string());

    return res;
  }

  inline Int operator/(const Int &other) const {
    if (!this->flag && !other.flag) {
      long long r = safe_div(this->x, other.x);

      if (r != LONG_LONG_MAX) {
				assert(Int(r).to_string() == ( Int2(this->x) / Int2(other.x)).to_string());

				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *quo = new bint<30>;
      bint<30> *rem = new bint<30>;

      bint<30>::div(a, b, quo, rem);

      delete a;
      delete b;
      delete rem;

      if (quo->size == 2) {
        bint<30>::to_long(quo, &r);

        delete quo;

				assert(Int(r).to_string() == ( Int2(this->x) / Int2(other.x)).to_string());

        return r;
      }

			assert(Int(quo->copy()).to_string() == ( Int2(this->x) / Int2(other.x)).to_string());

			return quo;
    }

    bint<30> *quo = new bint<30>();
    bint<30> *rem = new bint<30>();

    if (this->flag && !other.flag) {
      long long r;

      bint<30> *k = bint<30>::from(other.x);
      bint<30>::div(this->val, k, quo, rem);
      delete k;

      delete rem;

      if (quo->size == 2) {
        bint<30>::to_long(quo, &r);

        delete quo;

				assert(Int(r).to_string() == ( Int2(this->val->copy()) / Int2(other.x)).to_string());

				return r;
      }

			assert(Int(quo->copy()).to_string() == ( Int2(this->val) / Int2(other.x)).to_string());

			return quo;
    }

    if (!this->flag && other.flag) {

      long long r;

      bint<30> *k = bint<30>::from(this->x);

      bint<30>::div(k, other.val, quo, rem);

      delete k;

      delete rem;

      if (quo->size == 2) {
        bint<30>::to_long(quo, &r);

        delete quo;

				assert(Int(r).to_string() == ( Int2(this->x) / Int2(other.val->copy())).to_string());

				return r;
      }

      return quo;
    }

    bint<30>::div(this->val, other.val, quo, rem);

    delete rem;

		assert(Int(quo->copy()).to_string() == ( Int2(this->val->copy()) / Int2(other.val->copy())).to_string());

    return quo;
  }

  inline Int operator/(const Int &&other) const {
   if (!this->flag && !other.flag) {
      long long r = safe_div(this->x, other.x);

      if (r != LONG_LONG_MAX) {
				assert(Int(r).to_string() == ( Int2(this->x) / Int2(other.x)).to_string());

				return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *quo = new bint<30>;
      bint<30> *rem = new bint<30>;

      bint<30>::div(a, b, quo, rem);

      delete a;
      delete b;
      delete rem;

      if (quo->size == 2) {
        bint<30>::to_long(quo, &r);

        delete quo;

				assert(Int(r).to_string() == ( Int2(this->x) / Int2(other.x)).to_string());

        return r;
      }

			assert(Int(quo->copy()).to_string() == ( Int2(this->x) / Int2(other.x)).to_string());

			return quo;
    }

    bint<30> *quo = new bint<30>();
    bint<30> *rem = new bint<30>();

    if (this->flag && !other.flag) {
      long long r;

      bint<30> *k = bint<30>::from(other.x);
      bint<30>::div(this->val, k, quo, rem);
      delete k;

      delete rem;

      if (quo->size == 2) {
        bint<30>::to_long(quo, &r);

        delete quo;

				assert(Int(r).to_string() == ( Int2(this->val->copy()) / Int2(other.x)).to_string());

				return r;
      }

			assert(Int(quo->copy()).to_string() == ( Int2(this->val->copy()) / Int2(other.x)).to_string());

			return quo;
    }

    if (!this->flag && other.flag) {

      long long r;

      bint<30> *k = bint<30>::from(this->x);

      bint<30>::div(k, other.val, quo, rem);

      delete k;

      delete rem;

      if (quo->size == 2) {
        bint<30>::to_long(quo, &r);

        delete quo;

				assert(Int(r).to_string() == ( Int2(this->x) / Int2(other.val->copy())).to_string());

				return r;
      }

      return quo;
    }

    bint<30>::div(this->val, other.val, quo, rem);

    delete rem;

		assert(Int(quo->copy()).to_string() == ( Int2(this->val->copy()) / Int2(other.val->copy())).to_string());

    return quo;
  }

  inline Int operator/(const int z) const {

    if (!this->flag) {
      long long r = safe_div(this->x, z);

      if (r != LONG_LONG_MAX) {

				assert(Int(r).to_string() == ( Int2(this->x) / Int2(z)).to_string());
				return r;
      }

      bint<30> *res = new bint<30>();
      bint<30> *rem = new bint<30>();

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(z);

      bint<30>::div(a, b, res, rem);

      delete a;
      delete b;
      delete rem;

      if (res->size == 2) {
        long long r;

        bint<30>::to_long(res, &r);

        delete res;

				assert(Int(r).to_string() == ( Int2(this->x) / Int2(z)).to_string());
        return r;
      }

      return res;
    }

    bint<30> *res = new bint<30>();
    bint<30> *rem = new bint<30>();

    bint<30> *tmp = bint<30>::from(z);

    bint<30>::div(this->val, tmp, res, rem);

    delete tmp;
    delete rem;

    if (res->size == 2) {
      long long r;

      bint<30>::to_long(res, &r);

      delete res;

			assert(Int(r).to_string() == ( Int2(this->val->copy()) / Int2(z)).to_string());
      return r;
    }

		assert(Int(res->copy()).to_string() == ( Int2(this->val->copy()) / Int2(z)).to_string());

		return res;
  }

  inline Int operator%(const Int &other) const {

    if (!this->flag && !other.flag) {
      long long r = this->x % other.x;

      if (r != LONG_LONG_MAX) {
				assert(Int(r).to_string() == ( Int2(this->x) % Int2(other.x)).to_string());
        return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *quo = new bint<30>;
      bint<30> *rem = new bint<30>;

      bint<30>::div(a, b, quo, rem);

      delete a;
      delete b;
      delete quo;

      if (rem->size == 2) {
        bint<30>::to_long(rem, &r);

        delete rem;

				assert(Int(r).to_string() == ( Int2(this->x) % Int2(other.x)).to_string());
        return r;
      }

			assert(Int(rem->copy()).to_string() == ( Int2(this->x) % Int2(other.x)).to_string());
      return rem;
    }

    bint<30> *quo = new bint<30>();
    bint<30> *rem = new bint<30>();

    if (this->flag && !other.flag) {
      long long r;

      bint<30> *k = bint<30>::from(other.x);

      bint<30>::div(this->val, k, quo, rem);

      delete k;
      delete quo;

      if (rem->size == 2) {
        bint<30>::to_long(rem, &r);

        delete rem;

				assert(Int(r).to_string() == ( Int2(this->val->copy()) % Int2(other.x)).to_string());
        return r;
      }

			assert(Int(rem->copy()).to_string() == ( Int2(this->val->copy()) % Int2(other.x)).to_string());
      return rem;
    }

    if (!this->flag && other.flag) {

      long long r;

      bint<30> *k = bint<30>::from(this->x);

      bint<30>::div(k, other.val, quo, rem);

      delete k;
      delete quo;

      if (rem->size == 2) {
        bint<30>::to_long(rem, &r);

        delete rem;

				assert(Int(r).to_string() == ( Int2(this->x) % Int2(other.val->copy())).to_string());
        return r;
      }

			assert(Int(rem->copy()).to_string() == ( Int2(this->x) % Int2(other.val->copy())).to_string());
      return rem;
    }

    bint<30>::div(this->val, other.val, quo, rem);

    delete quo;

		assert(Int(rem->copy()).to_string() == ( Int2(this->val->copy()) % Int2(other.val->copy())).to_string());

		return rem;
  }

  inline Int operator%(const Int &&other) const {
    if (!this->flag && !other.flag) {
      long long r = this->x % other.x;

      if (r != LONG_LONG_MAX) {
				assert(Int(r).to_string() == ( Int2(this->x) % Int2(other.x)).to_string());
        return r;
      }

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(other.x);

      bint<30> *quo = new bint<30>;
      bint<30> *rem = new bint<30>;

      bint<30>::div(a, b, quo, rem);

      delete a;
      delete b;
      delete quo;

      if (rem->size == 2) {
        bint<30>::to_long(rem, &r);

        delete rem;

				assert(Int(r).to_string() == ( Int2(this->x) % Int2(other.x)).to_string());
        return r;
      }

			assert(Int(rem->copy()).to_string() == ( Int2(this->x) % Int2(other.x)).to_string());
      return rem;
    }

    bint<30> *quo = new bint<30>();
    bint<30> *rem = new bint<30>();

    if (this->flag && !other.flag) {
      long long r;

      bint<30> *k = bint<30>::from(other.x);

      bint<30>::div(this->val, k, quo, rem);

      delete k;
      delete quo;

      if (rem->size == 2) {
        bint<30>::to_long(rem, &r);

        delete rem;

				assert(Int(r).to_string() == ( Int2(this->val->copy()) % Int2(other.x)).to_string());
        return r;
      }

			assert(Int(rem->copy()).to_string() == ( Int2(this->val->copy()) % Int2(other.x)).to_string());
      return rem;
    }

    if (!this->flag && other.flag) {

      long long r;

      bint<30> *k = bint<30>::from(this->x);

      bint<30>::div(k, other.val, quo, rem);

      delete k;
      delete quo;

      if (rem->size == 2) {
        bint<30>::to_long(rem, &r);

        delete rem;

				assert(Int(r).to_string() == ( Int2(this->x) % Int2(other.val->copy())).to_string());
        return r;
      }

			assert(Int(rem->copy()).to_string() == ( Int2(this->x) % Int2(other.val->copy())).to_string());
      return rem;
    }

    bint<30>::div(this->val, other.val, quo, rem);

    delete quo;

		assert(Int(rem->copy()).to_string() == ( Int2(this->val->copy()) % Int2(other.val->copy())).to_string());

		return rem;
  }

  inline Int operator%(const int z) const {
    if (!this->flag) {
      long long r = this->x % z;

      if (r != LONG_LONG_MAX) {
				assert(Int(r).to_string() == ( Int2(this->x) % Int2(z)).to_string());
        return r;
      }

      bint<30> *res = new bint<30>();
      bint<30> *rem = new bint<30>();

      bint<30> *a = bint<30>::from(this->x);
      bint<30> *b = bint<30>::from(z);

      bint<30>::div(a, b, res, rem);

      delete a;
      delete b;
      delete res;

      if (rem->size == 2) {
        long long r;

        bint<30>::to_long(rem, &r);

        delete rem;

				assert(Int(r).to_string() == ( Int2(this->x) % Int2(z)).to_string());
        return r;
      }

			assert(Int(rem->copy()).to_string() == ( Int2(this->x) % Int2(z)).to_string());
      return rem;
    }

    bint<30> *res = new bint<30>();
    bint<30> *rem = new bint<30>();

    bint<30> *tmp = bint<30>::from(z);

    bint<30>::div(this->val, tmp, res, rem);

    delete tmp;
    delete res;

    if (rem->size == 2) {
      long long r;

      bint<30>::to_long(rem, &r);

      delete rem;

			assert(Int(r).to_string() == ( Int2(this->val->copy()) % Int2(z)).to_string());
      return r;
    }

		assert(Int(rem->copy()).to_string() == ( Int2(this->val->copy()) % Int2(z)).to_string());

		return rem;

  }

  inline Int operator+() {
		return *this;
	}

  inline Int operator-() {
		return *this * -1;
	}

  inline Int operator++() {

		*this = *this + 1;

		//assert(this->to_string() == ( Int2(this->val->copy()) + 1).to_string());

		return *this;
  }

  inline Int operator--() {
    *this = *this - 1;
		//assert(this->to_string() == ( Int2(this->val->copy()) - 1).to_string());
    return *this;
  }

  inline Int operator++(int) {
    Int tmp = *this;

		*this = *this + 1;

		//assert(this->to_string() == ( Int2(this->val->copy()) + 1).to_string());

		return tmp;
  }

  inline Int operator--(int) {
    Int tmp = *this;

		*this = *this - 1;

		//assert(this->to_string() == (Int2(this->val->copy()) - 1).to_string());

		return tmp;
  }

  inline bool operator==(const Int &other) const {
    if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x == other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) == 0;

      delete tmp;

			assert(r == (Int2(this->x) == Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) == 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) == Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) == 0;

		assert(r == (Int2(this->val->copy()) == Int2(other.val->copy())));

		return r;
  }

  inline bool operator==(const Int &&other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x == other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) == 0;

      delete tmp;

			assert(r == (Int2(this->x) == Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) == 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) == Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) == 0;

		assert(r == (Int2(this->val->copy()) == Int2(other.val->copy())));

		return r;
  }

  inline bool operator<(const Int &other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x < other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) < 0;

      delete tmp;

			assert(r == (Int2(this->x) < Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) < 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) < Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) < 0;

		assert(r == (Int2(this->val->copy()) < Int2(other.val->copy())));

		return r;
  }

  inline bool operator<(const Int &&other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x < other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) < 0;

      delete tmp;

			assert(r == (Int2(this->x) < Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) < 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) < Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) < 0;

		assert(r == (Int2(this->val->copy()) < Int2(other.val->copy())));

		return r;
  }


	inline bool operator<=(const Int &other) const {

		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x <= other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) <= 0;

      delete tmp;

			assert(r == (Int2(this->x) <= Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) <= 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) <= Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) <= 0;

		assert(r == (Int2(this->val->copy()) <= Int2(other.val->copy())));

		return r;
  }

  inline bool operator<=(const Int &&other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x <= other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) <= 0;

      delete tmp;

			assert(r == (Int2(this->x) <= Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) <= 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) <= Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) <= 0;

		assert(r == (Int2(this->val->copy()) <= Int2(other.val->copy())));

		return r;
  }

  inline bool operator>(const Int &other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x > other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) > 0;

      delete tmp;

			assert(r == (Int2(this->x) > Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) > 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) > Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) > 0;

		assert(r == (Int2(this->val->copy()) > Int2(other.val->copy())));

		return r;
  }

  inline bool operator>(const Int &&other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x > other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) > 0;

      delete tmp;

			assert(r == (Int2(this->x) > Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) > 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) > Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) > 0;

		assert(r == (Int2(this->val->copy()) > Int2(other.val->copy())));

		return r;
  }

  inline bool operator>=(const Int &other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x >= other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) >= 0;

      delete tmp;

			assert(r == (Int2(this->x) >= Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) >= 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) >= Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) >= 0;

		assert(r == (Int2(this->val->copy()) >= Int2(other.val->copy())));

		return r;
  }

  inline bool operator>=(const Int &&other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x >= other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) >= 0;

      delete tmp;

			assert(r == (Int2(this->x) >= Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) >= 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) >= Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) >= 0;

		assert(r == (Int2(this->val->copy()) >= Int2(other.val->copy())));

		return r;
  }

  inline bool operator!=(const Int &other) const {
		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x != other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) != 0;

      delete tmp;

			assert(r == (Int2(this->x) != Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) != 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) != Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) != 0;

		assert(r == (Int2(this->val->copy()) != Int2(other.val->copy())));

		return r;
  };

  inline bool operator!=(const Int &&other) const {

		if (!this->flag && !other.flag) {
      // printf("a) %lli  ==  %lli\n", this->x, other.x);
      return this->x != other.x;
    }

    if (!this->flag && other.flag) {
      // printf("b) %lli  ==  %s\n", this->x, other.val->to_string().c_str());

      bint<30> *tmp = bint<30>::from(this->x);

      bool r = bint<30>::compare(tmp, other.val) != 0;

      delete tmp;

			assert(r == (Int2(this->x) != Int2(other.val->copy())));

      return r;
    }

    if (this->flag && !other.flag) {
      // printf("c) %s  ==  %lli\n", this->val->to_string().c_str(), other.x);

      bint<30> *tmp = bint<30>::from(other.x);

      bool r = bint<30>::compare(this->val, tmp) != 0;

      delete tmp;

			assert(r == (Int2(this->val->copy()) != Int2(other.x)));
      return r;
    }

    // printf("d) %s  ==  %s\n", this->val->to_string().c_str(),
    // other.val->to_string().c_str());

		bool r = bint<30>::compare(this->val, other.val) != 0;

		assert(r == (Int2(this->val->copy()) != Int2(other.val->copy())));

		return r;
  };

  inline Int operator=(const Int &other) {

    if (this->flag && this->val)
      delete this->val;

    this->flag = other.flag;

    if (other.flag) {
      this->val = other.val->copy();
    } else {
      this->x = other.x;
    }

    return *this;
  }

  inline Int operator=(Int &&other) {

    if (this->flag && this->val)
      delete this->val;

    this->flag = other.flag;

    if (other.flag) {
      this->val = other.val->copy();
    } else {
      this->x = other.x;
    }

    return *this;

    // if(this->val) delete this->val;

    // this->val = other.val;
    // other.val = nullptr;

    // return *this;
  }

  std::string to_string();

  inline Int ceil_log2() {
    if (!flag) {
			assert((Int2(this->val->copy()).ceil_log2() == std::ceil(std::log2(x))));

			return std::ceil(std::log2(x));
		}
    return bint<30>::ceil_log2(this->val);
  }

  inline long long longValue() {
    if (!flag)
      return x;

    long long v = 0;

    short error = bint<30>::to_long(this->val, &v);

    if (error == -1) {
      // TODO: better error handling
      printf("long long overflow\n");
      exit(1);
    }

    return v;
  }

  inline double doubleValue() {
    if (!flag)
      return x;

    double v = 0.0;
    // TODO: if to_long returned -1, that means an overflow
    bint<30>::to_double(this->val, &v);
    return v;
  }

  friend Int gcd(const Int &a, const Int &b) {
    if (!a.flag && !b.flag) {
      long long z = safe_gcd(a.x, b.x);

			// assert(gcd(Int2(a.x), Int2(b.x)) == z);

			return z;
    }

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);

			bint<30> *r = bint<30>::gcd(k, b.val);

			// TODO: convert r to long long if it have size 2

			delete k;

			// assert(gcd(Int2(a.x), Int2(b.val->copy())) == Int2(r->copy()));

			return r;
    }
    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::gcd(a.val, k);
      delete k;

			// assert(gcd(Int2(a.val->copy()), Int2(b.x)) == Int2(r->copy()));

			return r;
    }

		bint<30>* g = bint<30>::gcd(a.val, b.val);

		// assert(gcd(Int2(a.val->copy()), Int2(b.val->copy())) == Int2(g->copy()));

		return g;
	}

  friend Int gcd(const Int &&a, const Int &&b) {
    if (!a.flag && !b.flag) {

			long long z = safe_gcd(a.x, b.x);
			// printf("%s   %lli\n", gcd(Int2(a.x), Int2(b.x)).to_string().c_str(), z);
			// assert(gcd(Int2(a.x), Int2(b.x)) == z);

			return z;
    }

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);

			bint<30> *r = bint<30>::gcd(k, b.val);

			// TODO: convert r to long long if it have size 2

			delete k;

			// assert(gcd(Int2(a.x), Int2(b.val->copy())) == Int2(r->copy()));

			return r;
    }
    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::gcd(a.val, k);
      delete k;

			// assert(gcd(Int2(a.val->copy()), Int2(b.x)) == Int2(r->copy()));

			return r;
    }

		bint<30>* g = bint<30>::gcd(a.val, b.val);

		// assert(gcd(Int2(a.val->copy()), Int2(b.val->copy())) == Int2(g->copy()));

		return g;
  }

  friend Int lcm(const Int &a, const Int &b) {
    if (!a.flag && !b.flag) {
      long long z = safe_lcm(a.x, b.x);

			// assert(lcm(Int2(a.x), Int2(b.x)) == z);

			return z;
    }

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);

			bint<30> *r = bint<30>::lcm(k, b.val);

			// TODO: convert r to long long if it have size 2

			delete k;

			// assert(lcm(Int2(a.x), Int2(b.val->copy())) == Int2(r->copy()));

			return r;
    }
    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::lcm(a.val, k);
      delete k;

			// assert(lcm(Int2(a.val->copy()), Int2(b.x)) == Int2(r->copy()));

			return r;
    }

		bint<30>* g = bint<30>::lcm(a.val, b.val);

		// assert(lcm(Int2(a.val->copy()), Int2(b.val->copy())) == Int2(g->copy()));

		return g;
  }

  friend Int lcm(const Int &&a, const Int &&b) {
    if (!a.flag && !b.flag) {
      long long z = safe_lcm(a.x, b.x);

			// assert(lcm(Int2(a.x), Int2(b.x)) == z);

			return z;
    }

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);

			bint<30> *r = bint<30>::lcm(k, b.val->copy());

			// TODO: convert r to long long if it have size 2

			delete k;

			// assert(lcm(Int2(a.x), Int2(b.val->copy())) == Int2(r->copy()));

			return r;
    }
    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::lcm(a.val, k);
      delete k;

			// assert(lcm(Int2(a.val->copy()), Int2(b.x)) == Int2(r->copy()));

			return r;
    }

		bint<30>* g = bint<30>::lcm(a.val, b.val);

		// assert(lcm(Int2(a.val->copy()), Int2(b.val->copy())) == Int2(g->copy()));

		return g;
  }

  inline void operator/=(Int v) { *this = *this / v; }
  inline void operator+=(Int v) { *this = *this + v; }
  inline void operator-=(Int v) { *this = *this - v; }
  inline void operator*=(Int v) { *this = *this * v; }

  explicit operator bool() const {
    if (!this->flag)
      return this->x != 0;
    return this->val->size > 0;
  }

  friend Int abs(const Int &&a) {
    if (!a.flag) return std::abs(a.x);
    return bint<30>::abs(a.val);
  }

  friend Int abs(const Int &a) {
    if (!a.flag) return std::abs(a.x);
    return bint<30>::abs(a.val);
  }

  friend Int fact(const Int &&a) {
    if (!a.flag) {
			long long z = safe_fact(a.x);

			if(z == LONG_LONG_MAX) {
				bint<30>* x = bint<30>::from(a.x);
				bint<30>* y = bint<30>::fact(x);

				delete x;

				return y;
			}

			assert(fact(Int2(a.x)) == Int2(z));

			return z;
		}

		bint<30> * g = bint<30>::fact(a.val);

		assert(fact(Int2(a.val->copy())) == Int2(g->copy()));

		return g;
	}

  friend Int fact(const Int &a) {
    if (!a.flag) {
			long long z = safe_fact(a.x);

			if(z == LONG_LONG_MAX) {
				bint<30>* x = bint<30>::from(a.x);
				bint<30>* y = bint<30>::fact(x);

				delete x;

				return y;
			}

			assert(fact(Int2(a.x)) == Int2(z));

			return z;
		}

		bint<30> * g = bint<30>::fact(a.val);

		assert(fact(Int2(a.val->copy())) == Int2(g->copy()));

		return g;
  }

  friend Int max(const Int &&a, Int &&b) {
    if (!a.flag && !b.flag)
      return std::max(a.x, b.x);

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);
      bint<30> *r = bint<30>::max(k, b.val);
      delete k;
      return r;
    }

    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::max(a.val, k);
      delete k;
      return r;
    }

    return bint<30>::max(a.val, b.val);
  }

  friend Int max(const Int &a, Int &b) {
    if (!a.flag && !b.flag)
      return std::max(a.x, b.x);

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);
      bint<30> *r = bint<30>::max(k, b.val);
      delete k;
      return r;
    }

    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::max(a.val, k);
      delete k;
      return r;
    }

    return bint<30>::max(a.val, b.val);
  }

  friend Int max(const Int &a, Int &&b) {

    if (!a.flag && !b.flag)
      return std::max(a.x, b.x);

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);
      bint<30> *r = bint<30>::max(k, b.val);
      delete k;
      return r;
    }

    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::max(a.val, k);
      delete k;
      return r;
    }

    return bint<30>::max(a.val, b.val);
  }

  friend Int min(const Int &&a, const Int &&b) {
    if (!a.flag && !b.flag)
      return std::min(a.x, b.x);

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);
      bint<30> *r = bint<30>::min(k, b.val);
      delete k;
      return r;
    }

    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::min(a.val, k);
      delete k;
      return r;
    }

    return bint<30>::min(a.val, b.val);
  }

  friend Int min(const Int &a, const Int &b) {
    if (!a.flag && !b.flag)
      return std::min(a.x, b.x);

    if (!a.flag && b.flag) {
      bint<30> *k = bint<30>::from(a.x);
      bint<30> *r = bint<30>::min(k, b.val);
      delete k;
      return r;
    }

    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::min(a.val, k);
      delete k;
      return r;
    }

    return bint<30>::min(a.val, b.val);
  }

  friend Int pow(const Int &&a, const Int &&b) {
    if (!a.flag && !b.flag) {
			long long z = safe_pow(a.x, b.x);

			if(z == LONG_LONG_MAX) {
				bint<30>* x = bint<30>::from(a.x);
				bint<30>* y = bint<30>::from(b.x);
				bint<30>* w = bint<30>::pow(x, y);

				delete x;
				delete y;

				assert(pow(Int2(a.x), Int2(b.x)) == Int2(w->copy()));

				return w;
			}

			assert(pow(Int2(a.x), Int2(b.x)) == Int2(z));

			return z;
		}
    if (!a.flag && b.flag) {

      bint<30> *k = bint<30>::from(a.x);
      bint<30> *r = bint<30>::pow(k, b.val);

			delete k;

			assert(pow(Int2(a.x), Int2(b.val->copy())) == Int2(r->copy()));

			return r;
    }

    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::pow(a.val, k);

			delete k;

			assert(pow(Int2(a.val->copy()), Int2(b.x)) == Int2(r->copy()));

			return r;
    }

		bint<30>* r = bint<30>::pow(a.val, b.val);

		assert(pow(Int2(a.val->copy()), Int2(b.val->copy())) == Int2(r->copy()));

    return r;
  }

  friend Int pow(const Int &a, const Int &b) {
    if (!a.flag && !b.flag) {
			long long z = safe_pow(a.x, b.x);

			if(z == LONG_LONG_MAX) {
				bint<30>* x = bint<30>::from(a.x);
				bint<30>* y = bint<30>::from(b.x);
				bint<30>* w = bint<30>::pow(x, y);

				delete x;
				delete y;

				assert(pow(Int2(a.x), Int2(b.x)) == Int2(w->copy()));

				return w;
			}

			assert(pow(Int2(a.x), Int2(b.x)) == Int2(z));

			return z;
		}
    if (!a.flag && b.flag) {

      bint<30> *k = bint<30>::from(a.x);
      bint<30> *r = bint<30>::pow(k, b.val);

			delete k;

			assert(pow(Int2(a.x), Int2(b.val->copy())) == Int2(r->copy()));

			return r;
    }

    if (a.flag && !b.flag) {
      bint<30> *k = bint<30>::from(b.x);
      bint<30> *r = bint<30>::pow(a.val, k);

			delete k;

			assert(pow(Int2(a.val->copy()), Int2(b.x)) == Int2(r->copy()));

			return r;
    }

		bint<30>* r = bint<30>::pow(a.val, b.val);

		assert(pow(Int2(a.val->copy()), Int2(b.val->copy())) == Int2(r->copy()));

    return r;
  }

  friend double pow(const Int &&a, const double b) {
    if (!a.flag) {
			double r = std::pow(a.x, b);
			assert(pow(Int2(a.x), b) == r);
			return r;
		}

	  double r = bint<30>::pow(a.val, b);

		assert(pow(Int2(a.val->copy()), b) == r);

		return r;
  }

  friend double pow(const Int &a, const double b) {
    if (!a.flag) {
			double r = std::pow(a.x, b);
			assert(pow(Int2(a.x), b) == r);
			return r;
		}

	  double r = bint<30>::pow(a.val, b);

		assert(pow(Int2(a.val->copy()), b) == r);

		return r;
  }

  friend Int isqrt(const Int &a) {
    if (!a.flag) {
			long long r = safe_isqrt(a.x);

			assert(isqrt(Int2(a.x)) == r);

			return r;
		}

    bint<30> *res = new bint<30>();

		bint<30>::isqrt(a.val, res, nullptr);

		assert(isqrt(Int2(a.val->copy())) == Int2(res->copy()));

		// TODO: convert to long long if the size of r is <= 2

    return res;
  }

  friend bool operator<(const unsigned int &a, const Int &v) {
    if (!v.flag)
      return a < v.x;

    bint<30> *tmp = bint<30>::from(a);

    bool res = bint<30>::compare(tmp, v.val) < 0;

    delete tmp;

    return res;
  }

  friend Int operator*(const int &a, const Int &v) {
    if (!v.flag)
      return a * v.x;

    bint<30> *tmp = bint<30>::from(a);
    bint<30> *res = new bint<30>();
    bint<30>::mul(tmp, v.val, res);
    delete tmp;
    return res;
  }

  friend bool operator>(const unsigned int &a, const Int &v) {
    if (!v.flag)
      return a > v.x;

    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) > 0;
    delete tmp;
    return res;
  }

  friend bool operator<=(const unsigned int &a, const Int &v) {
    if (!v.flag)
      return a <= v.x;

    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) <= 0;

    delete tmp;
    return res;
  }

  friend bool operator>=(const unsigned int &a, const Int &v) {
    if (!v.flag)
      return a >= v.x;

    bint<30> *tmp = bint<30>::from(a);
    bool res = bint<30>::compare(tmp, v.val) >= 0;
    delete tmp;
    return res;
  }

  friend Int operator+(const long long a, const Int &v) {
    if (!v.flag)
      return a + v.x;

    bint<30> *tmp = bint<30>::from(a);
    bint<30> *res = new bint<30>();

    bint<30>::add(tmp, v.val, res);

    delete tmp;

    return res;
  }

  friend Int operator-(const long long a, const Int &v) {
    if (!v.flag)
      return a - v.x;

    bint<30> *tmp = bint<30>::from(a);
    bint<30> *res = new bint<30>();

    bint<30>::sub(tmp, v.val, res);

    delete tmp;

    return res;
  }
};

inline Int2 Int2::operator+(const Int2 &other) const {
  bint<30> *tmp = new bint<30>();
  bint<30>::add(this->val, other.val, tmp);
  return Int2(tmp);
}

inline Int2 Int2::operator+(const Int2 &&other) const {
  bint<30> *tmp = new bint<30>();
  bint<30>::add(this->val, other.val, tmp);
  return Int2(tmp);
}

inline Int2 Int2::operator+(const int a) const {
  bint<30> *tmp = bint<30>::from(a);
  bint<30> *res = new bint<30>();
  bint<30>::add(this->val, tmp, res);

  delete tmp;

  return Int2(res);
}

inline Int2 Int2::operator-(const Int2 &other) const {
  bint<30> *tmp = new bint<30>();
  bint<30>::sub(this->val, other.val, tmp);
  return Int2(tmp);
}

inline Int2 Int2::operator-(const Int2 &&other) const {
  bint<30> *tmp = new bint<30>();
  bint<30>::sub(this->val, other.val, tmp);
  return Int2(tmp);
}

inline Int2 Int2::operator-(const int a) const {
  bint<30> *tmp = bint<30>::from(a);
  bint<30> *res = new bint<30>();
  bint<30>::sub(this->val, tmp, res);

  delete tmp;

  return Int2(res);
}

inline Int2 Int2::operator*(const Int2 &other) const {
  bint<30> *tmp = new bint<30>();
  bint<30>::mul(this->val, other.val, tmp);
  return Int2(tmp);
}

inline Int2 Int2::operator*(const Int2 &&other) const {
  bint<30> *tmp = new bint<30>();
  bint<30>::mul(this->val, other.val, tmp);
  return Int2(tmp);
}

inline Int2 Int2::operator*(const int a) const {
  bint<30> *tmp = bint<30>::from(a);
  bint<30> *res = new bint<30>();
  bint<30>::mul(this->val, tmp, res);

  delete tmp;

  return Int2(res);
}

inline Int2 Int2::operator/(const Int2 &other) const {
  bint<30> *quo = new bint<30>();
  bint<30> *rem = new bint<30>();

  bint<30>::div(this->val, other.val, quo, rem);

  delete rem;

  return Int2(quo);
}

inline Int2 Int2::operator/(const Int2 &&other) const {
  bint<30> *quo = new bint<30>();
  bint<30> *rem = new bint<30>();

  bint<30>::div(this->val, other.val, quo, rem);

  delete rem;

  return Int2(quo);
}

inline Int2 Int2::operator/(const int a) const {
  bint<30> *tmp = bint<30>::from(a);
  bint<30> *quo = new bint<30>();
  bint<30> *rem = new bint<30>();
  bint<30>::div(this->val, tmp, quo, rem);

  delete rem;
  delete tmp;

  return Int2(quo);
}

inline Int2 Int2::operator%(const Int2 &other) const {
  bint<30> *quo = new bint<30>();
  bint<30> *rem = new bint<30>();

  bint<30>::div(this->val, other.val, quo, rem);

  delete quo;

  return Int2(rem);
}

inline Int2 Int2::operator%(const Int2 &&other) const {
  bint<30> *quo = new bint<30>();
  bint<30> *rem = new bint<30>();

  bint<30>::div(this->val, other.val, quo, rem);

  delete quo;

  return Int2(rem);
}

inline Int2 Int2::operator%(const int a) const {
  bint<30> *tmp = bint<30>::from(a);
  bint<30> *quo = new bint<30>();
  bint<30> *rem = new bint<30>();
  bint<30>::div(this->val, tmp, quo, rem);

  delete quo;
  delete tmp;

  return Int2(rem);
}

inline Int2 Int2::operator+() { return *this; }
inline Int2 Int2::operator-() { return this->val->sign * -1; }

inline Int2 Int2::operator++() {
  *this = *this + 1;
  return *this;
}

inline Int2 Int2::operator--() {
  *this = *this - 1;
  return *this;
}

inline Int2 Int2::operator++(int) {
  Int2 tmp = *this;
  *this = *this + 1;
  return tmp;
}

inline Int2 Int2::operator--(int) {
  Int2 tmp = *this;
  *this = *this - 1;
  return tmp;
}

inline bool Int2::operator==(const Int2 &other) const {
  return bint<30>::compare(this->val, other.val) == 0;
}

inline bool Int2::operator==(const Int2 &&other) const {
  return bint<30>::compare(this->val, other.val) == 0;
}

inline bool Int2::operator<(const Int2 &other) const {
  return bint<30>::compare(this->val, other.val) < 0;
}
inline bool Int2::operator<(const Int2 &&other) const {
  return bint<30>::compare(this->val, other.val) < 0;
}
inline bool Int2::operator<=(const Int2 &other) const {
  return bint<30>::compare(this->val, other.val) <= 0;
}
inline bool Int2::operator<=(const Int2 &&other) const {
  return bint<30>::compare(this->val, other.val) <= 0;
}
inline bool Int2::operator>(const Int2 &other) const {
  return bint<30>::compare(this->val, other.val) > 0;
}
inline bool Int2::operator>(const Int2 &&other) const {
  return bint<30>::compare(this->val, other.val) > 0;
}
inline bool Int2::operator>=(const Int2 &other) const {
  return bint<30>::compare(this->val, other.val) >= 0;
}
inline bool Int2::operator>=(const Int2 &&other) const {
  return bint<30>::compare(this->val, other.val) >= 0;
}
inline bool Int2::operator!=(const Int2 &other) const {
  return bint<30>::compare(this->val, other.val) != 0;
};
inline bool Int2::operator!=(const Int2 &&other) const {
  return bint<30>::compare(this->val, other.val) != 0;
};

inline Int2 Int2::operator=(const Int2 &other) {
  if (this->val)
    delete this->val;

  this->val = other.val->copy();

  return *this;
}

inline Int2 Int2::operator=(Int2 &&other) {
  if (this->val)
    delete this->val;

  this->val = other.val;
  other.val = nullptr;

  return *this;
}

inline Int2 Int2::ceil_log2() { return bint<30>::ceil_log2(this->val); }

inline long long Int2::longValue() {
  long long v = 0;

  short error = bint<30>::to_long(this->val, &v);

  if (error == -1) {
    // TODO: better error handling
    printf("long long overflow\n");
    exit(1);
  }

  return v;
}

inline double Int2::doubleValue() {
  double v = 0.0;
  // TODO: if to_long returned -1, that means an overflow
  bint<30>::to_double(this->val, &v);
  return v;
}

#endif
