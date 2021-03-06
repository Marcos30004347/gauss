#include "Integer.hpp"
#include "gauss/Algebra/Int.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>

inline long long safe_mul(long long a, long long b, long long *c) {
  if (a == 0 || b == 0) {
    *c = 0;
    return LONG_LONG_OK;
  }

  long sign = (a < 0 ? -1 : +1) * (b < 0 ? -1 : +1);

  unsigned long long A = std::abs(a);
  unsigned long long B = std::abs(b);
  unsigned long long C = A * B;

  if (A == C / B) {
    if (C >= ((unsigned long long)LLONG_MAX) - 1) {
      return LONG_LONG_OVERFLOW;
    }

    *c = sign * ((long long)C);

    return LONG_LONG_OK;
  }

  return LONG_LONG_OVERFLOW;
}

inline long long safe_add(long long a, long long b, long long *c) {
  if (a > 0 && b > LLONG_MAX - a) {
    return LONG_LONG_OVERFLOW;
  } else if (a < 0 && b < LLONG_MIN - a) {
    return LONG_LONG_OVERFLOW;
  }

  *c = a + b;

  return LONG_LONG_OK;
}

inline long long safe_pow(long long a, long long b, long long *c) {
  if (b < 0) {

    if (safe_pow(a, -b, c) == LONG_LONG_OVERFLOW) {
      return LONG_LONG_OVERFLOW;
    }

    *c = 1 / *c;

    return LONG_LONG_OK;
  }

  long long r = 1;
  long long t = 1;

  while (b) {
    if (b % 2 == 1) {

      if (safe_mul(r, a, &t) == LONG_LONG_OVERFLOW) {
        return LONG_LONG_OVERFLOW;
      }

      r = t;
    }

    b = b >> 1;

    if (safe_mul(a, a, &t) == LONG_LONG_OVERFLOW) {
      return LONG_LONG_OVERFLOW;
    }

    a = t;
  }

  *c = r;

  return LONG_LONG_OK;
}

inline long long safe_sub(long long a, long long b, long long *c) {
  if ((b > 0) && (a < LLONG_MIN + b)) {
    return LONG_LONG_OVERFLOW;
  }
  if ((b < 0) && (a > LLONG_MAX + b)) {
    return LONG_LONG_OVERFLOW;
  }

  *c = a - b;
  return LONG_LONG_OK;
}

inline long long safe_div(long long a, long long b, long long *c) {
  *c = a / b;
  return LONG_LONG_OK;
}

inline long long safe_mod(long long a, long long b, long long *c) {
  *c = (b + (a % b)) % b;
  return LONG_LONG_OK;
}

inline long long safe_gcd(long long x, long long y, long long *c) {
  if (y == 0) {
    *c = x;
    return LONG_LONG_OK;
  }

  long long s = x < 0 ? -1 : +1;

  long long a = std::abs(x);
  long long b = std::abs(y);

  return safe_gcd(y, s * (a % b), c);
}

inline long long safe_lcm(long long a, long long b, long long *c) {
  safe_gcd(a, b, c);
  *c = (a / *c) * b;

  return LONG_LONG_OK;
}

inline long long safe_fact(long long a, long long *c) {
  if (a == 0 || a == 1) {
    *c = 1;

    return LONG_LONG_OK;
  }

  long long t;

  if (safe_fact(a - 1, &t) == LONG_LONG_OVERFLOW) {
    return LONG_LONG_OVERFLOW;
  }

  *c = a * t;

  return LONG_LONG_OK;
}

inline long long safe_isqrt(long long x, long long *c) {
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

  *c = r;

  return LONG_LONG_OK;
}

Int::Int() {
  this->flag = 0;
  this->x = 0;
}

Int Int::fromString(const char* a) {
	unsigned digits = strlen(a);

	Int x = 0;

	Int e = 1;

	for(unsigned i = 0; i < digits; i++) {
		x += (a[i] - '0') * e;
		e *= 10;
	}

	return x;

}

Int::Int(const Int &a) {
  if (!a.flag) {
    this->flag = 0;
    this->x = a.x;
  } else {
    this->flag = 1;
    this->val = a.val->copy();
  }
}

Int::Int(Int &&a) {
  if (!a.flag) {
    this->flag = 0;
    this->x = a.x;
  } else {
    this->flag = 1;
    this->val = a.val;

    a.flag = 0;
    a.val = 0;
  }
}

Int::Int(long int v) {
  if (v < LLONG_MAX) {
    this->flag = 0;
    this->x = v;
  } else {
    this->flag = 1;
    this->val = bint<30>::from<long int>(v);
  }
}

Int::Int(long long v) {
  if (v < LLONG_MAX) {
    this->flag = 0;
    this->x = v;
  } else {
    this->flag = 1;
    this->val = bint<30>::from<long long>(v);
  }
}

Int::Int(unsigned long long v) {
  if (v < (unsigned long long)LLONG_MAX) {
    this->flag = 0;
    this->x = v;
  } else {
    this->flag = 1;
    this->val = bint<30>::from<unsigned long long>(v);
  }
}

Int::Int(unsigned long int v) {
  if (v < (unsigned long int)LLONG_MAX) {
    this->flag = 0;
    this->x = v;
  } else {
    this->flag = 1;
    this->val = bint<30>::from<unsigned long long>(v);
  }
}

Int::Int(unsigned int v) {
  if (v < (unsigned int)LLONG_MAX) {
    this->flag = 0;
    this->x = v;
  } else {
    this->flag = 1;
    this->val = bint<30>::from<unsigned int>(v);
  }
}

Int::Int(int v) {
  this->flag = 0;
  this->x = v;
}

Int::Int(double b) {
  if (b < (double)LLONG_MAX) {
    this->flag = 0;
    this->x = b;
  } else {
    this->flag = 1;
    this->val = bint<30>::from(b);
  }
}

Int::~Int() {
  if (this->flag && this->val)
    delete this->val;
}

Int::Int(bint<30> *v) {
  this->flag = 1;
  this->val = v;
}

std::string Int::to_string() {
  if (!this->flag) {
		return std::to_string(x);
	}

  return this->val->to_string();
}

Int Int::operator+(const Int &other) const {
  switch ((this->flag << 1) | other.flag) {
  case 0: {
    long long r = 0;

    if (safe_add(this->x, other.x, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    bint<30> *c = new bint<30>;

    bint<30>::add(this->x, other.x, c);

    return c;
  }

  case 1: {
    bint<30> *t = new bint<30>();
    bint<30>::add(this->x, other.val, t);
    return t;
  }

  case 2: {
    bint<30> *t = new bint<30>();
		bint<30>::add(this->val, other.x, t);
    return t;
  }

  case 3: {
    bint<30> *t = new bint<30>();
		bint<30>::add(this->val, other.val, t);
		return t;
  }
  }

  assert(false);
  return Int();
}

Int Int::operator+(const Int &&other) const {
  switch ((this->flag << 1) | other.flag) {
  case 0: {
    long long r = 0;

    if (safe_add(this->x, other.x, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    bint<30> *c = new bint<30>;

    bint<30>::add(this->x, other.x, c);

    return c;
  }

  case 1: {
    bint<30> *t = new bint<30>();
    bint<30>::add(this->x, other.val, t);
    return t;
  }

  case 2: {
    bint<30> *t = new bint<30>();
    bint<30>::add(this->val, other.x, t);
    return t;
  }

  case 3: {
    bint<30> *t = new bint<30>();
    bint<30>::add(this->val, other.val, t);
    return t;
  }
  }

  assert(false);
  return Int();
}

Int Int::operator+(const int z) const {
  switch (flag) {
  case 0: {
    long long r = 0;

    if (safe_add(this->x, z, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    bint<30> *c = new bint<30>;

    bint<30>::add(this->x, z, c);

    return c;
  }

  case 1: {
    bint<30> *t = new bint<30>();
    bint<30>::add(this->val, z, t);
    return t;
  }
  }

  assert(false);

  return Int();
}

Int Int::operator-(const Int &other) const {
  switch ((this->flag << 1) | other.flag) {
  case 0: {
    long long r = 0;

    if (safe_sub(this->x, other.x, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    bint<30> *c = new bint<30>;

    bint<30>::sub(this->x, other.x, c);

    return c;
  }

  case 1: {
    bint<30> *t = new bint<30>();
    bint<30>::sub(this->x, other.val, t);
    return t;
  }

  case 2: {
    bint<30> *t = new bint<30>();
    bint<30>::sub(this->val, other.x, t);
    return t;
  }

  case 3: {
    bint<30> *t = new bint<30>();
    bint<30>::sub(this->val, other.val, t);
    return t;
  }
  }

  assert(false);
  return Int();
}

Int Int::operator-(const Int &&other) const {
 switch ((this->flag << 1) | other.flag) {
  case 0: {
    long long r = 0;

    if (safe_sub(this->x, other.x, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    bint<30> *c = new bint<30>;

    bint<30>::sub(this->x, other.x, c);

    return c;
  }

  case 1: {
    bint<30> *t = new bint<30>();
    bint<30>::sub(this->x, other.val, t);
    return t;
  }

  case 2: {
    bint<30> *t = new bint<30>();
    bint<30>::sub(this->val, other.x, t);
    return t;
  }

  case 3: {
    bint<30> *t = new bint<30>();
    bint<30>::sub(this->val, other.val, t);
    return t;
  }
  }

  assert(false);
  return Int();
}

Int Int::operator-(const int z) const {
  switch (flag) {
  case 0: {
    long long r = 0;

    if (safe_sub(this->x, z, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    bint<30> *res = new bint<30>();
    bint<30>::sub(this->x, z, res);

    if (res->size == 2) {
      long long r;

      bint<30>::to_long(res, &r);

      delete res;

      return r;
    }

    return res;
  }
  case 1: {

    bint<30> *res = new bint<30>();

		bint<30>::sub(this->val, z, res);

    if (res->size == 2) {
      long long r;

      bint<30>::to_long(res, &r);

      delete res;

      return r;
    }

    return res;
  }
  }

  assert(false);

  return Int();
}

Int Int::operator*(const Int &other) const {
	switch((this->flag << 1) | other.flag) {
	case 0: {
    long long r = 0;

    if (safe_mul(this->x, other.x, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    bint<30> *c = new bint<30>();

    bint<30>::mul(this->x, other.x, c);

    return c;
	}

	case 1: {
		bint<30> *t = new bint<30>();
		bint<30>::mul(this->x, other.val, t);
    return t;
	}

	case 2: {
		bint<30> *t = new bint<30>();
		bint<30>::mul(this->val, other.x, t);
    return t;
	}
	case 3: {
		bint<30> *t = new bint<30>();
		bint<30>::mul(this->val, other.val, t);
		return t;
	}
	}

	assert(false);
  return Int();
}

Int Int::operator*(const Int &&other) const {
	switch((this->flag << 1) | other.flag) {
	case 0: {
    long long r = 0;

    if (safe_mul(this->x, other.x, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    bint<30> *a = bint<30>::from(this->x);
    bint<30> *b = bint<30>::from(other.x);
    bint<30> *c = new bint<30>();

    bint<30>::mul(a, b, c);

    delete a;
    delete b;

    return c;
	}

	case 1: {
		bint<30> *t = new bint<30>();
    bint<30>::mul(this->x, other.val, t);
    return t;
	}

	case 2: {
		bint<30> *t = new bint<30>();
		bint<30>::mul(this->val, other.x, t);
    return t;
	}
	case 3: {
		bint<30> *t = new bint<30>();
		bint<30>::mul(this->val, other.val, t);
		return t;
	}
	}

	assert(false);
  return Int();
}

Int Int::operator*(const int z) const {
	switch (flag) {
	case 0: {
    long long r = 0;
    if (safe_mul(this->x, z, &r) != LONG_LONG_OVERFLOW) {

      return r;
    }

		bint<30> *t = new bint<30>();
    bint<30>::mul(this->x, z, t);
    return t;
	}
	case 1: {
		bint<30> *res = new bint<30>();
		bint<30>::mul(this->val, z, res);
		return res;
	}
	}

	assert(false);
	return Int();
}

Int Int::operator/(const Int &other) const {
  if (!this->flag && this->x == 0) {
    return 0;
  }

  if (this->flag && this->val->size == 0) {
    return 0;
  }

  if (!this->flag && !other.flag) {
    long long r = 0;

    if (safe_div(this->x, other.x, &r) != LONG_LONG_OVERFLOW) {

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

      return r;
    }

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

      return r;
    }

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

      return r;
    }

    return quo;
  }

  bint<30>::div(this->val, other.val, quo, rem);

  delete rem;

  return quo;
}

Int Int::operator/(const Int &&other) const {
  if (!this->flag && !other.flag) {
    long long r = 0;

    if (safe_div(this->x, other.x, &r) != LONG_LONG_OVERFLOW) {
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

      return r;
    }

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

      return r;
    }

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

      return r;
    }

    return quo;
  }

  bint<30>::div(this->val, other.val, quo, rem);

  delete rem;

  return quo;
}

Int Int::operator/(const int z) const {
  if (!this->flag) {
    long long r = 0;

    if (safe_div(this->x, z, &r) != LONG_LONG_OVERFLOW) {
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

    return r;
  }

  return res;
}

Int Int::operator%(const Int &other) const {

  if (!this->flag && !other.flag) {

    // assert(Int(this->x % other.x).to_string() == ( Int2(this->x) %
    // Int2(other.x)).to_string());
		return  this->x % other.x;

    // bint<30> *a = bint<30>::from(this->x);
    // bint<30> *b = bint<30>::from(other.x);

    // bint<30> *quo = new bint<30>;
    // bint<30> *rem = new bint<30>;

    // bint<30>::div(a, b, quo, rem);

    // delete a;
    // delete b;
    // delete quo;

    // if (rem->size == 2) {
    //   bint<30>::to_long(rem, &r);

    //   delete rem;

    //   // assert(Int(r).to_string() == ( Int2(this->x) %
    //   // Int2(other.x)).to_string());
    //   return r;
    // }

    // return rem;
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

      return r;
    }

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

      return r;
    }

    return rem;
  }

  bint<30>::div(this->val, other.val, quo, rem);

  delete quo;

  return rem;
}

Int Int::operator%(const Int &&other) const {
  if (!this->flag && !other.flag) {
		return  this->x % other.x;
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

      return r;
    }

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

      return r;
    }

    return rem;
  }

  bint<30>::div(this->val, other.val, quo, rem);

  delete quo;

  return rem;
}

Int Int::operator%(const int z) const {
  if (!this->flag) {
    return this->x % z;
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
    return r;
  }

  return rem;
}

Int Int::operator+() { return *this; }

Int Int::operator-() { return *this * -1; }

Int Int::operator++() {

  *this = *this + 1;

  // assert(this->to_string() == ( Int2(this->val->copy()) + 1).to_string());

  return *this;
}

Int Int::operator--() {
  *this = *this - 1;
  return *this;
}

Int Int::operator++(int) {
  Int tmp = *this;

  *this = *this + 1;

  return tmp;
}

Int Int::operator--(int) {
  Int tmp = *this;

  *this = *this - 1;

  return tmp;
}

bool Int::operator==(const Int &other) const {
  if (!this->flag && !other.flag) {
    return this->x == other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) == 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) == 0;
  }

	return bint<30>::compare(this->val, other.val) == 0;
}

bool Int::operator==(const Int &&other) const {
  if (!this->flag && !other.flag) {
    return this->x == other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) == 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) == 0;
  }

  return bint<30>::compare(this->val, other.val) == 0;
}

bool Int::operator<(const Int &other) const {
  if (!this->flag && !other.flag) {
    return this->x < other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) < 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) < 0;
  }

  return bint<30>::compare(this->val, other.val) < 0;
}

bool Int::operator<(const Int &&other) const {
  if (!this->flag && !other.flag) {
    return this->x < other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) < 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) < 0;
  }

  return bint<30>::compare(this->val, other.val) < 0;
}

bool Int::operator<=(const Int &other) const {
  if (!this->flag && !other.flag) {
    return this->x <= other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) <= 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) <= 0;
  }

  return bint<30>::compare(this->val, other.val) <= 0;
}

bool Int::operator<=(const Int &&other) const {
  if (!this->flag && !other.flag) {
    return this->x <= other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) <= 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) <= 0;
  }

  return bint<30>::compare(this->val, other.val) <= 0;
}

bool Int::operator>(const Int &other) const {
  if (!this->flag && !other.flag) {
    return this->x > other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) > 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) > 0;
  }

  return bint<30>::compare(this->val, other.val) > 0;
}

bool Int::operator>(const Int &&other) const {
  if (!this->flag && !other.flag) {
    return this->x > other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) > 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) > 0;
  }

  return bint<30>::compare(this->val, other.val) > 0;
}

bool Int::operator>=(const Int &other) const {
  if (!this->flag && !other.flag) {
    return this->x >= other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) >= 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) >= 0;
  }

  return bint<30>::compare(this->val, other.val) >= 0;
}

bool Int::operator>=(const Int &&other) const {
  if (!this->flag && !other.flag) {
    return this->x >= other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) >= 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) >= 0;
  }

  return bint<30>::compare(this->val, other.val) >= 0;
}

bool Int::operator!=(const Int &other) const {
  if (!this->flag && !other.flag) {
    return this->x != other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) != 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) != 0;
  }

  return bint<30>::compare(this->val, other.val) != 0;
};

bool Int::operator!=(const Int &&other) const {

  if (!this->flag && !other.flag) {
    return this->x != other.x;
  }

  if (!this->flag && other.flag) {
    return bint<30>::compare(this->x, other.val) != 0;
  }

  if (this->flag && !other.flag) {
    return bint<30>::compare(this->val, other.x) != 0;
  }

  return bint<30>::compare(this->val, other.val) != 0;
};

Int Int::operator=(const Int &other) {

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

Int Int::operator=(Int &&other) {
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

Int Int::ceil_log2() {
  if (!flag) {
    return std::ceil(std::log2(x));
  }
  return bint<30>::ceil_log2(this->val);
}

long long Int::longValue() {
  if (!flag)
    return x;

  long long v = 0;

  if (bint<30>::to_long(this->val, &v) == -1) {
		raise(error(ErrorCode::LONG_LONG_OVERFLOW, 0));
  }

  return v;
}

double Int::doubleValue() {
  if (!flag)
    return x;

  double v = 0.0;

	if(bint<30>::to_double(this->val, &v) == -1) {
		raise(error(ErrorCode::DOUBLE_OVERFLOW, 0));
	}

	return v;
}

Int gcd(const Int &a, const Int &b) {
  if (!a.flag && !b.flag) {
    long long z = 0;

    if (safe_gcd(a.x, b.x, &z) != LONG_LONG_OVERFLOW) {
      return z;
    }
  }

  if (!a.flag && b.flag) {
    bint<30> *k = bint<30>::from(a.x);

    bint<30> *r = bint<30>::gcd(k, b.val);

    // TODO: convert r to long long if it have size 2

    delete k;

    return r;
  }
  if (a.flag && !b.flag) {
    bint<30> *k = bint<30>::from(b.x);
    bint<30> *r = bint<30>::gcd(a.val, k);
    delete k;

    return r;
  }

  bint<30> *g = bint<30>::gcd(a.val, b.val);

  return g;
}

Int gcd(const Int &&a, const Int &&b) {
  if (!a.flag && !b.flag) {

    long long z = 0;
    if (safe_gcd(a.x, b.x, &z) != LONG_LONG_OVERFLOW) {
      return z;
    }
  }

  if (!a.flag && b.flag) {
    bint<30> *k = bint<30>::from(a.x);

    bint<30> *r = bint<30>::gcd(k, b.val);

    // TODO: convert r to long long if it have size 2

    delete k;

    return r;
  }
  if (a.flag && !b.flag) {
    bint<30> *k = bint<30>::from(b.x);
    bint<30> *r = bint<30>::gcd(a.val, k);
    delete k;

    return r;
  }

  bint<30> *g = bint<30>::gcd(a.val, b.val);

  return g;
}

Int lcm(const Int &a, const Int &b) {
  if (!a.flag && !b.flag) {
    long long z = 0;

    if (safe_lcm(a.x, b.x, &z) != LONG_LONG_OVERFLOW) {

      return z;
    }
  }

  if (!a.flag && b.flag) {
    bint<30> *k = bint<30>::from(a.x);

    bint<30> *r = bint<30>::lcm(k, b.val);

    // TODO: convert r to long long if it have size 2

    delete k;

    return r;
  }
  if (a.flag && !b.flag) {
    bint<30> *k = bint<30>::from(b.x);
    bint<30> *r = bint<30>::lcm(a.val, k);
    delete k;

    return r;
  }

  bint<30> *g = bint<30>::lcm(a.val, b.val);

  return g;
}

Int lcm(const Int &&a, const Int &&b) {
  if (!a.flag && !b.flag) {
    long long z = 0;
    if (safe_lcm(a.x, b.x, &z) != LONG_LONG_OVERFLOW) {
      return z;
    }
  }

  if (!a.flag && b.flag) {
    bint<30> *k = bint<30>::from(a.x);

    bint<30> *r = bint<30>::lcm(k, b.val->copy());

    // TODO: convert r to long long if it have size 2

    delete k;

    return r;
  }
  if (a.flag && !b.flag) {
    bint<30> *k = bint<30>::from(b.x);
    bint<30> *r = bint<30>::lcm(a.val, k);
    delete k;

    return r;
  }

  bint<30> *g = bint<30>::lcm(a.val, b.val);

  return g;
}

void Int::operator/=(Int v) { *this = *this / v; }
void Int::operator+=(Int v) { *this = *this + v; }
void Int::operator-=(Int v) { *this = *this - v; }
void Int::operator*=(Int v) { *this = *this * v; }

Int abs(const Int &&a) {
  if (!a.flag)
    return std::abs(a.x);
  return bint<30>::abs(a.val);
}

Int abs(const Int &a) {
  if (!a.flag)
    return std::abs(a.x);
  return bint<30>::abs(a.val);
}

Int fact(const Int &&a) {
  if (!a.flag) {
    long long z = 0;

    if (safe_fact(a.x, &z) != LONG_LONG_OVERFLOW) {
      return z;
    }

    bint<30> *x = bint<30>::from(a.x);

    return x;
  }

  bint<30> *g = bint<30>::fact(a.val);

  return g;
}

Int fact(const Int &a) {
  if (!a.flag) {
    long long z = 0;

    if (safe_fact(a.x, &z) != LONG_LONG_OVERFLOW) {
      return z;
    }

    bint<30> *x = bint<30>::from(a.x);

    return x;
  }

  bint<30> *g = bint<30>::fact(a.val);

  return g;
}

Int max(const Int &&a, Int &&b) {
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

Int max(const Int &a, Int &b) {
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

Int max(const Int &a, Int &&b) {

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

Int min(const Int &&a, const Int &&b) {
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

Int min(const Int &a, const Int &b) {
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

Int pow(const Int &&a, const Int &&b) {
  if (!a.flag && !b.flag) {
    long long z = 0;
    if (safe_pow(a.x, b.x, &z) != LONG_LONG_OVERFLOW) {
      return z;
    }

    bint<30> *x = bint<30>::from(a.x);
    bint<30> *y = bint<30>::from(b.x);
    bint<30> *w = bint<30>::pow(x, y);

    delete x;
    delete y;

    return w;
  }

  if (!a.flag && b.flag) {
    bint<30> *k = bint<30>::from(a.x);
    bint<30> *r = bint<30>::pow(k, b.val);

    delete k;

    return r;
  }

  if (a.flag && !b.flag) {
    bint<30> *k = bint<30>::from(b.x);
    bint<30> *r = bint<30>::pow(a.val, k);

    delete k;

    return r;
  }

  bint<30> *r = bint<30>::pow(a.val, b.val);

  return r;
}

Int pow(const Int &a, const Int &b) {
  if (!a.flag && !b.flag) {
    long long z = 0;

    if (safe_pow(a.x, b.x, &z) != LONG_LONG_OVERFLOW) {
      return z;
    }

    bint<30> *x = bint<30>::from(a.x);
    bint<30> *y = bint<30>::from(b.x);
    bint<30> *w = bint<30>::pow(x, y);

    delete x;
    delete y;

    return w;
  }
  if (!a.flag && b.flag) {

    bint<30> *k = bint<30>::from(a.x);
    bint<30> *r = bint<30>::pow(k, b.val);

    delete k;

    return r;
  }

  if (a.flag && !b.flag) {
    bint<30> *k = bint<30>::from(b.x);
    bint<30> *r = bint<30>::pow(a.val, k);

    delete k;

    return r;
  }

  bint<30> *r = bint<30>::pow(a.val, b.val);

  return r;
}

double pow(const Int &&a, const double b) {
  if (!a.flag) {
    return std::pow(a.x, b);
  }

  return bint<30>::pow(a.val, b);
}

double pow(const Int &a, const double b) {
  if (!a.flag) {
    return std::pow(a.x, b);
  }

  return bint<30>::pow(a.val, b);
}

Int isqrt(const Int &a) {
  if (!a.flag) {
    long long r = 0;

    if (safe_isqrt(a.x, &r) != LONG_LONG_OVERFLOW) {
      return r;
    }

    return r;
  }

  bint<30> *res = new bint<30>();

  bint<30>::isqrt(a.val, res, 0);

  // TODO: convert to long long if the size of r is <= 2

  return res;
}

bool operator<(const unsigned int a, const Int &v) {
  if (!v.flag)
    return a < v.x;

  bool res = bint<30>::compare(a, v.val) < 0;

  return res;
}

Int operator*(const int a, const Int &v) {
	return v * a;
}

bool operator>(const unsigned int a, const Int &v) {
  if (!v.flag)
    return a > v.x;

  bool res = bint<30>::compare(a, v.val) > 0;

  return res;
}

bool operator<=(const unsigned int a, const Int &v) {
  if (!v.flag)
    return a <= v.x;

  bool res = bint<30>::compare(a, v.val) <= 0;

  return res;
}

bool operator>=(const unsigned int a, const Int &v) {
  if (!v.flag)
    return a >= v.x;

  bool res = bint<30>::compare(a, v.val) >= 0;

  return res;
}

Int operator+(const long long a, const Int &v) {
  if (!v.flag)
    return a + v.x;

  bint<30> *res = new bint<30>();

  bint<30>::add(a, v.val, res);

  return res;
}

Int operator-(const long long a, const Int &v) {
  if (!v.flag)
    return a - v.x;

  bint<30> *res = new bint<30>();

  bint<30>::sub(a, v.val, res);

  return res;
}
