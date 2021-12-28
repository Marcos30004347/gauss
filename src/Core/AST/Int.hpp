#ifndef INT_HPP
#define INT_HPP

// references:
// https://cacr.uwaterloo.ca/hac/about/chap14.pdf
// https://github.com/python/cpython/blob/main/Objects/longobject.c
// The Art of Computer Programming Vol 2 by Donald E. Knuth

// TODO: better error handling

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <float.h>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <limits>
#include <math.h>
#include <sstream>
#include <vector>

#define pow2(e) (1 << e)

// computes n mod 2^e
#define modPow2(n, e) (n & ((1 << e) - 1))

// computes n / 2^e
#define quoPow2(n, e) (n >> e)

#define mulPow2(d, exp) (d << exp)

#define SUB(A, B, C)                                                           \
  {                                                                            \
    bint_t *tmp = new bint_t();                                                \
    sub(A, B, tmp);                                                            \
    delete C;                                                                  \
    C = tmp;                                                                   \
  }
#define MUL(x, y, z)                                                           \
  {                                                                            \
    bint_t *tmp = new bint_t();                                                \
    mul(x, y, tmp);                                                            \
    delete z;                                                                  \
    z = tmp;                                                                   \
  }

#if DBL_MANT_DIG == 53
#define EXP2_DBL_MANT_DIG 9007199254740992.0
#else
#define EXP2_DBL_MANT_DIG (ldexp(1.0, DBL_MANT_DIG))
#endif

inline int high_bit(uint32_t x) {
  const int blen[32] = {0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
                     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

  int msb = 0;

  while (x >= 32) {
    msb += 6;
    x >>= 6;
  }

  return msb + blen[x];
}


inline int ull_ceil_log2(unsigned long long x)
{
  static const unsigned long long t[6] = {
    0xFFFFFFFF00000000ull,
    0x00000000FFFF0000ull,
    0x000000000000FF00ull,
    0x00000000000000F0ull,
    0x000000000000000Cull,
    0x0000000000000002ull
  };

  int y = (((x & (x - 1)) == 0) ? 0 : 1);
  int j = 32;
  int i;

  for (i = 0; i < 6; i++) {
    int k = (((x & t[i]) == 0) ? 0 : j);
    y += k;
    x >>= k;
    j >>= 1;
  }

  return y;
}

inline int msb(unsigned int v) {
  static const int pos[32] = {0, 1, 28, 2, 29, 14, 24, 3,
    30, 22, 20, 15, 25, 17, 4, 8, 31, 27, 13, 23, 21, 19,
    16, 7, 26, 12, 18, 6, 11, 5, 10, 9};
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v = (v >> 1) + 1;
  return pos[(v * 0x077CB531UL) >> 27];
}

#define arith_rshift(I, J) ((I) < 0 ? -1 - ((-1 - (I)) >> (J)) : (I) >> (J))

template <
	int exp = 30,
	typename single_type = uint32_t,
	typename double_type = uint64_t,
	typename ssingle_type = int32_t,
	typename sdouble_type = int64_t
	>
class bint {
	static_assert(exp <= sizeof(single_type) * CHAR_BIT - 1, "base expoent needs to be smaller than the number of bits in the single_type");
	static_assert(sizeof(single_type) <= 32/CHAR_BIT, "single type can have at most 32 bits");
	static_assert(sizeof(double_type) >= 2 * sizeof(single_type),
                "sizeof(double_type) needs to be at least twice"
                "as big as sizeof(single_type)");
  // digit2_t is a type capable of holding
  // at least two elements of type single_type

  static const single_type base = ((single_type)1 << exp);
	static const single_type mask = ((single_type)(base - 1));

public:
  // digit one is the type used to store every digit of the number base 2^exp
  using digit_t = single_type;
  using digit2_t = double_type;
	using sdigit_t = ssingle_type;
	using sdigit2_t = sdouble_type;

	using bint_t = bint<exp, digit_t, digit2_t, sdigit_t, sdigit2_t>;

  digit_t *digit;
  size_t size;
  short sign;

  bint(digit_t *d, size_t s, short sign = 1) : digit{d}, size{s}, sign{sign} {}
  bint() : digit{nullptr}, size{0}, sign{1} {}

  ~bint() {
    if (digit)
      delete[] digit;
  }

  bint_t *copy() {
    bint_t *t = new bint_t();

    t->size = this->size;
    t->digit = new digit_t[this->size];
    t->sign = this->sign;

    memcpy(t->digit, this->digit, this->size * sizeof(digit_t));

    return t;
  }

  void resize(uint64_t s) {
    size = s;

    if (digit)
      delete digit;

    if (s == 0) {
      digit = nullptr;
      return;
    }

    digit = new digit_t[s];

		memset(digit, 0, sizeof(digit_t)*s);
  }

  // shift the bits in a to the left m times and save the
  // result on z
  static digit_t digits_lshift(digit_t *a, size_t l, int d,
                               digit_t *z) {
    digit_t carry = 0;

    assert(0 <= d && (long unsigned int)d < CHAR_BIT*sizeof(digit_t));

    for (size_t i = 0; i < l; i++) {
      // shift digits and combine them
      // with the carry
      digit2_t acc = (digit2_t)a[i] << d | carry;
      // save the exp bits
      z[i] = (digit_t)acc & mask;
      // take the remaining bits and
      // save it to the carry
      carry = (digit_t)(acc >> exp);
    }

    return carry;
  }

  // shift the bits in a to the left m times and save the
  // result on z
  static digit_t digits_rshift(digit_t *a, size_t l, int d,
                               digit_t *z) {
    digit_t carry = 0;

    // mask with last m bits set
    digit_t msk = ((digit_t)1 << d) - 1U;

    assert(0 <= d && (long unsigned int)d < CHAR_BIT*sizeof(digit_t));

    for (size_t i = l; i-- > 0;) {
      digit2_t acc = (digit2_t)carry << exp | a[i];
      carry = (digit_t)acc & msk;
      z[i] = (digit_t)(acc >> d);
    }

    return carry;
  }

  void trim() {
    if (size == 0)
      return;

    size_t k = size - 1;

    while (k > 0 && !digit[k])
      k--;

    if (!digit[k]) {
      delete digit;

      size = 0;
      digit = nullptr;
    } else {
      size = k + 1;
      digit = (digit_t *)realloc(digit, sizeof(digit_t) * size);
    }
  }

  // convert x to base 2^exp
  template <typename T> static bint_t* from(T x) {
    short sign = 1;

    if(x == 0) return new bint_t();

    if (x < 0) {
      sign = -1;

      x = -x;
    }

    size_t s = 10;

    digit_t *v = new digit_t[s];

    size_t i = 0;

    T q = quoPow2(x, exp);
    v[0] = modPow2(x, exp);

    while (q > 0) {
      i = i + 1;

      if (i >= s) {
        s = s * 10;
        v = (digit_t *)realloc(v, sizeof(digit_t) * s);
      }

      x = q;

      q = quoPow2(x, exp);
      v[i] = modPow2(x, exp);
    }

		v = (digit_t*)realloc(v, sizeof(digit_t) * i + 1);

    return new bint(v, i + 1, sign);
  }

  static bint_t* from(double y) {
		short sign = 1;
		double x = 0;

		y = std::modf(y, &x);

		if(x == 0) {
			return new bint_t();
		}

		double b = std::pow(2, exp);

		if (x < 0) {
      sign = -1;
      x = -x;
    }

    size_t s = 10;

    digit_t *v = new digit_t[s];

    size_t i = 0;

    v[0] = fmod(x, b);

    double q = (x - v[0]) / b;

    while (q > 0) {
      i = i + 1;

      if (i >= s) {
        s = s * 10;
        v = (digit_t *)realloc(v, sizeof(digit_t) * s);
      }

      x = q;

      v[i] = fmod(x, b);
      q = (x - v[i]) / b;
    }

    return new bint(v, i + 1, sign);
	}
	// sum the absolute values of the big integers with digits x[0...a]
  // and y[0...b] and save in z[0...a + 1]. It's assumed that a >= b.
  // All the memory should be pre-allocated before execution.
  static void abs_add_digits(bint_t *x, bint_t *y, bint_t *z) {
    digit_t carry = 0;

    size_t i = 0;

    size_t a = x->size;
    size_t b = y->size;

    if (a < b) {
      bint_t *t = x;
      x = y;
      y = t;

      size_t c = a;
      a = b;
      b = c;
    }

    z->resize(x->size + 1);

    for (; i < b; ++i) {
      carry += x->digit[i] + y->digit[i];
      z->digit[i] = (digit_t)(carry & mask);
      carry >>= exp;
    }

    for (; i < a; ++i) {
      carry += x->digit[i];
      z->digit[i] = (digit_t)(carry & mask);
      carry >>= exp;
    }

    z->digit[i] = carry;

    return z->trim();
  }

  // subtract the absolute values of the big integers with digits x[0...a]
  // and y[0...b] and save in z[0...a]. It's assumed that a >= b.
  // All the memory should be pre-allocated before execution.
  static void abs_sub_digits(bint_t *x, bint_t *y, bint_t *z) {
    short sign = 1;

		size_t a = x->size;
    size_t b = y->size;

		// if(a == 0) {
		// 	z->resize(b);
		// 	z->sign = -y->sign;
		// 	memcpy(z->digit, y->digit, sizeof(digit_t)*b);
		// 	return;
		// }

		// if(b == 0) {
		// 	z->resize(a);
		// 	z->sign = x->sign;
		// 	memcpy(z->digit, x->digit, sizeof(digit_t)*a);
		// 	return;
		// }
    size_t i;

    if (b > a) {
      sign = -1;

      bint_t *t = x;
      x = y;
      y = t;

      size_t c = a;
      a = b;
      b = c;
    } else if (a == b) {
      // Get the index of the digit that x and y differ
      long long i = (long long)a;

      while (--i >= 0 && x->digit[i] == y->digit[i]){}

			if(i < 0) {
				z->resize(0);
				z->sign = 1;
				return;
			}

      // if (i == 0 && x->digit[i] == y->digit[i]) {
      //   z->resize(0);
      //   z->sign = 1;
			// 	printf("B\n");
      //   return;
      // }

      if (x->digit[i] < y->digit[i]) {
        bint_t *t = x;
        x = y;
        y = t;

        size_t c = a;
        a = b;
        b = c;

        sign = -1;
      }

      a = b = i + 1;
    }

		digit_t borrow = 0;

    z->resize(x->size + 1);

		for (i = 0; i < b; ++i) {
      borrow = x->digit[i] - y->digit[i] - borrow;
      z->digit[i] = borrow & mask;
      borrow >>= exp;
      borrow &= 1;
    }

    for (; i < a; ++i) {
      borrow = x->digit[i] - borrow;
      z->digit[i] = borrow & mask;
      borrow >>= exp;
      borrow &= 1;
    }

    z->sign = sign;

    return z->trim();
  }

  // Square the bint of digits x[0...a] and store the
  // result in w[0...2*a]
  static void abs_square_digits(bint_t *x, bint_t *w) {
    size_t i = 0;

    w->resize(2 * x->size);

    for (; i < 2 * x->size; i++) {
      w->digit[i] = 0;
    }

    for (i = 0; i < x->size; i++) {
      digit2_t carry;
      digit2_t xi = x->digit[i];

      digit_t *pw = w->digit + (i << 1);
      digit_t *px = x->digit + (i + 1);
      digit_t *pe = x->digit + x->size;

      carry = *pw + xi * xi;
      *pw++ = (digit_t)(carry & mask);
      carry = quoPow2(carry, exp);

      xi <<= 1;

      while (px != pe) {
        carry += *pw + *px++ * xi;
        *pw++ = (digit_t)(carry & mask);
        carry = quoPow2(carry, exp);
        assert(carry <= (mask << 1));
      }

      if (carry) {
        carry += *pw;
        *pw++ = (digit_t)(carry & mask);
        carry = quoPow2(carry, exp);
      }

      if (carry) {
        *pw += (digit_t)(carry & mask);
      }

      assert((carry >> exp) == 0);
    }

    return w->trim();
  }

  // multiply the big integer with digits in x[0...a] with y[0...b]
  // and store the result in z[a + b]. It is assumed that a >= b.
  // All the space should be pre-alocated before execution
  static void abs_mul_digits(bint_t *x, bint_t *y, bint_t *z) {
    if (x->size == 0 || y->size == 0) {
      return z->resize(0);
    }

    z->resize(x->size + y->size);

    for (size_t i = 0; i < x->size; i++) {
      digit2_t carry = 0;
      digit2_t xi = x->digit[i];

      digit_t *pz = z->digit + i;
      digit_t *py = y->digit;

      digit_t *pe = y->digit + y->size;

      while (py < pe) {
        carry += *pz + *py++ * xi;
        *pz++ = (digit_t)(carry & mask);
        carry >>= exp;
        assert(carry <= (mask << 1));
      }
      if (carry) {
        *pz += (digit_t)(carry & mask);
      }

      assert((carry >> exp) == 0);
    }

    return z->trim();
  }

  static void add(bint_t *x, bint_t *y, bint_t *z) {
		if(x->size <= 1 && y->size <= 1) {
			digit_t a = x->size > 0 ? x->digit[0] : 0;
			digit_t b = y->size > 0 ? y->digit[0] : 0;

			int64_t c = (int64_t)a * x->sign + (int64_t)b * y->sign;

			short s = 1;

			if(c < 0) {
				s = -1;
				c = -c;
			}

			digit2_t v = (digit2_t)c;
			digit_t A = (digit_t)v & mask;
			digit_t B = (digit_t)(v >> exp) & mask;

			if(B) {
				z->resize(2);
				z->sign = s;
				z->digit[0] = A;
				z->digit[1] = B;
			} else if(A) {
				z->resize(1);
				z->sign = s;
				z->digit[0] = A;
			} else {
				z->resize(0);
			}

			return;
		}

		if (x->sign < 0) {
      if (y->sign < 0) {
        z->sign = -1 * z->sign;
        return abs_add_digits(x, y, z);
      }

      return abs_sub_digits(y, x, z);
    }

    if (y->sign < 0) {
      return abs_sub_digits(x, y, z);
    }

    return abs_add_digits(x, y, z);
  }

  static void sub(bint_t *x, bint_t *y, bint_t *z) {
		if (x->sign < 0) {
      if (y->sign < 0) {
        return abs_sub_digits(y, x, z);
      }

			abs_add_digits(x, y, z);

			if (z->size) {
				z->sign = -1 * z->sign;
			}

			return;
    }

    if (y->sign < 0) {
      return abs_add_digits(x, y, z);
    }

    return abs_sub_digits(x, y, z);
  }

  static void mul(bint_t *x, bint_t *y, bint_t *z) {
    z->sign = x->sign * y->sign;
    abs_mul_digits(x, y, z);

    return z->trim();
  }

  static short fast_mod(bint_t *x, bint_t *y, bint_t *rem) {
    digit_t left = x->digit[0];
    digit_t righ = y->digit[0];

    rem->resize(1);
    rem->sign = x->sign * y->sign;

		rem->digit[0] = left % righ;

    return 1;
  }

  static short fast_div(bint_t *x, bint_t *y, bint_t *quo) {
		digit_t left = x->digit[0];
    digit_t righ = y->digit[0];

    quo->resize(1);
    quo->sign = x->sign * y->sign;
		quo->digit[0] = left / righ;

    return 1;
  }

  static short div(bint_t *x, bint_t *y, bint_t *quo, bint_t *rem) {
    // Following The Art of Computer Programming, Vol.2, section 4.3.1,
    // Algorithm D.
    size_t m = x->size;
    size_t n = y->size;
		//printf("HAHAHA\n");

		//printf("%s\n", x->to_string().c_str());
		//printf("%s\n", y->to_string().c_str());

		if (m == 0) {
      quo->resize(0);
      rem->resize(0);

			//printf("return\n");
      return 1;
    }

    if (n == 0) {
			// TODO: throw division by zero error
			//printf("return\n");
      return 0;
    }

		if(m < n || (m == n && x->digit[m - 1] < y->digit[n - 1]) ) {
			quo->resize(0);
			quo->sign = 1;

			rem->resize(x->size);
			rem->sign = x->sign;

			for(size_t i = 0; i < x->size; i++)
				rem->digit[i] = x->digit[i];

			//printf("return\n");
			return 1;
		}

    if (m == 1 && n == 1) {
      fast_div(x, y, quo);
			quo->trim();

      if(rem) {
				fast_mod(x, y, rem);
				rem->trim();
			}

			//printf("return\n");
			return 1;
    }

		quo->sign = x->sign * y->sign;
		rem->sign = x->sign * y->sign;

		if (n < 2) {
      digit2_t r = 0;

      digit_t k = y->digit[0];

      quo->resize(m);

      digit_t *in = x->digit + m;
      digit_t *ou = quo->digit + m;
      long long t = m;

      while (--t >= 0) {
        digit_t hi;
        r = (r << exp) | *--in;
        *--ou = hi = (digit_t)(r / k);
        r -= (digit2_t)hi * k;
      }

      if(rem) {
				rem->resize(1);
				rem->digit[0] = r;
			}

      quo->trim();

			if(rem) rem->trim();

			//printf("return\n");
      return 1;
    }

    digit_t carry = 0;

    // D1: Normalization

    // because we are using binary base 2^exp
    // this operatios return the expoent of
    // 2^k such that v*2^k have v[n - 1] >= floor(base / 2)
    int d = exp - high_bit(y->digit[n - 1]);

    bint_t u;
    bint_t v;

    u.resize(m + 1);
    v.resize(n);

    // sinse we are using binary base, shifting
    // the value left is the same as multiplying
    // by 2^d
    // Multiply y by 2^d
    carry = digits_lshift(y->digit, y->size, d, v.digit);

    assert(carry == 0);

    // Multiply x by 2^d
    carry = digits_lshift(x->digit, x->size, d, u.digit);

		if (carry != 0 || u.digit[m - 1] >= v.digit[n - 1]) {
      u.digit[m] = carry;
      m = m + 1;
    }

    // D2 - loop from D3 to D7
    digit_t v1 = v.digit[n - 1];
    digit_t v2 = v.digit[n - 2];

    long long j = m - n;

		//printf("---> %lli\n", j);

		assert(j >= 0);
    quo->resize(j);

    digit_t *u0 = u.digit;
    digit_t *v0 = v.digit;

    digit_t *uj = u0 + j;
    digit_t *qj = quo->digit + j;

    while (uj-- > u0) {
      // D3
      digit_t ut = uj[n];

			assert(ut <= v1);
      // uu = (u[j + n]*b + u[j + n - 1])
      digit2_t uu = ((digit2_t)ut << exp) | uj[n - 1];

      digit_t q = (digit_t)(uu / v1);
      digit_t r = (digit_t)(uu - (digit2_t)v1 * q);

      // test q >= b or q*v[n - 2] > b*r + u[j + n - 2]
      // r << exp | uj[n - 2] = r*b + u[j + n - 2]
			//q >= base
			while ((digit2_t)v2 * q > (((digit2_t)r << exp) | uj[n - 2])) {
        q = q - 1;
        r = r + v1;

        if (r >= base)
          break;
      }

			assert(q <= base);

      // D4 replace (u[j + n], u[j + n - 1],...,u[j])b
      // by (u[j + n], u[j + n - 1],...,u[j])b - q*(0,v[n-1],...,v[1],v[0])b
      sdigit_t borrow = 0;

      for (size_t i = 0; i < n; ++i) {
        sdigit2_t z = (sdigit_t)uj[i] + borrow - (sdigit2_t)q * (sdigit2_t)v0[i];
        uj[i] = (digit_t)z & mask;
        borrow = (sdigit_t)arith_rshift(z, exp);
      }

      // D5
      assert((sdigit_t)ut + borrow == -1 || (sdigit_t)ut + borrow == 0);

      // D6
      if ((sdigit_t)ut + borrow < 0) {
				digit_t carry = 0;

				for (size_t i = 0; i < n; ++i) {
          carry += uj[i] + v0[i];
          uj[i] = carry & mask;
          carry >>= exp;
        }

        q = q - 1;
      }

      assert(q < base);
      *--qj = q;
    }

    rem->resize(n);

    carry = digits_rshift(u.digit, n, d, rem->digit);

    assert(carry == 0);

    quo->trim();
    rem->trim();

		//printf("return\n");
    return 1;
  }

  // TODO: to string of the number in base 10
  // TODO: add construction from strings and const char* types

  static short compare(bint_t *v0, bint_t *v1) {
		if(v0->size == 0 && v1->size == 0) return 0;

		if (v0->sign != v1->sign)
      return v0->sign > v1->sign ? 1 : -1;
    if (v0->size != v1->size)
      return v0->size > v1->size ? 1 : -1;

    size_t i = v0->size - 1;

    while (i > 0 && v0->digit[i] == v1->digit[i])
      --i;

    if (i == 0 && v0->digit[0] == v1->digit[0])
      return 0;

    return v0->digit[i] > v1->digit[i] ? 1 : -1;
  }

  static double frexp(bint_t *a, size_t *e, bool *overflow) {
    // from https://github.com/python/cpython/blob/main/Objects/longobject.c
    size_t m = a->size;

    if (m == 0) {
      *e = 0;
      return 0;
    }

    double dx;

    digit_t *x_digits = new digit_t[2 + (DBL_MANT_DIG + 1) / exp]{0};

    static const int half_even_correction[8] = {0, -1, -2, 1, 0, -1, 2, 1};

    size_t al = high_bit(a->digit[m - 1]);
    size_t mx = std::numeric_limits<size_t>::max();

    if (m >= (mx - 1) / exp + 1 &&
        (m > (mx - 1) / exp + 1 || al > (mx - 1) % exp + 1)) {
      *overflow = true;
      return -1;
    }

    al = (m - 1) * exp + al;

    size_t x_size = 0;

    if (al <= DBL_MANT_DIG + 2) {
      size_t shift_d = (DBL_MANT_DIG + 2 - al) / exp;
      size_t shift_b = (DBL_MANT_DIG + 2 - al) % exp;

      x_size = shift_d;
      digit_t rem = digits_lshift(a->digit, m, shift_b, x_digits + x_size);

      x_size += m;

      x_digits[x_size++] = rem;
    } else {
      size_t shift_d = (al - DBL_MANT_DIG - 2) / exp;
      size_t shift_b = (al - DBL_MANT_DIG - 2) % exp;
      digit_t rem =
          digits_rshift(a->digit + shift_d, m - shift_d, shift_b, x_digits);

      x_size = m - shift_d;

      if (rem) {
        x_digits[0] |= 1;
      } else {
        while (shift_d > 0) {
          if (a->digit[--shift_d]) {
            x_digits[0] |= 1;
            break;
          }
        }
      }
    }

    x_digits[0] += half_even_correction[x_digits[0] & 7];

    dx = x_digits[--x_size];

    while (x_size > 0) {
      dx = dx * base + x_digits[--x_size];
    }

    dx /= 4.0 * EXP2_DBL_MANT_DIG;

    if (dx == 1.0) {
      if (al == mx) {
        *overflow = true;
        return -1;
      }
      dx = 0.5;
      al += 1;
    }

    *e = al;

    *overflow = false;
    return a->sign < 0 ? -dx : dx;
  }

  static short to_double(bint_t *b, double* x) {
    if (b->size == 0) {
			*x = 0;
			return 1;
		}

    if (b->size == 1) {
			*x = (double)b->digit[0] * b->sign;
      return 1;
		}

    if (exp * b->size <= CHAR_BIT * sizeof(unsigned long long)) {
      digit_t *z = b->digit;

      unsigned long long v = 0;

			for (size_t i = 0; i < b->size - 1; i += 2) {
        digit2_t x = (digit2_t)z[i + 1] << exp;
        digit2_t y = (digit2_t)z[i];
        v |= (unsigned long long)((digit2_t)(x | y)) << (i * exp);
      }

			*x = ((double)v) * b->sign;
      return 1;
    }

    bool overflow = false;

    size_t e = 0.0;

    *x = frexp(b, &e, &overflow);

		if ((*x == -1 && overflow) || e > DBL_MAX_EXP) {
			return -1;
    }

		*x = std::ldexp(*x, e) * b->sign;

		return 1;
  }

  static short to_long(bint_t *b, long long* v) {
    if (b->size == 0) {
			*v = 0;
			return 1;
		}

    if (b->size == 1) {
			*v = (long long)b->digit[0] * b->sign;
			return 1;
		}
    if (exp * b->size <= CHAR_BIT * sizeof(long long)) {
      digit_t *z = b->digit;

      *v = 0;

      for (size_t i = 0; i < b->size - 1; i += 2) {
        digit2_t x = (digit2_t)z[i + 1] << exp;
        digit2_t y = (digit2_t)z[i];
        *v |= (long long)((digit2_t)(x | y)) << (i * exp);
      }

			*v = *v * b->sign;

			return 1;
    }

    // overflow
    return -1;
  }

	static bint_t *ceil_log2(bint_t* a) {
		bint_t* z = bint_t::from(0);
		bint_t* x = bint_t::from(0);

		bint_t* e = bint_t::from(exp);
		bint_t* b = bint_t::from(a->size - 1);

		bint_t* w = bint::from(ull_ceil_log2(a->digit[a->size - 1]));

		mul(e, b, z);
		add(z, w, x);

		delete e;
		delete b;
		delete z;
		delete w;

		return x;
	}

  static bint_t *abs(bint_t *a) {
    bint_t *b = a->copy();
    b->sign = 1;
    return b;
  }

  static bint_t *max(bint_t *a, bint_t *b) {
    return compare(a, b) > 0 ? a->copy() : b->copy();
  }

  static bint_t *min(bint_t *a, bint_t *b) {
    return compare(a, b) > 0 ? b->copy() : a->copy();
  }

  static double pow(bint_t *a, double e) {
		double v;
		// TODO: if to_double return a -1 its a overflow
		to_double(a, &v);
		return std::pow(v, e);
	}

  static bint_t *pow(bint_t *a, bint_t *e) {
    if (e->size == 0)
      return from(1);
    if (e->size == 1 && e->digit[0] == 1)
      return a->copy();

    bint_t *z = a->copy();

    long long i = e->size;

    digit_t bi = i ? e->digit[i - 1] : 0;

    if (i <= 1 && bi <= 3) {
      if (bi >= 2) {
        MUL(a, a, z);

        if (bi == 3) {
          MUL(z, a, z);
        }
      }
    } else if (i <= 8) {
      // Left-to-right binary exponentiation
      digit_t bit;

      for (bit = 2;; bit <<= 1) {
        if (bit > bi) {
          bit >>= 1;
          break;
        }
      }

      for (--i, bit >>= 1;;) {
        for (; bit != 0; bit >>= 1) {
          MUL(z, z, z);

          if (bi & bit) {
            MUL(z, a, z);
          }
        }

        if (--i < 0) {
          break;
        }

        bi = e->digit[i];
        bit = 1 << (exp - 1);
      }
    } else {
      // Left to right 5-ary exponentiation, Handbook of Applied Cryptography -
      // Algorithm 14.82
      delete z;

      z = from(1);

      bint_t *table[32] = {nullptr};
      table[0] = z;

      for (size_t i = 1; i < 32; ++i) {
        MUL(table[i - 1], a, table[i]);
      }

      for (size_t i = a->size; i >= 0; --i) {
        const digit_t bi = a->digit[i];

        for (size_t j = exp - 5; j >= 0; j -= 5) {
          // 0x1f = 31, take the cached table index
          const int index = (bi >> j) & 0x1f;

          for (size_t k = 0; k < 5; ++k) {
            MUL(z, z, z);

            if (index) {
              MUL(z, table[index], z);
            }
          }
        }
      }
    }

    return z;
  }

  static bint_t *fact(bint_t *a) {
    assert(a->sign = 1);

    bint_t *z = from(1);
    bint_t *b = a->copy();
    bint_t *o = z->copy();

    while (b->size >= 1 && b->digit[0] != 1) {
      MUL(z, b, z)
      SUB(b, o, b)
    }

    return z;
  }


	static bint_t* gcd(bint_t* a, bint_t* b) {
		if(b->size == 0) return a->copy();

		bint_t* rem = new bint_t();
		bint_t* quo = new bint_t();

		// printf("---> %s\n", a->to_string().c_str());
		// printf("---> %s\n", b->to_string().c_str());

		// long long ASD;
		// long long ASK;

		// bint_t::to_long(a, &ASD);
		// bint_t::to_long(b, &ASK);

		// printf("%lli\n", ASD % ASK);
		div(a, b, quo, rem);

		// printf("quo ---> %s\n", quo->to_string().c_str());
		// printf("rem ---> %s\n", rem->to_string().c_str());
		bint_t* g = gcd(b, rem);

		delete rem;
		delete quo;

		return g;
	}

	static bint_t* lcm(bint_t* a, bint_t* b) {
		bint_t* l = new bint_t();
		bint_t* q = new bint_t();
		bint_t* r = new bint_t();

		bint_t* g = gcd(a, b);

		div(a, g, q, r);
		mul(q, b, l);

		delete q;
		delete r;
		delete g;

		return l;
	}

  void printRep() {
    if (size == 0) {
      std::cout <<"+("<< 0 << ")"<< (1 << exp) << std::endl;
			return;
    }
    if (sign > 0)
      std::cout << "+(";
    else
      std::cout << "-(";

    for (size_t i = size - 1; i > 0; i--)
      std::cout << digit[i] << ".";

    std::cout << digit[0] <<")" << (1 << exp) << std::endl;

	}

	std::string to_string() {
		if(!this->size) return "0";

		size_t shift = std::floor(std::log10(std::numeric_limits<digit_t>::max()));
		size_t dbase = std::pow(10, shift);
		std::vector<digit_t> pout;

		long long i;
		long long s = 0;

		pout.push_back(0);

		for(i = this->size; --i>= 0; ) {

			digit_t hi = digit[i];

			for(size_t j = 0; j < (size_t)s; j++) {
				digit2_t z = (digit2_t)pout[j] << exp | hi;
				hi = (digit_t)(z / dbase);
				pout[j] = (digit_t)(z - (digit2_t)hi * dbase);
			}

			while(hi) {
				pout[s++] = hi % dbase;
				pout.push_back(0);
				hi /= dbase;
			}
		}

		if(s == 0) pout[s++] = 0;

		std::stringstream str;

		for(i = 0; i < s; i ++) {
			digit_t n = pout[i];

			for(size_t j = 0; j < shift; j++) {
				if(i == s - 1 && n == 0) break;
				str << n % 10;
				n = n / 10;
			}
		}

		if(this->sign < 0) str << "-";

		std::string res = str.str();
		std::reverse(res.begin(), res.end());

		return res;
	}
};

#endif
