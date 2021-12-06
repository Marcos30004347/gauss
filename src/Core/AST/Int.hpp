#ifndef INT_HPP
#define INT_HPP

// references:
// https://cacr.uwaterloo.ca/hac/about/chap14.pdf
// https://github.com/python/cpython/blob/main/Objects/longobject.c
// The Art of Computer Programming Vol 2 by Donald E. Knuth

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <limits>
#include <float.h>
#include <math.h>


#define pow2(e) (1 << e)

// computes n mod 2^e
#define modPow2(n, e) (n & ((1 << e) - 1))

// computes n / 2^e
#define quoPow2(n, e) (n >> e)

#define mulPow2(d, exp) (d << exp)

// Code taken from pycore_bitutils.h from the
// python programming language
static inline int high_bit(uint32_t x) {
  const int K[32] = {0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
                     5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

  int k = 0;

  while (x >= 32) {
    k += 6;
    x >>= 6;
  }

  return k + K[x];
}

#define arith_rshift(I, J) ((I) < 0 ? -1 - ((-1 - (I)) >> (J)) : (I) >> (J))

template <char exp = 30, typename single_type = uint32_t,
          typename double_type = uint64_t>
class bint {
  static_assert(sizeof(double_type) >= 2 * sizeof(single_type),
                "sizeof(double_type) needs to be at least twice"
                "as big as sizeof(single_type)");

  // digit2_t is a type capable of holding
  // at least two elements of type single_type

  static const single_type base = pow2(exp);
  static const single_type mask = base - (single_type)1;

public:
  // digit one is the type used to store every digit of the number base 2^exp
  using digit_t = single_type;
  using digit2_t = double_type;
  using bint_t = bint<exp, digit_t, digit2_t>;

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
    t->sing = this->sign;
    memcpy(this->digit, t->digit, this->size * sizeof(digit_t));
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
  }

  // shift the bits in a to the left m times and save the
  // result on z
  static digit_t digits_lshift(digit_t *a, size_t length, size_t d, digit_t *z) {
    digit_t carry = 0;

    assert(0 <= d && d < sizeof(digit_t));

    for (size_t i = 0; i < length; i++) {
      // shift digits and combine them
      // with the carry
      digit2_t acc = (digit2_t)a << d | carry;
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
  static digit_t digits_rshift(digit_t *a, size_t length, size_t d, digit_t *z) {
    digit_t carry = 0;

    // mask with last m bits set
    digit_t mask = pow2(d) - 1;

    assert(0 <= d && d < sizeof(digit_t));


    for (size_t i = length; i-- > 0;) {
      digit2_t acc = (digit2_t)carry << exp | a[i];
      carry = (digit_t)acc & mask;
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
  template <typename T> static bint<exp, digit_t, digit2_t> from(T x) {
    short sign = 1;

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

    return bint(v, i + 1, sign);
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
      i = a;

      while (i > 0 && x->digit[i] == y->digit[i])
        i = i - 1;

      if (i == 0 && x->digit[i] == y->digit[i]) {
        z->resize(0);
        z->size = 0;
        z->sign = 1;
        return;
      }

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
        carry = quoPow2(carry, exp);
        assert(carry <= mask);
      }

      if (carry) {
        *pz += (digit_t)(carry & mask);
      }

      assert((carry >> exp) == 0);
    }

    return z->trim();
  }

  static void add(bint_t *x, bint_t *y, bint_t *z) {
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
        return abs_sub_digits(y, z, z);
      }
      abs_add_digits(x, y, z);
      if (z->size) {
        z->sign = -1 * z->sign;
      }
    }

    if (y->sign < 0) {
      return abs_add_digits(x, y, z);
    }

    return abs_sub_digits(x, y, z);
  }

  static short fast_mod(bint_t *x, bint_t *y, bint_t *rem) {
    int32_t left = x->digit[0];
    int32_t righ = y->digit[0];

    assert(x->size == 1);
    assert(y->size == 1);

    rem->resize(1);
    rem->sign = y->sign;

    if (x->sign == y->sign) {
      rem->digit[0] = left % righ;
    } else {
      rem->digit[0] = righ - 1 - (left - 1) % righ;
    }

    return 1;
  }

  static short fast_div(bint_t *x, bint_t *y, bint_t *quo) {
    int32_t left = x->digit[0];
    int32_t righ = y->digit[0];

    assert(x->size == 1);
    assert(y->size == 1);

    quo->resize(1);
    quo->sign = y->sign;

    if (x->sign == y->sign) {
      quo->digit[0] = left / righ;
    } else {
      quo->digit[0] = -1 - (left - 1) / righ;
    }

    return 1;
  }
  static short div(bint_t *x, bint_t *y, bint_t *quo, bint_t *rem) {
    // Following The Art of Computer Programming, Vol.2, section 4.3.1,
    // Algorithm D.
    size_t m = x->size;
    size_t n = y->size;

    if (quo)
      quo->sign = x->sign * y->sign;
    if (rem)
      rem->sign = x->sign * y->sign;

    if (m == 0) {
      quo->resize(0);
      rem->resize(0);
      return 1;
    }

    if (n == 0) {
      return 0;
    }

    if (m == 1 && n == 1) {
      fast_div(x, y, quo);
      fast_mod(x, y, rem);

      quo->trim();
      rem->trim();

      return 1;
    }

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

      rem->resize(1);
      rem->digit[0] = r;

      quo->trim();
      rem->trim();

      return 1;
    }

    digit_t carry = 0;

    // D1: Normalization

    // because we are using binary base 2^exp
    // this operatios return the expoent of
    // 2^x = floor(2^exp/(v[n - 1] + 1))
    digit_t d = exp - high_bit(y->digit[n - 1]);

    bint_t u;
    bint_t v;

    u.resize(m + 1);
    u.digit[m] = 0;

    v.resize(n);

    // sinse we are using binary base, shifting
    // the value left is the same as multiplying
    // by 2^d

    // Multiply y by 2^d
    carry = digits_lshift(y->digit, y->size, d, v->digit);
    assert(carry == 0);

    // Multiply x by 2^d
    carry = digits_lshift(x->digit, x->size, d, u->digit);
    if (carry != 0 || u.digit[m - 1] >= v.digit[n - 1]) {
      u.digit[m] = carry;
      m = m + 1;
    }

    // D2 - loop from D3 to D7
    digit_t v1 = v.digit[n - 1];
    digit_t v2 = v.digit[n - 2];

    size_t j = m - n;

    assert(j >= 0);

    quo->resize(j);

    digit_t *u0 = u.digit;
    digit_t *v0 = v.digit;

    digit_t *uj = u0 + j;
    digit_t *qj = quo->digit + j;

    while (uj-- > u0) {
      // D3
      digit_t ut = uj[n];

      // uu = (u[j + n]*b + u[j + n - 1])
      digit2_t uu = ((digit2_t)ut << exp) | uj[n - 1];

      digit_t q = (digit_t)(uu / v1);
      digit_t r = (digit_t)(uu - (digit2_t)v1 * q);

      // test q >= b or q*v[n - 2] > b*r + u[j + n - 2]
      // r << exp | uj[n - 2] = r*b + u[j + n - 2]
      while (q > base || (digit2_t)v2 * q > ((digit2_t)r << exp | uj[n - 2])) {
        q = q - 1;
        r = r + v1;

        if (r >= base)
          break;
      }

      // D4 replace (u[j + n], u[j + n - 1],...,u[j])b
      // by (u[j + n], u[j + n - 1],...,u[j])b - q*(0,v[n-1],...,v[1],v[0])b
      int32_t borrow = 0;
      for (size_t i = 0; i < n; ++i) {
        int64_t z = (int32_t)uj[i] + borrow - (int64_t)q * (int64_t)v0[i];
        uj[i] = (digit_t)z & mask;
        borrow = (int32_t)arith_rshift(z, exp);
      }

      // D5
      assert((int32_t)ut + borrow == -1 || (int32_t)ut + borrow == 0);

      // D6
      if ((int32_t)ut + borrow < 0) {
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

    carry = rshift(&u, d, rem);

    assert(carry == 0);

    rem->trim();
    quo->trim();

    return 1;
  }

  // TODO: to string of the number in base 10
  // TODO: add pow methods
  // TODO: add fact method
  // TODO: add gcd method
  // TODO: add lcm method
  // TODO: add construction from strings and const char* types

  static short compare(bint_t *v0, bint_t *v1) {
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

  void print() {
    for (size_t i = size - 1; i > 0; i--) {
      digit[i] * pow2(exp);
    }
  }

	static double frexp(bint_t* a, size_t* e, bool* overflow) {
		// from https://github.com/python/cpython/blob/main/Objects/longobject.c
		size_t m = a->size;

		if(m == 0) {
			*e = 0;
			return 0;
		}

		double dx;

		digit_t* x_digits = new digit_t[2 + (DBL_MANT_DIG + 1) / exp];
		bint_t x = bint_t(x_digits, 2 + (DBL_MANT_DIG + 1) / exp);

		static const int half_even_correction[8] = {0, -1, -2, 1, 0, -1, 2, 1};

		size_t al = high_bit(a->digit[m - 1]);
		size_t mx = std::numeric_limits<size_t>::max();

		if(m >= (m - 1) / exp + 1 && (m > (mx - 1) / exp + 1 || al > (mx - 1) % exp + 1)) {
			*overflow = true;
			return -1;
		}

		al = (m - 1) * exp + al;

		size_t x_size = 0;

		if(al <= DBL_MANT_DIG + 2) {
			size_t shift_d = (DBL_MANT_DIG + 2 - al) / exp;
			size_t shift_b = (DBL_MANT_DIG + 2 - al) % exp;

			x_size = shift_d;
			digit_t rem = digits_lshift(a->digit, m, shift_b, x_digits + x_size);

			x_size += m;

			x.digits[x_size++] = rem;
		} else {
			size_t shift_d = (al + 2 - al) / exp;
			size_t shift_b = (al + 2 - al) % exp;
			digit_t rem = digits_rshift(a->digit + shift_d, m - shift_d, shift_b, x_digits);

			x_size = m - shift_d;

			if(rem) {
				x_digits[0] |= 1;
			} else {
				while(shift_d > 0) {
					if(a->digit[--shift_d]) {
						x_digits[0] |= 1;
						break;
					}
				}
			}
		}

		x_digits[0] += half_even_correction[x_digits[0] & 7];

		dx = x_digits[--x_size];

		while(x_size > 0) {
			dx = dx * base + x_digits[--x_size];
		}

		dx /= 4.0 * std::ldexp(1.0, DBL_MANT_DIG);

		if(dx == 1.0) {
			if(al == mx) {
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

  static double to_double(bint_t *b) {
    if (b->size == 0)
      return 0;

    if (b->size == 1)
      return (double)b->digit[0] * b->sign;

    if (b->size == 2) {
			size_t margin = std::numeric_limits<digit_t>::digits - exp;
			digit_t carry = b->digits[1] & ((1 << margin) - 1);

			digit_t v0 = (carry << exp) & (b->digit[0] & mask);
			digit_t v1 = b->digit[1] >> margin;

			return (double)((digit2_t)v0 & (digit2_t)v1) * b->sign;
		}

		bool overflow = false;

		size_t e = 0.0;

		double x = frexp(b, &e, &overflow);

		if((x == -1 && overflow)|| e > FLT_MAX_EXP) {
			//TODO: ERROR OVERFLOW
			return -1;
		}

		return std::ldexp(x, e);
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

  void printRep() {
    if (size == 0) {
      std::cout << 0 << std::endl;
    }
    if (sign > 0)
      std::cout << "+";
    else
      std::cout << "-";

    for (size_t i = size - 1; i > 0; i--)
      std::cout << digit[i] << ".";
    std::cout << digit[0] << std::endl;
  }
};

#endif
