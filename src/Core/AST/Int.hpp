#ifndef INT_HPP
#define INT_HPP

// references:
// https://cacr.uwaterloo.ca/hac/about/chap14.pdf
// https://github.com/python/cpython/blob/main/Objects/longobject.c
// The Art of Computer Programming Vol 2 by Donald E. Knuth

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cassert>

#define pow2(e) (1 << e)

// computes n mod 2^e
#define modPow2(n, e) (n & ((1 << e) - 1))

// computes n / 2^e
#define quoPow2(n, e) (n >> e)

template <char exp = 30, typename single_type = uint32_t, typename double_type = uint64_t>
class bint {
	using digit2_t = double_type;

	static const single_type base = pow2(exp);
	static const single_type mask = base - 1;
public:

	using digit_t = single_type;

  digit_t *digit;
  size_t size;

  bint(digit_t *d, size_t s) : digit{d}, size{s} {}
  ~bint() { delete[] digit; }

  void trim() {
    while (!digit[size - 1])
      size--;
		digit = realloc(digit, sizeof(digit_t)*size);
  }

  // convert x to base 2^exp
  template <typename T> static bint<exp, digit_t, digit2_t> from(T x) {
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

    return bint(v, i + 1);
  }


	// sum the absolute values of the big integers with digits x[0...a]
	// and y[0...b] and save in z[0...a + 1]. It's assumed that a >= b.
	// All the memory should be pre-allocated before execution.
	static short add(digit_t* x, size_t a, digit_t* y, size_t b, digit_t* z) {
		digit_t carry = 0;

		size_t i = 0;

		for(; i < b; ++i) {
			carry += x[i] + y[i];
			z[i] = (digit_t)(carry & mask);
			carry >>= exp;
		}

		for(; i < a; ++i) {
			carry += x[i];
			z[i] = (digit_t)(carry & mask);
			carry >>= exp;
		}

		z[i] = carry;
		return 1;
	}

	// subtract the absolute values of the big integers with digits x[0...a]
	// and y[0...b] and save in z[0...a]. It's assumed that a >= b.
	// All the memory should be pre-allocated before execution.
	static short sub(digit_t* x, size_t a, digit_t* y, size_t b, digit_t* z) {
		short sign = 1;

		size_t i;


		// Get the index of the digit that x and y differ
		if(a == b) {
			i = a;
			while(--i >= 0 && x[i] == y[i]);

			if(i < 0) {
				z[0] = 0;
				return 1;
			}

			if(x[i] < y[i]) {
				sign = -1;

				digit_t *t = x;
				size_t   c = a;

				x = y;
				y = t;
				a = b;
				b = c;
			}

			a = i + 1;
			b = a;
		}

		digit_t borrow = 0;

		for(i = 0; i < b; ++i) {
			borrow = x[i] - y[i] - borrow;
			z[i] = borrow & mask;
			borrow >>= exp;
			borrow &= 1;
		}

		for(; i < a; ++i) {
			borrow = x[i] - borrow;
			z[i] = borrow & mask;
			borrow >>= exp;
			borrow &= 1;
		}

		return sign;
	}

	// Square the bint of digits x[0...a] and store the
	// result in w[0...2*a]
	static short square(digit_t* x, size_t a, digit_t* w) {
		size_t i = 0;

		for(; i < 2*a; i++) {
			w[i] = 0;
		}

		for(i = 0; i < a; i++) {
			digit2_t carry;
			digit2_t xi = x[i];

			digit_t *pw = w + (i << 1);
			digit_t *px = x + (i + 1);
			digit_t *pe = x + a;

			carry = *pw + xi * xi;
			*pw++ = (digit_t)(carry & mask);
			carry = quoPow2(carry, exp);

			xi <<= 1;

			while(px != pe) {
				carry += *pw + *px++ * xi;
				*pw++ = (digit_t)(carry & mask);
				carry = quoPow2(carry, exp);
				assert(carry <= (mask << 1));
			}

			if(carry) {
				carry += *pw;
				*pw++ = (digit_t)(carry & mask);
				carry = quoPow2(carry, exp);
			}

			if(carry) {
				*pw += (digit_t)(carry & mask);
			}

			assert((carry >> exp) == 0);
		}

		return 1;
	}

	// multiply the big integer with digits in x[0...a] with y[0...b]
	// and store the result in z[a + b]. It is assumed that a >= b.
	// All the space should be pre-alocated before execution
	static short mul(digit_t* x, size_t a, digit_t* y, size_t b, digit_t* z) {
		for(size_t i = 0; i < a; i++) {
			digit2_t carry = 0;
			digit2_t xi = x[i];

			digit_t *pz = z + i;
			digit_t *py = y;

			digit_t *pe = y + b;

			while(py < pe) {
				carry += *pz + *py++ * xi;
				*pz++ = (digit_t)(carry & mask);
				carry = quoPow2(carry, exp);
				assert(carry <= mask);
			}

			if(carry) {
				*pz += (digit_t)(carry & mask);
			}

			assert((carry >> exp) == 0);
		}

		return 1;
	}

	//TODO: currently, the binary operations only work
	//      if the size of the first operand is bigger
	//      or equal than the size of the second, fix
	//      that by considering adding a sign flag.

 	//TODO: fast division for bints of size 1

	//TODO: long division

	//TODO: to string of the number in base 10

	void print() {
		for(size_t i = size - 1; i > 0; i--) {
			 digit[i] * pow2(exp);
		}
	}

  void printRep() {
    for (size_t i = size - 1; i > 0; i--)
      std::cout << digit[i] << ".";
    std::cout << digit[0] << std::endl;
  }
};

#endif
