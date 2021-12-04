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
#include <cstring>

#define pow2(e) (1 << e)

// computes n mod 2^e
#define modPow2(n, e) (n & ((1 << e) - 1))

// computes n / 2^e
#define quoPow2(n, e) (n >> e)

template <char exp = 30, typename single_type = uint32_t, typename double_type = uint64_t>
class bint {
	static_assert(sizeof(double_type) >= 2*sizeof(single_type),
								"sizeof(double_type) needs to be at least twice"
								"as big as sizeof(single_type)");

	// digit2_t is a type capable of holding
	// at least two elements of type single_type

	static const single_type base = pow2(exp);
	static const single_type mask = base - 1;
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
	~bint() { if(digit) delete[] digit; }

	void resize(uint64_t s) {
		size = s;

		if(digit) delete digit;

		if(s == 0) {
			digit = nullptr;
			return;
		}

		digit = new digit_t[s];
	}

  void trim() {
		if(size == 0) return;

    size_t k = size - 1;

    while (k > 0 && !digit[k])
      k--;

		if(!digit[k]) {
			delete digit;

			size = 0;
			digit = nullptr;
		} else {
			size = k + 1;
			digit = (digit_t*)realloc(digit, sizeof(digit_t)*size);
		}
  }

  // convert x to base 2^exp
  template <typename T> static bint<exp, digit_t, digit2_t> from(T x) {
    short sign = 1;

		if(x < 0) {
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
	static void abs_add_digits(bint_t* x, bint_t* y, bint_t* z) {
		digit_t carry = 0;

		size_t i = 0;

		size_t a = x->size;
		size_t b = y->size;

		if(a < b) {
			bint_t* t = x;
			x = y;
			y = t;

			size_t c = a;
			a = b;
			b = c;
		}

		z->resize(x->size + 1);

		for(; i < b; ++i) {
			carry += x->digit[i] + y->digit[i];
			z->digit[i] = (digit_t)(carry & mask);
			carry >>= exp;
		}

		for(; i < a; ++i) {
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
	static void abs_sub_digits(bint_t* x, bint_t* y, bint_t* z) {
		short sign = 1;

		size_t a = x->size;
		size_t b = y->size;

		size_t i;

		if(b > a) {
			sign = -1;

			bint_t* t = x;
			x = y;
			y = t;

			size_t c = a;
			a = b;
			b = c;
		} else if(a == b) {
			// Get the index of the digit that x and y differ
			i = a;

			while(i > 0 && x->digit[i] == y->digit[i]) i = i - 1;

			if(i == 0 && x->digit[i] == y->digit[i]) {
				z->resize(0);
				z->size = 0;
			  z->sign = 1;
				return;
			}

			if(x->digit[i] < y->digit[i]) {
				bint_t* t = x;
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

		for(i = 0; i < b; ++i) {
			borrow = x->digit[i] - y->digit[i] - borrow;
			z->digit[i] = borrow & mask;
			borrow >>= exp;
			borrow &= 1;
		}

		for(; i < a; ++i) {
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
	static void abs_square_digits(bint_t* x, bint_t* w) {
		size_t i = 0;

		w->resize(2*x->size);

		for(; i < 2*x->size; i++) {
			w->digit[i] = 0;
		}

		for(i = 0; i < x->size; i++) {
			digit2_t carry;
			digit2_t xi = x->digit[i];

			digit_t *pw = w->digit + (i << 1);
			digit_t *px = x->digit + (i + 1);
			digit_t *pe = x->digit + x->size;

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

		return w->trim();
	}

	// multiply the big integer with digits in x[0...a] with y[0...b]
	// and store the result in z[a + b]. It is assumed that a >= b.
	// All the space should be pre-alocated before execution
	static void abs_mul_digits(bint_t* x, bint_t* y, bint_t* z) {
		if(x->size == 0 || y->size == 0) {
			return z->resize(0);
		}

		z->resize(x->size + y->size);

		for(size_t i = 0; i < x->size; i++) {
			digit2_t carry = 0;
			digit2_t xi = x->digit[i];

			digit_t *pz = z->digit + i;
			digit_t *py = y->digit;

			digit_t *pe = y->digit + y->size;

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

		return z->trim();
	}


	static void add(bint_t* x, bint_t* y, bint_t* z) {
		if(x->sign < 0) {
			if(y->sign < 0) {
				z->sign = -1 * z->sign;
				return abs_add_digits(x, y, z);
			}

			return abs_sub_digits(y, x, z);
		}

		if(y->sign < 0) {
			return abs_sub_digits(x, y, z);
		}

		return abs_add_digits(x, y, z);
		/*
			if(y->size == 0) {
			z->resize(x->size);
			z->digit = (digit_t*)memcpy(x->digit, z->digit,
			sizeof(digit_t)*x->size); z->sign = x->sign; z->size = x->size; return;
			}

			if(x->size == 0) {
			z->resize(y->size);
			z->digit = (digit_t*)memcpy(y->digit, z->digit,
			sizeof(digit_t)*y->size); z->sign = y->sign; z->size = y->size; return;
			}


			// -x + -y = -z or x + y = z
			if(x->sign == y->sign) {

			if (y->size > x->size) {
			t = x;
			x = y;
			y = t;
			}

			z->resize(x->size + 1);
			z->sign = x->sign;

			abs_add_digits(x->digit, x->size, y->digit, y->size,
			z->digit);

			return z->trim();
			}

			// -x + y = z
			if(x->sign < 0) {
			if(y->size > x->size) {
			z->sign = -1;
			t = x;
			x = y;
			y = t;
			} else {
			z->sign = 1;
			}

			x->printRep();
			y->printRep();

			z->resize(x->size + 1);

			z->sign *= abs_sub_digits(x->digit, x->size, y->digit,
			y->size, z->digit); z->printRep(); return z->trim();
			}

			// x + -y

			if(y->size > x->size) {
			z->sign = -1;
			t = x;
			x = y;
			y = t;
			} else {
			z->sign = 1;
			}

			z->resize(x->size + 1);

			z->sign *= abs_sub_digits(x->digit, x->size, y->digit,
			y->size, z->digit);

			return z->trim();
		*/
	}

	static void sub(bint_t* x, bint_t* y, bint_t * z) {
		if(x->sign < 0) {
			if(y->sign < 0) {
				return abs_sub_digits(y, z, z);
			}
			abs_add_digits(x, y, z);
			if(z->size) {
				z->sign = -1 * z->sign;
			}
		}

		if(y->sign < 0) {
			return abs_add_digits(x, y, z);
		}

		return abs_sub_digits(x, y, z);
		// if(x->size == 0 && y->size == 0) {
		// 	z->resize(0);
		// 	return;
		// }

		// if(y->size == 0) {
		// 	z->resize(x->size);

		// 	z = memcpy(x->digit, z->digit, sizeof(digit_t)*x->size);
		// 	z->sign = x->sign;
		// 	z->size = x->size;
		// 	return;
		// }

		// if(x->size == 0) {
		// 	z->resize(y->size);

		// 	z = memcpy(y->digit, z->digit, sizeof(digit_t)*y->size);
		// 	z->sign = -y->sign;
		// 	z->size = y->size;
		// 	return;
		// }

		// bint_t* t = nullptr;

		// // a - -b = a + b
		// if(x->sign > 0 && y->sign < 0) {
		// 	if(y->size > x->size) {
		// 		t = x;
		// 		x = y;
		// 		y = t;
		// 	}

		// 	z->resize(x->size + 1);

		// 	z->sign = 1;

		// 	z->sign *= abs_add_digits(x->digit, x->size, y->digit, y->size, z->digit);

		// 	return z->trim();
		// }

		// // -a - b = -(a + b)
		// if(x->sign < 0 && y->sign > 0) {
		// 	if(y->size > x->size) {
		// 		t = x;
		// 		x = y;
		// 		y = t;
		// 	}

		// 	z->resize(x->size + 1);

		// 	z->sign = 1;

		// 	z->sign *= abs_add_digits(x->digit, x->size, y->digit, y->size, z->digit);

		// 	return z->trim();
		// }

		// // -a - -b = -a + b
		// if(x->sign < 0 && y->sign < 0) {
		// 	if(y->size > x->size) {
		// 		z->sign = 1;
		// 		t = x;
		// 		x = y;
		// 		y= t;
		// 	} else if(x->size > y->size) {
		// 		z->sign = -1;
		// 	} else {
		// 		z->sign = 1;
		// 	}

		// 	z->resize(x->size + 1);

		// 	z->sign *= abs_sub_digits(x->digit, x->size, y->digit, y->size, z->digit);

		// 	return z->trim();
		// }

		// // a - b
		// if(y->size > x->size) {
		// 	z->sign = -1;
		// 	t = x;
		// 	x = y;
		// 	y = t;
		// } else {
		// 	z->sign = 1;
		// }

		// z->resize(x->size + 1);

		// z->sign *= abs_sub_digits(x->digit, x->size, y->digit, y->size, z->digit);

		// return z->trim();
	}



 	//TODO: fast division for bints of size 1
	//TODO: long division
	//TODO: to string of the number in base 10
	//TODO: add pow methods
	//TODO: add max/min methods
	//TODO: add abs method
	//TODO: add fact method
	//TODO: add gcd method
	//TODO: add lcm method
	//TODO: add construction from strings and const char* types


	static short compare(bint_t* v0, bint_t* v1) {
		if(v0->sign != v1->sign) return v0->sign > v1->sign ? 1 : -1;
		if(v0->size != v1->size) return v0->size > v1->size ? 1 : -1;

		size_t i = v0->size - 1;

		while(i > 0 && v0->digit[i] == v1->digit[i]) --i;

		if(i == 0 && v0->digit[0] == v1->digit[0]) return 0;

		return v0->digit[i] > v1->digit[i] ? 1 : -1;
	}

	void print() {
		for(size_t i = size - 1; i > 0; i--) {
			digit[i] * pow2(exp);
		}
	}

  void printRep() {
		if(size == 0) {
			std::cout << 0 << std::endl;
		}
		if(sign > 0) std::cout << "+";
		else std::cout << "-";

    for (size_t i = size - 1; i > 0; i--)
      std::cout << digit[i] << ".";
    std::cout << digit[0] << std::endl;
  }
};

#endif
