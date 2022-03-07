#ifndef INTEGER_H
#define INTEGER_H

#include "Int.hpp"

#include <climits>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include <string>

#define LONG_LONG_OK 1
#define LONG_LONG_OVERFLOW LLONG_MAX

struct Int {
  Int(bint<30> *v);

  char flag = 0;

	union {
		long long x = 0;
		bint<30> *val;
	};

  void to_long_if_small() {
    if (!this->flag || this->val->size > 2)
      return;

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
  Int(unsigned long int);
  Int(unsigned int);
  Int(int);
  Int(double);
  Int(const Int &);
  Int(Int &&);
  ~Int();
	std::string to_string();
	static Int fromString(const char*);

	Int operator+(const Int &other) const ;
	Int operator+(const Int &&other) const ;
	Int operator+(const int z) const;
	Int operator-(const Int &other) const ;
	Int operator-(const Int &&other) const;
	Int operator-(const int z) const ;
	Int operator*(const Int &other) const;
	Int operator*(const Int &&other) const ;
	Int operator*(const int z) const ;
	Int operator/(const Int &other) const ;
	Int operator/(const Int &&other) const;
	Int operator/(const int z) const ;
	Int operator%(const Int &other) const ;
	Int operator%(const Int &&other) const;
	Int operator%(const int z) const ;
	Int operator+() ;
	Int operator-();
	Int operator++();
	Int operator--() ;
	Int operator++(int) ;
	Int operator--(int) ;
	bool operator==(const Int &other) const ;
	bool operator==(const Int &&other) const ;
	bool operator<(const Int &other) const;
	bool operator<(const Int &&other) const;
	bool operator<=(const Int &other) const ;
	bool operator<=(const Int &&other) const ;
	bool operator>(const Int &other) const;
	bool operator>(const Int &&other) const;
	bool operator>=(const Int &other) const ;
	bool operator>=(const Int &&other) const ;
	bool operator!=(const Int &other) const ;
	bool operator!=(const Int &&other) const ;
	Int operator=(const Int &other) ;
	Int operator=(Int &&other) ;
	Int ceil_log2() ;
	long long longValue() ;
	double doubleValue() ;
  friend Int gcd(const Int &a, const Int &b);
  friend Int gcd(const Int &&a, const Int &&b) ;
  friend Int lcm(const Int &a, const Int &b);
  friend Int lcm(const Int &&a, const Int &&b);
	void operator/=(Int v);
	void operator+=(Int v);
	void operator-=(Int v);
	void operator*=(Int v);

  explicit operator bool()  const {
		if (!this->flag)
			return this->x != 0;
		return this->val->size > 0;
	}

	// operator long long() const {
	// 	if(!flag) return x;

	// 	long long r;

	// 	bint<30>::to_long(val, &r);

	// 	return r;
	// }

	// operator unsigned long() const {
	// 	if(!flag) return x;

	// 	long long r;

	// 	bint<30>::to_long(val, &r);

	// 	return r;
	// }


  friend Int abs(const Int &&a);

  friend Int abs(const Int &a);

  friend Int fact(const Int &&a);
  friend Int fact(const Int &a);
  friend Int max(const Int &&a, Int &&b);
  friend Int max(const Int &a, Int &b);
  friend Int max(const Int &a, Int &&b);
  friend Int min(const Int &&a, const Int &&b) ;
  friend Int min(const Int &a, const Int &b);
  friend Int pow(const Int &&a, const Int &&b);
  friend Int pow(const Int &a, const Int &b);
  friend Int isqrt(const Int &a);
  friend Int operator*(const int a, const Int &v) ;
  friend Int operator+(const long long a, const Int &v) ;
  friend Int operator-(const long long a, const Int &v);
  friend double pow(const Int &&a, const double b);
  friend double pow(const Int &a, const double b);
  friend bool operator<(const unsigned int a, const Int &v);
  friend bool operator>(const unsigned int a, const Int &v);
  friend bool operator<=(const unsigned int a, const Int &v);
  friend bool operator>=(const unsigned int a, const Int &v);
};

#endif
