#ifndef INTEGER_H
#define INTEGER_H

#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip> 

struct Int
{
	std::vector<int> a;

	static const int base = 1000000000;
	static const int base_digits = 9; 

	int sign;

	int size();

	Int operator ^(const Int &v);
	
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
	Int operator-(const Int &v) const;  
	Int operator-() const;
	Int operator/(const Int &v) const; 
	Int operator%(const Int &v) const;
	Int operator*(const Int &v) const;
	Int operator/(int v) const;
	Int operator*(int v) const ;
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

	friend std::pair<Int, Int> divmod(const Int &a1, const Int &b1) 
	{
		int norm = 	Int::base / (b1.a.back() + 1);
	
		Int a = a1.abs() * norm;
		Int b = b1.abs() * norm;
		Int q, r;
	
		q.a.resize(a.a.size());
 
		for (int i = a.a.size() - 1; i >= 0; i--) 
		{
			r *= 	Int::base;
			r += a.a[i];
			int s1 = r.a.size() <= b.a.size() ? 0 : r.a[b.a.size()];
			int s2 = r.a.size() <= b.a.size() - 1 ? 0 : r.a[b.a.size() - 1];
			int d = ((long long) 	Int::base * s1 + s2) / b.a.back();
			r -= b * d;

			while (r < 0)
			{
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
 
	friend Int gcd(const Int &a, const Int &b) 
	{
		return b.isZero() ? a : gcd(b, a % b);
	}
	
	friend Int lcm(const Int &a, const Int &b) 
	{
		return a / gcd(a, b) * b;
	}
 
	void read(const std::string &s);
 
	friend std::istream& operator>>(std::istream &stream, Int &v) 
	{
		std::string s;
		stream >> s;
		v.read(s);
		return stream;
	}
 
	friend std::ostream& operator<<(std::ostream &stream, const Int &v) {
		if (v.sign == -1)
		{
			stream << '-';
		}

		stream << (v.a.empty() ? 0 : v.a.back());

		for (int i = (int) v.a.size() - 2; i >= 0; --i)
		{
			stream << std::setw(Int::base_digits) << std::setfill('0') << v.a[i];
		}

		return stream;
	}
 
	static std::vector<int> convert_base(const std::vector<int> &a, int old_digits, int new_digits);
	static std::vector<long long> karatsubaMultiply(const std::vector<long long> &a, const std::vector<long long> &b);
};

Int abs(Int a);
Int fact(Int i);
Int pow(Int base,Int exp);
Int max(Int a, Int b);
Int min(Int a, Int b);

bool operator<(const unsigned int& a, const Int &v);
Int operator*(const int& a, const Int &v);
bool operator>(const unsigned int& a, const Int &v);
bool operator<=(const unsigned int& a, const Int &v);
bool operator>=(const unsigned int& a, const Int &v);



#endif
