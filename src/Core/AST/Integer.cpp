#include "Integer.hpp"

int Int::size()
{
	if(a.empty())return 0;
	int ans=(a.size()-1)*Int::base_digits;
	int ca=a.back();

	while(ca)
	{
		ans++,ca/=10;
	}

	return ans;
}

Int Int::operator ^(const Int &v)
{
	Int ans=1, a=*this,b=v;

	while(!b.isZero())
	{
		if(b%2)
		{
			ans*=a;
		}

		a*=a,b/=2;
	}
	return ans;
}

std::string Int::to_string()
{
	std::stringstream ss;
	ss << *this;
	std::string s;
	ss >> s;
	return s;
}

int Int::sumof()
{
	std::string s = to_string();
	int ans = 0;
	for(auto c : s)  ans += c - '0';
	return ans;
}

Int::Int():sign(1) {}

Int::Int(long long v) 
{
	sign = 1;

	a.clear();

	if (v < 0)
	{
		sign = -1, v = -v;
	}

	for (; v > 0; v = v / 	Int::base)
	{
		a.push_back(v % 	Int::base);
	}
}

Int::Int(const std::string &s) 
{
	read(s);
}

Int Int::operator+(const Int &v) const 
{
	if (sign == v.sign) 
	{
		Int res = v;

		for (int i = 0, carry = 0; i < (int) std::max(a.size(), v.a.size()) || carry; ++i) 
		{
			if (i == (int) res.a.size())
			{
				res.a.push_back(0);
			}

			res.a[i] += carry + (i < (int) a.size() ? a[i] : 0);
	
			carry = res.a[i] >= 	Int::base;
	
			if (carry)
			{
				res.a[i] -= 	Int::base;
			}
		}
		return res;
	}
	return *this - (-v);
}

Int Int::operator-(const Int &v) const 
{
	if (sign == v.sign) 
	{
		if (abs() >= v.abs()) 
		{
			Int res = *this;
			
			for (int i = 0, carry = 0; i < (int) v.a.size() || carry; ++i) 
			{
				res.a[i] -= carry + (i < (int) v.a.size() ? v.a[i] : 0);
				carry = res.a[i] < 0;
				if (carry)
					res.a[i] += 	Int::base;
			}
			
			res.trim();
			return res;
		}

		return -(v - *this);
	}

	return *this + (-v);
}

void Int::operator*=(int v) 
{
	if (v < 0)
	{
		sign = -sign, v = -v;
	}

	for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) 
	{
		if (i == (int) a.size())
		{
			a.push_back(0);
		}

		long long cur = a[i] * (long long) v + carry;
		carry = (int) (cur / 	Int::base);
		a[i] = (int) (cur % 	Int::base);
	}

	trim();
}

Int Int::operator*(int v) const 
{
	Int res = *this;
	res *= v;
	return res;
}

void Int::operator++()
{
	*this = *this + 1;
}

void Int::operator--()
{
	*this = *this - 1;
}

void Int::operator++(int)
{
	*this = *this + 1;
}

void Int::operator--(int)
{
	*this = *this - 1;
}

void Int::operator*=(long long v) 
{
	if (v < 0)
	{
		sign = -sign, v = -v;
	}

	if(v > 	Int::base)
	{
		*this = *this * (v / 	Int::base) * 	Int::base + *this * (v % 	Int::base);
		return ;
	}

	for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) 
	{
		if (i == (int) a.size())
		{
			a.push_back(0);
		}

		long long cur = a[i] * (long long) v + carry;
		carry = (int) (cur / 	Int::base);
		a[i] = (int) (cur % 	Int::base);
	}

	trim();
}

Int Int::operator*(long long v) const 
{
	Int res = *this;
	res *= v;
	return res;
}

Int Int::operator/(const Int &v) const 
{
	return divmod(*this, v).first;
}

Int Int::operator%(const Int &v) const 
{
	return divmod(*this, v).second;
}

void Int::operator/=(int v) 
{
	if (v < 0)
		sign = -sign, v = -v;
	for (int i = (int) a.size() - 1, rem = 0; i >= 0; --i) {
		long long cur = a[i] + rem * (long long) 	Int::base;
		a[i] = (int) (cur / v);
		rem = (int) (cur % v);
	}
	trim();
}

Int Int::operator/(int v) const 
{
	Int res = *this;
	res /= v;
	return res;
}

int Int::operator%(int v) const 
{
	if (v < 0)
	{
		v = -v;
	}

	int m = 0;

	for (int i = a.size() - 1; i >= 0; --i)
	{
		m = (a[i] + m * (long long) 	Int::base) % v;
	}

	return m * sign;
}

void Int::operator+=(const Int &v) 
{
	*this = *this + v;
}

void Int::operator-=(const Int &v) 
{
	*this = *this - v;
}

void Int::operator*=(const Int &v) 
{
	*this = *this * v;
}

void Int::operator/=(const Int &v) 
{
	*this = *this / v;
}

bool Int::operator<(const Int &v) const 
{
	if (sign != v.sign)
	{
		return sign < v.sign;
	}

	if (a.size() != v.a.size())
	{
		return a.size() * sign < v.a.size() * v.sign;
	}

	for (int i = a.size() - 1; i >= 0; i--)
	{
		if (a[i] != v.a[i])
		{
			return a[i] * sign < v.a[i] * sign;
		}
	}

	return false;
}

bool Int::operator>(const Int &v) const 
{
	return v < *this;
}

bool Int::operator<=(const Int &v) const 
{
	return !(v < *this);
}

bool Int::operator>=(const Int &v) const 
{
	return !(*this < v);
}

bool Int::operator==(const Int &v) const 
{
	return !(*this < v) && !(v < *this);
}

bool Int::operator!=(const Int &v) const 
{
	return *this < v || v < *this;
}

void Int::trim() 
{
	while (!a.empty() && !a.back())
	{
		a.pop_back();
	}

	if (a.empty())
	{
		sign = 1;
	}
}

bool Int::isZero() const 
{
	return a.empty() || (a.size() == 1 && !a[0]);
}

Int Int::operator-() const 
{
	Int res = *this;
	res.sign = -sign;
	return res;
}

Int Int::abs() const 
{
	Int res = *this;
	res.sign *= res.sign;
	return res;
}

long long Int::longValue() const 
{
	long long res = 0;
	for (int i = a.size() - 1; i >= 0; i--)
		res = res * 	Int::base + a[i];
	return res * sign;
}

void Int::read(const std::string &s) 
{
	sign = 1;
	a.clear();
	int pos = 0;
	
	while (pos < (int) s.size() && (s[pos] == '-' || s[pos] == '+')) 
	{
		if (s[pos] == '-')
			sign = -sign;
		++pos;
	}
	
	for (int i = s.size() - 1; i >= pos; i -= Int::base_digits) 
	{
		int x = 0;
		for (int j = std::max(pos, i - Int::base_digits + 1); j <= i; j++)
			x = x * 10 + s[j] - '0';
		a.push_back(x);
	}
	
	trim();
}


std::vector<int> Int::convert_base(const std::vector<int> &a, int old_digits, int new_digits) {
	std::vector<long long> p(std::max(old_digits, new_digits) + 1);
	p[0] = 1;

	for (int i = 1; i < (int) p.size(); i++)
	{
		p[i] = p[i - 1] * 10;
	}
	
	std::vector<int> res;
	
	long long cur = 0;
	int cur_digits = 0;
	
	for (int i = 0; i < (int) a.size(); i++) 
	{
		cur += a[i] * p[cur_digits];
		cur_digits += old_digits;
		
		while (cur_digits >= new_digits) 
		{
			res.push_back(int(cur % p[new_digits]));
			cur /= p[new_digits];
			cur_digits -= new_digits;
		}
	}

	res.push_back((int) cur);

	while (!res.empty() && !res.back())
	{
		res.pop_back();
	}
	return res;
}

std::vector<long long> Int::karatsubaMultiply(const std::vector<long long> &a, const std::vector<long long> &b) {
	int n = a.size();

	std::vector<long long> res(n + n);

	if (n <= 32) 
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				res[i + j] += a[i] * b[j];
			}
		}
		return res;
	}

	int k = n >> 1;

	std::vector<long long> a1(a.begin(), a.begin() + k);
	std::vector<long long> a2(a.begin() + k, a.end());
	std::vector<long long> b1(b.begin(), b.begin() + k);
	std::vector<long long> b2(b.begin() + k, b.end());
	std::vector<long long> a1b1 = karatsubaMultiply(a1, b1);
	std::vector<long long> a2b2 = karatsubaMultiply(a2, b2);

	for (int i = 0; i < k; i++)
	{
		a2[i] += a1[i];
	}

	for (int i = 0; i < k; i++)
	{
		b2[i] += b1[i];
	}

	std::vector<long long> r = karatsubaMultiply(a2, b2);

	for (int i = 0; i < (int) a1b1.size(); i++)
	{
		r[i] -= a1b1[i];
	}

	for (int i = 0; i < (int) a2b2.size(); i++)
	{
		r[i] -= a2b2[i];
	}

	for (int i = 0; i < (int) r.size(); i++)
	{
		res[i + k] += r[i];
	}

	for (int i = 0; i < (int) a1b1.size(); i++)
	{
		res[i] += a1b1[i];
	}

	for (int i = 0; i < (int) a2b2.size(); i++)
	{
		res[i + n] += a2b2[i];
	}

	return res;
}

Int Int::operator*(const Int &v) const 
{
	std::vector<int> a6 = convert_base(this->a, Int::base_digits, 6);
	std::vector<int> b6 = convert_base(v.a, Int::base_digits, 6);

	std::vector<long long> a(a6.begin(), a6.end());
	std::vector<long long> b(b6.begin(), b6.end());

	while (a.size() < b.size())
	{
		a.push_back(0);
	}

	while (b.size() < a.size())
	{
		b.push_back(0);
	}

	while (a.size() & (a.size() - 1))
	{
		a.push_back(0), b.push_back(0);
	}

	std::vector<long long> c = karatsubaMultiply(a, b);
	
	Int res;
	
	res.sign = sign * v.sign;

	for (int i = 0, carry = 0; i < (int) c.size(); i++) 
	{
		long long cur = c[i] + carry;
		res.a.push_back((int) (cur % 1000000));
		carry = (int) (cur / 1000000);
	}

	res.a = convert_base(res.a, 6, Int::base_digits);
	res.trim();
	return res;
}

Int abs(Int a)
{
	return a.abs();
}

Int fact(Int i)
{
	if(i == 1 || i == 0)
	{
		return Int(1);
	}
	
	return i * fact(i-1);
}

Int pow(Int base, Int exp){
    if (exp == 0){
        return 1;
    }
	
    if (exp == 1){
        return base;
    }
	
    Int oneHalf = pow(base, exp / 2);
	
    return oneHalf * oneHalf * pow(base,exp % 2);
}

Int max(Int a, Int b)
{
	return a > b ? a : b;
}

Int min(Int a, Int b)
{
	return a < b ? a : b;
}

bool operator<(const unsigned int& a, const Int &v)
{
	return Int(a) < v;
}

Int operator*(const int& a, const Int &v)
{
	return Int(a) * v;
}

bool operator>(const unsigned int& a, const Int &v)
{
	return Int(a) > v;
}

bool operator<=(const unsigned int& a, const Int &v)
{
	return Int(a) <= v;
}

bool operator>=(const unsigned int& a, const Int &v)
{
	return Int(a) >= v;
}
