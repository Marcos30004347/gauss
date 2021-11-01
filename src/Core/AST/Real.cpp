#include "Real.hpp"

#include <cmath>
#include <limits>
#include <iostream>

Real::Real(double v)
{
	double i;

	Int e = 0;

	while(modf(v, &i) >= std::numeric_limits<double>::epsilon())
	{
		v *= 10;
		e++;
	}

	Int b = 10;

	Int n = v;
	Int d = pow(b, e);

	Int u = gcd(n, d);

	num 	= n / u;
	den = d / u;
}



Real::Real(std::string v)
{
	std::string delimiter = ".";

	size_t pos = 0;

	std::string token;

	pos = v.find(delimiter);

	if(pos != std::string::npos)
	{
			token = v.substr(0, pos);
			v.erase(0, pos + delimiter.length());

			num = Int(token);
			den = Int(v);
	}
	else
	{
		num = Int(v);
		den = 1;
	}
}

Real::Real(Int v)
{
	num = v;
	den = 1;
}

Real::Real(Int n, Int d)
{
	Int u = gcd(n, d);

	num 	= n / u;
	den = d / u;
}

Real::Real()
{
	this->num = 0;
	this->den = 1;
}

Real Real::operator+(const Real& b) const
{
	Int k = this->num * b.den + this->den * b.num;
	Int j = b.den * this->den;

	Int u = gcd(k, j);

	return Real(k/u, j/u);
}

Real Real::operator-(const Real& b) const
{
	Int k = this->num * b.den - this->den * b.num;
	Int j = b.den * this->den;

	Int u = gcd(k, j);

	return Real(k/u, j/u);
}

Real Real::operator*(const Real& b) const
{
	Int k = this->num * b.num;
	Int j = b.den * this->den;

	Int u = gcd(k, j);

	return Real(k/u, j/u);
}

Real Real::operator/(const Real& b) const
{
	Int k = this->num * b.den;
	Int j = b.num * this->den;

	Int u = gcd(k, j);

	return Real(k/u, j/u);
}

Real Real::operator-() const
{
	return *this * Real(Int(-1), Int(1));
}

bool Real::operator>(const Real& other) const
{
	return this->num * other.den > other.num * this->den;
}

bool Real::operator>=(const Real& other) const
{
	return this->num * other.den >= other.num * this->den;
}

bool Real::operator<(const Real& other) const
{
	return this->num * other.den < other.num * this->den;
}

bool Real::operator<=(const Real& other) const
{
	return this->num * other.den <= other.num * this->den;
}

bool Real::operator==(const Real& other) const
{
	return this->num == other.num && other.den == this->den;
}

bool Real::operator!=(const Real& other) const
{
	return this->num != other.num || other.den != this->den;
}

Real Real::nan()
{
	Real v;
	v.not_a_number = true;
	return v;
}

Int Real::toInt() const
{
	return this->num/this->den;
}

bool Real::isNan() const
{
	return this->not_a_number;
}

Int Real::numerator() const
{
	return this->num;
}

Int Real::denominator() const
{
	return this->den;
}

Real abs(Real v)
{
	if(v.numerator() < 0 && v.denominator() > 0)
	{
		return Real(-1 * v.numerator(), v.denominator());
	}

	if(v.numerator() > 0 && v.denominator() < 0)
	{
		return Real(v.numerator(), -1 * v.denominator());
	}

	if(v.numerator() < 0 && v.denominator() < 0)
	{
		return Real(-1 * v.numerator(), -1 * v.denominator());
	}

	return v;
}

Real pow(Real a, Int b)
{
	Int t = pow(a.numerator(), b);
	Int j = pow(a.denominator(), b);

	Int u = gcd(t, j);

	return Real(t / u, j / u);
}

Real nthRoot(Real num, Int n, Real precision)
{
	Real x = num / 2;

	Real dx = (num / pow(x, n - 1) - x) / n;

	while(dx >= precision || dx <= -precision)
	{
		x = x + dx;
		dx = (num / pow(x, n - 1) - x ) / n;
	}

	return x;
}

Real sqrt(Real v, Real precision)
{
	return nthRoot(v, 2, precision);
}


std::string Real::to_string()
{
	if(this->den == 1)
	{
		return this->num.to_string();
	}

	return this->num.to_string() + "/" + this->den.to_string(); 
}

Real Real::eps = Real(Int(1), pow(Int(10), Int(22)));

Real Real::computeEulerConstant(Real eps)
{
	Real e = Real(2.0);
	Int c = 2;
	Real t = e;
	
	do 
	{
		t = e;
		e = e + Real(1, fact(c));
		c = c + 1;
	} while(abs(e - t) > eps);
	
	return e;
}

Real Real::computePiConstant(Real eps)
{
	Real pi = Real(Int(0));

	Real t = pi;
	Int n  = 0;
	
	do 
	{
		t = pi;
		pi = pi + Real(pow(Int(-1), n), 2*n + 1);
		n = n + 1;
		std::cout << pi.to_string() << std::endl;
	} while(abs(pi - t) > eps);

	return pi*Int(4);
}
