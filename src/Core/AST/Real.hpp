#ifndef REAL_H
#define REAL_H

#include <string>
#include "Integer.hpp"

class Real
{
	private:

	Int num;
	Int den;

	bool not_a_number = false;

	public:

	static Real eps;

	static Real computeEulerConstant(Real eps);
	static Real computePiConstant(Real eps);

	bool isNan() const;

	static Real nan();

	Real();
	Real(double);
	Real(std::string);
	Real(Int v);
	Real(Int n, Int d);

	Int numerator() const;
	Int denominator() const;

	Real operator+(const Real&) const;
	Real operator-(const Real&) const;
	Real operator*(const Real&) const;
	Real operator/(const Real&) const;
	Real operator-() const;

	Int toInt() const;

	bool operator>(const Real&) const;
	bool operator<(const Real&) const;
	bool operator<=(const Real&) const;
	bool operator>=(const Real&) const;
	bool operator==(const Real&) const;
	bool operator!=(const Real&) const;

	std::string to_string();
};

Real computeEulerConstant(Real eps = Real::eps);

Real pow(Real, Int);
Real abs(Real);
Real nthRoot(Real, Int, Real eps = Real::eps);
Real sqrt(Real, Real eps = Real::eps);
#endif
