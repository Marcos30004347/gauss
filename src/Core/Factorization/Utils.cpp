#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Debug/Assert.hpp"

#include <cmath>
#include <limits>
#include <random>
#include <vector>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace factorization {

Int comb(Int n, Int k)
{
	return fact(n) / (fact(k) * fact(n - k));
}

Int landauMignotteBound(ast::Expr u, ast::Expr x)
{
	double P;
	Int d, cn;
	Expr p, lc, n, t1, t2;

	p = algebraicExpand(u);

	n = degree(p, x);

	lc = leadCoeff(p, x);

	assert(
		lc.kind() == Kind::Integer,
		"Landau Mignote Bound works on polynomial"
		"with integer coefficients only"
	);

	P = 0;

	d = std::floor(n.value().longValue() * 0.5);

	cn = lc.value();

	// iterate over all factors of u(x)
	while(p != 0)
	{
		lc = leadCoeff(p, x);

		assert(
			lc.kind() == Kind::Integer,
			"Landau Mignote Bound works on polynomial"
			"with integer coefficients only"
		);

		P = P + lc.value().longValue() * lc.value().longValue();

		t1 = power(x, n);

		t2 = mulPoly(lc, t1);

		t1 = subPoly(p, t2);

		p = t1;
	}

	P = std::sqrt(P);

	double B = 0.0;

	B = B + comb(d - 1, std::floor(d.longValue() * 0.5) - 1).longValue() * P;
	B = B + comb(d - 1, std::floor(d.longValue() * 0.5)).longValue() * cn.longValue();

	return std::ceil(B);
}

Int norm(Expr u, Expr L, Expr K, long long i)
{
	if(i == L.size())
	{
		assert(
			u.kind() == Kind::Integer,
			"Polynomial needs to have"
			"integer coefficients in K[L...]"
		);

		return u.value();
	}

	Int k = 0;

	Expr q, p, t, c, n;

	n = degree(u, L[i]);

	p = algebraicExpand(u);

	for(Int e = n.value(); e >= 0; e--)
	{
		c = coeff(u, L[i], e);

		k = max(abs(norm(c, L, K, i + 1)), abs(k));

		t = c * power(L[i], e);

		q = subPoly(p, t);

		p = algebraicExpand(q);
	}

	return k;
}



Int norm(Expr u, Expr x)
{
	Int k = 0;

	Expr q, p, t, c, n;

	n = degree(u, x);

	p = algebraicExpand(u);

	for(Int e = n.value(); e >= 0; e--)
	{
		c = coeff(u, x, e);

		assert(c.kind() == Kind::Integer, "coeffs needs to be integers");

		k = max(abs(c.value()), abs(k));

		t = c * power(x, e);

		q = subPoly(p, t);

		p = algebraicExpand(q);

	}

	return k;
}

Int l1norm(Expr u, Expr L, Expr K, long long i)
{
	if(i == L.size())
	{
		assert(
			u.kind() == Kind::Integer,
			"Polynomial needs to have"
			"integer coefficients in K[L...]"
		);

		return abs(u.value());
	}

	Int k = 0;

	Expr q, p, t, c, n;

	n = degree(u, L[i]);

	p = algebraicExpand(u);

	for(Int e = n.value(); e >= 0; e--)
	{
		c = coeff(u, L[i], e);

		k = abs(norm(c, L, K, i + 1)) + k;

		t = c * power(L[i], e);

		q = subPoly(p, t);

		p = algebraicExpand(q);
	}

	return k;
}


Int l1norm(Expr u, Expr x)
{
	Int k = 0;

	Expr q, p, t, c, n;

	n = degree(u, x);

	p = algebraicExpand(u);

	for(Int e = n.value(); e >= 0; e--)
	{
		c = coeff(u, x, e);

		assert(c.kind() == Kind::Integer, "coeffs needs to be integers");

		k = abs(c.value()) + k;

		t = c*power(x, e);

		q = subPoly(p, t);

		p = algebraicExpand(q);
	}

	return k;
}

Int random(long long min, long long max)
{
	std::random_device dev;

	std::mt19937 rng(dev());

	std::uniform_int_distribution<std::mt19937::result_type> dist(min, max);

	return Int((long long)dist(rng));
}


Expr sortTerms(Expr& F) {
	std::vector<Expr> O = F.operands();

	sort(O);

	return Expr(F.kind(), O);
}

}
