#include "MathSystem/Factorization/Utils.hpp"
#include "MathSystem/Polynomial/Polynomial.hpp"

#include <climits>
#include <cmath>
#include <cstddef>
#include <limits>

//TODO: move random stuff to a single file
#ifndef WASM_BUILD
#include <random>
#endif

#include <vector>


using namespace alg;
using namespace polynomial;

namespace factorization {

Int comb(Int n, Int k)
{
	return fact(n) / (fact(k) * fact(n - k));
}

// Int landauMignotteBound(expr u, expr x)
// {
// 	double P;
// 	Int d, cn;
// 	expr p, lc, n, t1, t2;

// 	p = expand(u);

// 	n = degree(p, x);

// 	lc = leadCoeff(p, x);

// 	assert(lc.kind() == kind::INT);

// 	P = 0;

// 	d = std::floor(n.value().longValue() * 0.5);

// 	cn = lc.value();

// 	// iterate over all factors of u(x)
// 	while(p != 0)
// 	{
// 		lc = leadCoeff(p, x);

// 		assert(
// 			lc.kind() == kind::INT		);

// 		P = P + lc.value().longValue() * lc.value().longValue();

// 		t1 = pow(x, n);

// 		t2 = mulPoly(lc, t1);

// 		t1 = subPoly(p, t2);

// 		p = t1;
// 	}

// 	P = std::sqrt(P);

// 	double B = 0.0;

// 	B = B + comb(d - 1, std::floor(d.longValue() * 0.5) - 1).longValue() * P;
// 	B = B + comb(d - 1, std::floor(d.longValue() * 0.5)).longValue() * cn.longValue();

// 	return std::ceil(B);
// }

// Int norm(expr u, expr L, expr K, size_t i)
// {
// 	if(i == L.size())
// 	{
// 		assert(
// 			u.kind() == kind::INT		);

// 		return u.value();
// 	}

// 	Int k = 0;

// 	expr q, p, t, c, n;

// 	n = degree(u, L[i]);

// 	p = expand(u);

// 	for(Int e = n.value(); e >= 0; e--)
// 	{
// 		c = coeff(u, L[i], e);

// 		k = max(abs(norm(c, L, K, i + 1)), abs(k));

// 		t = c * pow(L[i], e);

// 		q = subPoly(p, t);

// 		p = expand(q);
// 	}

// 	return k;
// }



// Int norm(expr u, expr x)
// {
// 	Int k = 0;

// 	expr q, p, t, c, n;

// 	n = degree(u, x);

// 	p = expand(u);

// 	for(Int e = n.value(); e >= 0; e--)
// 	{
// 		c = coeff(u, x, e);

// 		assert(c.kind() == kind::INT);

// 		k = max(abs(c.value()), abs(k));

// 		t = c * pow(x, e);

// 		q = subPoly(p, t);

// 		p = expand(q);

// 	}

// 	return k;
// }

Int normPolyExpr(expr u)
{
	if(u.kind() == kind::INT) {
		return u.value();
	}

	Int k = 0;

	expr q, p, t, c, n;

	for(Int i = 0; i < u.size(); i++) {
		assert(u[i].kind() == kind::MUL && u[i].size() == 2);
		assert(u[i][0].kind() == kind::INT);
		k = max(abs(u[i][0].value()), abs(k));
	}

	return k;
}



Int normPolyExpr(expr u, expr L, expr K)
{
	assert(K.identifier() == "Z");

	if(L.size() == 0) {
		assert(u.kind() == kind::INT);
		return abs(u.value());
	}

	if(u.kind() == kind::INT) {
		return abs(u.value());
	}

	Int k = 0;

	expr R = rest(L);

	for(Int i = 0; i < u.size(); i++) {
		k = max(abs(normPolyExpr(u[i][0], R, K)), k);
	}

	return k;
}

// Int l1norm(expr u, expr L, expr K, size_t i)
// {
// 	if(i == L.size())
// 	{
// 		assert(
// 			u.kind() == kind::INT		);

// 		return abs(u.value());
// 	}

// 	Int k = 0;

// 	expr q, p, t, c, n;

// 	n = degree(u, L[i]);

// 	p = expand(u);

// 	for(Int e = n.value(); e >= 0; e--)
// 	{
// 		c = coeff(u, L[i], e);

// 		k = abs(norm(c, L, K, i + 1)) + k;

// 		t = c * pow(L[i], e);

// 		q = subPoly(p, t);

// 		p = expand(q);
// 	}

// 	return k;
// }


// Int l1norm(expr u, expr x)
// {
// 	Int k = 0;

// 	expr q, p, t, c, n;

// 	n = degree(u, x);

// 	p = expand(u);

// 	for(Int e = n.value(); e >= 0; e--)
// 	{
// 		c = coeff(u, x, e);

// 		assert(c.kind() == kind::INT);

// 		k = abs(c.value()) + k;

// 		t = c*pow(x, e);

// 		q = subPoly(p, t);

// 		p = expand(q);
// 	}

// 	return k;
// }



Int l1normPolyExpr(expr u)
{
	Int k = 0;

	expr q, p, t, c, n;

	for(Int i = 0; i < u.size(); i++) {
		assert(u[i].kind() == kind::MUL && u[i].size() == 2);
		assert(u[i][0].kind() == kind::INT);
		k = abs(u[i][0].value()) + k;
	}

	return k;
}


Int random(long long min, long long max)
{
	std::random_device dev;

	std::mt19937 rng(dev());

	std::uniform_int_distribution<std::mt19937::result_type> dist(min, max);

	return (long long)dist(rng);
}

expr sortTerms(expr& F) {
	expr k = F;

	if(is(&k, kind::LIST)) {
		k.expr_list->sortMembers();
		return k;
	}

	if(is(&k, kind::SET)) {
		return k;
	}

	sort(&k);

	return k;
}

}
