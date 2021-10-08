#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Zp.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Debug/Assert.hpp"

#include <cmath>
#include <limits>
#include <random>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace factorization {

long fact(long p)
{
	long f = 1;

	while(p > 0)
	{
		f = f * p;
		p = p - 1;
	}

	return f;
}

long comb(long n, long k)
{
	return fact(n) / (fact(k) * fact(n - k));
}

long landauMignotteBound(ast::AST* u, ast::AST* x)
{
	double P;
	long d, cn;
	AST *p, *lc, *n, *t1, *t2;

	p = algebraicExpand(u);

	n = degree(p, x);

	lc = leadCoeff(p, x);

	assert(
		lc->kind() == Kind::Integer, 
		"Landau Mignote Bound works on polynomial"
		"with integer coefficients only"
	);

	P = 0;

	d = std::floor(n->value() * 0.5);

	cn = lc->value();
	
	delete n;
	delete lc;

	// iterate over all factors of u(x)
	while(p->isNot(0))
	{
		lc = leadCoeff(p, x);

		assert(
			lc->kind() == Kind::Integer, 
			"Landau Mignote Bound works on polynomial"
			"with integer coefficients only"
		);

		P = P + lc->value() * lc->value();

		t1 = power(x->copy(), n->copy());
	
		t2 = mulPoly(lc, t1);
		
		delete t1;

		t1 = subPoly(p, t2);

		delete p;

		p = t1;

		delete lc;
	}	

	delete p;

	P = std::sqrt(P);

	double B = 0.0;

	B += comb(d - 1, std::floor(d * 0.5) - 1) * P;
	B += comb(d - 1, std::floor(d * 0.5)) * cn;

	return std::ceil(B);
}

// AST* repeatSquaring(AST* a, AST* x, long n, long m)
// {
// 	AST *t1, *t2, *b[64];

// 	long v = n;

// 	long k = 0;

// 	while (v >>= 1) k++;

// 	b[k] = a->copy();

// 	for (long i = k - 1; i>= 0; i--)
// 	{
// 		t1 = mulPoly(b[i + 1], b[i + 1]);

// 		t2 = reduceAST(t1);

// 		delete t1;

// 		t1 = sZp(t2, x, m);
	
// 		delete t2;

// 		if(n & (1 << i))
// 		{
// 			t2 = mulPoly(t1, a);

// 			b[i] = sZp(t2, x, m);

// 			delete t2;
// 		}
// 		else
// 		{
// 			b[i] = sZp(t1, x, m);
// 		}

// 		delete t1;
// 	}

// 	for(int i = 1; i<=k; i++) delete b[i];
	
// 	return b[0];
// }



long norm(AST* u, AST* L, AST* K, long i)
{
	if(i == L->numberOfOperands())
	{
		assert(
			u->kind() == Kind::Integer, 
			"Polynomial needs to have"
			"integer coefficients in K[L...]"
		);

		return u->value();
	}

	long k = 0;

	AST *q, *p, *t, *e, *c, *n;
	
	n = degree(u, L->operand(i));
	
	p = algebraicExpand(u);

	for(long j = n->value(); j >= 0; j--)
	{
		e = integer(j);
	
		c = coeff(u, L->operand(i), e);

		k = std::max(std::abs(norm(c, L, K, i + 1)), std::abs(k));
	
		t = mul({c, power(L->operand(i)->copy(), e)});
	
		q = subPoly(p, t);

		delete p;
	
		p = algebraicExpand(q);	
	
		delete t;
	
		delete q;
	}

	delete p;

	delete n;

	return k;
}


long l1norm(AST* u, AST* L, AST* K, long i)
{
	if(i == L->numberOfOperands())
	{
		assert(
			u->kind() == Kind::Integer, 
			"Polynomial needs to have"
			"integer coefficients in K[L...]"
		);

		return std::abs(u->value());
	}

	long k = 0;

	AST *q, *p, *t, *e, *c, *n;
	
	n = degree(u, L->operand(i));
	
	p = algebraicExpand(u);

	for(long j = n->value(); j >= 0; j--)
	{
		e = integer(j);
	
		c = coeff(u, L->operand(i), n);
	
		k = std::abs(norm(c, L, K, i + 1)) + k;
	
		t = mul({c, power(L->operand(i)->copy(), e)});
	
		q = subPoly(p, t);

		delete p;
	
		p = algebraicExpand(q);	
	
		delete t;
	
		delete q;
	}

	delete p;

	delete n;

	return k;
}


long random(long min, long max)
{
	std::random_device dev;
	
	std::mt19937 rng(dev());
	
	std::uniform_int_distribution<std::mt19937::result_type> dist(min, max);
	
	return dist(rng);
}

}
