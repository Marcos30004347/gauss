#include "Wang.hpp"

#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include <limits>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace factorization {

long norm(AST* u, AST* L, AST* K, long i = 0)
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
	
		c = coeff(u, L->operand(i), n);
	
		k = std::max(std::abs(norm(c, L, K, i + 1)), k);
	
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

long gcd(long  a, long  b) {
	if (a == 0)
	{
		return b;
	}

	return gcd(b % a, a);
}

bool nondivisors(long G, AST* F, long d, AST* L, AST* K, long* p)
{
	assert(p != nullptr, "Array d[i] cant be null!");
	assert(G != 0, "G needs to be different from zero!");
	assert(d != 0, "c needs to be different from zero!");

	long i, j, k, q, r;

	AST *Fi;
	
	k = F->numberOfOperands();

	long* x = new long[k + 1];

	x[0] = G * d;

	for(i = 1; i <= k; i++)
	{
		Fi = F->operand(i - 1);
	
		q = norm(Fi, L, K);
	
		for(j = i - 1; j >= 0; j--)
		{
			while(r != 1)
			{
				r = x[j];

				r = gcd(r, q);
				q = q / r;
			}

			if(q == 1)
			{
				return false;
			}
		}
	
		x[i] = q;
	}	

	for(i = 1; i <= k; i++)
	{
		p[i - 1] = x[i];
	}

	delete x;
	
	return true;
}


}
