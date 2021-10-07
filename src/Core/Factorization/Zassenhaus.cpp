#include "Zassenhaus.hpp"
#include "Berlekamp.hpp"
#include "Utils.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Polynomial/Zp.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Calculus/Calculus.hpp"
#include <cmath>

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace simplification;

namespace factorization {

long log(double x, double base) 
{
	return (long)(std::log(x) / std::log(base));
}

// From modern computer algebra by Gathen
AST* zassenhaus(AST* f, AST* L, AST* K)
{
	assert(L->numberOfOperands() == 1, "L should have one symbol");
	assert(K->identifier() == "Z", "Zassenhaus work only on the integers");

	bool stop = false;

	long i, l, p, A, B, C, gamma, bound, gcd;

	AST *x, *n, *c, *z, *b, *a, *F, *D, *E, *q, *H;

	x = L->operand(0);

	n = degree(f, x);

	if(n->is(1))
	{
		delete n;
	
		return list({ f->copy() });
	}

	q = integer(p);

	A = norm(f, L, K);

	b = leadCoeff(f, x);

	B = std::abs(std::sqrt(n->value() + 1) * std::pow(2, n->value()) * A * b->value());

	C = std::pow(n->value() + 1, 2 * n->value()) * std::pow(A, 2 * n->value() - 1);

	gamma = std::ceil(2 * log2(C));

	a = list({});

	// choose a prime number p such that f be square free in Zp[x]
	// and such that p dont divide lc(f)
	for(i = 1; primes[i] <= 2 * gamma * std::log(gamma); i++)
	{
		p = primes[i];
		
		if(b->value() % p == 0)
		{
			continue;
		}

		F = Zp(f, x, p);
	
		D = derivate(F, x);

		E = Zp(D, x, p);
		
		delete D;

		D = gcdGPE_Zp(F, E, x, p);

		gcd = D->value();
	
		delete E;		
		delete F;
		delete D;
			
		if(b->value() % p && gcd == 1)
		{
			break;
		}	
	}

	l = log(2*B + 1, p);

	H = berlekampFactors(f, x, q);

	delete z;
	delete c;
	delete b;
}


}
