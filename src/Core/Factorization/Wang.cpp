#include "Wang.hpp"
#include "Utils.hpp"
#include "SquareFree.hpp"
#include "Berlekamp.hpp"

#include "Core/Algebra/Set.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Polynomial/Zp.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Calculus/Calculus.hpp"
#include <limits>
#include <random>

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace simplification;

namespace factorization {

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

AST* groundLeadCoeff(AST* f, AST* L)
{
	long i = 0;
	
	AST* p = f->copy();

	AST* t = nullptr;
	
	for(i = 0; i < L->numberOfOperands(); i++)
	{
		t = leadCoeff(p, L->operand(i));
		
		delete p;
		
		p = t;
	}

	return p;	
}

ast::AST* trialDivision(ast::AST* f, ast::AST* F, ast::AST* L, ast::AST* K)
{
	AST *v, *q, *r, *d, *t = list({});

	bool stop = false;

	f = f->copy();

	long i, k;

	for(i = 0; i < F->numberOfOperands(); i++)
	{
		k = 0;
		v = F->operand(i);

		while(!stop)
		{
			d =	recPolyDiv(f, v, L, K);

			q = d->operand(0);
			r = d->operand(1);

			if(r->is(0))
			{
				delete f;
			
				f = q;
			
				k = k + 1;
			}

			stop = !r->is(0);

			delete d;
		}

		t->includeOperand(list({ v->copy(), integer(k) }));
	}

	delete f;

	return t;
}


AST* univariateFactors(AST* f, AST* L, AST* K)
{
	assert(L->numberOfOperands() == 1, "L needs to have just one variable");

	AST *ct, *pr, *n, *x, *lc, *t1, *T;

	x = L->operand(0);

	ct = cont(f, L, K);

	pr = pp(f, ct, L, K);

	n = degree(pr, x);
	
	lc = leadCoeff(f, x);

	if(lc->value() < 0)
	{
		t1 = mul({integer(-1), ct});
		ct = reduceAST(t1);
		
		delete t1;
		
		t1 = mul({integer(-1), pr});
		pr = reduceAST(t1);
		
		delete t1;
	}

	if(n->value() <= 0)
	{
		return list({ ct, list({ integer(1), integer(1) }) });
	}

	if(n->value() == 1)
	{
		return list({ ct, list({ pr, integer(1) }) });
	}

	t1 = pr;

	pr = squareFreePart(t1, L, K);

	delete t1;

	T = list({});
}

AST* factors(AST* f, AST* L, AST* K)
{
	long i;

	AST *t1, *v, *n, *F, *T, *t, *R, *u, *ct, *pr, *lc;

	if(L->numberOfOperands() == 1)
	{
		ct = cont(f, L, K);

		pr = pp(f, ct, L, K);

		n = degree(pr, L->operand(0));
	

	}

	ct = cont(f, L, K);

	pr = pp(f, ct, L, K);

	lc = groundLeadCoeff(pr, L);

	if(lc->kind() == Kind::Integer && lc->value() < 0)
	{
		t1 = mul({integer(-1), ct});
		ct = reduceAST(t1);
		
		delete t1;
		
		t1 = mul({integer(-1), pr});
		pr = reduceAST(t1);
	
		delete t1;
	}

	v = cont(pr, L, K);
	
	t1 = pr;

	pr = pp(t1, v, L, K);

	delete t1;

	F = list({});

	n = degree(pr, L->operand(0));
	
	if(n->value() > 0)
	{
		t = squareFreePart(pr, L, K);

		T = factorsWang(t, L, K);
	
		t1 = trialDivision(f, T, L, K);

		for(i = 0; i < t1->numberOfOperands(); i++)
		{
			F->includeOperand(t1->operand(i));
			
			t1->removeOperand(i);
		}	
	
		delete t1;
	
		delete T;
		delete t;
	}

	R = rest(L);

	T = factors(v, R, K);

	u = T->operand(1);

	while(u->numberOfOperands())
	{
		F->includeOperand(list({ u->operand(0)->operand(0), u->operand(0)->operand(1) }),0);
		
		u->operand(0)->removeOperand(0L);
		u->operand(0)->removeOperand(1L);

		u->deleteOperand(0);
	}

	delete T;

	delete R;

	delete n;

	delete pr;
	delete lc;

	return list({ ct, F });
}


AST* factorsWang(AST* f, AST* L, AST* K)
{
	return nullptr;
}


// bool getEvaluationPoint(ast::AST* F, ast::AST* L, ast::AST* K, long mod)
// {
// 	long c, i, t, *a;

// 	AST* s = set({});

// 	t = L->numberOfOperands();
// 	a = new long[t];

// 	while(true)
// 	{
// 		c = 0;

// 		while(c < 3)
// 		{
// 			for(i = 0; i < t; i++)
// 			{
// 				a[i] = random() % mod;
// 			}


			

// 		}

// 		mod = mod + 1;
// 	}

// 	delete a;
// }

}
