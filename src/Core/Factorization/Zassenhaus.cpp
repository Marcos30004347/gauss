#include "Utils.hpp"
#include "Hensel.hpp"
#include "Berlekamp.hpp"
#include "Zassenhaus.hpp"
#include "SquareFree.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/GaloisField/GaloisField.hpp"

#include <cmath>

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace galoisField;
using namespace simplification;

namespace factorization {

void subsetsRec(AST* arr, AST* data, AST* s, int start, int end, int index, int r)
{
	AST* c;

	long i, j;

	if (index == r)
	{
		s->includeOperand(set({}));

		c = s->operand(s->numberOfOperands() - 1);

		for(j = 0; j < r; j++)
		{
			c->includeOperand(data->operand(j));
		}
	}
	else
	{
		for (i = start; i <= end && end - i + 1 >= r - index; i++)
		{
			data->includeOperand(arr->operand(i)->copy());
			subsetsRec(arr, data, s, i+1, end, index+1, r);
			data->removeOperand(data->numberOfOperands() - 1);
		}
	}
}

AST* subset(AST* s, long r)
{
	long n = s->numberOfOperands();
	
	AST* d, *res;

	d = set({});

	res = set({});

	subsetsRec(s, d, res, 0, n - 1, 0, r);

	delete d;
	
	return res;
}

// from Algorithms for Computer Algebra Geddes
AST* cantorZassenhausDDF(AST* v, AST* x, long p)
{
	long i;

	AST *h, *f, *t, *G, *g, *n;
	
	i = 1;

	h = x->copy();

	f = v->copy();

	G = list({});

	g = integer(1);

	n = degree(f, x);

	while(n->value() >= 2*i)
	{
		t = powModPolyGf(h, f, x, p, p, true);
	
		delete h;

		h = t;
		
		t = subPolyGf(h, x, x, p, true);

		delete g;

		g = gcdPolyGf(t, f, x, p, true);
	
		delete t;
	
		if(g->isNot(1))
		{
			G->includeOperand(list({ g->copy(), integer(i) }));

			t = quoPolyGf(f, g, x, p, true);

			delete f;

			f = t;

			t = remPolyGf(h, f, x, p, true);
		
			delete h;
			
			h = t;
		}

		delete n;
	
		n = degree(f, x);

		i = i + 1;
	};

	if(f->isNot(1))
	{
		G->includeOperand(list({ f->copy(), degree(f, x) }));
	}

	delete g;
	delete h;
	delete f;
	delete n;

	return G;
}

// from Algorithms for Computer Algebra Geddes
AST* cantorZassenhausEDF(AST* a, AST* x, long n, long p)
{
	long m, i;

	AST *g, *da, *F, *v, *h, *k, *f1, *f2, *t;

	da = degree(a, x);
	
	if(da->value() <= n)
	{
		delete da;

		return list({ a->copy() });
	}

	m = da->value() / n;

	F = list({ a->copy() });

	while(F->numberOfOperands() < m)
	{
		v = randPolyGf(2*n - 1, x, p);

		if(p == 2)
		{
			t = v->copy();

			for(i = 0; i < std::pow(2, n * m - 1); i++)
			{
				h = powModPolyGf(t, a, x, 2, p, true);
				
				k = addPolyGf(v, h, x, p, true);
				
				delete v;

				delete t;

				t = h;
				v = k;
			}

			delete t;
		}
		else
		{
			h = powModPolyGf(v, a, x, (std::pow(p, n) - 1) / 2, p, true);

			delete v;
		
			v = h;

			k = integer(1);

			h = subPolyGf(v, k, x, p, true);
		
			delete v;
			
			delete k;
		
			v = h;
		}
	
		g = gcdPolyGf(a, v, x, p, true);

		if(g->isNot(1) && !g->match(a))
		{
			k = quoPolyGf(a, g, x, p, true);
			
			f1 = cantorZassenhausEDF(g, x, n, p);
			f2 = cantorZassenhausEDF(k, x, n, p);
			
			delete k;
			
			delete F;
		
			F = append(f1, f2);

			delete f1;
			delete f2;
		}

		delete v;
		delete g;
	}

	delete da;

	return F;
}

// From modern computer algebra by Gathen
AST* zassenhaus(AST* f, AST* x)
{
	bool stop = false;

	long s, i, j, l, p, A, B, C, gamma, gcd;

	AST *g, *n, *b, *F, *D, *E, *q, *H, *Z, *G, *T, *S, *M, *u, *v, *gi, *L, *K;

	L = list({ x->copy() });

	K = symbol("Z");

	f = f->copy();

	// zero = integer(0);


	n = degree(f, x);

	if(n->is(1))
	{
		delete n;
	
		return list({ f->copy() });
	}

	A = norm(f, x);
	printf("f = %s\n", f->toString().c_str());
	b = leadCoeff(f, x);

	B = std::sqrt(n->value() + 1) * std::pow(2, n->value()) * A * b->value();

	C = std::pow(n->value() + 1, 2 * n->value()) * std::pow(A, 2 * n->value() - 1);

	gamma = std::ceil(2 * log2(C));
	
	printf("B = %li\n", B);
	printf("A = %li\n", A);

	printf("b = %s\n", b->toString().c_str());

	// printf("%li\n",  2 * gamma * std::log(gamma));

	// choose a prime number p such that f be square free in Zp[x]
	// and such that p dont divide lc(f)
	for(i = 1; primes[i] <= 2 * gamma * std::log(gamma); i++)
	{
		p = primes[i];

		if(b->value() % p == 0)
		{
			continue;
		}

		F = gf(f, x, p, true);
	
		D = derivate(F, x);

		E = gf(D, x, p, true);
		
		delete D;

		D = gcdPolyGf(F, E, x, p, false);

		gcd = D->value();
	
		delete E;		
		delete F;
		delete D;
			
		if(b->value() % p && gcd == 1)
		{
			break;
		}	
	}

	q = integer(p);

	printf("p = %li\n", p);

	l = std::ceil(std::log(2*B + 1) / std::log(p));

	printf("l = %li\n", l);

	v = pp(f, L, K);

	u = gf(v, x, p, true);

	Z = berlekampFactors(u, x, q);

	printf("berlekamp = %s\n", Z->toString().c_str());

	printf("FACTORING  = %s\n",  u->toString().c_str());
	printf("FACTORING  = %s\n",  v->toString().c_str());

	g = multifactorHenselLifting(f, Z, x, p, l);

	printf("hensel = %s\n", g->toString().c_str());

	T = set({});
	
	for(i = 0; i < g->numberOfOperands(); i++)
	{
		T->includeOperand(integer(i));
	}

	F = list({});

	s = 1;

	printf("T = %s\n", T->toString().c_str());

	while(2*s <= T->numberOfOperands())
	{
		stop = false;

		M = subset(T, s);
	
		for(j = 0; j < M->numberOfOperands(); j++)
		{
			S = M->operand(j);
		
			H = mul({b->copy()});
			G = mul({b->copy()});

			Z = difference(T, S);
		
			for(i = 0; i < S->numberOfOperands(); i++)
			{
				gi = g->operand(S->operand(i)->value());
				
				G->includeOperand(gi->copy());
			}
		
			for(i = 0; i < Z->numberOfOperands(); i++)
			{
				gi = g->operand(Z->operand(i)->value());
				
				H->includeOperand(gi->copy());
			}
		
			u = gf(G, x, std::pow(p, l), true);
			v = gf(H, x, std::pow(p, l), true);
		
			delete G;
			delete H;
			
			G = u;
			H = v;

			if(norm(G, x) > std::pow(p, l) / 2)
			{
				delete Z;
				delete G;
				delete H;
			
				continue;
			}

			if(norm(H, x) > std::pow(p, l) / 2)
			{
				delete Z;
				delete G;
				delete H;
			
				continue;
			}

			if(l1norm(G, x) * l1norm(H, x) <= B)
			{
				delete T;

				T = Z;

				F->includeOperand(pp(G, L, K));

				delete f;

				f = pp(H, L, K);

				delete b;

				b = leadCoeff(f, x);

				stop = true;

				break;
			}
		}
	
		delete M;
	
		if(!stop) s = s + 1;
	}

	delete L;
	delete K;

	return F;
}


}
