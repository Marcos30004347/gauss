#include "Utils.hpp"
#include "Hensel.hpp"
#include "Berlekamp.hpp"
#include "Zassenhaus.hpp"
#include "SquareFree.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
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

double log(double x, double base) 
{
	return std::log(x) / std::log(base);
}

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

// AST* gcd_Zp(AST* f, AST* g, AST* L, AST* p)
// {
// 	AST *t, *a, *b;

// 	a = f->copy();
// 	b = g->copy();

// 	while(b->isNot(0))
// 	{
// 		t = a;
		
// 		a = b;
		
// 		b = remainderGPE_Zp(t, b, L->operand(0), p->value());
		
// 		delete t;
// 	}

// 	t = monic_Zp(a, L->operand(0), p->value());
	
// 	delete a;
// 	delete b;

// 	return t;
// }

// AST* gcd_sZp(AST* f, AST* g, AST* L, AST* p)
// {
// 	AST *t, *a, *b;

// 	a = f->copy();
// 	b = g->copy();

// 	while(b->isNot(0))
// 	{
// 		t = a;
		
// 		a = b;
		
// 		b = remainderGPE_sZp(t, b, L->operand(0), p->value());
		
// 		delete t;
// 	}

// 	t = monic_sZp(a, L->operand(0), p->value());
	
// 	delete a;
// 	delete b;

// 	return t;
// }


AST* divGf(AST* a, AST* b, AST* x, long p, bool symmetric = false)
{
	AST* da = degree(a, x);
	AST* db = degree(b, x);

	if(da->value() < db->value())
	{
		delete da;
		delete db;

		return list({ integer(0), a->copy() });
	}

	long k, j, s, e, lb, d;

	AST *dq, *dr, *q, *r;
	AST *t1, *t2, *t3, *ex;
	
	// a = a->copy();
	// b = b->copy();

	AST** A = new AST*[da->value() + 1];
	AST** B = new AST*[db->value() + 1];

	for(k = da->value(); k >= 0; k--)
	{
		ex = integer(k);

		A[k] = coeff(a, x, ex);
	
		delete ex;
	}

	for(k = db->value(); k >= 0; k--)
	{
		ex = integer(k);
	
		B[k] = coeff(b, x, ex);
	
		delete ex;
	}


	dq = integer(da->value() - db->value());
	dr = integer(db->value() - 1);

	t1 = leadCoeff(b, x);	

	lb = modInverse(t1->value(), p);

	delete t1;

	for(k = da->value(); k >= 0; k--)
	{
		// ex = integer(k);
		// t1 = coeff(a, x, ex);
		t1 = A[k]->copy();
	
		s = std::max(0L, k - dq->value());
		e = std::min(dr->value(), k);
	
		// delete ex;
		
		for(j = s; j <= e; j++)
		{
			// ex = integer(j);
			// cb = coeff(b, x, ex);

			// delete ex;
		
			// ex = integer(k - j + db->value());

			// ca = coeff(a, x, ex);
	
			// delete ex;
		
			// t2 = mulPoly(cb, ca);

			t2 = mulPoly(B[j], A[k - j + db->value()]);
			
			// delete ca;
			// delete cb;

			t3 = subPoly(t1, t2);

			delete t1;
			delete t2;

			t1 = t3;
		}

		t3 = reduceAST(t1);
		
		delete t1;
		
		t1 = t3;	

		if(symmetric)
		{
			t2 = integer(sZp(t1->value(), p));
		}
		else 
		{
			t2 = integer(Zp(t1->value(), p));
		}

		delete t1;
		
		t1 = t2;
		
		if(t1->value() < 0)
		{
			t2 = integer(t1->value() + p);

			delete t1;

			t1 = t2;
		}

		if(da->value() - k <= dq->value())
		{

			t3 = integer(lb);
		
			t2 = mulPoly(t1, t3);
		
			delete t1;
			delete t3;
		
			t1 = reduceAST(t2);

			if(symmetric)
			{
				t2 = integer(sZp(t1->value(), p));
			}
			else 
			{
				t2 = integer(Zp(t1->value(), p));
			}
		
			delete t1;
		
			t1 = t2;
		}
	
		
		// ex = integer(k);
		
		// ca = coeff(a, x, ex);
		
		// t3 = power(x->copy(), ex->copy());
		
		// t4 = mulPoly(ca, t3);

		// delete t3;

		// t3 = subPoly(a, t4);

		// delete a;

		// a = algebraicExpand(t3);
		
		// delete t3;

	
		// t3 = power(x->copy(), ex->copy());
		
		// t4 = mulPoly(t1, t3);
	
		// delete t3;
	
		// t3 = addPoly(a, t4);

		// delete a;

		// a = algebraicExpand(t3);
		
		// delete ex;
		// delete ca;
		delete A[k];
		
		A[k] = t1;
	
	}

	q = add({integer(0)});
	r = add({integer(0)});

	d = 0;
	for(k = da->value() - dq->value(); k <= da->value(); k++)
	{
		// ex = integer(k);
	
		// q->includeOperand(
		// 	mul({
		// 		coeff(a, x, ex),
		// 		power(x->copy(), integer(d))
		// 	})
		// );
		q->includeOperand(
			mul({
				A[k]->copy(),
				power(x->copy(), integer(d))
			})
		);

		// delete ex;

		d = d + 1;
	}

	d = 0;
	for(k = 0; k <= dr->value(); k++)
	{

		// ex = integer(k);

		// r->includeOperand(
		// 	mul({
		// 		coeff(a, x, ex),
		// 		power(x->copy(), integer(d))
		// 	})
		// );

		// delete ex;
		r->includeOperand(
			mul({
				A[k]->copy(),
				power(x->copy(), integer(d))
			})
		);

		d = d + 1;
	}

	delete da;
	delete db;

	delete dq;
	delete dr;
	
	if(symmetric)
	{
		t1 = sZp(q, x, p);
		t2 = sZp(r, x, p);
	}
	else 
	{
		t1 = Zp(q, x, p);
		t2 = Zp(r, x, p);	
	}

	delete q;
	delete r;

	return list({ t1, t2 });
}

AST* remGf(AST* a, AST* b, AST* x, long p, bool symmetric = false)
{
	AST* d = divGf(a, b, x, p, symmetric);

	AST* r = d->operand(1L);

	d->removeOperand(1L);

	delete d;

	return r;
}

AST* quoGf(AST* a, AST* b, AST* x, long p, bool symmetric = false)
{
	AST* d = divGf(a, b, x, p, symmetric);

	AST* q = d->operand(0L);

	d->removeOperand(0L);

	delete d;

	return q;
}

AST* monicGf(AST* f, AST* x, long p, bool symmetric = false)
{
	if(f->is(0))
	{
		return integer(0);
	}

	AST* lc = leadCoeff(f, x);

	AST* F = quoGf(f, lc, x, p, symmetric);

	return list({ lc, F });
}

AST* gcdGf(AST* a, AST* b, AST* x, long p, bool symmetric = false)
{
	AST* da = degree(a, x);
	AST* db = degree(b, x);
	
	if(da->kind() == Kind::MinusInfinity || db->value() > da->value())
	{
		delete da;
		delete db;
		
		return gcdGf(b, a, x, p, symmetric);
	}

	AST *t1;

	a = a->copy();
	b = b->copy();

	while(b->isNot(0) && db->kind() != Kind::MinusInfinity && db->value() >= 0)
	{
		t1 = a;
		
		a = b;

		b = remGf(t1, b, x, p, symmetric);

		delete t1;

		delete db;

		db = degree(b, x);
	}

	delete b;

	delete da;
	delete db;

	// printf("RETURNING %s\n", a->toString().c_str());
	b = monicGf(a, x, p, symmetric);
	// printf("RETURNING %s\n", b->toString().c_str());

	delete a;

	a = b->operand(1L);
	
	b->removeOperand(1L);
	
	delete b;
	
	return a;
}

AST* addGf(AST* f, AST* g, AST* x, long p, bool symmetric = false)
{
	AST *t, *u;

	u = addPoly(f, g);

	if(symmetric)
	{
		t = sZp(u, x, p);
	}
	else
	{
		t = Zp(u, x, p);
	}

	delete u;

	return t;
}

AST* subGf(AST* f, AST* g, AST* x, long p, bool symmetric = false)
{
	AST *t, *u;

	u = subPoly(f, g);

	if(symmetric)
	{
		t = sZp(u, x, p);
	}
	else
	{
		t = Zp(u, x, p);
	}

	delete u;

	return t;
}

AST* mulGf(AST* f, AST* g, AST* x, long p, bool symmetric = false)
{
	AST *t, *u;

	u = mulPoly(f, g);

	if(symmetric)
	{
		t = sZp(u, x, p);
	}
	else
	{
		t = Zp(u, x, p);
	}

	delete u;

	return t;
}

AST* powGf(AST* f, AST* g, AST* x, long n, long p, bool symmetric = false)
{
	AST *a, *b, *t;

	if(n == 0) return integer(1);
	
	a = f->copy();

	b = integer(1);

	while(n > 1)
	{
		if(n % 2 == 0)
		{
			t = mulGf(a, a, x, p, symmetric);

			a = remGf(t, g, x, p, symmetric);
			
			delete t;
		
			n = n / 2;
		}
		else
		{
			t = mulGf(a, b, x, p, symmetric);
		
			b = remGf(t, g, x, p, symmetric);

			delete t;

			t = mulGf(a, a, x, p, symmetric);

			a = remGf(t, g, x, p, symmetric);

			delete t;
		
			n = (n - 1) / 2;
		}
	}


	t = mulGf(a, b, x, p, symmetric);

	delete a;
	delete b;

	a = remGf(t, g, x, p, symmetric);

	delete t;

	return a;
}

AST* randGf(long d, AST* x, long p)
{
	if(d == 0 || d == 1)
	{
		return integer(random(0, p));
	} 
	AST* r = add({});

	for(long i = d; i >= 2; i--)
	{
		r->includeOperand(mul({
			integer(random(0, p)),
			power(x->copy(), integer(i))
		}));
	}
	long k;

	k = random(0, p);
	
	if(k != 0)
	{
		r->includeOperand(mul({
			integer(k),
			x->copy()
		}));
	}

	k = random(0, p);

	if(k != 0)
	{
		r->includeOperand(integer(k));
	}

	return r;
}

// from Algorithms for Computer Algebra Geddes
AST* distinctDegreeFactorization(AST* v, AST* L, AST* K, AST* q)
{
	assert(K->identifier() == "Z", "distinct degree only works on Zp[x]");
	
	long i, p;

	AST *x, *h, *f, *t, *G, *g, *n;
	
	x = L->operand(0);

	i = 1;

	h = x->copy();

	f = v->copy();

	G = list({});

	g = integer(1);

	p = q->value();

	n = degree(f, x);

	while(n->value() >= 2*i)
	{
		t = powGf(h, f, x, p, p, true);
	
		delete h;

		h = t;
		
		t = subGf(h, x, x, p, true);

		delete g;

		g = gcdGf(t, f, x, p, true);
	
		if(g->isNot(1))
		{
			G->includeOperand(list({ g->copy(), integer(i) }));

			t = quoGf(f, g, x, p, true);

			delete f;

			f = t;

			t = remGf(h, f, x, p, true);
		
			delete h;
			
			h = t;
		}
	
		delete n;
	
		n = degree(f, x);

		i = i + 1;
	};

	if(f->isNot(1))
	{
		G->includeOperand(list({f->copy(), degree(f, x)}));
	}

	delete h;
	delete f;
	delete n;

	return G;
}

// from Algorithms for Computer Algebra Geddes
AST* equalDegreeFactorization(AST* a, AST* x, long n, long p)
{
	long m, i;

	AST *g, *da, *F, *v, *h, *k, *f1, *f2;

	da = degree(a, x);
	
	if(da->value() <= n)
	{
		return list({ a->copy() });
	}

	m = da->value() / n;

	F = list({ a->copy() });

	while(F->numberOfOperands() < m)
	{
		v = randGf(2*n - 1, x, p);

		if(p == 2)
		{
			for(i = 0; i < std::pow(2, n * m - 1); i++)
			{
				h = powGf(v, a, x, 2, p, true);
				
				k = addGf(v, h, x, p, true);
				
				delete v;
				
				v = k;
			}
		}
		else
		{
			h = powGf(v, a, x, (std::pow(p, n) - 1) / 2, p, true);
		
			delete v;
		
			v = h;
			k = integer(1);
			h = subGf(v, k, x, p, true);
		
			delete v;
			delete k;
		
			v = h;
		}

		g = gcdGf(a, v, x, p, true);

		if(g->isNot(1) && !g->match(a))
		{
			k = quoGf(a, g, x, p, true);
			
			f1 = equalDegreeFactorization(g, x, n, p);
			f2 = equalDegreeFactorization(k, x, n, p);
			
			delete k;
		
			while(f1->numberOfOperands() > 0)
			{
				F->includeOperand(f1->operand(0L));
				f1->removeOperand(0L);
			}
		
			while(f2->numberOfOperands() > 0)
			{
				F->includeOperand(f2->operand(0L));
				f2->removeOperand(0L);
			}
		
			delete f1;
			delete f2;
		}
	}

	return F;



	// F = list({ f->copy() });

	// a = randGf(2*n->value() - 1, x, p);

	// g = gcdGf(a, f, x, p, true);

	// if(g->isNot(1))
	// {
	// 	return list({ g });
	// }

	// b = powGf(a, f, x, std::pow(p, d) -1 , p, true);

	return nullptr;
}

// // assert that g* and h* in Z[x] will have max-norm at most
// // p^l/2 satisfying g* = b * g[i] | i in S 
// // and h* = b * h[i] | i not in S and
// bool isInvalidSet(AST* g, AST* S, long z, long b, long p, long l,  AST* L,  AST* K)
// {
// 	long i, q;

// 	AST *G, *gi, *c, *zero, *t1, *t2;
	
// 	if(b == 1)
// 	{
// 		q = 1;
	
// 		zero = integer(0);
	
// 		for(i = 0; i < S->numberOfOperands(); i++)
// 		{
// 			gi = g->operand(S->operand(i)->value());

// 			c = coeff(gi, L->operand(0), zero);

// 			q = c->value() * q;
// 		}
	
// 		q = q % (long)std::pow(p, l);

// 		if(q > std::pow(p, l) / 2)
// 		{
// 			q = q - std::pow(p, l);
// 		}
	
// 		if(q == 0)
// 		{
// 			return false;
// 		}
	
// 		return z % q != 0;
// 	}

// 	zero = integer(0);

// 	t1 = mul({ integer(b) });
	
// 	for(i = 0; i < S->numberOfOperands(); i++)
// 	{
// 		gi = g->operand(S->operand(i)->value());
// 		t1->includeOperand(gi->copy());
// 	}

// 	t2 = sZp(t1, L->operand(0), std::pow(p, l));
	
// 	delete t1;

// 	G = pp(t2, L, K);
	
// 	delete t2;

// 	c = coeff(G, L->operand(0), zero);

// 	q = c->value();
	
// 	delete c;
	
// 	delete G;

// 	return q != 0 && z % q != 0;
// }

// From modern computer algebra by Gathen
AST* zassenhaus(AST* f, AST* L, AST* K)
{
	assert(L->numberOfOperands() == 1, "L should have one symbol");
	assert(K->identifier() == "Z", "Zassenhaus work only on the integers");

	bool stop = false;

	long s, i, j, l, p, A, B, C, gamma, gcd;

	AST *g, *x, *n, *b, *F, *D, *E, *q, *H, *Z, *G, *T, *S, *M, *u, *v, *gi;

	f = f->copy();

	// zero = integer(0);

	x = L->operand(0);

	n = degree(f, x);

	if(n->is(1))
	{
		delete n;
	
		return list({ f->copy() });
	}

	A = norm(f, L, K);
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

	q = integer(p);

	printf("p = %li\n", p);

	l = std::ceil(std::log(2*B + 1) / std::log(p));

	printf("l = %li\n", l);

	v = pp(f, L, K);

	u = sZp(v, x, p);

	Z = berlekampFactors(u, x, q);

	printf("berlekamp = %s\n", Z->toString().c_str());

	AST* T0 = mulPoly(Z->operand(0), Z->operand(1));
	AST* T1 = mulPoly(Z->operand(2), Z->operand(3));

	printf("FACTORING  = %s\n",  sZp(mulPoly(T0, T1), x, p)->toString().c_str());
	printf("FACTORING  = %s\n",  u->toString().c_str());
	printf("FACTORING  = %s\n",  v->toString().c_str());

	g = multifactorHenselLifting(f, Z, L, K, p, l);

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
		
			u = sZp(G, x, std::pow(p, l));
			v = sZp(H, x, std::pow(p, l));
		
			delete G;
			delete H;
			
			G = u;
			H = v;

			if(norm(G, L ,K) > std::pow(p, l) / 2)
			{
				delete Z;
				delete G;
				delete H;
			
				continue;
			}

			if(norm(H, L ,K) > std::pow(p, l) / 2)
			{
				delete Z;
				delete G;
				delete H;
			
				continue;
			}

			if(l1norm(G, L, K) * l1norm(H, L, K) <= B)
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

	return F;
}


}
