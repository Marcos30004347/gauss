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


AST* divGf(AST* a, AST* b, AST* x, long p, bool symetric = true)
{
	printf("OTRA = %s\n", divideGPE_sZp(a, b, x, p)->toString().c_str());
	printf("OTRA = %s\n", divideGPE_Zp(a, b, x, p)->toString().c_str());

	printf("\nnum = %s\n", a->toString().c_str());
	printf("den = %s\n", b->toString().c_str());

	AST* da = degree(a, x);
	AST* db = degree(b, x);

	if(da->value() < db->value())
	{
		delete da;
		delete db;

		return list({ integer(0), a->copy() });
	}

	long k, j, s, e, lb, d;

	a = a->copy();
	b = b->copy();

	AST *dq, *dr, *q, *r;
	AST *t1, *t2, *t3, *t4, *t5, *ex, *ca, *cb;

	dq = integer(da->value() - db->value());
	dr = integer(db->value() - 1);

	t1 = leadCoeff(b, x);	

	lb = modInverse_p(t1->value(), p);


	// if(symetric)
	// {
	// 	lb = sZp(lb, p);
	// }
	// else 
	// {
	// 	lb = Zp(lb, p);
	// }

	delete t1;

	for(k = da->value(); k >= 0; k--)
	{
		ex = integer(k);
		t1 = coeff(a, x, ex);
		
		s = std::max(0L, k - dq->value());
		e = std::min(dr->value(), k);
		
		delete ex;
		
		for(j = s; j <= e; j++)
		{
			ex = integer(j);
			cb = coeff(b, x, ex);

			delete ex;
		
			ex = integer(k - j + db->value());
			ca = coeff(a, x, ex);
			
			delete ex;
		
			t2 = mulPoly(cb, ca);
			
			delete ca;
			delete cb;

			t3 = subPoly(t1, t2);

			delete t1;
			delete t2;

			t1 = t3;
		}

		// if(symetric)
		// {
		// 	t2 = integer(sZp(t1->value(), p));
		// }
		// else 
		// {
		// 	t2 = integer(Zp(t1->value(), p));
		// }
		t2 = integer((t1->value() % p));

		delete t1;
		
		t1 = t2;
		
		if(t1->value() < 0)
		{
			printf("AAAAAAAAAAAAAAAAAAAA\n");
			t2 = integer(t1->value() + p);

			delete t1;

			t1 = t2;
		}

		if(k >= db->value())
		{
			t3 = integer(lb);
		
			t2 = mulPoly(t1, t3);
		
			delete t1;
			delete t3;
		
			t1 = t2;

			// if(symetric)
			// {
			// 	t2 = integer(sZp(t1->value(), p));
			// }
			// else 
			// {
			// 	t2 = integer(Zp(t1->value(), p));
			// }
			t2 = integer((t1->value() % p));
		
			delete t1;
		
			t1 = t2;
		}
		
		ex = integer(k);
		
		ca = coeff(a, x, ex);
		
		t3 = power(x->copy(), ex->copy());
		
		t4 = mulPoly(ca, t3);

		delete t3;

		t3 = subPoly(a, t4);

		delete a;

		a = algebraicExpand(t3);
		
		delete t3;

		t3 = power(x->copy(), ex->copy());
		
		t4 = mulPoly(t1, t3);
	
		delete t3;
	
		t3 = addPoly(a, t4);

		delete a;

		a = algebraicExpand(t3);
	}

	q = add({});
	r = add({});

	printf("==========> %s\n", a->toString().c_str());

	d = 0;
	for(k = da->value() - dq->value(); k <= da->value(); k++)
	{
		ex = integer(k);
	
		q->includeOperand(
			mul({
				coeff(a, x, ex),
				power(x->copy(), integer(d))
			})
		);

		delete ex;
		d = d + 1;
	}

	d = 0;
	for(k = 0; k <= dr->value(); k++)
	{
		ex = integer(k);
	
		r->includeOperand(
			mul({
				coeff(a, x, ex),
				power(x->copy(), integer(d))
			})
		);

		delete ex;
		d = d + 1;
	}

	printf("==========> %s\n", q->toString().c_str());
	printf("==========> %s\n", r->toString().c_str());

	delete da;
	delete db;

	delete dq;
	delete dr;
	
	if(symetric)
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
	printf("res = (%s * %s) + %s\n", b->toString().c_str(), t1->toString().c_str(), t2->toString().c_str());
	printf("res = %s\n\n", sZp(addPoly(t2, sZp(mulPoly(b, t1),x,p)),x,p)->toString().c_str() );

	return list({ t1, t2 });
}

AST* remGf(AST* a, AST* b, AST* x, long p, bool symetric = true)
{

	AST* d = divGf(a, b, x, p, symetric);

	AST* r = d->operand(1L);

	d->removeOperand(1L);



	delete d;

	return r;
}

AST* quoGf(AST* a, AST* b, AST* x, long p, bool symetric = true)
{
	AST* d = divGf(a, b, x, p, symetric);
	AST* r = d->operand(0L);

	d->removeOperand(0L);

	delete d;

	return r;
}


AST* gcdGf(AST* a, AST* b, AST* x, long p, bool symetric = true)
{
	AST* da = degree(a, x);
	AST* db = degree(b, x);
	
	if(b->is(0) || db == integer(-1))
	{
		return b->copy();
	}
	
	if(db->value() > da->value())
	{
		delete da;
		delete db;
		
		return gcdGf(b, a, x, p, symetric);
	}

	AST *t1;

	a = a->copy();
	b = b->copy();


	while(db->kind() != Kind::MinusInfinity && db->value() >= 0)
	{
		t1 = a;
		
		a = b;

		b = remGf(t1, b, x, p, symetric);

		printf("%s\n", b->toString().c_str());

		delete t1;

		delete db;

		db = degree(b, x);
	}

	delete b;

	delete da;
	delete db;

	return a;
}

AST* powGf(AST* f, AST* g, AST* x, long n, long p)
{
	AST *a, *b, *t, *u;

	if(n == 0) return integer(1);
	
	a = f->copy();

	b = integer(1);

	// AST* t1 = add({
	// 	mul({ integer(5), power(symbol("x"), integer(7)) }),
	// 	mul({ integer(8), power(symbol("x"), integer(5)) }),
	// 	power(symbol("x"), integer(2)),
	// 	symbol("x"),
	// 	integer(1)
	// });

	// AST* t2 = add({
	// 	power(symbol("x"), integer(2)),
	// 	integer(1)
	// });
	// printf("-> %s\n", divGf(t1, t2, x, p)->toString().c_str());
	// printf("-> %s\n", divideGPE_sZp(t1, t2, x, p)->toString().c_str());


	while(n > 1)
	{
		if(n % 2 == 0)
		{
			t = mulPoly(a, a);
		
			delete a;
		
			u = t->copy(); //sZp(t, x, p);
			
			delete t;

			// a = remainderGPE_sZp(u, g, x, p);
			a = remGf(u, g, x, p);
			
			// printf("-> %s\n", divGf(u, g, x, p)->toString().c_str());
			// printf("-> %s\n", divideGPE_sZp(u, g, x, p)->toString().c_str());

			delete u;
		
			n = n / 2;
		}
		else
		{
			t = mulPoly(a, b);
		
			delete b;
		
			u = sZp(t, x, p);
		
			delete t;
					
			// b = remainderGPE_sZp(u, g, x, p);
			b = remGf(u, g, x, p);

			// printf("-> %s\n", divGf(u, g, x, p)->toString().c_str());
			// printf("-> %s\n", divideGPE_sZp(u, g, x, p)->toString().c_str());

			delete u;

			t = mulPoly(a, a);
		
			delete a;
			
			u = sZp(t, x, p);

			delete t;

			// a = remainderGPE_sZp(u, g, x, p);
			a = remGf(u, g, x, p);
			// printf("-> %s\n", divGf(u, g, x, p)->toString().c_str());
			// printf("-> %s\n", divideGPE_sZp(u, g, x, p)->toString().c_str());

			delete u;
		
			n = (n - 1) / 2;
		}
	}


	t = mulPoly(a, b);

	u = sZp(t, x, p);

	delete a;
	delete b;

	delete t;

	a = remainderGPE_sZp(u, g, x, p);

	delete u;

	return a;
}


AST* distinctDegreeFactorization(AST* v, AST* L, AST* K, AST* q)
{
	long i, p;

	AST *x, *h, *f, *t, *G, *g, *k;
	
	x = L->operand(0);

	i = 0;
	h = x->copy();
	f = v->copy();

	G = list({});

	g = integer(0);

	// printf("%s\n", remainderGPE_sZp(power(symbol("x"), integer(9)), v, symbol("x"), 3)->toString().c_str());

	// return nullptr;
	p = q->value();

	do
	{
		i = i + 1;
	
		printf("aaa\n");

		printf("h = %s\n", h->toString().c_str());
		printf("v = %s\n", v->toString().c_str());
		printf("q = %s\n", q->toString().c_str());
	
		printf("*****\n");

		t = powGf(h, v, x, p, p);

		delete h;

		h = t;

		printf("h = %s\n", h->toString().c_str());

		t = subPoly(h, x);

		delete g;
		
		g = gcdGf(t, f, x, p, false);

		printf("REM = %s\n", divGf(f, t, x, p)->toString().c_str() );
		printf("REM = %s\n", sZp(mulPoly(t, quoGf(f, t, x, p)), x, p)->toString().c_str() );

		printf("g = %s\n", g->toString().c_str());

		delete t;

		t = quoGf(f, g, x, p);

		delete f;

		f = t;

		printf("f = %s\n", f->toString().c_str());

		G->includeOperand(g);
	
	} while (f->isNot(1));
	
	return G;
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
