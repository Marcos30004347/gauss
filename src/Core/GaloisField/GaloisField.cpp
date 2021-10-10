#include "GaloisField.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Simplification/Simplification.hpp"

#include <random>
#include <limits>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace galoisField {

long mod(long a, long b, bool symmetric) {
	
	long n = (b + (a % b)) % b;

	if(symmetric)
	{
		if(0 <= n && n <= b/2) 
		{
			return n;
		}

		return n - b;
	}

	return n;
}

long randomGf(long p, bool symmetric)
{
	std::random_device dev;
	
	std::mt19937 rng(dev());
	
	std::uniform_int_distribution<std::mt19937::result_type> dist(
		std::numeric_limits<long>::min(), 
		std::numeric_limits<long>::max() 
	);
	
	return mod(dist(rng), p, symmetric);
}

long inverseGf(long a, long b, bool symmetric) {
	long t, nt, r, nr, q, tmp;
	
	if (b < 0) b = -b;
	if (a < 0) a = b - (-a % b);
	
	t   = 0;  
	nt  = 1;  
	r 	= b;  
	nr  = a % b;
	
	while (nr != 0) 
	{
		q 	= r / nr;
		tmp = nt;  
		nt 	= t - q * nt;  
		t 	= tmp;
		tmp = nr;  
		nr 	= r - q * nr;  
		r 	= tmp;
	}

	if (r > 1)
	{
		printf("%li have no inverse mod %li\n", a, b);
		exit(1);
	}

	if (t < 0) t += b;

	return mod(t, b, symmetric);
}

long quoGf(long s, long t, long p, bool symmetric) {
	return mod((s * inverseGf(t,p,symmetric)), p, symmetric);
}

AST* gf(AST* u, AST* x, int s, bool symmetric) 
{
	AST* k = algebraicExpand(u);

	AST* p = new AST(Kind::Addition);

	AST* d = degree(k, x);

	for(int i = 0; i <= d->value(); i++) {
		AST* n = integer(i);

		AST* c = coeff(k, x, n);

		p->includeOperand(
			mul({
				integer(mod(c->value(), s, symmetric)),
				power(x->copy(), n)
			})
		);

		delete c;
	}
	
	AST* r = algebraicExpand(p);

	delete d;
	delete k;
	delete p;

	return r;
}

AST* divPolyGf(AST* a, AST* b, AST* x, long p, bool symmetric)
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

	lb = inverseGf(t1->value(), p, symmetric);

	delete t1;

	for(k = da->value(); k >= 0; k--)
	{
		t1 = A[k]->copy();
	
		s = std::max(0L, k - dq->value());
		e = std::min(dr->value(), k);
		
		for(j = s; j <= e; j++)
		{
			t2 = mulPoly(B[j], A[k - j + db->value()]);
			
			t3 = subPoly(t1, t2);

			delete t1;
			delete t2;

			t1 = t3;
		}

		t3 = reduceAST(t1);
		
		delete t1;
		
		t1 = t3;	

		t2 = integer(mod(t1->value(), p, symmetric));

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
		
			t2 = integer(mod(t1->value(), p, symmetric));

			delete t1;
		
			t1 = t2;
		}

		delete A[k];
		
		A[k] = t1;
	}

	q = add({integer(0)});
	r = add({integer(0)});

	d = 0;
	for(k = da->value() - dq->value(); k <= da->value(); k++)
	{
		q->includeOperand(
			mul({
				A[k]->copy(),
				power(x->copy(), integer(d))
			})
		);

		d = d + 1;
	}

	d = 0;
	for(k = 0; k <= dr->value(); k++)
	{
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
	
	t1 = gf(q, x, p, symmetric);
	t2 = gf(r, x, p, symmetric);

	delete q;
	delete r;

	return list({ t1, t2 });
}

AST* remPolyGf(AST* a, AST* b, AST* x, long p, bool symmetric)
{
	AST* d = divPolyGf(a, b, x, p, symmetric);

	AST* r = d->operand(1L);

	d->removeOperand(1L);

	delete d;

	return r;
}

AST* quoPolyGf(AST* a, AST* b, AST* x, long p, bool symmetric)
{
	AST* d = divPolyGf(a, b, x, p, symmetric);

	AST* q = d->operand(0L);

	d->removeOperand(0L);

	delete d;

	return q;
}

AST* monicPolyGf(AST* f, AST* x, long p, bool symmetric)
{
	if(f->is(0))
	{
		return integer(0);
	}

	AST* lc = leadCoeff(f, x);

	AST* F = quoPolyGf(f, lc, x, p, symmetric);

	return list({ lc, F });
}

AST* gcdPolyGf(AST* a, AST* b, AST* x, long p, bool symmetric)
{
	AST* da = degree(a, x);
	AST* db = degree(b, x);
	
	if(da->kind() == Kind::MinusInfinity || db->value() > da->value())
	{
		delete da;
		delete db;
		
		return gcdPolyGf(b, a, x, p, symmetric);
	}

	AST *t1;

	a = a->copy();
	b = b->copy();

	while(b->isNot(0) && db->kind() != Kind::MinusInfinity && db->value() >= 0)
	{
		t1 = a;
		
		a = b;

		b = remPolyGf(t1, b, x, p, symmetric);

		delete t1;

		delete db;

		db = degree(b, x);
	}

	delete b;

	delete da;
	delete db;

	b = monicPolyGf(a, x, p, symmetric);

	delete a;

	a = b->operand(1L);
	
	b->removeOperand(1L);
	
	delete b;
	
	return a;
}

AST* addPolyGf(AST* f, AST* g, AST* x, long p, bool symmetric)
{
	AST *t, *u;

	u = addPoly(f, g);

	t = gf(u, x, p, symmetric);

	delete u;

	return t;
}

AST* subPolyGf(AST* f, AST* g, AST* x, long p, bool symmetric)
{
	AST *t, *u;

	u = subPoly(f, g);

	t = gf(u, x, p, symmetric);

	delete u;

	return t;
}

AST* mulPolyGf(AST* f, AST* g, AST* x, long p, bool symmetric)
{
	AST *t, *u;

	u = mulPoly(f, g);

	t = gf(u, x, p, symmetric);

	delete u;

	return t;
}

AST* powModPolyGf(AST* f, AST* g, AST* x, long n, long p, bool symmetric)
{
	AST *a, *b, *t;

	if(n == 0) return integer(1);
	
	a = f->copy();

	b = integer(1);

	while(n > 1)
	{
		if(n % 2 == 0)
		{
			t = mulPolyGf(a, a, x, p, symmetric);

			a = remPolyGf(t, g, x, p, symmetric);
	
			delete t;
		
			n = n / 2;
		}
		else
		{
			t = mulPolyGf(a, b, x, p, symmetric);
		
			b = remPolyGf(t, g, x, p, symmetric);

			delete t;

			t = mulPolyGf(a, a, x, p, symmetric);

			a = remPolyGf(t, g, x, p, symmetric);

			delete t;
		
			n = (n - 1) / 2;
		}
	}

	t = mulPolyGf(a, b, x, p, symmetric);

	delete a;
	delete b;

	a = remPolyGf(t, g, x, p, symmetric);

	delete t;

	return a;
}

AST* randPolyGf(long d, AST* x, long p, bool symmetric)
{
	long k = 0;

	if(d == 0)
	{
		return integer(randomGf(p, symmetric));
	} 

	if(d == 1)
	{
		k = randomGf(p, symmetric);

		if(k == 0) return x->copy();
	
		return add({x->copy(), integer(k)});
	}

	AST* r = add({ power(symbol("x"), integer(d)) });

	for(long i = d - 1; i >= 2; i--)
	{
		k = randomGf(p, symmetric);
		
		if(k != 0)
		{
			r->includeOperand(mul({
				integer(k),
				power(x->copy(), integer(i))
			}));
		}
	}

	k = randomGf(p, symmetric);
	
	if(k != 0)
	{
		r->includeOperand(mul({
			integer(k),
			x->copy()
		}));
	}

	k = randomGf(p, symmetric);

	if(k != 0)
	{
		r->includeOperand(integer(k));
	}

	return r;
}

AST* extendedEuclidGf(AST* f, AST* g, AST* x, long p, bool sym)
{
	if(f->is(0) || g->is(0))
	{
		return list({integer(1), integer(0), integer(0)});
	}

	AST *t, *s, *p0, *i, *lc, *k1, *r0, *p1, *r1, *t0, *t1, *t2, *t3, *s0, *s1, *Q, *R, *T;

	t1 = monicPolyGf(f, x, p, sym); 
	t2 = monicPolyGf(g, x, p, sym); 

	p0 = t1->operand(0);
	r0 = t1->operand(1);

	p1 = t2->operand(0);
	r1 = t2->operand(1);

	t1->removeOperand(0L);
	t1->removeOperand(0L);
	t2->removeOperand(0L);
	t2->removeOperand(0L);
	
	delete t1;
	delete t2;

	if(f->is(0))
	{
		t1 = integer(0);
	
		t2 = integer(inverseGf(p1->value(), p, sym));
	
		t3 = r1;

		delete p0;
		delete p1;
		delete r0;
	
		return list({t1, t2, t3});
	}

	if(g->is(0))
	{
		t1 = integer(inverseGf(p0->value(), p, sym));
		t2 = integer(0);
		t3 = r0;

		delete p0;
		delete p1;
		delete r1;
	
		return list({t1, t2, t3});
	}

	s0 = integer(inverseGf(p0->value(), p, sym));
	s1 = integer(0);
	
	t0 = integer(0);
	t1 = integer(inverseGf(p1->value(), p, sym));

	while(true)
	{
		T = divPolyGf(r0, r1, x, p, sym);
		
		Q = T->operand(0L);
		R = T->operand(1L);
		
		T->removeOperand(0L);
		T->removeOperand(0L);
		
		delete T;
	
		if(R->is(0))
		{
			delete Q;
			delete R;
			
			break;
		}

		T = monicPolyGf(R, x, p, sym);
		
		delete r0;
	
		r0 = r1;
	
		lc = T->operand(0L);
		r1 = T->operand(1L);
		
		T->removeOperand(0L);
		T->removeOperand(0L);
		
		delete T;
	
		i = integer(inverseGf(lc->value(), p, sym));
	
		k1 = mulPolyGf(s1, Q, x, p, sym);
		s  = subPolyGf(s0, k1, x, p, sym);

		delete k1;
	
		k1 = mulPolyGf(t1, Q, x, p, sym);
		t  = subPolyGf(t0, k1, x, p, sym);

		delete k1;

		delete s0;
		delete t0;
	
		s0 = s1;
		t0 = t1;

		s1 = mulPolyGf(s, i, x, p, sym);
		t1 = mulPolyGf(t, i, x, p, sym);

		delete Q;
		delete R;

		delete lc;
		delete i;
	}

	delete s0;
	delete t0;
	delete r0;

	return list({ r1, s1, t1 });
}



}
