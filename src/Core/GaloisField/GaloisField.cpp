/**
 * @file GaloisField.cpp
 * @author Marcos Vincius Moreira Santos (marcos30004347@gmail.com)
 * @brief This file implement some of the Finite field methods used in the library
 * @version 0.1
 * @date 2021-10-11
 * 
 * @ref Michael Monagan - In-place Arithmetic for Polynomials over Zn
 * 
 * @copyright Copyright (c) 2021
 */

#include "GaloisField.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Debug/Assert.hpp"

#include "Core/Simplification/Simplification.hpp"

#include <random>
#include <limits>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace galoisField {

long long mod(long long a, long long b, bool symmetric) {
	
	long long n = (b + (a % b)) % b;

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

long long randomGf(long long p, bool symmetric)
{
	std::random_device dev;
	
	std::mt19937 rng(dev());
	
	std::uniform_int_distribution<std::mt19937::result_type> dist(
		std::numeric_limits<long long>::min(), 
		std::numeric_limits<long long>::max() 
	);
	
	return mod(dist(rng), p, symmetric);
}

long long inverseGf(long long a, long long b, bool symmetric) {
	long long t, nt, r, nr, q, tmp;
	
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

long long quoGf(long long s, long long t, long long p, bool symmetric) {
	return mod((s * inverseGf(t,p,symmetric)), p, symmetric);
}

AST* gf(AST* u, AST* x, long long s, bool symmetric) 
{
	AST* k = algebraicExpand(u);

	if(k->kind() == Kind::Fail || k->kind() == Kind::Undefined)
	{
		return k;
	}

	if(k->kind() == Kind::MinusInfinity || k->kind() == Kind::Infinity)
	{
		delete k;
	
		return undefined();
	}

	if(k->kind() == Kind::Integer)
	{
		long p = k->value();
	
		delete k;
	
		return integer(mod(p, s, symmetric));
	}

	if(k->kind() == Kind::Symbol)
	{
		return k;
	}

	if(k->kind() == Kind::Fraction)
	{
		assert(
			k->operand(0)->kind() == Kind::Integer, 
			"numerator of a fraction needs to be a integer"
		);
		
		assert(
			k->operand(1)->kind() == Kind::Integer, 
			"denominator of a fraction needs to be a integer"
		);
		
		long n = k->operand(0)->value();
		long d = k->operand(1)->value();
	
		delete k;
	
		return integer(mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric));
	}

	if(k->kind() == Kind::Derivative)
	{
		AST* p =  new AST(Kind::Derivative,{
			gf(k->operand(0), x, s, symmetric),
			k->operand(1)->copy()
		});
	
		delete k;
	
		return p;
	}

	if(k->kind() == Kind::Integral)
	{
		AST* p = new AST(Kind::Integral,{
			gf(k->operand(0), x, s, symmetric),
			k->operand(1)->copy()
		});

		delete k;
	
		return p;
	}

	if(k->kind() == Kind::Factorial)
	{
		if(k->operand(0)->kind() == Kind::Integer)
		{
			AST* f = reduceAST(k);
	
			AST* p = gf(f, x, s, symmetric);
		
			delete f;
			
			delete k;
		
			return p;
		}

		return k;
	}

	if(k->kind() == Kind::Division)
	{
		AST* p = div(gf(k->operand(0), x, s, symmetric), gf(k->operand(1), x, s, symmetric));
		AST* t = reduceAST(p);
		AST* r = gf(t, x, s, symmetric);
	
		delete p;
		delete t;
		delete k;
	
		return r;
	}

	if(k->kind() == Kind::Power)
	{
		AST* p = power(gf(k->operand(0), x, s, symmetric), k->operand(1)->copy());
	
		delete k;
	
		return p;
	}

	if(k->kind() == Kind::FunctionCall)
	{
		return k;
	}

	if(k->kind() == Kind::Multiplication)
	{
		AST* p = new AST(Kind::Multiplication);
	
		for(long i = 0; i < k->numberOfOperands(); i++) {
			p->includeOperand(gf(k->operand(i), x, s, symmetric));
		}
	
		delete k;
	
		return p;
	}

	if(k->kind() == Kind::Addition || k->kind() == Kind::Subtraction)
	{
		AST* p = new AST(k->kind());

		AST* d = degree(k, x);

		for(long i = 0; i <= d->value(); i++) {
			AST* n = integer(i);

			AST* c = coeff(k, x, n);

			p->includeOperand(
				mul({
					gf(c, x, s, symmetric),
					power(x->copy(), n)
				})
			);

			delete c;
		}
		
		AST* r = reduceAST(p);

		delete d;
		delete k;
		delete p;

		return r;
	}

	return k;

	// AST* k = algebraicExpand(u);

	// AST* p = new AST(Kind::Addition);

	// AST* d = degree(k, x);

	// for(int i = 0; i <= d->value(); i++) {
	// 	AST* n = integer(i);

	// 	AST* c = coeff(k, x, n);

	// 	p->includeOperand(
	// 		mul({
	// 			integer(mod(c->value(), s, symmetric)),
	// 			power(x->copy(), n)
	// 		})
	// 	);

	// 	delete c;
	// }
	
	// AST* r = algebraicExpand(p);

	// delete d;
	// delete k;
	// delete p;

	// return r;
}



AST* gf(AST* u, long long s, bool symmetric) 
{
	AST* k = algebraicExpand(u);
	if(k->kind() == Kind::Fail || k->kind() == Kind::Undefined)
	{
		return k;
	}

	if(k->kind() == Kind::MinusInfinity || k->kind() == Kind::Infinity)
	{
		delete k;
	
		return undefined();
	}

	if(k->kind() == Kind::Integer)
	{

		long p = k->value();
	
		delete k;
	
		return integer(mod(p, s, symmetric));
	}

	if(k->kind() == Kind::Symbol)
	{
		return k;
	}

	if(k->kind() == Kind::Fraction)
	{
		assert(
			k->operand(0)->kind() == Kind::Integer, 
			"numerator of a fraction needs to be a integer"
		);
		
		assert(
			k->operand(1)->kind() == Kind::Integer, 
			"denominator of a fraction needs to be a integer"
		);
		
		long n = k->operand(0)->value();
		long d = k->operand(1)->value();
	
		delete k;
	
		return integer(mod(mod(n, s, symmetric) * inverseGf(d, s, symmetric), s, symmetric));
	}

	if(k->kind() == Kind::Derivative)
	{
		AST* p =  new AST(Kind::Derivative,{
			gf(k->operand(0), s, symmetric),
			k->operand(1)->copy()
		});
	
		delete k;
	
		return p;
	}

	if(k->kind() == Kind::Integral)
	{
		AST* p = new AST(Kind::Integral,{
			gf(k->operand(0), s, symmetric),
			k->operand(1)->copy()
		});

		delete k;
	
		return p;
	}

	if(k->kind() == Kind::Factorial)
	{
		if(k->operand(0)->kind() == Kind::Integer)
		{
			AST* f = reduceAST(k);
	
			AST* p = gf(f, s, symmetric);
		
			delete f;
			
			delete k;
		
			return p;
		}

		return k;
	}

	if(k->kind() == Kind::Division)
	{
		AST* p = div(gf(k->operand(0), s, symmetric), gf(k->operand(1), s, symmetric));
		AST* t = reduceAST(p);
		AST* r = gf(t, s, symmetric);
	
		delete p;
		delete t;
		delete k;
	
		return r;
	}

	if(k->kind() == Kind::Power)
	{
		AST* p = power(gf(k->operand(0), s, symmetric), k->operand(1)->copy());
	
		delete k;
	
		return p;
	}

	if(k->kind() == Kind::FunctionCall)
	{
		return k;
	}

	if(k->kind() == Kind::Multiplication)
	{
		AST* p = new AST(Kind::Multiplication);
	
		for(long i = 0; i < k->numberOfOperands(); i++) {
			p->includeOperand(gf(k->operand(i), s, symmetric));
		}

		delete k;
	
		return p;
	}

	if(k->kind() == Kind::Addition || k->kind() == Kind::Subtraction)
	{
		AST* p = new AST(k->kind());

		for(long i = 0; i < k->numberOfOperands(); i++) {
			p->includeOperand(
				mul({
					gf(k->operand(i), s, symmetric),
				})
			);
		}
		
		AST* r = reduceAST(p);

		delete k;
		delete p;

		return r;
	}

	return k;
}



AST* divPolyGf(AST* a, AST* b, AST* x, long long p, bool symmetric)
{
	AST* da = degree(a, x);
	AST* db = degree(b, x);

	if(da->value() < db->value())
	{
		delete da;
		delete db;

		return list({ integer(0), a->copy() });
	}

	long long k, j, s, e, lb, d;

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
	
		s = std::max(0LL, k - dq->value());
		e = std::min(dr->value(), k);

		for(j = s; j <= e; j++)
		{
			t2 = mulPoly(B[j], A[k - j + db->value()]);
			
			t3 = subPoly(t1, t2);

			delete t2;

			delete t1;

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
		
			delete t2;
		
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

	for(k = 0; k <= da->value(); k++)
	{
		delete A[k];
	}

	delete[] A;

	for(k = 0; k <= db->value(); k++)
	{
		delete B[k];
	}

	delete[] B;

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

AST* remPolyGf(AST* a, AST* b, AST* x, long long p, bool symmetric)
{
	AST* d = divPolyGf(a, b, x, p, symmetric);

	AST* r = d->operand(1L);

	d->removeOperand(1L);

	delete d;

	return r;
}

AST* quoPolyGf(AST* a, AST* b, AST* x, long long p, bool symmetric)
{
	AST* d = divPolyGf(a, b, x, p, symmetric);

	AST* q = d->operand(0L);

	d->removeOperand(0L);

	delete d;

	return q;
}

AST* monicPolyGf(AST* f, AST* x, long long p, bool symmetric)
{
	if(f->is(0))
	{
		return integer(0);
	}

	AST* lc = leadCoeff(f, x);

	AST* F = quoPolyGf(f, lc, x, p, symmetric);

	return list({ lc, F });
}

AST* gcdPolyGf(AST* a, AST* b, AST* x, long long p, bool symmetric)
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


	delete da;
	delete db;

	delete b;

	b = monicPolyGf(a, x, p, symmetric);

	delete a;

	a = b->operand(1L);
	
	b->removeOperand(1L);
	
	delete b;
	
	return a;
}

AST* addPolyGf(AST* f, AST* g, AST* x, long long p, bool symmetric)
{
	AST *t, *u;

	u = addPoly(f, g);

	t = gf(u, x, p, symmetric);

	delete u;

	return t;
}

AST* subPolyGf(AST* f, AST* g, AST* x, long long p, bool symmetric)
{
	AST *t, *u;

	u = subPoly(f, g);

	t = gf(u, x, p, symmetric);

	delete u;

	return t;
}

AST* mulPolyGf(AST* f, AST* g, AST* x, long long p, bool symmetric)
{
	AST *t, *u;

	u = mulPoly(f, g);

	t = gf(u, x, p, symmetric);

	delete u;

	return t;
}

AST* powModPolyGf(AST* f, AST* g, AST* x, long long n, long long p, bool symmetric)
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

			delete a;

			a = remPolyGf(t, g, x, p, symmetric);
	
			delete t;
		
			n = n / 2;
		}
		else
		{
			t = mulPolyGf(a, b, x, p, symmetric);

			delete b;
			
			b = remPolyGf(t, g, x, p, symmetric);

			delete t;

			t = mulPolyGf(a, a, x, p, symmetric);

			delete a;
	
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

AST* randPolyGf(long long d, AST* x, long long p, bool symmetric)
{
	long long k = 0;

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

AST* extendedEuclidGf(AST* f, AST* g, AST* x, long long p, bool sym)
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

		delete s;
		delete t;
	}

	delete s0;
	delete t0;
	delete r0;
	delete p0;
	delete p1;

	return list({ r1, s1, t1 });
}



}
