#include "Zp.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace simplification;

namespace polynomial {

int pow(int x, unsigned int y, unsigned int m) {
	if (y == 0)
		return 1;
	int p = (pow(x, y / 2, m) % m);
	p = ((p * p) % m);

	return (y % 2 == 0) ? p : (x * p) % m;
}
 
int gcd(int a, int b) {
	if (a == 0)
		return b;
	return gcd(b % a, a);
}

int extended_euclidean(int a, int b, int& x, int& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
  
	  int x1, y1;
	
    int d = extended_euclidean(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return d;
}

long modInverse(long a, long b) {
	int t, nt, r, nr, q, tmp;
	
	if (b < 0) b = -b;
	if (a < 0) a = b - (-a % b);
	
	t = 0;  nt = 1;  r = b;  nr = a % b;
	
	while (nr != 0) 
	{
		q = r/nr;
		tmp = nt;  
		nt = t - q*nt;  
		t = tmp;
		tmp = nr;  
		nr = r - q*nr;  
		r = tmp;
	}

	if (r > 1)
	{
		printf("%li have no inverse mod %li\n", a, b);
		exit(1);
	}

	if (t < 0) t += b;

	return t;

}

 
long gcd_sZp(long a, long b) {
	if (a == 0)
	{
		return b;
	}

	return gcd(sZp(b, a), a);
}

long modInverse_sZp(long a, long p) {
	long g = gcd_sZp(a, p);

	if (g != 1 && g!= -1) 
	{
		printf("Inverse of %li in sZ%li doesn't exist!\n", a, p);
		abort();
	} else 
	{
		return pow(a, p - 2, p);
	}
}

long sZp(long b, long m) 
{
	b = mod(b, m);

	if(0 <= b && b <= m/2) 
	{
		return b;
	}

	return b - m;
}

long Zp(long b, long m) 
{
	return mod(b, m);
}

long division_Zp(long s, long t, long p) {
	return mod((s * modInverse(t,p)), p);
}

long division_sZp(long s, long t, long p) {
	return sZp(mod(s * modInverse(t, p), p), p);
}

long mul_Zp(long s, long t, long p) {
	return mod((s * t), p);
}

long mul_sZp(long s, long t, long p) {
	return sZp(mod(s * t, p), p);
}

// Tnn in Symbolic Algebra, this function
// gives the non negative of u(x) 
// defined in Z[x] projected on Zp[x]
AST* Zp(AST* u, AST* x, int s) {
	AST* u_ = algebraicExpand(u);

	AST* Tnn_u = new AST(Kind::Addition);

	AST* d = degree(u_, x);

	for(int i=0; i<=d->value(); i++) {
		AST* d = integer(i);

		AST* c = coeff(u_, x, d);
		
		Tnn_u->includeOperand(mul({ integer(mod(c->value(), s)), power(x->copy(), d) }));
		
		delete c;
	}


	AST* r = expandAST(Tnn_u);

	delete d;
	delete u_;
	delete Tnn_u;

	return r;
}

// Ts in Symbolic Algebra, this function
// gives the symetric projection of u(x) 
// defined in Z[x] projected on Zp[x]
AST* sZp(AST* u, AST* x, int s) {
	// TODO: assume that u(x) is already expanded
	AST* u_ = algebraicExpand(u);

	AST* Tnn_u = new AST(Kind::Addition);

	AST* d = degree(u_, x);

	for(int i=0; i <= d->value(); i++) {

		AST* e = integer(i);
		AST* c = coeff(u_, x, e);
		// AST* c  = expandAST(c_);

		if(i > 0)
		{
			Tnn_u->includeOperand(
				mul({
					integer(sZp(c->value(), s)),
					power(x->copy(), e)
				})
			);
		}
		else
		{
			Tnn_u->includeOperand(integer(sZp(c->value(), s)));
		}

		// delete c_;
		delete c;
	}

	AST* r = reduceAST(Tnn_u);
	
	delete d;
	delete u_;
	delete Tnn_u;

	return r;
}

AST* divideGPE_Zp(AST* u, AST* v, AST* x, int p) 
{
	AST *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *t9, *t10;

	AST* q = integer(0);
	AST* r = u->copy();

	AST* m = degree(r, x);
	AST* n = degree(v, x);

	AST* lcv = leadCoeff(v, x);

	while(
		m->kind() != Kind::MinusInfinity &&
		(m->kind() == Kind::Integer && n->kind() == Kind::Integer &&
		m->value() >= n->value())
	) {
	
		AST* lcr = leadCoeff(r, x);
	
		AST* s = integer(division_Zp(lcr->value(), lcv->value(), p));

		t1 = power(x->copy(), sub({m->copy(), n->copy()}));
		t2 = mulPoly(s, t1);
		t3 = addPoly(q, t2);

		delete q;
	
		q = reduceAST(t3);
	
		delete t1;
		delete t2;
		delete t3;

		t1 = power(x->copy(), m->copy());
		t2 = mulPoly(lcr, t1);
		t3 = subPoly(r, t2);

		t4 = power(x->copy(), n->copy());
		t5 = mulPoly(lcv, t4);
		t6 = subPoly(v, t5);

		t7 = mulPoly(t6, s);
		t8 = power(x->copy(), sub({m->copy(), n->copy()}));
		t9 = mulPoly(t7, t8);
		t10 = subPoly(t3, t9);
	
		delete r;
	
		r = reduceAST(t10);
	
		delete t1;
		delete t2;
		delete t3;
		delete t4;
		delete t5;
		delete t6;
		delete t7;
		delete t8;
		delete t9;
		delete t10;
	
		delete m;
		delete lcr;
		delete s;
	
		m = degree(r, x);
	}

	AST* res = list({ Zp(q, x, p), Zp(r, x, p) });
	
	delete q;
	delete r;

	delete m;
	delete n;
	delete lcv;

	return res;
}


AST* divideGPE_sZp(AST* u, AST* v, AST* x, int p) 
{
	AST *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *t9, *t10;

	AST* q = integer(0);
	AST* r = u->copy();

	AST* m = degree(r, x);
	AST* n = degree(v, x);

	AST* lcv = leadCoeff(v, x);

	while(
		m->kind() != Kind::MinusInfinity &&
		(m->kind() == Kind::Integer && n->kind() == Kind::Integer &&
		m->value() >= n->value())
	) {
	
		AST* lcr = leadCoeff(r, x);
	
		AST* s = integer(division_sZp(lcr->value(), lcv->value(), p));

		t1 = power(x->copy(), sub({m->copy(), n->copy()}));
		t2 = mulPoly(s, t1);
		t3 = addPoly(q, t2);

		delete q;
	
		q = reduceAST(t3);

		delete t1;
		delete t2;
		delete t3;

		t1 = power(x->copy(), m->copy());
		t2 = mulPoly(lcr, t1);
		t3 = subPoly(r, t2);

		t4 = power(x->copy(), n->copy());
		t5 = mulPoly(lcv, t4);
		t6 = subPoly(v, t5);

		t7 = mulPoly(t6, s);
		t8 = power(x->copy(), sub({m->copy(), n->copy()}));
		t9 = mulPoly(t7, t8);
		t10 = subPoly(t3, t9);
	
		delete r;
	
		r = reduceAST(t10);
	
		delete m;
		delete lcr;
		delete s;
	
		m = degree(r, x);
	}

	AST* res = list({ sZp(q,x,p), sZp(r,x,p) });

	delete q;
	delete r;
	delete m;
	delete n;
	delete lcv;

	return res;
}

AST* remainderGPE_Zp(AST* u, AST* v, AST* x, int p) {
	AST* k = divideGPE_Zp(u,v,x,p);
	AST* r = k->operand(1)->copy();
	delete k;
	return r;
}

AST* quotientGPE_Zp(AST* u, AST* v, AST* x, int p) {
	AST* k = divideGPE_Zp(u,v,x,p);
	AST* q = k->operand(0)->copy();
	delete k;
	return q;
}

AST* remainderGPE_sZp(AST* u, AST* v, AST* x, int p) {
	AST* k = divideGPE_sZp(u,v,x,p);

	AST* r = k->operand(1)->copy();

	delete k;

	return r;
}

AST* quotientGPE_sZp(AST* u, AST* v, AST* x, int p) {
	AST* k = divideGPE_sZp(u,v,x,p);
	AST* q = k->operand(0)->copy();
	delete k;
	return q;
}


AST* gcdGPE_Zp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return integer(0);
	}

	AST* U = u->copy();
	AST* V = v->copy();

	while (V->kind() != Kind::Integer || V->value() != 0) {
		AST* R = remainderGPE_Zp(U, V, x, p);
		delete U;
		U = V->copy();
		delete V;
		V = R->copy();
		delete R;
	}

	AST* lco = leadCoeff(U, x);
	
	AST* e = mul({
		integer(division_Zp(1, lco->value(), p)),
		U->copy()
	});
	
	delete lco;

	AST* res = Zp(e, x, p);

	delete U;
	delete V;
	delete e;

	return res;
}


AST* gcdGPE_sZp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return integer(0);
	}

	AST* U = u->copy();
	AST* V = v->copy();

	while (V->kind() != Kind::Integer ||(V->kind() == Kind::Integer && V->value() != 0)) {
		AST* R = remainderGPE_sZp(U, V, x, p);

		delete U;
		U = V->copy();
		
		delete V;
		V = R->copy();
		
		delete R;
	}

	AST* lco = leadCoeff(U,x);
	
	AST* e = mul({ integer(division_sZp(1, lco->value(), p)), U->copy() });

	AST* res = sZp(e, x, p);

	delete U;
	delete V;
	delete e;
	delete lco;

	return res;
}


AST* extendedEuclideanAlgGPE_Zp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return list({ integer(0), integer(0), integer(0) });
	}

	AST* U 		= u->copy();
	AST* V 		= v->copy();
	AST* App 	= integer(1);
	AST* Ap 	= integer(0);
	AST* Bpp 	= integer(0);
	AST* Bp 	= integer(1);

	while (
		V->kind() != Kind::Integer ||
		(V->kind() == Kind::Integer && V->value() != 0)
	) {
		AST* d = divideGPE_Zp(U,V,x,p);
	
		AST* q = d->operand(0);
		AST* r = d->operand(1);

		AST* A_ = sub({ App->copy(), mul({q->copy(), Ap->copy()}) });
		AST* B_ = sub({ Bpp->copy(), mul({q->copy(), Bp->copy()}) });

		AST* A = Zp(A_,x, p);
		AST* B = Zp(B_,x, p);
	
		delete A_;
		delete B_;
	
		delete App;
		App = Ap->copy();

		delete Ap;
		Ap 	= A->copy();

		delete Bpp;
		Bpp = Bp->copy();

		delete Bp;
		Bp 	= B->copy();

		delete U;
		U = V->copy();

		delete V;
		V = r->copy();

		delete A;
		delete B;
	
		delete d;
	}

	AST* c = leadCoeff(U, x);

	AST* App__ = mul({ App->copy(), integer(modInverse(c->value(), p)) });
	AST* App_ = Zp(App__, x, p);
	delete App;
	delete App__;
	App = App_;

	AST* Bpp__ = mul({ Bpp->copy(), integer(modInverse(c->value(), p)) });
	AST* Bpp_ = Zp(Bpp__, x, p);
	delete Bpp;
	delete Bpp__;
	Bpp = Bpp_;
	
	AST* U__ = mul({U->copy(), integer(modInverse(c->value(), p))});
	AST* U_ = Zp(U__, x, p);
	delete U;
	delete U__;
	U = U_;

	delete Ap;
	delete Bp;
	delete V;
	delete c;
	
	return list({ U, App, Bpp });
}


ast::AST* extendedEuclideanAlgGPE_sZp(AST* u, AST* v, AST* x, int p) {


	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) return list({ integer(0), integer(0), integer(0) });

	AST* U 		= u->copy();
	AST* V 		= v->copy();
	AST* App 	= integer(1);
	AST* Ap 	= integer(0);
	AST* Bpp 	= integer(0);
	AST* Bp 	= integer(1);

	while (V->kind() != Kind::Integer || V->value() != 0) {

		AST* d = divideGPE_sZp(U,V,x,p);
		
		AST* q = d->operand(0);
		AST* r = d->operand(1);

		AST* A_ = sub({ App->copy(), mul({q->copy(), Ap->copy()}) });
		AST* B_ = sub({ Bpp->copy(), mul({q->copy(), Bp->copy()}) });
		
		AST* A = sZp(A_,x, p);
		AST* B = sZp(B_,x, p);
	
		delete A_;
		delete B_;
	
		delete App;
		App = Ap->copy();

		delete Ap;
		Ap 	= A->copy();

		delete Bpp;
		Bpp = Bp->copy();

		delete Bp;
		Bp 	= B->copy();

		delete A;
		delete B;

		delete U;
		U = V->copy();

		delete V;
		V = r->copy();

		delete d;
	}

	AST* c = leadCoeff(U, x);

	AST* App__ = mul({ App->copy(), integer(modInverse(mod(c->value(),p), p)) });
	delete App;
	App = sZp(App__, x, p);
	delete App__;

	AST* Bpp__ = mul({ Bpp->copy(), integer(modInverse(mod(c->value(),p), p)) });
	delete Bpp;
	Bpp = sZp(Bpp__, x, p);
	delete Bpp__;
	
	AST* U__ = mul({U->copy(), integer(modInverse(mod(c->value(),p), p)) });
	delete U;
	U = sZp(U__, x, p);
	delete U__;

	delete Ap;
	delete Bp;
	delete V;
	delete c;
	
	return list({ U, App, Bpp });
}


bool isRowOfZeros(AST* M, int n, int j)
{
	for(int i=0; i<n; i++)
	{
		if(M->operand(j)->operand(i)->isNot(0))
			return false;
	}

	return true;
}

AST* nullSpace_sZp(AST* M, signed long q)
{
	assert(
		M->numberOfOperands() >= 1,
		"The matrix should have at least one row"
	);

	assert(
		M->numberOfOperands() == M->operand(0)->numberOfOperands(),
		"The matrix should be square"
	);

	M = M->copy();


	int k, i, j, n = M->numberOfOperands();

	for(k=0; k < n; k++)
	{
		for(i = k; i < n && M->operand(k)->operand(i)->is(0); i++)
		{}

		if(i < n)
		{
			// Normalized column i 
			
			signed long d = M->operand(k)->operand(i)->value();
		
			for(j = 0; j < n; j++)
			{
				signed long Mji = M->operand(j)->operand(i)->value();
				signed long p = division_sZp(Mji, d, q);
			
				M->operand(j)->deleteOperand(i);
				M->operand(j)->includeOperand(integer(p), i);
			}

			// Switch column i with column k
			for(j = 0; j < n; j++)
			{
				AST* Mji = M->operand(j)->operand(i)->copy();
				AST* Mjk = M->operand(j)->operand(k)->copy();

				M->operand(j)->deleteOperand(i);
				M->operand(j)->includeOperand(Mjk, i);
			
				M->operand(j)->deleteOperand(k);
				M->operand(j)->includeOperand(Mji, k);
			}

			// Eliminate rest of row k via column operations
			for(i = 0; i < n; i++)
			{
				if(i != k)
				{
					signed long Mki = M->operand(k)->operand(i)->value();

					for(j = 0; j < n; j++)
					{
						signed long col_i = M->operand(j)->operand(i)->value();
						signed long col_k = M->operand(j)->operand(k)->value();

						signed long tt = sZp(col_i - col_k * Mki, q);

						M->operand(j)->deleteOperand(i);
						M->operand(j)->includeOperand(integer(tt) ,i);
					}
				}
				
			}
		}
	}

	// Convert M to M-I
	for(i = 0; i < n; i++)
	{
		signed long tt = M->operand(i)->operand(i)->value() - 1;

		M->operand(i)->deleteOperand(i);
		M->operand(i)->includeOperand(integer(tt), i);
	}

	i = 0;
	j = 0;

	AST* v = list({});

	while(j < n)
	{
		while(isRowOfZeros(M, n, j) && j < n)
		{
			// AST* r = list({});
			// for(int k = 0; k < n; k++)
			// {
			// 	r->includeOperand(integer(0));
			// }
			// v->includeOperand(r);

			j = j + 1;
		}

		if(j < n)
		{
			i = i + 1;

			AST* r = list({});

			for(int k = 0; k < n; k++)
			{
				r->includeOperand(integer(-1 * M->operand(j)->operand(k)->value()));
			}

			v->includeOperand(r);
		}

		j = j + 1;
	}

	delete M;

	return v;
}

AST* monic_Zp(AST* f, AST* x, long p)
{
	if(f->is(0))
	{
		return integer(0);
	}

	AST* lc = leadCoeff(f, x);

	AST* F = quotientGPE_Zp(f, lc, x, p);

	return list({ lc, F });
}

AST* monic_sZp(AST* f, AST* x, long p)
{
	if(f->is(0))
	{
		return integer(0);
	}

	AST* lc = leadCoeff(f, x);

	AST* F = quotientGPE_sZp(f, lc, x, p);

	return list({ lc, F });
}


AST* extendedGCDGf(AST* f, AST* g, AST* x, long p)
{
	if(f->is(0) || g->is(0))
	{
		return list({integer(1), integer(0), integer(0)});
	}

	AST *t, *s, *p0, *i, *lc, *k0, *k1, *r0, *p1, *r1, *t0, *t1, *t2, *t3, *s0, *s1, *Q, *R, *T;

	t1 = monic_Zp(f, x, p); 
	t2 = monic_Zp(g, x, p); 

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
	
		t2 = integer(modInverse(p1->value(), p));
	
		t3 = r1;

		delete p0;
		delete p1;
		delete r0;
	
		return list({t1, t2, t3});
	}

	if(g->is(0))
	{
		t1 = integer(modInverse(p0->value(), p));
		t2 = integer(0);
		t3 = r0;

		delete p0;
		delete p1;
		delete r1;
	
		return list({t1, t2, t3});
	}

	s0 = integer(modInverse(p0->value(), p));
	s1 = integer(0);
	
	t0 = integer(0);
	t1 = integer(modInverse(p1->value(), p));

	while(true)
	{
		T = divideGPE_Zp(r0, r1, x, p);
		
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

		T = monic_Zp(R, x, p);
		
		delete r0;
	
		r0 = r1;
	
		lc = T->operand(0L);
		r1 = T->operand(1L);
		
		T->removeOperand(0L);
		T->removeOperand(0L);
		
		delete T;
	
		i = integer(modInverse(lc->value(), p));

		k0 = mulPoly(s1, Q);
	
		k1 = Zp(k0, x, p);
	
		delete k0;
		
		k0 = subPoly(s0, k1);

		delete k1;
	
		s = Zp(k0, x, p);

		delete k0;

		k0 = mulPoly(t1, Q);
	
		k1 = Zp(k0, x, p);
	
		delete k0;
		
		k0 = subPoly(t0, k1);
	
		delete k1;
	
		t = Zp(k0, x, p);

		delete k0;

		delete s0;
		delete t0;
	
		s0 = s1;
		t0 = t1;

		k0 = mulPoly(s, i);
		k1 = mulPoly(t, i);

		s1 = Zp(k0, x, p);
		t1 = Zp(k1, x, p);
		
		delete k0;
		delete k1;

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
