#include "Hensel.hpp"

#include "Core/GaloisField/GaloisField.hpp"

#include <cmath>

using namespace ast;
using namespace algebra;
using namespace galoisField;
using namespace polynomial;

namespace factorization {

AST* leadCoeffReplace(AST* ux, AST* x, AST* c)
{
	AST* lc = leadCoeff(ux, x);
	
	AST* de = degree(ux, x);

	AST* px = mul({ lc, power(x->copy(), de->copy()) });

	AST* rx = sub({ux->copy(), px });

	AST* kx = algebraicExpand(rx);

	AST* ox = add({kx, mul({ c->copy(), power(x->copy(), de)})});

	AST* zx = algebraicExpand(ox);

	delete rx;
	delete ox;

	return zx;
}

AST* normalize(AST* ux, AST* x)
{
	AST* lc = leadCoeff(ux, x);

	AST* px = mul({ power(lc, integer(-1)), ux->copy() });

	AST* kx = algebraicExpand(px);

	delete px;

	return kx;
}

// AST* univariateHensel(AST* ax, AST* x, AST* p, AST* ux_1, AST* wx_1, AST* B, AST* zeta, bool symmetric)
// {

// 	AST* tmp = nullptr;
// 	AST* gam = nullptr;
// 	AST* tal = nullptr;
// 	AST* qx  = nullptr;
// 	AST* rx  = nullptr;
// 	AST* cx  = nullptr;

// 	AST* Z = symbol("Z");

// 	AST* L = list({ x->copy() });

// 	// 1. Define polynomial and its modulo p factors
// 	AST* alph = leadCoeff(ax, x);

// 	if(zeta->kind() == Kind::Undefined)
// 	{
// 		zeta = alph->copy();
// 	}
// 	else
// 	{
// 		zeta = zeta->copy();
// 	}

// 	tmp = mul({ zeta->copy(), ax->copy() });
	
// 	ax = algebraicExpand(tmp);

// 	delete tmp;

// 	// Normalization maybe wrong 
// 	AST* kx = normalize(ux_1, x);
// 	AST* zx = normalize(wx_1, x);

// 	ux_1 = mulPolyGf(zeta, kx, x, p->value(), symmetric);
// 	wx_1 = mulPolyGf(alph, zx, x, p->value(), symmetric);

// 	delete kx;
// 	delete zx;

// 	// printf("u[1](x) = %s\n", ux_1->toString().c_str());
// 	// printf("w[1](x) = %s\n", wx_1->toString().c_str());

// 	// 2. Apply extended Euclidean algorithm to ux_1, wx_2 defined in Zp[x]
// 	AST* l = extendedEuclidGf(ux_1, wx_1, x, p->value(), symmetric); //extendedEuclideanAlgGPE_sZp(ux_1, wx_1, x, p->value());
	
// 	AST* sx = l->operand(1)->copy();
// 	AST* tx = l->operand(2)->copy();

// 	// printf("s(x) = %s\n", sx->toString().c_str());
// 	// printf("t(x) = %s\n", tx->toString().c_str());

// 	delete l;

// 	// 3. Initialization for iteration
// 	AST* ux = leadCoeffReplace(ux_1, x, zeta);
// 	AST* wx = leadCoeffReplace(wx_1, x, alph);

// 	AST* fx = sub({ ax->copy(), mul({ ux->copy(), wx->copy() }) });
	
// 	AST* ex = algebraicExpand(fx);

// 	delete fx;

// 	long modulus = p->value();

// 	// 4. Iterate until either the factorization in Z[x] is obtained 
// 	// or else the bound on modulus is reached
// 	while(ex->isNot(0) && modulus < 2 * B->value() * zeta->value())
// 	{
// 		// 4.1 Solve in the domain Zp[x] the polynomial equation
// 		tmp = div(ex->copy(), integer(modulus)); // MAYBE DIV IN sZp[x]

// 		cx = algebraicExpand(tmp);

// 		delete tmp;

// 		// tmp = mul({ sx->copy(), cx->copy() });

// 		delete gam;

// 		gam = mulPolyGf(sx, cx, x, p->value(), symmetric);// sZp(tmp, x, p->value()); 

// 		// delete tmp;

// 		// tmp = mul({ tx->copy(), cx->copy() });

// 		delete tal;
	
// 		tal = mulPolyGf(tx, cx, x, p->value(), symmetric);// sZp(tmp, x, p->value()); 

// 		delete tmp;

// 		tmp = quoPolyGf(gam, wx_1, x, p->value(), symmetric);

// 		qx = tmp->operand(0)->copy();
// 		rx = tmp->operand(1)->copy();

// 		delete tmp;

// 		delete gam;

// 		gam = rx;

// 		tmp = mulPolyGf(qx, ux_1, x, p->value(), symmetric); //add({ tal->copy(), mul({ qx, ux_1->copy() }) });

// 		delete tal;

// 		tal = addPolyGf(tal, tmp, x, p->value(), symmetric); //sZp(tmp, x, p->value());

// 		delete tmp;

// 		// 4.2 Update factors and compute the error
// 		ux = add({ ux, mul({ tal->copy(), integer(modulus) })});
// 		wx = add({ wx, mul({ gam->copy(), integer(modulus) })});

// 		tmp = sub({ ax->copy(), mul({ ux->copy(), wx->copy() }) });
	
// 		delete ex;

// 		ex = algebraicExpand(tmp);

// 		delete tmp;
	
// 		tmp = algebraicExpand(ux);
		
// 		delete ux;
		
// 		ux = tmp;

// 		tmp = algebraicExpand(wx);

// 		delete wx;

// 		wx = tmp;

// 		modulus = modulus * p->value();

// 		delete cx;
// 	}
	

// 	AST* lf = nullptr;	
	
// 	// 5. Check termination status
// 	if(ex->is(0))
// 	{
// 		// Factorization obtained - remove contents
// 		AST* d = cont(ux, x);

// 		AST* px = quotientGPE(ux, d, x);
	
// 		AST* k = integer(zeta->value() / d->value());

// 		AST* kx = quotientGPE(wx, k, x);

// 		delete d;
// 		delete k;

// 		// Note: a(x) <- a(x)/y would restore a(x) to its input value
// 		lf = list({ px, kx });
// 	}
// 	else
// 	{
// 		lf = new AST(Kind::Fail);
// 	}
// 	delete ux;
// 	delete wx;

// 	delete gam;
// 	delete tal;
// 	delete sx;
// 	delete tx;

// 	delete ex;
// 	delete Z;
// 	delete L;
// 	delete ax;
// 	delete zeta;
// 	delete alph;
// 	delete ux_1;
// 	delete wx_1;

// 	return lf;
// }

AST* henselSep(AST* f, AST* g, AST* h, AST* s, AST* t, AST* x, long m, bool symmetric)
{
	AST *one, *e, *q, *r, *G, *H, *b, *c, *d, *S, *T, *t1, *t2, *t3, *t4;

	one = integer(1);

	t2 = mulPolyGf(g, h, x, m * m, symmetric);
	e  = subPolyGf(f, t2, x, m * m, symmetric);

	delete t2;

	t2 = mulPolyGf(s, e, x, m*m, symmetric);
	t1 = divPolyGf(t2, h, x, m * m, symmetric);

	delete t2;

	q = t1->operand(0);
	r = t1->operand(1);

	t1->removeOperand(0L);
	t1->removeOperand(0L);

	delete t1;

	t2 = mulPolyGf(t, e, x, m * m, symmetric); //sZp(t1, x, m * m);
	t3 = mulPolyGf(q, g, x, m * m, symmetric);//sZp(t1, x, m * m);

	t4 = addPolyGf(t2, t3, x, m * m, symmetric);

	delete t2;
	delete t3;

	G = addPolyGf(g, t4, x, m * m, symmetric); //sZp(t1, x, m * m);
	H = addPolyGf(h, r, x, m * m, symmetric); // sZp(t2, x, m * m);

	delete t4;

	t2 = mulPolyGf(s, G, x, m*m, symmetric); //sZp(t1, x, m * m);
	t3 = mulPolyGf(t, H, x, m * m, symmetric);// sZp(t1, x, m * m);
	t4 = addPolyGf(t2, t3, x, m * m, symmetric);
	b  = subPolyGf(t4, one, x, m * m, symmetric); //sZp(t1, x, m * m);

	delete t2;
	delete t3;
	delete t4;

	delete one;

	t2 = mulPolyGf(s, b, x, m * m, symmetric);// sZp(t1, x, m * m);
	t3 = divPolyGf(t2, H, x, m * m, symmetric); //divideGPE_sZp(t2, H, x, m * m);

	c = t3->operand(0);
	d = t3->operand(1);

	t3->removeOperand(0L);
	t3->removeOperand(0L);

	delete t2;
	delete t3;

	S  = subPolyGf(s, d, x, m * m, symmetric); //sZp(t1, x, m * m);
	t2 = mulPolyGf(t, b, x, m * m, symmetric);// sZp(t1, x, m * m);
	t3 = mulPolyGf(c, G, x, m * m, symmetric);// sZp(t1, x, m * m);
	t1 = addPolyGf(t2, t3, x, m * m, symmetric);
	T  = subPolyGf(t, t1, x, m * m, symmetric); // sZp(t2, x, m * m);
	
	delete t2;
	delete t3;
	delete t1;

	delete e;
	delete q;
	delete r;
	delete b;
	delete c;
	delete d;

	return list({G, H, S, T});
}

// long euclidExtended(long a, long b, long *x, long *y)
// {
//     if (a == 0)
//     {
//         *x = 0;
//         *y = 1;
//         return b;
//     }
 
//     long x1, y1; 
//     long gcd = euclidExtended(b%a, a, &x1, &y1);

//     *x = y1 - (b/a) * x1;
//     *y = x1;
 
//     return gcd;
// }

AST* multifactorHenselLifting(AST* v, AST* H, AST* x, long p, long l, bool symmetric)
{
	long i, j, r, k, d, a;

	AST *f, *fi, *lc, *t1, *t2, *g, *h, *s;
	AST *t, *e, *T, *H0, *H1, *F0, *F1, *F;

	lc = leadCoeff(v, x);

	f = quotientGPE(v, lc, x);

	r = H->numberOfOperands();

	if(r == 1)
	{
		// using extended euclid to compute 1 = lc mod(p^l)
		// euclidExtended(lc->value(), std::pow(p, l), &a, &b);
		
		// IN CASE OF ERROR, CHANGE SYMMETRIC TO ALWAYS BE TRUE
		a = inverseGf(lc->value(), std::pow(p, l), symmetric);

		t1 = integer(mod(a, std::pow(p, l), symmetric));
	
		fi = mulPolyGf(f, t1, x, std::pow(p, l), symmetric);
	
		delete t1;
	
		delete lc;
	
		delete f;
	
		return list({ fi });
	}

	k = std::floor(r / 2.0);
	d = std::ceil(log2(l));
	
	g = mul({ lc->copy() });
	h = mul({});

	for(i=0; i<k; i++)
	{
		g->includeOperand(H->operand(i)->copy());
	}

	for(i=k; i<r; i++)
	{
		h->includeOperand(H->operand(i)->copy());
	}

	t1 = gf(g, x, p, symmetric);
	t2 = gf(h, x, p, symmetric);

	delete g;
	delete h;

	g = t1;
	h = t2;

	e = extendedEuclidGf(g, h, x, p, symmetric);
	
	s = e->operand(1);
	t = e->operand(2);

	e->removeOperand(2);
	e->removeOperand(1);

	delete e;

	for(j = 1; j <= d; j++)
	{
		T = henselSep(f, g, h, s, t, x, std::pow(p, std::pow(2, j - 1)), symmetric);
		
		delete g;
		delete h;
		delete s;
		delete t;
		
		g = T->operand(0);
		h = T->operand(1);
		s = T->operand(2);
		t = T->operand(3);
		
		T->removeOperand(0L);
		T->removeOperand(0L);
		T->removeOperand(0L);
		T->removeOperand(0L);
		
		delete T;
	}

	delete s;
	delete t;

	H0 = list({});
	H1 = list({});

	for(i = 0; i < k; i++)
	{
		H0->includeOperand(H->operand(i)->copy());
	}

	for(i = k; i < r; i++)
	{
		H1->includeOperand(H->operand(i)->copy());
	}

	delete f;

	F0 = multifactorHenselLifting(g, H0, x, p, l, symmetric);
	F1 = multifactorHenselLifting(h, H1, x, p, l, symmetric);

	delete H0;
	delete H1;

	delete g;
	delete h;
	printf("%s\n", lc->toString().c_str());
	delete lc;

	F = list({});

	while(F0->numberOfOperands() > 0)
	{
		F->includeOperand(F0->operand(0));
		F0->removeOperand(0L);
	}

	while(F1->numberOfOperands() > 0)
	{
		F->includeOperand(F1->operand(0));
		F1->removeOperand(0L);
	}

	delete F0;
	delete F1;

	return F;
}


}
