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

AST* henselSep(AST* f, AST* g, AST* h, AST* s, AST* t, AST* x, Int m, bool symmetric)
{
	AST *one, *e, *q, *r, *G, *H, *b, *c, *d, *S, *T, *t1, *t2, *t3, *t4;

	printf("m =%li\n", m);
	printf("f = %s\n", f->toString().c_str());
	printf("g = %s\n", g->toString().c_str());
	printf("h = %s\n", h->toString().c_str());
	printf("s = %s\n", s->toString().c_str());
	printf("t = %s\n", t->toString().c_str());

	one = integer(1);

	t2 = mulPolyGf(g, h, x, m * m, symmetric);

	// t2 = mulPoly(g, h);
	printf("t2 = %s\n", t2->toString().c_str());

	t3 = subPoly(f, t2);

	e  = subPolyGf(f, t2, x, m * m, symmetric);

	printf("e = %s\n", e->toString().c_str());
	delete t2;

	t2 = mulPolyGf(s, e, x, m * m, symmetric);
	t1 = divPolyGf(t2, h, x, m * m, symmetric);

	delete t2;

	q = t1->operand(0);
	r = t1->operand(1);

	t1->removeOperand(0L);
	t1->removeOperand(0L);

	delete t1;

	t2 = mulPoly(t, e);
	t3 = mulPoly(q, g);

	t4 = addPoly(t2, t3);
	printf("t2 = %Lf\n", 89420766.0L * -45566718075060.0L);
	printf("t2 = %s\n", t->toString().c_str());
	printf("t2 = %s\n", e->toString().c_str());
	printf("t2 = %s\n", t2->toString().c_str());
	printf("t3 = %s\n", t3->toString().c_str());
	printf("t4 = %s\n", t4->toString().c_str());

	delete t2;
	delete t3;

	G = addPolyGf(g, t4, x, m * m, symmetric);
	H = addPolyGf(h, r, x, m * m, symmetric);

	printf("G = %s\n", G->toString().c_str());
	printf("H = %s\n", H->toString().c_str());

	delete t4;

	t2 = mulPolyGf(s, G, x, m*m, symmetric);
	t3 = mulPolyGf(t, H, x, m * m, symmetric);
	t4 = addPolyGf(t2, t3, x, m * m, symmetric);
	b  = subPolyGf(t4, one, x, m * m, symmetric);

	delete t2;
	delete t3;
	delete t4;

	delete one;

	t2 = mulPolyGf(s, b, x, m * m, symmetric);
	t3 = divPolyGf(t2, H, x, m * m, symmetric);

	c = t3->operand(0);
	d = t3->operand(1);

	t3->removeOperand(0L);
	t3->removeOperand(0L);

	delete t2;
	delete t3;

	S  = subPolyGf(s, d, x, m * m, symmetric); 
	t2 = mulPolyGf(t, b, x, m * m, symmetric);
	t3 = mulPolyGf(c, G, x, m * m, symmetric);
	t1 = addPolyGf(t2, t3, x, m * m, symmetric);
	T  = subPolyGf(t, t1, x, m * m, symmetric);
	
	delete t2;
	delete t3;
	delete t1;

	delete e;
	delete q;
	delete r;
	delete b;
	delete c;
	delete d;

	return list({ G, H, S, T });
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

AST* multifactorHenselLifting(AST* v, AST* H, AST* x, Int p, Int l, bool symmetric)
{
	printf("*****\n");

	printf("v = %s\n", v->toString().c_str());
	printf("H = %s\n", H->toString().c_str());
	printf("p^l = %s^%s\n", p.to_string(), l.to_string());
	
	Int i, j, r, k, d;

	Int a;

	AST *f, *fi, *lc, *t1, *g, *h, *s;
	AST *t, *e, *T, *H0, *H1, *F0, *F1, *F;

	lc = leadCoeff(v, x);

	f = v->copy();

	r = H->numberOfOperands();
	printf("A\n");

	if(r == 1)
	{
		a = inverseGf(lc->value(), pow(p, l), symmetric);

		t1 = integer(mod(a, pow(p, l), symmetric));
	
		fi = mulPolyGf(f, t1, x, pow(p, l), symmetric);
	
		delete t1;
	
		delete lc;
	
		delete f;
	
		return list({ fi });
	}
	printf("B\n");

	k = r / 2;
	d = std::ceil(log2(l.longValue()));
	
	g = lc->copy();
	h = gf(H->operand(k), p, symmetric);  //integer(1);

	for(i=0; i<k; i++)
	{
		t1 = mulPolyGf(g, H->operand(i), x, p, symmetric);
		delete g;
		g = t1;
	}

	printf("h1 = %s\n", h->toString().c_str());
	printf("%s\n", H->toString().c_str());
	for(i=k + 1; i<r; i++)
	{
		printf("fi = %s\n", H->operand(i)->toString().c_str());

		t1 = mulPolyGf(h, H->operand(i), x, p, symmetric);
		delete h;
		h = t1;
	}

	printf("---> g = %s\n", g->toString().c_str());
	printf("---> h = %s\n", h->toString().c_str());

	e = extendedEuclidGf(g, h, x, p, symmetric);
	
	s = e->operand(1);
	t = e->operand(2);

	e->removeOperand(2);
	e->removeOperand(1);

	delete e;
	printf("C\n");

	for(j = 1; j <= d; j++)
	{
		printf("aaa\n");
		T = henselSep(f, g, h, s, t, x, pow(p, pow(2, j - 1)), symmetric);
		
		delete g;
		delete h;
		delete s;
		delete t;
		
		g = T->operand(0);
		h = T->operand(1);
		s = T->operand(2);
		t = T->operand(3);
		printf("===> g = %s\n", g->toString().c_str());
		printf("===> h = %s\n", h->toString().c_str());
		T->removeOperand(0L);
		T->removeOperand(0L);
		T->removeOperand(0L);
		T->removeOperand(0L);
		
		delete T;
	}
	printf("D\n");

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
