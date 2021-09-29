#include "Core/Simplification/Simplification.hpp"
#include "Core/Debug/Assert.hpp"
#include "Resultant.hpp"

using namespace ast;
using namespace algebra;
using namespace simplification;

namespace polynomial {

AST* univariateResultant(AST* ux, AST* vx, AST* x)
{
	AST* m = degreeGPE(ux, x);
	AST* n = degreeGPE(vx, x);

	if(n->is(0))
	{
		delete n;
	
		return power(vx->copy(), m);
	}


	AST* r = remainderGPE(ux, vx, x);
	
	if(r->is(0))
	{
		delete m;
		delete n;
		delete r;

		return integer(0);
	}

	AST* s = degreeGPE(r, x);
	AST* l = coefficientGPE(vx, x, n);

	AST* k = mul({
		power(integer(-1), mul({m->copy(), n->copy()})),
		power(l->copy(), sub({m->copy(), s->copy()})),
		univariateResultant(vx, r, x)
	});

	AST* e = algebraicExpand(k);

	delete m;
	delete n;
	delete l;
	delete s;
	delete k;
	delete r;

	return e;
}

AST* multivariateResultant(AST* u, AST* v, AST* L, AST* K)
{
	AST* x = L->operand(0)->copy();

	AST* m = degreeGPE(u, x);
	AST* n = degreeGPE(v, x);

	// if(m->isLessThan(n))
	if(m->value() < n->value())
	{
		AST* k = mul({
			power(integer(-1), mul({m->copy(), n->copy()})),
			multivariateResultant(v, u, L, K)
		});

		delete m;
		delete n;

		AST* t = algebraicExpand(k);

		delete k;
		delete x;

		return t;
	}

	if(n->is(0))
	{
		AST* k =  power(v->copy(), m->copy());
	
		delete m;
		delete n;
		delete x;
	
		return k;
	}

	AST* delta = add({
		m->copy(),
		mul({ integer(-1), n->copy() }),
		integer(1)
	});

	AST* r = pseudoRemainder(u, v, x);

	if(r->is(0))
	{
		delete r;
		delete delta;
	
		delete m;
		delete n;
		delete x;

		return integer(0);
	}

	AST* s = degreeGPE(r, x);

	AST* e = mul({
		power(integer(-1), mul({m->copy(), n->copy()})),
		multivariateResultant(v, r, L, K)
	});

	AST* w = algebraicExpand(e);
	
	delete e;

	AST* l = leadingCoefficientGPE(v, x);

	AST* k = add({ mul({ delta->copy(), n->copy() }), mul({integer(-1), m->copy()}), s});

	AST* z = power(l, k);

	AST* f = algebraicExpand(z);

	delete z;

	AST* g = recQuotient(w, f, L, K);

	delete x;
	delete m;
	delete n;

	delete f;
	delete delta;
	delete r;
	delete w;

	return g;
}


AST* srPolynomialResultantRec(AST* u, AST* v, AST* L, AST* K, AST* i, AST* delta_prev, AST* gamma_prev)
{
	assert(u->isNot(0), "Polynomial should be non-zero");
	assert(v->isNot(0), "Polynomial should be non-zero");

	AST* x = L->operand(0)->copy();
	AST* m = degreeGPE(u, x);
	AST* n = degreeGPE(v, x);


	if(m->value() < n->value())
	{
		AST* e = mul({
			power(integer(-1), mul({m->copy(), n->copy()})),
			srPolynomialResultantRec(v, u, L, K, i, delta_prev, gamma_prev)
		});
	
	
		AST* k = reduceAST(e);

		delete x;
		delete m;
		delete n;
		delete e;
	
		return k;
	}

	if(n->is(0))
	{
		AST* e = power(v->copy(), m->copy());

		AST* k = reduceAST(e);

		delete x;
		delete m;
		delete n;
		delete e;

		return k;
	}

	AST* r = pseudoRemainder(u, v, x);

	if(r->is(0))
	{
	
		delete r;
		delete x;
		delete m;
		delete n;
	
		return integer(0);
	}

	AST* delta = add({
		m->copy(),
		mul({integer(-1), n->copy()}),
		integer(1)
	});

	AST* R = rest(L);

	AST* gamma 	= nullptr;
	AST* beta 	= nullptr;

	if(i->is(1))
	{
		gamma = integer(-1);

		AST* tmp = power(integer(-1), delta->copy());

		beta = reduceAST(tmp);

		delete tmp;
	}
	else
	{
		AST* f = leadingCoefficientGPE(u, x);
	
		AST* tmp1 = power(mul({integer(-1), f->copy()}), sub({delta_prev->copy(), integer(1)}));
		AST* tmp2 = power(gamma_prev->copy(), sub({delta_prev->copy(), integer(2)}));
	
		AST* tmp3 = algebraicExpand(tmp1);
		AST* tmp4 = algebraicExpand(tmp2);
		
		delete tmp1;
		delete tmp2;
		
		gamma = recQuotient(tmp3, tmp4, R, K);
	
		delete tmp3;
		delete tmp4;
	
		AST* tmp5 = mul({
			integer(-1),
			f->copy(),
			power(gamma->copy(), sub({ delta->copy(), integer(1) }))
		});
	
		beta = algebraicExpand(tmp5);

		delete tmp5;

		delete f;
	}

	AST* t = recQuotient(r, beta, L, K);

	delete r;
	r = t;

	// Note: original algorithm didnt have this here,
	// but on testing there is a case where r = 0, 
	// so it will not be possible to run the 
	// srPolynomialResultantRec method with r
	// in the following lines.
	if(r->is(0))
	{
		delete gamma;
		delete beta;
		delete m;
		delete n;
		delete x;
		delete delta;
		delete R;
		delete r;

		return integer(0);
	}

	AST* tmp1 = add({i->copy(), integer(1)});
	AST* tmp2 = reduceAST(tmp1);
	
	delete tmp1;

	AST* tmp3 = mul({
		power(integer(-1), mul({ m->copy(), n->copy() })),
		power(beta->copy(), n->copy()),
		srPolynomialResultantRec(v, r, L, K, tmp2, delta, gamma)
	});

	delete tmp2;

	AST* w = algebraicExpand(tmp3);

	delete tmp3;

	AST* l = coefficientGPE(v, x, n);

	AST* s = degreeGPE(r, x);

	delete r;

	AST* k = add({
	 mul({ delta->copy(), n->copy() }),
	 mul({integer(-1), m->copy()}),
	 s->copy()
	});

	AST* tmp4 = power(l->copy(), k->copy());

	delete l;
	delete k;

	AST* f = algebraicExpand(tmp4);

	delete tmp4;

	AST* o = recQuotient(w, f, L, K);
	
	delete s;

	delete delta;

	delete w;
	delete f;

	delete m;
	delete n;

	delete gamma;
	delete beta;

	delete x;

	delete R;

	return o;
}

AST* srPolynomialResultant(AST*	u, AST* v, AST* L, AST* K)
{
	AST* x = L->operand(0)->copy();
	AST* R = rest(L);

	AST* m = degreeGPE(u, x);
	AST* n = degreeGPE(v, x);

	AST* cont_u = polynomialContentSubResultant(u, x, R, K);
	AST* pp_u = recQuotient(u, cont_u, L, K);
	AST* cont_v = polynomialContentSubResultant(v, x, R, K);
	AST* pp_v = recQuotient(v, cont_v, L, K);
	
	AST* i = integer(1);
	AST* delta = integer(0);
	AST* g = integer(0);

	AST* s = srPolynomialResultantRec(pp_u, pp_v, L, K, i, delta, g);

	delete i;
	delete delta;
	delete g;

	AST* t = mul({
		power(cont_u, n),
		power(cont_v, m),
		s
	});

	AST* k = reduceAST(t);

	delete t;
	delete x;
	delete R;
	delete pp_u;
	delete pp_v;

	return k;
}


AST* polyRemSeqRec(AST* Gi2, AST* Gi1, AST* L, AST* hi2, AST* K)
{
	AST *Gi, *hi1, *di2, *t1, *t2, *t3, *t4, *t5, *t6, *nk, *cnt, *ppk, *R, *r, *x;

	x = L->operand(0);

	if(Gi1->is(0))
	{
		return list({ integer(1), integer(0) });
	}

	t4 = pseudoRemainder(Gi2, Gi1, x);
	printf("r = %s\n", t4->toString().c_str());
	
	if(t4->is(0))
	{
		delete t4;

		nk = degreeGPE(Gi1, x);

		if(nk->value() > 0)
		{
			R = rest(L);
			
			cnt = polynomialContentSubResultant(Gi1, x, R, K);
		
			ppk = recQuotient(Gi1, cnt, L, K);
			
			r = list({ ppk, integer(0) });

			delete cnt;
			delete R;
		}
		else
		{
			r = list({integer(1), Gi1->copy()});
		}		

		delete nk;
	
		return r;
	}

	di2 = sub({ degreeGPE(Gi2, x), degreeGPE(Gi1, x) });

	t1 = add({ di2->copy(), integer(1) });
	t2 = power(integer(-1), t1);
	t3 = mul({t2, t4});
	t4  = algebraicExpand(t3);

	delete t3;
	
	t1 = leadingCoefficientGPE(Gi2, x);
	t2 = power(hi2->copy(), di2->copy());
	t3 = mul({t1, t2});

	t5 = algebraicExpand(t3);

	delete t3;
	printf("A\n");
	Gi = recQuotient(t4, t5, L, K);
	printf("A\n");

	printf("Gi %s\n", Gi->toString().c_str());
	

	delete t4;
	delete t5;

	t1 = leadingCoefficientGPE(Gi1, x);
	t2 = power(t1, di2->copy());

	t3 = hi2->copy();
	t4 = sub({integer(1), di2->copy()});
	t5 = power(t3, t4);
	t6 = mul({t2, t5});

	hi1 = algebraicExpand(t6); // h4

	delete t6;

	// AST* hi, *di1, *gi, *Hi;
	// AST* _Gi = op->operand(1)->copy();

	// di1 = sub({ degreeGPE(Gi1, x), degreeGPE(_Gi, x) });
	
	// t1 = leadingCoefficientGPE(_Gi, x);
	// t2 = power(t1, di1->copy()); // gi^d[i-1]

	// t3 = sub({integer(1), di1->copy()});
	// t4 = power(hi1->copy(), t3);
	// t5 = mul({t2, t4});

	// hi = algebraicExpand(t5);

	// delete t5;

	// gi = leadingCoefficientGPE(_Gi, x);
	
	// t1 = mul({hi->copy(), _Gi->copy()});
	// t2 = algebraicExpand(t1);


	// Hi = recQuotient(t2, gi, L, K);

	// printf("Hi = %s\n", Hi->toString().c_str());

	// delete t1;
	// delete t2;

	t3 = polyRemSeqRec(Gi1, Gi, L, hi1, K);

	delete Gi;
	delete hi1;
	delete di2;

	return t3;
}

AST* polyRemSeq(AST* F1, AST* F2, AST* L, AST* K)
{
	if(F1->kind() == Kind::Integer && F2->kind() == Kind::Integer)
	{
		return integerGCD(F1, F2);
	}

	AST* x = L->operand(0);
	AST* m = degreeGPE(F1, x);
	AST* n = degreeGPE(F2, x);
	
	if(m->value() < n->value())
	{
		delete m;
		delete n;

		return polyRemSeq(F2, F1, L, K);
	}
	
	printf("f = %s\n", F1->toString().c_str());
	printf("g = %s\n", F2->toString().c_str());
	delete m;
	delete n;

	AST *t1, *t2, *t3, *t4, *t5;
	AST *G1, *G2, *G3, *h2, *nk, *ppk, *cnt, *R, *r;

	G1 = F1;
	G2 = F2;

	// compute G[3]
	t1 = sub({ degreeGPE(F1, x), degreeGPE(F2, x) });
	t2 = add({ t1, integer(1) });

	t3 = power(integer(-1), t2);

	t4 = pdiv(G1, G2, x)->operand(1)->copy();
	printf("r = %s\n", t4->toString().c_str());

	t5 = mul({t3, t4});

	printf("A\n");
	G3 = algebraicExpand(t5);
	printf("A\n");


	delete t5;

	if(G3->is(0))
	{
		nk = degreeGPE(G2, x);

		if(nk->value() > 0)
		{
			R = rest(L);

			cnt = polynomialContent(G2, x, R, K);
			ppk = recQuotient(G2, cnt, L, K);
		
			r = list({ppk, integer(0)});

			delete cnt;
			delete R;
		}
		else
		{
			r = list({integer(1), G2->copy()});
		}		

		delete nk;
	
		return r;
	}

	// compute h[2]
	t1 = leadingCoefficientGPE(G2, x);
	t2 = sub({ degreeGPE(F1, x), degreeGPE(F2, x) });
	t3 = power(t1, t2);
	h2 = reduceAST(t3);

	delete t3;
	printf("Gi %s\n", G3->toString().c_str());
	t3 = polyRemSeqRec(G2, G3, L, h2, K);
	
	delete G3;
	delete h2;

	return t3;
}

}
