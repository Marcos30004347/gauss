#include "Core/Simplification/Simplification.hpp"
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

	AST* d = add({
		m->copy(),
		mul({integer(-1), n->copy()}),
		integer(1)
	});

	AST* r = remainderGPE(u, v, x);

	if(r->is(0))
	{
		delete r;
		delete d;
	
		return integer(0);
	}

	AST* s = degreeGPE(r, x);

	AST* e = mul({
		power(integer(-1), mul({m->copy(), n->copy()})),
		multivariateResultant(v, r, L, K)
	});

	AST* w = algebraicExpand(e);

	AST* l = leadingCoefficientGPE(v, x);

	AST* k = add({
	 mul({ d->copy(), n->copy() }),
	 mul({integer(-1), m->copy()}),
	 s->copy()
	});

	AST* z = power(l, k);

	AST* f = algebraicExpand(z);

	AST* g = recQuotient(w, f, L, K);

	delete x;
	delete m;
	delete n;
	delete s;
	delete e;
	delete z;
	delete f;
	delete d;
	delete r;
	delete w;

	return g;
}


AST* srPolynomialResultantRec(AST* u, AST* v, AST* L, AST* K, AST* i, AST* dp, AST* gp)
{
	AST* x = L->operand(0);
	AST* m = degreeGPE(u, x);
	AST* n = degreeGPE(v, x);


	if(m->value() < n->value())
	{
		AST* e = mul({
			power(integer(-1), mul({m->copy(), n->copy()})),
			srPolynomialResultantRec(v, u, L, K, i, dp, gp)
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

	AST* d = add({
		m->copy(),
		mul({integer(-1), n->copy()}),
		integer(1)
	});

	AST* R = rest(L);

	AST* gamma 	= nullptr;
	AST* beta 	= nullptr;

	AST* tmp 		= nullptr;
	
	if(i->is(1))
	{
		gamma = integer(-1);
		tmp = power(integer(-1), d->copy());
		beta = reduceAST(tmp);
		
		delete tmp;
	}
	else
	{
		AST* f = leadingCoefficientGPE(u, x);
	
		AST* tmp1 = power(mul({integer(-1), f->copy()}), sub({dp->copy(), integer(1)}));
		AST* tmp2 = power(gamma, sub({dp->copy(), integer(2)}));
	
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
			power(gamma->copy(), sub({d->copy(), integer(1)}))
		});
	
		beta = algebraicExpand(tmp5);

		delete tmp5;

		delete f;
	}


	AST* t = recQuotient(r, beta, L, K);

	delete r;
	r = t;

	AST* tmp1 = add({i->copy(), integer(1)});
	AST* tmp2 = reduceAST(tmp1);
	
	delete tmp1;

	AST* tmp3 = mul({
		power(integer(-1), mul({ m->copy(), n->copy() })),
		power(beta->copy(), n->copy()),
		srPolynomialResultantRec(v, r, L, K, tmp2, d, gamma)
	});
	
	delete tmp2;

	AST* w = algebraicExpand(tmp3);

	delete tmp3;

	AST* l = coefficientGPE(v, x, n);

	AST* s = degreeGPE(r, x);

	delete r;

	AST* k = add({
	 mul({ d->copy(), n->copy() }),
	 mul({integer(-1), m->copy()}),
	 s->copy()
	});

	AST* tmp4 = power(l->copy(), k->copy());

	delete l;
	delete k;

	AST* f = algebraicExpand(tmp4);

	delete tmp4;

	AST* o = recQuotient(w, f, L, K);

	delete w;
	delete f;
	delete L;
	delete K;

	delete m;
	delete n;

	delete gamma;
	delete beta;

	return o;
}

AST* srPolynomialResultant(AST*	u, AST* v, AST* L, AST* K)
{

	AST* x = L->operand(0);
	AST* R = rest(L);

	AST* m = degreeGPE(u, x);
	AST* n = degreeGPE(v, x);

	AST* cont_u = polynomialContent(u, x, R, K);
	AST* pp_u = recQuotient(u, cont_u, L, K);
	AST* cont_v = polynomialContent(v, x, R, K);
	AST* pp_v = recQuotient(v, cont_v, L, K);

	AST* i = integer(1);
	AST* d = integer(0);
	AST* g = integer(0);

	printf("***************\n");

	AST* s = srPolynomialResultantRec(pp_u, pp_v, L, K, i, d, g);

	printf("***************\n");
	delete i;
	delete d;
	delete g;

	AST* t = mul({
		power(cont_u, n),
		power(cont_v, m),
		s
	});

	AST* k = reduceAST(t);

	delete x;
	delete R;
	delete pp_u;
	delete pp_v;

	return k;
}

}
