#include "SquareFree.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/GaloisField/GaloisField.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace galoisField;
using namespace simplification;

namespace factorization {

AST* squareFreeFactorization(AST* ax, AST* x)
{
	AST *ox, *bx, *cx, *wx, *yx, *zx, *qx, *tx;

	long i = 1;

	ox = integer(1);

	bx = derivate(ax, x);
	cx = gcdGPE(ax, bx, x);
	wx = quotientGPE(ax, cx, x);

	while(cx->isNot(1))
	{
		yx = gcdGPE(wx, cx, x);
		zx = quotientGPE(wx, yx, x);

		ox = mul({ ox, power(zx, integer(i)) });

		i = i + 1;

		delete wx;

		wx = yx;

		qx = quotientGPE(cx, yx, x);

		delete cx;

		cx = qx;
	}

	delete bx;
	delete cx;

	ox = mul({ ox , power(wx, integer(i)) });

	tx = reduceAST(ox);

	delete ox;

	return tx;
}

AST* squareFreeFactorization2(AST* ax, AST* x)
{
	AST *ox, *bx, *cx, *wx, *yx, *kx, *zx, *gx, *tx, *rx, *ux;
	
	unsigned int i = 1;

	ox = integer(1);
	
	bx = derivate(ax, x);
	cx = gcdGPE(ax, bx, x);

	if(cx->is(1))
	{
		wx = ax->copy();
	}
	else
	{
		wx = quotientGPE(ax, cx, x);
		yx = quotientGPE(bx, cx, x);

		kx = derivate(wx, x);
		zx = subPoly(yx, kx);

		delete kx;

		while(zx->isNot(0))
		{
			gx = gcdGPE(wx, zx, x);
	
			ox = mul({ ox, power(gx->copy(), integer(i))});

			i = i + 1;

			tx = quotientGPE(wx, gx, x);

			delete wx;
	
			wx = tx;

			yx = quotientGPE(zx, gx, x);

			rx = derivate(wx, x);

			delete zx;
			zx = subPoly(yx, rx);

			delete rx;
			delete gx;
		}

		delete zx;
	}

	ox = mul({ox, power(wx, integer(i))});

	ux = reduceAST(ox);

	delete cx;
	delete bx;
	delete ox;

	return ux;
}

AST* squareFreeFactorizationFiniteField(AST* ax, AST* x, AST* q, bool symmetric)
{
	AST* p = q->copy();

	unsigned int i = 1;

	AST* ox = integer(1);
	AST* ux = derivate(ax, x);

	AST* bx = gf(ux, x, p->value());

	delete ux;

	if(bx->isNot(0))
	{
		AST* cx = gcdPolyGf(ax, bx, x, p->value(), symmetric);
		AST* wx = quoPolyGf(ax, cx, x, p->value(), symmetric);

		while(wx->isNot(1))
		{
			AST* yx = gcdPolyGf(wx, cx, x, p->value(), symmetric);
			AST* zx = quoPolyGf(wx, yx, x, p->value(), symmetric);

			ox = mul({ ox, power(zx, integer(i))});

			i = i + 1;

			delete wx;
			wx = yx;

			AST* kx = quoPolyGf(cx, yx, x, p->value(), symmetric);

			delete cx;
			cx = kx;
		}

		if(cx->isNot(1))
		{
			AST* kx = add({});
			AST* deg = degree(cx, x);

			for(Int i = 0; i <= deg->value(); i++)
			{
				AST* j = integer(i);

				kx->includeOperand(mul({
					coeff(cx, x, j),
					power(x->copy(), integer(i/p->value()))
				}));

				delete j;
			}

			delete cx;
			cx = reduceAST(kx);

			delete deg;
			delete kx;

			AST* sx = squareFreeFactorizationFiniteField(cx, x, q);

			delete cx;
			cx = sx;

			ox = mul({ ox, power(cx->copy(), integer(p->value())) });
		}

		delete cx;
		delete wx;
	}
	else
	{
		AST* deg = degree(ax, x);
		AST* kx = add({});

		for(Int i = 0; i <= deg->value(); i++)
		{
			AST* j = integer(i);
			
			kx->includeOperand(
				mul({
					coeff(ax, x, j),
					power(x->copy(), integer(i/p->value()))
				})
			);

			delete j;
		}

		delete deg;
		delete ax;
	
		ax = kx;

		AST* sx = squareFreeFactorizationFiniteField(ax, x, q);

		delete ox;
	
		ox = power(sx, integer(p->value()));
	}

	AST* tx = reduceAST(ox);

	delete ox;
	delete bx;

	return tx;
}

bool isSquareFreeInZp(AST* f, AST* x, long p, bool symmetric)
{
	bool r = false;

	AST *lc, *t, *k, *v, *g;

	if(f->is(0))
	{
		return true;
	}
 	
	lc = leadCoeff(f, x);
	
	v = quoPolyGf(f, lc, x, p, symmetric);

	delete lc;

	k = derivate(v, x);
	
	t = gf(k, x, p, symmetric);
	
	delete k;

	g = gcdPolyGf(v, t, x, p, symmetric);

	delete t;
	delete v;

	r = g->is(1);

	delete g;

	return r;	
}

AST* squareFreePart(AST* f, AST* L, AST* K)
{
	AST *g, *u, *v, *s;
	
	long i;

	g = f->copy();

	for(i = 0; i < L->numberOfOperands(); i++)
	{
		printf("f %s\n", f->toString().c_str());
		u = reduceAST(derivate(f, L->operand(i)));
		printf("diff %s\n", u->toString().c_str());
		v = mvPolyGCD(g, u, L, K);
		printf("gcd %s\n", v->toString().c_str());
		
		delete g;
		
		g = v;
		
		delete u;
	}
	printf("quo\n");

	s = recQuotient(f, g, L, K);

	delete g;

	g = pp(s, L, K);

	delete s;
	
	AST* R = list({});

	for(i = 0; i < L->numberOfOperands(); i++)
	{
		if(!g->freeOf(L->operand(i)))
		{
			R->includeOperand(L->operand(i)->copy());
		}
	}

	return list({g, R});
}

bool isSquareFree(ast::AST* f, ast::AST* x, ast::AST* K)
{
	long e = 1;

	AST *k, *n, *g;

	if(f->is(0))
	{
		return true;
	}

	k = derivate(f, x);
	g = gcdGPE(f, k, x);
	n = degree(g, x);
	
	e = n->value().longValue();
	
	delete k;
	delete g;
	delete n;

	return !e;
}


}
