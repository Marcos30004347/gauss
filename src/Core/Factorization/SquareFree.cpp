#include "SquareFree.hpp"

#include "Core/Debug/Assert.hpp"
#include "Core/Polynomial/Zp.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
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

AST* squareFreeFactorizationFiniteField(AST* ax, AST* x, AST* q)
{
	assert(q->kind() == Kind::Power, "p is not a order of a Galois Field, should be q = p^m");

	AST* p = q->operand(0);

	unsigned int i = 1;

	AST* ox = integer(1);
	AST* ux = derivate(ax, x);

	AST* bx = Zp(ux, x, p->value());

	delete ux;

	if(bx->isNot(0))
	{
		AST* cx = gcdGPE_Zp(ax, bx, x, p->value());
		AST* wx = quotientGPE_Zp(ax, cx, x, p->value());

		while(wx->isNot(1))
		{
			AST* yx = gcdGPE_Zp(wx, cx, x, p->value());
			AST* zx = quotientGPE_Zp(wx, yx, x, p->value());

			ox = mul({ ox, power(zx, integer(i))});

			i = i + 1;

			delete wx;
			wx = yx;

			AST* kx = quotientGPE_Zp(cx, yx, x, p->value());

			delete cx;
			cx = kx;
		}

		if(cx->isNot(1))
		{
			AST* kx = add({});
			AST* deg = degree(cx, x);

			for(unsigned int i = 0; i <= deg->value(); i++)
			{
				AST* j = integer(i);

				kx->includeOperand(mul({
					coefficientGPE(cx, x, j),
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

		for(unsigned int i = 0; i <= deg->value(); i++)
		{
			AST* j = integer(i);
			
			kx->includeOperand(
				mul({
					coefficientGPE(ax, x, j),
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


}
