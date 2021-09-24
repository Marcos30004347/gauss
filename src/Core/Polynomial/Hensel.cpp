#include "Hensel.hpp"
#include "Zp.hpp"
#include "Core/Algebra/List.hpp"

using namespace ast;
using namespace algebra;

namespace polynomial {

AST* leadCoeffReplace(AST* ux, AST* x, AST* c)
{
	AST* lc = leadingCoefficientGPE(ux, x);
	
	AST* de = degreeGPE(ux, x);

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
	AST* lc = leadingCoefficientGPE(ux, x);

	AST* px = mul({ power(lc, integer(-1)), ux->copy() });

	AST* kx = algebraicExpand(px);

	delete px;

	return kx;
}

AST* univariateHensel(AST* ax, AST* x, AST* p, AST* ux_1, AST* wx_1, AST* B, AST* zeta)
{
	// printf("************\n");

	AST* tmp = nullptr;
	AST* gam = nullptr;
	AST* tal = nullptr;
	AST* qx  = nullptr;
	AST* rx  = nullptr;
	AST* cx  = nullptr;

	AST* Z = symbol("Z");

	AST* L = list({ x->copy() });

	// 1. Define polynomial and its modulo p factors
	AST* alph = leadingCoefficientGPE(ax, x);

	if(zeta->kind() == Kind::Undefined)
	{
		zeta = alph->copy();
	}
	else
	{
		zeta = zeta->copy();
	}

	tmp = mul({ zeta->copy(), ax->copy() });
	
	ax = algebraicExpand(tmp);

	delete tmp;

	// printf("a = %s\n", alph->toString().c_str());
	// printf("y = %s\n", zeta->toString().c_str());
	// printf("a(x) = %s\n", ax->toString().c_str());
	// printf("u[1](x) = %s\n", ux_1->toString().c_str());
	// printf("w[1](x) = %s\n", wx_1->toString().c_str());

	// Normalization maybe wrong 
	AST* kx = mul({ zeta->copy(), normalize(ux_1, x) });
	AST* zx = mul({ alph->copy(), normalize(wx_1, x) });

	ux_1 = sZp(kx, x, p->value());
	wx_1 = sZp(zx, x, p->value());

	delete kx;
	delete zx;

	// printf("u[1](x) = %s\n", ux_1->toString().c_str());
	// printf("w[1](x) = %s\n", wx_1->toString().c_str());

	// 2. Apply extended Euclidean algorithm to ux_1, wx_2 defined in Zp[x]
	AST* l = extendedEuclideanAlgGPE_sZp(ux_1, wx_1, x, p->value());
	
	AST* sx = l->operand(1)->copy();
	AST* tx = l->operand(2)->copy();

	// printf("s(x) = %s\n", sx->toString().c_str());
	// printf("t(x) = %s\n", tx->toString().c_str());

	delete l;

	// 3. Initialization for iteration
	AST* ux = leadCoeffReplace(ux_1, x, zeta);
	AST* wx = leadCoeffReplace(wx_1, x, alph);

	AST* fx = sub({ ax->copy(), mul({ ux->copy(), wx->copy() }) });
	
	AST* ex = algebraicExpand(fx);

 	// printf("u(x) = %s\n", ux->toString().c_str());
 	// printf("w(x) = %s\n", wx->toString().c_str());
	// printf("e(x) = %s\n", ex->toString().c_str());
	
	delete fx;

	long modulus = p->value();

	// 4. Iterate until either the factorization in Z[x] is obtained 
	// or else the bound on modulus is reached
	while(ex->isNot(0) && modulus < 2 * B->value() * zeta->value())
	{
		// 4.1 Solve in the domain Zp[x] the polynomial equation
		tmp = div(ex->copy(), integer(modulus)); // MAYBE DIV IN sZp[x]

		cx = algebraicExpand(tmp);

		delete tmp;

		tmp = mul({ sx->copy(), cx->copy() });

		delete gam;

		gam = sZp(tmp, x, p->value()); 

		delete tmp;

		tmp = mul({ tx->copy(), cx->copy() });

		delete tal;
	
		tal = sZp(tmp, x, p->value()); 

		delete tmp;

		tmp = divideGPE_sZp(gam, wx_1, x, p->value());

		qx = tmp->operand(0)->copy();
		rx = tmp->operand(1)->copy();

		delete tmp;

		delete gam;

		gam = rx;

		tmp = add({ tal->copy(), mul({ qx, ux_1->copy() }) });

		delete tal;

		tal = sZp(tmp, x, p->value());

		delete tmp;

		// 4.2 Update factors and compute the error
		ux = add({ ux, mul({ tal->copy(), integer(modulus) })});
		wx = add({ wx, mul({ gam->copy(), integer(modulus) })});

		tmp = sub({ ax->copy(), mul({ ux->copy(), wx->copy() }) });
	
		delete ex;

		ex = algebraicExpand(tmp);

		delete tmp;
	
		tmp = algebraicExpand(ux);
		
		delete ux;
		
		ux = tmp;

		tmp = algebraicExpand(wx);

		delete wx;

		wx = tmp;
	
		// printf("\n");
		// printf("%li\n", modulus);
		// printf("tal(x) = %s\n", tal->toString().c_str());
		// printf("gam(x) = %s\n", gam->toString().c_str());
		// printf("u(x) = %s\n",   ux->toString().c_str());
		// printf("w(x) = %s\n",   wx->toString().c_str());
		// printf("e(x) = %s\n",   ex->toString().c_str());
	
		modulus = modulus * p->value();

		delete cx;
	}
	

	AST* lf = nullptr;	
	
	// 5. Check termination status
	if(ex->is(0))
	{
		// Factorization obtained - remove contents
		AST* d = cont(ux, x);

		AST* px = quotientGPE(ux, d, x);
	
		AST* k = integer(zeta->value() / d->value());

		AST* kx = quotientGPE(wx, k, x);

		delete d;
		delete k;

		// Note: a(x) <- a(x)/y would restore a(x) to its input value
		lf = list({ px, kx });
	}
	else
	{
		lf = new AST(Kind::Fail);
	}
	delete ux;
	delete wx;

	delete gam;
	delete tal;
	delete sx;
	delete tx;

	delete ex;
	delete Z;
	delete L;
	delete ax;
	delete zeta;
	delete alph;
	delete ux_1;
	delete wx_1;

	return lf;
}

}
