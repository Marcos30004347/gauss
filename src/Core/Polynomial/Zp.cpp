#include "Zp.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Expand/Expand.hpp"

using namespace ast;
using namespace expand;
using namespace algebra;

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

int modInverse_p(int a, int p) {
	int g = gcd(a, p);

	if (g != 1 && g!= -1) {
		printf("Inverse of %i in Z%i doesn't exist!\n", a, p);
		abort();
	} else {
		return pow(a, p - 2, p);
	}
}

int sZp(int b, int m) 
{
	b = mod(b, m);

	if(0 <= b && b <= m/2) 
	{
		return b;
	}

	return b - m;
}

int division_Zp(int s, int t, int p) {
	return mod((s * modInverse_p(t,p)), p);
}

int division_sZp(int s, int t, int p) {
	return sZp(mod(s * modInverse_p(t, p), p), p);
}

int mul_Zp(int s, int t, int p) {
	return mod((s * t), p);
}

int mul_sZp(int s, int t, int p) {
	return sZp(mod(s * t, p), p);
}

// Tnn in Symbolic Algebra, this function
// gives the non negative of u(x) 
// defined in Z[x] projected on Zp[x]
AST* Zp(AST* u, AST* x, int s) {
	AST* u_ = expandAST(u);

	AST* Tnn_u = new AST(Kind::Addition);

	AST* d = degreeGPE(u_, x);

	for(int i=0; i<=d->value(); i++) {
		AST* d = integer(i);

		AST* c = coefficientGPE(u_, x, d);
		
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
	AST* u_ = algebraicExpand(u);

	AST* Tnn_u = new AST(Kind::Addition);

	AST* d = degreeGPE(u_, x);

	for(int i=0; i<=d->value(); i++) {

		AST* e = integer(i);
		AST* c_ = coefficientGPE(u_, x, e);
		AST* c = expandAST(c_);

		Tnn_u->includeOperand(
			mul({
				integer(sZp(mod(c->value(), s), s)),
				power(x->copy(), e)
			})
		);

		delete c_;
		delete c;
	}

	AST* r = algebraicExpand(Tnn_u);
	
	delete d;
	delete u_;
	delete Tnn_u;

	return r;
}

AST* divideGPE_Zp(AST* u, AST* v, AST* x, int p) {
	AST* q = integer(0);
	AST* r = u->copy();

	AST* m = degreeGPE(r, x);
	AST* n = degreeGPE(v, x);

	AST* lcv = leadingCoefficientGPE(v, x);

	while(
		m->kind() != Kind::MinusInfinity &&
		(m->kind() == Kind::Integer && n->kind() == Kind::Integer &&
		m->value() >= n->value())
	) {
	
		AST* lcr = leadingCoefficientGPE(r, x);
	
		AST* s = integer(division_Zp(lcr->value(), lcv->value(), p));
		
		AST* q_ = add({
			q->copy(),
			mul({
				s->copy(),
				power(
					x->copy(),
					sub({
						m->copy(),
						n->copy()
					})
				)
			})
		});
		
		delete q;
		q = algebraicExpand(q_);
		delete q_;

		AST* r_ = sub({
			sub({
				r->copy(),
				mul({
					lcr->copy(),
					power(x->copy(), m->copy())
				})
			}),
			mul({
				sub({
					v->copy(),
					mul({
						lcv->copy(),
						power(x->copy(), n->copy())
					}),
				}),
				s->copy(),
				power(
					x->copy(),
					sub({m->copy(), n->copy()})
				)
			})
		});

		delete r;
		r = algebraicExpand(r_);
		delete r_;

		delete m;
		delete lcr;
		delete s;
	
		m = degreeGPE(r, x);
	}

	AST* res = list({ Zp(q, x, p), Zp(r, x, p) });
	
	delete q;
	delete r;

	delete m;
	delete n;
	delete lcv;

	return res;
}


AST* divideGPE_sZp(AST* u, AST* v, AST* x, int p) {

	AST* q = integer(0);
	AST* r = u->copy();

	AST* m = degreeGPE(r, x);
	AST* n = degreeGPE(v, x);

	AST* lcv = leadingCoefficientGPE(v, x);

	while(
		m->kind() != Kind::MinusInfinity &&
		(m->kind() == Kind::Integer && n->kind() == Kind::Integer &&
		m->value() >= n->value())
	) {
	
		AST* lcr = leadingCoefficientGPE(r, x);
	
		AST* s = integer(division_sZp(lcr->value(), lcv->value(), p));

		AST* q_ = add({
			q->copy(),
			mul({
				s->copy(),
				power(
					x->copy(),
					sub({
						m->copy(),
						n->copy()
					})
				)
			})
		});
		
		
		delete q;
	
		q = algebraicExpand(q_);
	
		delete q_;

		AST* r_ = sub({
			sub({
				r->copy(),
				mul({
					lcr->copy(),
					power(x->copy(), m->copy())
				})
			}),
			mul({
				sub({
					v->copy(),
					mul({
						lcv->copy(),
						power(x->copy(), n->copy())
					}),
				}),
				s->copy(),
				power(
					x->copy(),
					sub({m->copy(), n->copy()})
				)
			})
		});

		delete r;

		r = algebraicExpand(r_);

		delete r_;
		delete m;
		delete lcr;
		delete s;
	
		m = degreeGPE(r, x);
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

	AST* lco = leadingCoefficientGPE(U, x);
	
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

	AST* lco = leadingCoefficientGPE(U,x);
	
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

	AST* c = leadingCoefficientGPE(U, x);

	AST* App__ = mul({ App->copy(), integer(modInverse_p(c->value(), p)) });
	AST* App_ = Zp(App__, x, p);
	delete App;
	delete App__;
	App = App_;

	AST* Bpp__ = mul({ Bpp->copy(), integer(modInverse_p(c->value(), p)) });
	AST* Bpp_ = Zp(Bpp__, x, p);
	delete Bpp;
	delete Bpp__;
	Bpp = Bpp_;
	
	AST* U__ = mul({U->copy(), integer(modInverse_p(c->value(), p))});
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

	AST* c = leadingCoefficientGPE(U, x);

	AST* App__ = mul({ App->copy(), integer(modInverse_p(mod(c->value(),p), p)) });
	delete App;
	App = sZp(App__, x, p);
	delete App__;

	AST* Bpp__ = mul({ Bpp->copy(), integer(modInverse_p(mod(c->value(),p), p)) });
	delete Bpp;
	Bpp = sZp(Bpp__, x, p);
	delete Bpp__;
	
	AST* U__ = mul({U->copy(), integer(modInverse_p(mod(c->value(),p), p)) });
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

}
