#include "Zp.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Expand/Expand.hpp"

using namespace ast;
using namespace expand;
using namespace algebra;

namespace polynomial {

int mod(int a, int b) {
	return (b + (a%b)) % b;
}

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
	if (g != 1) {
		printf("Inverse of %i in Z%i doesn't exist!\n", a, p);
		abort();
	} else {
		return pow(a, p - 2, p);
	}
}

int division_Zp(int s, int t, int p) {
	return mod((s*modInverse_p(t,p)), p);
}

int S(int b, int m) {
	if(0 <= b && b <= m/2) {
		return b;
	}

	return b - m;
}

int division_Sp(int s, int t, int p) {
	return S(mod((s*modInverse_p(t,p)), p), p);
}

AST* Tnn(AST* u, AST* x, int s) {
	AST* u_ = expandAST(u);

	AST* Tnn_u = new AST(Kind::Addition);

	AST* d = degreeGPE(u_, x);

	for(int i=0; i<=d->value(); i++) {
		AST* d = integer(i);

		AST* c_ = coefficientGPE(u_, x, d);
		AST* c = expandAST(c_);
		
		Tnn_u->includeOperand(mul({ integer(mod(c->value(), s)), power(x->deepCopy(), d) }));
		
		delete c_;
		delete c;
	}


	AST* r = expandAST(Tnn_u);

	delete d;
	delete u_;
	delete Tnn_u;

	return r;
}

AST* Ts(AST* u, AST* x, int s) {
	AST* u_ = expandAST(u);

	AST* Tnn_u = new AST(Kind::Addition);

	AST* d = degreeGPE(u_, x);

	for(int i=0; i<=d->value(); i++) {

		AST* d = integer(i);
		AST* c_ = coefficientGPE(u_, x, d);
		AST* c = expandAST(c_);

		Tnn_u->includeOperand(mul({integer(S(mod(c->value(), s), s)), power(x->deepCopy(), d)}));
		
		delete c;
		delete c_;
	}

	AST* r = expandAST(Tnn_u);
	
	delete d;
	delete u_;
	delete Tnn_u;

	return r;
}

AST* divideGPE_Zp(AST* u, AST* v, AST* x, int p) {

	AST* q = integer(0);
	AST* r = u->deepCopy();

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
			q->deepCopy(),
			mul({
				s->deepCopy(),
				power(
					x->deepCopy(),
					sub({
						m->deepCopy(),
						n->deepCopy()
					})
				)
			})
		});
		
		
		delete q;
		q = Tnn(q_,x,p);
		delete q_;

		AST* r_ = sub({
			sub({
				r->deepCopy(),
				mul({
					lcr->deepCopy(),
					power(x->deepCopy(), m->deepCopy())
				})
			}),
			mul({
				sub({
					v->deepCopy(),
					mul({
						lcv->deepCopy(),
						power(x->deepCopy(), n->deepCopy())
					}),
				}),
				s->deepCopy(),
				power(
					x->deepCopy(),
					sub({m->deepCopy(), n->deepCopy()})
				)
			})
		});

		delete r;
		r = Tnn(r_, x, p);
		delete r_;

		delete m;
		delete lcr;
		delete s;
	
		m = degreeGPE(r, x);
	}

	AST* res = list({ Tnn(q, x, p), Tnn(r, x, p) });
	
	delete q;
	delete r;

	delete m;
	delete n;
	delete lcv;

	return res;
}


AST* divideGPE_Sp(AST* u, AST* v, AST* x, int p) {

	AST* q = integer(0);
	AST* r = u->deepCopy();

	AST* m = degreeGPE(r, x);
	AST* n = degreeGPE(v, x);

	AST* lcv = leadingCoefficientGPE(v, x);

	while(
		m->kind() != Kind::MinusInfinity &&
		(m->kind() == Kind::Integer && n->kind() == Kind::Integer &&
		m->value() >= n->value())
	) {
	
		AST* lcr = leadingCoefficientGPE(r, x);
	
		AST* s = integer(division_Sp(lcr->value(), lcv->value(), p));

		AST* q_ = add({
			q->deepCopy(),
			mul({
				s->deepCopy(),
				power(
					x->deepCopy(),
					sub({
						m->deepCopy(),
						n->deepCopy()
					})
				)
			})
		});
		
		
		delete q;
	
		q = Ts(q_,x,p);
	
		delete q_;

		AST* r_ = sub({
			sub({
				r->deepCopy(),
				mul({
					lcr->deepCopy(),
					power(x->deepCopy(), m->deepCopy())
				})
			}),
			mul({
				sub({
					v->deepCopy(),
					mul({
						lcv->deepCopy(),
						power(x->deepCopy(), n->deepCopy())
					}),
				}),
				s->deepCopy(),
				power(
					x->deepCopy(),
					sub({m->deepCopy(), n->deepCopy()})
				)
			})
		});

		delete r;

		r = Ts(r_, x, p);

		delete r_, m, lcr, s;
	
		m = degreeGPE(r, x);
	}

	AST* res = list({Ts(q,x,p), Ts(r,x,p)});
	
	delete q;
	delete r;
	delete m;
	delete n;
	delete lcv;

	return res;
}

AST* remainderGPE_Zp(AST* u, AST* v, AST* x, int p) {
	AST* res = divideGPE_Zp(u,v,x,p);
	AST* r = res->operand(1)->deepCopy();
	delete res;
	return r;
}

AST* quotientGPE_Zp(AST* u, AST* v, AST* x, int p) {
	AST* res = divideGPE_Zp(u,v,x,p);
	AST* q = res->operand(0)->deepCopy();
	delete res;
	return q;
}

AST* remainderGPE_Sp(AST* u, AST* v, AST* x, int p) {
	AST* res = divideGPE_Sp(u,v,x,p);
	AST* r = res->operand(1)->deepCopy();
	delete res;
	return r;
}

AST* quotientGPE_Sp(AST* u, AST* v, AST* x, int p) {
	AST* res = divideGPE_Sp(u,v,x,p);
	AST* q = res->operand(0)->deepCopy();
	delete res;
	return q;
}


AST* gcdGPE_Zp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return integer(0);
	}

	AST* U = u->deepCopy();
	AST* V = v->deepCopy();

	while (V->kind() != Kind::Integer ||(V->kind() == Kind::Integer && V->value() != 0)) {
		AST* R = remainderGPE_Zp(U, V, x, p);
		delete U;
		U = V->deepCopy();
		delete V;
		V = R->deepCopy();
		delete R;
	}

	AST* lco = leadingCoefficientGPE(U,x);
	
	AST* e = mul({ integer(division_Zp(1, lco->value(), p)), U->deepCopy() });
	
	delete lco;

	AST* res = Tnn(e, x, p);

	delete U;
	delete V;
	delete e;

	return res;
}


AST* gcdGPE_Sp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return integer(0);
	}

	AST* U = u->deepCopy();
	AST* V = v->deepCopy();

	while (V->kind() != Kind::Integer ||(V->kind() == Kind::Integer && V->value() != 0)) {
		AST* R = remainderGPE_Sp(U, V, x, p);
		delete U;
		U = V->deepCopy();
		delete V;
		V = R->deepCopy();
		delete R;
	}

	AST* lco = leadingCoefficientGPE(U,x);
	
	AST* e = mul({ integer(division_Sp(1, lco->value(), p)), U->deepCopy() });

	AST* res = Ts(e, x, p);

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

	AST* U 		= u->deepCopy();
	AST* V 		= v->deepCopy();
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

		AST* A_ = sub({ App->deepCopy(), mul({q->deepCopy(), Ap->deepCopy()}) });
		AST* B_ = sub({ Bpp->deepCopy(), mul({q->deepCopy(), Bp->deepCopy()}) });

		AST* A = Tnn(A_,x, p);
		AST* B = Tnn(B_,x, p);
	
		delete A_;
		delete B_;
	
		delete App;
		App = Ap->deepCopy();

		delete Ap;
		Ap 	= A->deepCopy();

		delete Bpp;
		Bpp = Bp->deepCopy();

		delete Bp;
		Bp 	= B->deepCopy();

		delete U;
		U = V->deepCopy();

		delete V;
		V = r->deepCopy();

		delete A;
		delete B;
	
		delete d;
	}

	AST* c = leadingCoefficientGPE(U, x);

	AST* App__ = mul({ App->deepCopy(), integer(modInverse_p(c->value(), p)) });
	AST* App_ = Tnn(App__, x, p);
	delete App;
	delete App__;
	App = App_;

	AST* Bpp__ = mul({ Bpp->deepCopy(), integer(modInverse_p(c->value(), p)) });
	AST* Bpp_ = Tnn(Bpp__, x, p);
	delete Bpp;
	delete Bpp__;
	Bpp = Bpp_;
	
	AST* U__ = mul({U->deepCopy(), integer(modInverse_p(c->value(), p))});
	AST* U_ = Tnn(U__, x, p);
	delete U;
	delete U__;
	U = U_;

	delete Ap;
	delete Bp;
	delete V;
	delete c;
	
	return list({ U, App, Bpp });
}


ast::AST* extendedEuclideanAlgGPE_Sp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return list({ integer(0), integer(0), integer(0) });
	}

	AST* U 		= u->deepCopy();
	AST* V 		= v->deepCopy();
	AST* App 	= integer(1);
	AST* Ap 	= integer(0);
	AST* Bpp 	= integer(0);
	AST* Bp 	= integer(1);

	while (
		V->kind() != Kind::Integer ||
		(V->kind() == Kind::Integer && V->value() != 0)
	) {
		AST* d = divideGPE_Sp(U,V,x,p);
	
		AST* q = d->operand(0);
		AST* r = d->operand(1);

		AST* A_ = sub({ App->deepCopy(), mul({q->deepCopy(), Ap->deepCopy()}) });
		AST* B_ = sub({ Bpp->deepCopy(), mul({q->deepCopy(), Bp->deepCopy()}) });
		
		AST* A = Ts(A_,x, p);
		AST* B = Ts(B_,x, p);
	
		delete A_;
		delete B_;
	
		delete App;
		App = Ap->deepCopy();

		delete Ap;
		Ap 	= A->deepCopy();

		delete Bpp;
		Bpp = Bp->deepCopy();

		delete Bp;
		Bp 	= B->deepCopy();

		delete A;
		delete B;

		delete U;
		U = V->deepCopy();

		delete V;
		V = r->deepCopy();

		delete d;
	}

	AST* c = leadingCoefficientGPE(U, x);

	AST* App__ = mul({ App->deepCopy(), integer(modInverse_p(c->value(), p)) });
	AST* App_ = Ts(App__, x, p);
	delete App;
	delete App__;
	App = App_;

	AST* Bpp__ = mul({ Bpp->deepCopy(), integer(modInverse_p(c->value(), p)) });
	AST* Bpp_ = Ts(Bpp__, x, p);
	delete Bpp;
	delete Bpp__;
	Bpp = Bpp_;
	
	AST* U__ = mul({U->deepCopy(), integer(modInverse_p(c->value(), p))});
	AST* U_ = Ts(U__, x, p);
	delete U;
	delete U__;
	U = U_;

	delete Ap;
	delete Bp;
	delete V;
	delete c;
	
	return list({ U, App, Bpp });
}




}
