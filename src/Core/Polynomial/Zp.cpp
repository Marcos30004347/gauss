
/**
 * NOTTHING HERE WAS TESTED YET, EVERYTHING COULD BE WRONG!!!!!!
 */

#include "Zp.hpp"
#include "Core/Expand/Expand.hpp"
using namespace ast;
using namespace expand;
using namespace algebra;

namespace polynomial {

int mod(int a, int b) {
	return (b + (a%b)) % b;
}

int power(int x, unsigned int y, unsigned int m) {
	if (y == 0)
		return 1;
	int p = (power(x, y / 2, m) % m);
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
		return power(a, p - 2, p);
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
		AST* e = pow(x->deepCopy(), inte(i));
		AST* c_ = coefficientGPE(u_, e);
		AST* c = expandAST(c_);
		
		Tnn_u->includeOperand(mul({inte(mod(c->value(), s)), e}));
		
		delete c, c_;
	}


	AST* r = expandAST(Tnn_u);

	delete u_, Tnn_u;
	return r;
}

AST* Ts(AST* u, AST* x, int s) {
	AST* u_ = expandAST(u);

	AST* Tnn_u = new AST(Kind::Addition);

	AST* d = degreeGPE(u_, x);

	for(int i=0; i<=d->value(); i++) {

		AST* e = pow(x->deepCopy(), inte(i));
		AST* c_ = coefficientGPE(u_, e);
		AST* c = expandAST(c_);

		Tnn_u->includeOperand(mul({inte(S(mod(c->value(), s), s)), e}));
		
		delete c, c_;
	}


	AST* r = expandAST(Tnn_u);
	delete u_, Tnn_u;

	return r;
}

std::pair<AST*, AST*> divideGPE_Zp(AST* u, AST* v, AST* x, int p) {

	AST* q = inte(0);
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
	
		AST* s = inte(division_Zp(lcr->value(), lcv->value(), p));

		AST* q_ = add({
			q->deepCopy(),
			mul({
				s->deepCopy(),
				pow(
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
					pow(x->deepCopy(), m->deepCopy())
				})
			}),
			mul({
				sub({
					v->deepCopy(),
					mul({
						lcv->deepCopy(),
						pow(x->deepCopy(), n->deepCopy())
					}),
				}),
				s->deepCopy(),
				pow(
					x->deepCopy(),
					sub({m->deepCopy(), n->deepCopy()})
				)
			})
		});

		delete r;

		r = Tnn(r_, x, p);
		
		delete r_, m, lcr, s;
	
		m = degreeGPE(r, x);
	
	}

	std::pair<AST*, AST*> res = { Tnn(q,x,p), Tnn(r,x,p) };
	
	delete q, r, m, n, lcv;

	return res;
}


std::pair<AST*, AST*> divideGPE_Sp(AST* u, AST* v, AST* x, int p) {

	AST* q = inte(0);
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
	
		AST* s = inte(division_Sp(lcr->value(), lcv->value(), p));

		AST* q_ = add({
			q->deepCopy(),
			mul({
				s->deepCopy(),
				pow(
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
					pow(x->deepCopy(), m->deepCopy())
				})
			}),
			mul({
				sub({
					v->deepCopy(),
					mul({
						lcv->deepCopy(),
						pow(x->deepCopy(), n->deepCopy())
					}),
				}),
				s->deepCopy(),
				pow(
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

	std::pair<AST*, AST*> res = {Ts(q,x,p), Ts(r,x,p)};
	
	delete q, r, m, n, lcv;

	return res;
}

AST* remainderGPE_Zp(AST* u, AST* v, AST* x, int p) {
	std::pair<AST*, AST*> res = divideGPE_Zp(u,v,x,p);
	delete res.first;
	return res.second;
}

AST* quotientGPE_Zp(AST* u, AST* v, AST* x, int p) {
	std::pair<AST*, AST*> res = divideGPE_Zp(u,v,x,p);
	delete res.second;
	return res.first;
}

AST* remainderGPE_Sp(AST* u, AST* v, AST* x, int p) {
	std::pair<AST*, AST*> res = divideGPE_Sp(u,v,x,p);
	delete res.first;
	return res.second;
}

AST* quotientGPE_Sp(AST* u, AST* v, AST* x, int p) {
	std::pair<AST*, AST*> res = divideGPE_Sp(u,v,x,p);
	delete res.second;
	return res.first;
}


AST* gcdGPE_Zp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return inte(0);
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
	
	AST* e = mul({ inte(division_Zp(1, lco->value(), p)), U->deepCopy() });
	
	delete lco;

	AST* res = Tnn(e, x, p);

	delete U, V, e;

	return res;
}


AST* gcdGPE_Sp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return inte(0);
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
	
	AST* e = mul({ inte(division_Sp(1, lco->value(), p)), U->deepCopy() });

	AST* res = Ts(e, x, p);

	delete U, V, e, lco;

	return res;

}


std::vector<AST*> extendedEuclideanAlgGPE_Zp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return { inte(0), inte(0), inte(0) };
	}

	AST* U 		= u->deepCopy();
	AST* V 		= v->deepCopy();
	AST* App 	= inte(1);
	AST* Ap 	= inte(0);
	AST* Bpp 	= inte(0);
	AST* Bp 	= inte(1);

	while (
		V->kind() != Kind::Integer ||
		(V->kind() == Kind::Integer && V->value() != 0)
	) {
		std::pair<AST*, AST*> d = divideGPE_Zp(U,V,x,p);
	
		AST* q = d.first;
		AST* r = d.second;

		AST* A_ = sub({ App->deepCopy(), mul({q->deepCopy(), Ap->deepCopy()}) });
		AST* B_ = sub({ Bpp->deepCopy(), mul({q->deepCopy(), Bp->deepCopy()}) });
		
		AST* A = Tnn(A_,x, p);
		AST* B = Tnn(B_,x, p);
	
		delete A_, B_;
	
		delete App;
		App = Ap->deepCopy();

		delete Ap;
		Ap 	= A->deepCopy();

		delete Bpp;
		Bpp = Bp->deepCopy();

		delete Bp;
		Bp 	= B->deepCopy();

		delete A, B;

		delete U;
		U = V->deepCopy();

		delete V;
		V = r->deepCopy();

		delete q, r;
	}

	AST* c = leadingCoefficientGPE(U, x);

	AST* App__ = mul({ App->deepCopy(), inte(modInverse_p(c->value(), p)) });
	AST* App_ = Tnn(App__, x, p);
	delete App, App__;
	App = App_;

	AST* Bpp__ = mul({ Bpp->deepCopy(), inte(modInverse_p(c->value(), p)) });
	AST* Bpp_ = Tnn(Bpp__, x, p);
	delete Bpp, Bpp__;
	Bpp = Bpp_;
	
	AST* U__ = mul({U->deepCopy(), inte(modInverse_p(c->value(), p))});
	AST* U_ = Tnn(U__, x, p);
	delete U, U__;
	U = U_;

	delete Ap, Bp, V, c;
	
	return { U, App, Bpp };
}


std::vector<ast::AST*> extendedEuclideanAlgGPE_Sp(AST* u, AST* v, AST* x, int p) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return { inte(0), inte(0), inte(0) };
	}

	AST* U 		= u->deepCopy();
	AST* V 		= v->deepCopy();
	AST* App 	= inte(1);
	AST* Ap 	= inte(0);
	AST* Bpp 	= inte(0);
	AST* Bp 	= inte(1);

	while (
		V->kind() != Kind::Integer ||
		(V->kind() == Kind::Integer && V->value() != 0)
	) {
		std::pair<AST*, AST*> d = divideGPE_Sp(U,V,x,p);
	
		AST* q = d.first;
		AST* r = d.second;

		AST* A_ = sub({ App->deepCopy(), mul({q->deepCopy(), Ap->deepCopy()}) });
		AST* B_ = sub({ Bpp->deepCopy(), mul({q->deepCopy(), Bp->deepCopy()}) });
		
		AST* A = Ts(A_,x, p);
		AST* B = Ts(B_,x, p);
	
		delete A_, B_;
	
		delete App;
		App = Ap->deepCopy();

		delete Ap;
		Ap 	= A->deepCopy();

		delete Bpp;
		Bpp = Bp->deepCopy();

		delete Bp;
		Bp 	= B->deepCopy();

		delete A, B;

		delete U;
		U = V->deepCopy();

		delete V;
		V = r->deepCopy();

		delete q, r;
	}

	AST* c = leadingCoefficientGPE(U, x);

	AST* App__ = mul({ App->deepCopy(), inte(modInverse_p(c->value(), p)) });
	AST* App_ = Ts(App__, x, p);
	delete App, App__;
	App = App_;

	AST* Bpp__ = mul({ Bpp->deepCopy(), inte(modInverse_p(c->value(), p)) });
	AST* Bpp_ = Ts(Bpp__, x, p);
	delete Bpp, Bpp__;
	Bpp = Bpp_;
	
	AST* U__ = mul({U->deepCopy(), inte(modInverse_p(c->value(), p))});
	AST* U_ = Ts(U__, x, p);
	delete U, U__;
	U = U_;

	delete Ap, Bp, V, c;
	
	return { U, App, Bpp };
}



}
