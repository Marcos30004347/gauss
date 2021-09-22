#include "Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Calculus/Calculus.hpp"

#include <numeric>

using namespace ast;
using namespace expand;
using namespace simplification;
using namespace algebra;
using namespace calculus;


namespace polynomial {

void includeVariable(std::vector<AST*>& vars, AST* u) {
	bool included = false;

	for(AST* k : vars) {
		if(k->match(u)){
			included = true;
			break;
		}
	}

	if(!included) {
		vars.push_back(u->copy());
	}
}

bool isGeneralMonomial(AST* u, AST* v) {
	AST* S;
	if(v->kind() != Kind::Set) {
		S = set({v->copy()});
	} else {
		S = v->copy();
	}

	if(exists(S, u)) {
		delete S;
		return true;
	} else if(u->kind() == Kind::Power){
		AST* b = u->operand(0);
		AST* e = u->operand(1);

		if(exists(S, b) && e->kind() == Kind::Integer && e->value() > 1) {
			delete S;
			return true;
		}
	} else if(u->kind() == Kind::Multiplication) {
		for(unsigned int i=0; i<u->numberOfOperands(); i++) {
			if(isGeneralMonomial(u->operand(i), S) == false) {
				delete S;
				return false;
			}
		}
		delete S;
		return true;
	}

	bool r = u->freeOfElementsInSet(S);
	
	delete S;
	return r;
}

bool isGerenalPolynomial(AST* u, AST* v) {
	AST* S;

	if(v->kind() != Kind::Set) {
		S = set({ v->copy() });
	} else {
		S = v->copy();
	}

	if(u->kind() != Kind::Addition && u->kind() != Kind::Subtraction) {
		bool r = isGeneralMonomial(u, S);
		delete S;
		return r;
	}

	if(exists(S, u)) {
		delete S;
		return true;
	}

	for(unsigned int i=0; i<u->numberOfOperands(); i++) {
		if(isGeneralMonomial(u->operand(i), S) == false) {
			delete S;
			return false;
		}
	}

	delete S;
	return true;
}

AST* coeffVarMonomial(AST* u, AST* S) {
	if(!isGeneralMonomial(u, S))
		return undefined();

	if(isConstant(u))
		return list({ u->copy(), integer(1) });

	if(exists(S, u))
		return list({ integer(1), u->copy() });
	
	if(u->kind() == Kind::Power && exists(S, u->operand(0)))
		return list({ integer(1), u->copy() });

	if(u->kind() == Kind::Multiplication) {
		AST* C = list({});
		AST* V = list({});
	
		for(unsigned int i=0; i<u->numberOfOperands(); i++) {
			AST* L = coeffVarMonomial(u->operand(i), S);
		
			AST* CL = list({L->operand(0)->copy()});
			AST* VL = list({L->operand(1)->copy()});
	
			AST* C_ = join(C, CL);
			AST* V_ = join(V, VL);
			
			delete C;
			delete V;
			
			C = C_;
			V = V_;
		
			delete L;
			delete CL;
			delete VL;
		}

		AST* coefs = mul({});
		AST* vars  = mul({});
		
		for(unsigned int i=0; i<C->numberOfOperands(); i++) {
			if(
				C->operand(i)->kind() == Kind::Integer &&
				C->operand(i)->value() == 1
			) continue;
				
				coefs->includeOperand(C->operand(i)->copy());
		}
		
		for(unsigned int i=0; i<V->numberOfOperands(); i++) {
			if(
				V->operand(i)->kind() == Kind::Integer &&
				V->operand(i)->value() == 1
			) continue;
			vars->includeOperand(V->operand(i)->copy());
		}
		
		delete C;
		delete V;
	
		if(coefs->numberOfOperands() == 0) {
			delete coefs;
			coefs = integer(1);
		} else if(coefs->numberOfOperands() == 1) {
			AST* coefs_ = coefs->operand(0)->copy();
			delete coefs;
			coefs = coefs_;
		}
	
		if(vars->numberOfOperands() == 0) {
			delete vars;
			vars = integer(1);
		} else if(vars->numberOfOperands() == 1) {
			AST* vars_ = vars->operand(0)->copy();
			delete vars;
			vars = vars_;
		}
	
		return list({ coefs, vars });
	}

	return list({ u->copy(), integer(1) });
}

AST* collectTerms(AST* u, AST* S) {
	if(u->kind() != Kind::Addition) {
		AST* L = coeffVarMonomial(u, S);
		if(L->kind() == Kind::Undefined) {
			delete L;
			return undefined();
		}

		return u->copy();
	}

	if(exists(S, u)) {
		return u->copy();
	}

	int N = 0;

	AST* T = list({});

	for(unsigned int i=0; i<u->numberOfOperands(); i++) {
		AST* f = coeffVarMonomial(u->operand(i), S);
		
		if(f->kind() == Kind::Undefined) {
			delete f;
			return undefined();
		}
	
		int j = 1;
		bool combined = false;
		
		while(!combined && j <= N) {
			int j_ = j - 1;

			if(f->operand(1)->match(T->operand(j_)->operand(1))) {
				
				AST* Tj = list({
					add({ 
						T->operand(j_)->operand(0)->copy(),
						f->operand(0)->copy()
					}),
					f->operand(1)->copy()
				}); 

				AST* Tj_ = T->operand(j_);
				T->removeOperand(j_);
				
				delete Tj_;
				
				T->includeOperand(Tj, j_);
				
				combined = true;
			}

			j = j+1;
		}

		if(!combined) {
			T->includeOperand(f->copy(), N);
			N = N + 1;
		}
	
		delete f;
	}

	AST* v = add({});

	for(int j=0; j<N; j++) {
		if(
			T->operand(j)->operand(1)->kind() == Kind::Integer &&
			T->operand(j)->operand(1)->value() == 1
		) {
			v->includeOperand(T->operand(j)->operand(0)->copy());
		} else {
			v->includeOperand(mul({
				T->operand(j)->operand(0)->copy(),
				T->operand(j)->operand(1)->copy(),
			}));
		}
	}

	delete T;

	if(v->numberOfOperands() == 0) {
		delete v;
		return integer(0);
	}

	if(v->numberOfOperands() == 1) {
		AST* v_ = v->operand(0)->copy();
		delete v;
		v = v_;
	}

	return v;
}

AST* degreeGME(AST* u, AST* v) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return new AST(Kind::MinusInfinity);
	
	if(isConstant(u))
		return integer(0);

	AST* S;
	if(v->kind() != Kind::Set) {
		S = set({ v->copy() });
	} else {
		S = v->copy();
	}

	if(exists(S, u)) {

		delete S;
		return integer(1);
	} else if(u->kind() == Kind::Power){
		AST* b = u->operand(0);
		AST* e = u->operand(1);

		if(exists(S, b) && isConstant(e)) {
			delete S;
			return e->copy();
		}

	} else if(u->kind() == Kind::Multiplication) {
		AST* deg = integer(0);
		for(unsigned int i=0; i<u->numberOfOperands(); i++) {
			AST* deg_ = degreeGME(u->operand(i), S);
			if(deg_->value() > deg->value()) {
				delete deg;
				deg = deg_;
			} else {
				delete deg_;
			}
		}
		delete S;
		return deg;
	}
	
	delete S;
	return integer(0);
}

AST* degreeGPE(AST* u, AST* v) {
	AST* S;

	if(u->kind() == Kind::Integer && u->value() == 0) {
		return new AST(Kind::MinusInfinity);
	}
	
	if(v->kind() != Kind::Set) {
		S = set({v->copy()});
	} else {
		S = v->copy();
	}

	if(u->kind() != Kind::Addition && u->kind() != Kind::Subtraction) {
		AST* r = degreeGME(u, S);
		delete S;
		return r;
	}

	if(exists(S, u)) {
		delete S;
		return integer(1);
	}

	AST* deg = integer(0);

	for(unsigned int i=0; i<u->numberOfOperands(); i++) {
		AST* deg_ = degreeGME(u->operand(i), S);

		if(deg_->value() > deg->value()) {
			delete deg;
			deg = deg_;
		} else {
			delete deg_;
		}
	}

	delete S;

	return deg;
}

AST* variables(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction
	)	return set({});

	if(u->kind() == Kind::Power) {
		AST* b = u->operand(0);
		AST* e = u->operand(1);

		if(e->kind() == Kind::Integer && e->value() > 1)
			return set({ b->copy() });

		return set({ u->copy() });
	}

	if(
		u->kind() == Kind::Addition ||
		u->kind() == Kind::Subtraction ||
		u->kind() == Kind::Multiplication
	) {
		AST* S = set({});

		for(unsigned int i=0; i<u->numberOfOperands(); i++) {
			AST* S_ = variables(u->operand(i));
			AST* S__ = unification(S, S_);
			delete S;
			delete S_;
			S = S__;
		}

		return S;
	}

	return set({ u->copy() });
}

AST* coefficientGME(AST* u, AST* x) 
{
	if(u->match(x)) 
	{
		return list({ integer(1), integer(1) });
	}

	if(u->kind() == Kind::Power) 
	{
		AST* b = u->operand(0);
		AST* e = u->operand(1);

		if(b->match(x) && e->kind() == Kind::Integer && e->value() > 0) 
		{
			return list({ integer(1), e->copy() });
		}

	} 
	else if(u->kind() == Kind::Multiplication) 
	{
		AST* m = integer(0);
		AST* c = u->copy();

		for(unsigned int i=0; i<u->numberOfOperands(); i++) 
		{
			AST* f =	coefficientGME(u->operand(i), x);
			
			if(f->kind() == Kind::Undefined) 
			{
				delete m;
				delete c;
				delete f;
		
				return undefined();
			} 
			if(
				f->operand(1)->kind() != Kind::Integer ||
				f->operand(1)->value() != 0
			) 
			{
				delete m;
				delete c;
	
				m = f->operand(1)->copy();
				AST* c_ = div(u->copy(), power(x->copy(), m->copy()));
				c = algebraicExpand(c_);
				delete c_;
			}
	
			delete f;
		}

		return list({ c, m });
	}

	if(u->freeOf(x)) 
	{
		return list({ u->copy(), integer(0) });
	}

	return undefined();
}

AST* coefficientGPE(AST* u, AST* x, AST* j) 
{
	
	if(u->kind() != Kind::Addition && u->kind() != Kind::Subtraction) {
		AST* f = coefficientGME(u, x);
		
		if(f->kind() == Kind::Undefined) return f;
		
		if(j->match(f->operand(1))) {
			AST* k = f->operand(0)->copy();
			
			delete f;

			return k;
		}

		delete f;

		return integer(0);
	}

	if(x->match(u)) {
		if(j->kind() == Kind::Integer && j->value() == 1) {
			return integer(1);
		}

		return integer(0);
	}

	AST* c = integer(0);

	for(unsigned int i=0; i<u->numberOfOperands(); i++) 
	{
		AST* f = coefficientGME(u->operand(i), x);

		if(f->kind() == Kind::Undefined) return f;
		
		if(j->match(f->operand(1))) {
			AST* k = f->operand(0)->copy();

			if(c->kind() == Kind::Integer && c->value() == 0) {
				delete c;
				c = new AST(u->kind());
				c->includeOperand(k);
			} else {
				c->includeOperand(k);
			}
		}

		delete f;
	}

	if(c->kind()!= Kind::Integer && c->numberOfOperands() == 1) {
		AST* l = c->operand(0)->copy();
		delete c;
		return l;
	}
	
	return c;
}

AST* leadingCoefficientGPE(AST* u, AST* x) {
	// assert(
	// 	!isConstant(x),
	// 	"leadingCoefficientGPE: 'param(x)=%s' "
	// 	"cant be a constant expression",
	// 	x->toString().c_str()
	// );
	AST* d = degreeGPE(u, x);

	AST* lc = coefficientGPE(u, x, d);

	AST* r = algebraicExpand(lc);

	delete d;
	delete lc;

	return r;
}

AST* divideGPE(AST* u, AST* v, AST* x) {
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

		AST* s = div(lcr->copy(), lcv->copy());

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

	AST* res = list({ algebraicExpand(q), algebraicExpand(r) });
	
	delete q;
	delete r;
	delete m;
	delete n;
	delete lcv;

	return res;
}

AST* quotientGPE(AST* u, AST* v, AST* x) {
	AST* res = divideGPE(u,v,x);
	AST* r = res->operand(0)->copy();
	delete res;
	return r;
}

AST* remainderGPE(AST* u, AST* v, AST* x) {
	AST* res = divideGPE(u,v,x);
	AST* r = res->operand(1)->copy();
	delete res;
	return r;
}

AST* expandGPE(AST* u, AST* v, AST* x, AST* t) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return integer(0);

	AST* d = divideGPE(u, v, x);

	AST* q = d->operand(0);
	AST* r = d->operand(1);

	AST* expoent = add({
		mul({
			t->copy(),
			expandGPE(q, v, x, t)
		}),
		r->copy()
	});

	AST* res = algebraicExpand(expoent);

	delete expoent;
	delete d;

	return res;
}


AST* gcdGPE(AST* u, AST* v, AST* x) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return integer(0);
	}

	AST* U = u->copy();
	AST* V = v->copy();

	while (V->kind() != Kind::Integer || (V->kind() == Kind::Integer && V->value() != 0)) {
		AST* R = remainderGPE(U, V, x);

		delete U;

		U = V;
		V = R;
	}

	AST* e = mul({div(integer(1), leadingCoefficientGPE(U,x)), U});
	AST* res = algebraicExpand(e);

	delete e;
	delete V;

	return res;
}

AST* extendedEuclideanAlgGPE(AST* u, AST* v, AST* x) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return list({ integer(0), integer(0), integer(0) });
	}

	AST* U 		= u->copy();
	AST* V 		= v->copy();

	AST* App 	= integer(1), *Ap = integer(0), *Bpp = integer(0), *Bp = integer(1);

	while (V->kind() != Kind::Integer || V->value() != 0) {
		AST* d = divideGPE(U,V,x);

		AST* q = d->operand(0);
		AST* r = d->operand(1);

		AST* A_ = sub({
			App->copy(),
			mul({
				q->copy(),
				Ap->copy()
			})
		});

		AST* B_ = sub({
			Bpp->copy(),
			mul({
				q->copy(),
				Bp->copy()
			})
		});

		AST* A = algebraicExpand(A_);
		AST* B = algebraicExpand(B_);
		
		// printf("A %s\n", A->toString().c_str());
		// printf("B %s\n", B->toString().c_str());
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

	AST* App__ = div(App->copy(), c->copy());
	AST* App_ = algebraicExpand(App__);
	delete App__;

	delete App;
	App = App_;

	AST* Bpp__ = div(Bpp->copy(), c->copy());
	AST* Bpp_ = algebraicExpand(Bpp__);
	delete Bpp__;

	delete Bpp;
	Bpp = Bpp_;

	AST* U__ = div(U->copy(), c->copy());
	AST* U_ = algebraicExpand(U__);
	delete U__;
	
	delete U;
	U = U_;

	delete c;
	delete Ap;
	delete Bp;
	delete V;

	return list({ U, App, Bpp });
}

AST* algMulInverseAST(AST* v, AST* p, AST* a) {
	// TODO: assert following statements
	// 1. a is a symbol that represents an algebraic number
	// 2. p is a monoic irrecudctible polynomial in Q[a] with degree(p,a) >= 2
	// 3. v is a non zero polynomial in Q(a) with degree(v) < degree(p)
	

	AST* w = extendedEuclideanAlgGPE(v,p,a);
	AST* r = w->operand(1)->copy();
	delete w;
	return r;
}

AST* algDivideAST(AST* u, AST* v, AST* p, AST* a) {
	// TODO: assert following statements
	// a is a symbol that represents an algebraic number;
	// p is a monic, irrreducible polynomial in Q[α] with deg(p, α) ≥ 2;
	// u and v are both polynomials in Q(a) with degree < deg(p) and v != 0;
	AST* w = algMulInverseAST(v, p, a);

	AST* e = mul({u->copy(), w->copy()});

	AST* k = algebraicExpand(e);

	AST* r = remainderGPE(k, p, a);

	delete w;
	delete e;
	delete k;

	return r;
}

AST* algCoeffSimp(AST* u, AST* x, AST* p, AST* a) {
	// assert(
	// 	x->kind() == Kind::Symbol,
	// 	"algCoeffSimp: 'param(x)=%s' needs to be a symbol",
	// 	x->toString().c_str()
	// );


	AST* d = degreeGPE(u, x);

	if(d->value() == 0) {
		delete d;
		return remainderGPE(u, p, a);
	}

	AST* r = new AST(u->kind());

	for(int i=0; i <= d->value(); i++) {
		AST* d = integer(i);

		AST* coeff_ = coefficientGPE(u, x, d);
		AST* coeff = algebraicExpand(coeff_);
		AST* k = remainderGPE(coeff, p, a);
		
		delete coeff_;
		delete coeff;
		
		r->includeOperand(mul({k, power(x->copy(), d)}));
	}
	
	AST* res = algebraicExpand(r);
	
	delete d;
	delete r;
	
	return res;
}

std::vector<AST*> algPolynomialDivisionAST(AST* u, AST* v, AST* x, AST* p, AST* a) {
	// TODO: assert following statements
	// u, v : polynomials in Q(a)[x] with v != 0;
	// x : a symbol;
	// a : a symbol that represents an algebraic number;
	// p : a monic, irrreducible polynomial in Q[a] with degree ≥ 2;

	AST* q = integer(0);
	AST* r = u->copy();
	AST* m = degreeGPE(r, x);
	AST* n = degreeGPE(v, x);
	AST* lcv = leadingCoefficientGPE(v, x);

	
	AST* p_ = deepReplace(p, x, a);


	while(
		m->kind() != Kind::MinusInfinity && (
			m->kind() == Kind::Integer && n->kind() == Kind::Integer &&
			m->value() >= n->value()
		)
	) {
	
		AST* lcr = leadingCoefficientGPE(r, x);

		AST* s = algDivideAST(lcr, lcv, p_, a);
		
		AST* q_ = add({
			q->copy(),
			mul({
				s->copy(),
				power(x->copy(),
				sub({m->copy(), n->copy()}))
			})
		});

	
		delete q;
	
		q = algebraicExpand(q_);
	
		delete q_;

		AST* e = sub({
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
					})
				}),
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
	
		AST* r_ = algebraicExpand(e);
	
		delete e;
		delete r;
	
		r = algCoeffSimp(r_, x, p_, a);
	
		delete r_;
		
		delete m;

		m = degreeGPE(r, x);
		
		delete s;
		delete lcr;		
	}

	delete m;
	delete n;
	delete lcv;
	delete p_;

	return { q, r };
}

AST* algPolynomialRemainderAST(AST* u, AST* v, AST* x, AST* p, AST* a) {
	std::vector<AST*> res = algPolynomialDivisionAST(u,v,x,p,a);
	delete res[0];
	return res[1];
}
AST* algPolynomialQuotientAST(AST* u, AST* v, AST* x, AST* p, AST* a) {
	std::vector<AST*> res = algPolynomialDivisionAST(u,v,x,p,a);
	delete res[1];
	return res[0];
}

AST* algPolynomialGCDAST(AST* u, AST* v, AST* x, AST* p, AST* a) {
	AST* U = u->copy();
	AST* V = v->copy();

	while(
		V->kind() != Kind::Integer ||
		V->value() != 0
	) {
		AST* R = algPolynomialRemainderAST(U, V, x, p, a);

		delete U;

		U = V->copy();

		delete V;

		V = R->copy();

		delete R;
	}

	AST* r = algMonicAST(U, x, p, a);	
	
	delete U;
	delete V;
	
	return r;
}

AST* algMonicAST(AST* u,AST* x, AST* p,AST* a) {
	AST* lc = leadingCoefficientGPE(u, x);
	AST* k_ = algPolynomialQuotientAST(u, lc, x, p, a);
	delete lc;
	return k_;
}


AST* recPolyDiv(AST* u, AST* v, AST* L, AST* K) {
	if(L->numberOfOperands() == 0) {
		AST* d_ = div(u->copy(), v->copy());
		AST* d = algebraicExpand(d_);
		delete d_;
		
		if(K->identifier() == "Z") {
			if(d->kind() == Kind::Integer) {
				return list({ d, integer(0) });
			}
			
			delete d;
			
			return list({ integer(0), u->copy() });
		}

		if(K->identifier() == "Q") {
			AST* L_ = list({ d, integer(0) });
			return L_;
		}

		// Undefined integral domain
		return undefined();
	}

	AST* x = first(L);
	AST* r = u->copy();

	AST* m = degreeGPE(r, x);
	AST* n = degreeGPE(v, x);
	
	AST* q = integer(0);
	
	AST* lcv = leadingCoefficientGPE(v, x);

	while(
		m->kind() != Kind::MinusInfinity &&
		m->value() >= n->value()
	) {
		AST* lcr = leadingCoefficientGPE(r, x);
		AST* restL = rest(L);
		
		AST* d = recPolyDiv(lcr, lcv, restL, K);
		
		delete restL;
		
		if(d->operand(1)->kind() != Kind::Integer || d->operand(1)->value() != 0) {
			AST* result = algebraicExpand(q);
			
			delete x;
			delete m;
			delete n;
			delete q;
			delete d;
			delete lcv;
			delete lcr;
			
			return list({ result, r });
		} else {
			AST* c = d->operand(0)->copy();
			
			AST* q_ = add({
				q->copy(),
				mul({
					c->copy(),
					power(
						x->copy(),
						sub({ m->copy(), n->copy() })
					)
				})
			});

			delete q;

			q = algebraicExpand(q_);

			delete q_;
			
			AST* r_ = sub({
				r->copy(),
				mul({
					c->copy(),
					v->copy(),
					power(
						x->copy(),
						sub({ m->copy(), n->copy() })
					)
				})
			});

			delete r;
			r = algebraicExpand(r_);
			delete r_;

			delete m;


			m = degreeGPE(r, x);
		
			delete c;
		}
	
		delete lcr;
		delete d;
	}

	AST* result = algebraicExpand(q);

	delete x;
	delete m;
	delete n;
	delete q;
	delete lcv;

	return list({ result, r });
}

AST* recQuotient(AST* u, AST* v, AST* L, AST* K) {
	AST* r = recPolyDiv(u, v, L, K);
	AST* q = r->operand(0)->copy();
	delete r;
	return q;
}

AST* recRemainder(AST* u, AST* v, AST* L, AST* K) {
	AST* r = recPolyDiv(u, v, L, K);
	AST* q = r->operand(1)->copy();
	delete r;
	return q;
}


AST* pseudoDivision(AST* u, AST* v, AST* x) {
	AST* p = integer(0);
	AST* s = u->copy();
	AST* m = degreeGPE(s, x);
	AST* n = degreeGPE(v, x);


	AST* e_ = add({
		sub({
			m->copy(),
			n->copy()
		}),
		integer(1)
	});
	AST* e = algebraicExpand(e_);
	delete e_;

	AST* zero = integer(0);

	AST* delta = max(e, zero);

	delete e;
	delete zero;
	
	// AST* ex = power(x->copy(), n->copy());
	AST* lcv = coefficientGPE(v, x, n);
	// delete ex;

	int tal = 0;
	while(m->kind() != Kind::MinusInfinity && m->value() >= n->value()) {
		// AST* ex_ = power(x->copy(), m->copy());
		AST* lcs = coefficientGPE(s, x, m);
		// delete ex_;

		p = add({
			mul({lcv->copy(), p}),
			mul({lcs->copy(), power(x->copy(), sub({m->copy(), n->copy()}))}),
		});

		AST* s_ = sub({
			mul({lcv->copy(), s}),
			mul({
				lcs->copy(),
				v->copy(),
				power(
					x->copy(),
					sub({
						m->copy(),
						n->copy(),
					})
				)
			}),
		});
	
		s = algebraicExpand(s_);
		delete s_;
	
		tal = tal + 1;

		delete m;
		m = degreeGPE(s, x);
		delete lcs;

	}

	AST* resQ = mul({
		power(
			lcv->copy(),
			sub({delta->copy(), integer(tal)})
		),
		p->copy()
	});

	AST* resR = mul({
		power(
			lcv->copy(),
			sub({ delta->copy(), integer(tal) })
		),
		s->copy()
	});


	AST* res = list({
		algebraicExpand(resQ),
		algebraicExpand(resR)
	});
	
	delete p;
	delete s;
	delete m;
	delete n;
	delete resQ;
	delete resR;
	delete delta;
	delete lcv;
	
	return res;
}

AST* pseudoQuotient(AST* u, AST* v, AST* x) {
	AST* r = pseudoDivision(u, v, x);
	AST* q = r->operand(0)->copy();
	delete r;
	return q;
}

AST* pseudoRemainder(AST* u, AST* v, AST* x) {
	AST* r = pseudoDivision(u, v, x);
	AST* q = r->operand(1)->copy();
	delete r;
	return q;
}




AST* getNormalizationFactor(AST* u, AST* L, AST* K) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return integer(0);

	if(isConstant(u)) {
		if(isGreaterZero(u)) {
			if(K->identifier() == "Z") {
				return integer(1);
				// return u->copy();
			}
			if(K->identifier() == "Q") {
				return power(u->copy(), integer(-1));
			}

			// undefined integral domain
			return undefined();
		} else {
			if(K->identifier() == "Z") {
				return integer(-1);
				// return mul({ integer(-1), u->copy() });
			}
			if(K->identifier() == "Q") {
				return mul({ integer(-1), power(u->copy(), integer(-1)) });
			}

			// undefined integral domain
			return undefined();
		}
	}

	if(L->numberOfOperands() == 0)
		return undefined();
	
	AST* lc = leadingCoefficientGPE(u, L->operand(0));
	AST* restL = rest(L);

	AST* c = getNormalizationFactor(lc, restL, K);

	delete restL;
	delete lc;

	return c;
}

AST* normalizePoly(AST* u, AST* L, AST* K) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return integer(0);

	AST* u__ = mul({ getNormalizationFactor(u, L, K), u->copy() });

	AST* u_ = algebraicExpand(u__);

	delete u__;
	
	return u_;
}

// Finds the content of u with respect to x using
// the auxiliary variables R with coefficient domain K,
// with is Z or Q
AST* polynomialContent(AST* u, AST* x, AST* R, AST* K) {
	if(u->kind() == Kind::Integer && u->value() == 0) {
		return integer(0);
	}

	AST* deg = degreeGPE(u, x);

	AST* gcd = coefficientGPE(u, x, deg);

	for(int i=deg->value() - 1; i>= 0; i--) {
		AST* d = integer(i);
		AST* coef = coefficientGPE(u, x, d);
		
		delete d;

		AST* gcd_ = mvPolyGCD(gcd, coef, R, K);

		delete gcd;
		delete coef;
		gcd = gcd_;
	}

	delete deg;

	return gcd;
}

// AST* polynomialPrimitivePart(AST* u, AST* x, AST* R, AST* K)
// {
// 	return recQuotient(u, x, R, K);
// }

AST* mvPolyGCDRec(AST* u, AST* v, AST* L, AST* K) 
{
	if(L->numberOfOperands() == 0) {
		if(K->identifier() == "Z") { return integerGCD(u, v); }
		if(K->identifier() == "Q") { return integer(1); }
	}

	AST* x = first(L);
	AST* R = rest(L);

	AST* cont_u = polynomialContent(u, x, R, K);
	AST* cont_v = polynomialContent(v, x, R, K);

	AST* d = mvPolyGCDRec(cont_u, cont_v, R, K);

	AST* pp_u = recQuotient(u, cont_u, L, K);
	AST* pp_v = recQuotient(v, cont_v, L, K);


	while(pp_v->kind() != Kind::Integer || pp_v->value() != 0) {
		
		AST* r = pseudoRemainder(pp_u, pp_v, x);

		AST* pp_r;
		
		if(r->kind() == Kind::Integer && r->value() == 0) {
			pp_r = integer(0);
		} else {
			AST* cont_r = polynomialContent(r, x, R, K);
			pp_r = recQuotient(r, cont_r, L, K);
			
			delete cont_r;
		}

		delete r;
	
		delete pp_u;
		pp_u = pp_v;
		pp_v = pp_r;
	}
	
	AST* k = mul({ d->copy(), pp_u->copy() });
	AST* result = algebraicExpand(k);

	delete k;

	delete x;
	delete d;
	delete cont_u;
	delete cont_v;
	delete pp_u;
	delete pp_v;
	delete R;

	return result;
}


AST* mvPolyGCD(AST* u, AST* v, AST* L, AST* K) {
	if(u->kind() == Kind::Integer && u->value() == 0) {
		return normalizePoly(v, L, K);
	}

	if(v->kind() == Kind::Integer && v->value() == 0) {
		return normalizePoly(u, L, K);
	}

	AST* gcd = mvPolyGCDRec(u, v, L, K);

	AST* r = normalizePoly(gcd, L, K);

	delete gcd;
	
	return r;
}

AST* leadingMonomial(AST* u, AST* L) {
	if(L->numberOfOperands() == 0) {
		return u->copy();
	}

	AST* x = first(L);
	AST* m = degreeGPE(u, x);

	AST* c = coefficientGPE(u, x, m);

	AST* restL = rest(L);

	AST* r_ = mul({
		power(x->copy(), m->copy()),
		leadingMonomial(c, restL)
	});

	delete c;
	delete m;
	delete x;
	delete restL;

	AST* r = algebraicExpand(r_);

	delete r_;

	return r;
}

bool wasSimplified(AST* u) {
	if(u->kind() == Kind::Symbol)
		return true;

	if(u->kind() == Kind::Integer)
		return true;

	if(u->kind() == Kind::Division)
		return false;
	
	if(u->kind() == Kind::Multiplication) {
		
		for(unsigned int i=0; i<u->numberOfOperands(); i++) {
			if(u->operand(i)->kind() == Kind::Fraction)
			{
				return false;
			}

			if(
				u->operand(i)->kind() == Kind::Power &&
				u->operand(i)->operand(1)->kind() == Kind::Integer &&
				u->operand(i)->operand(1)->value() < 0
			)	
			{
				return false;
			}
		}

		return true;
	}

	return false;
}


/**
 * Return summation(u[i]/v) if v divides u[i]
 */
AST* G(AST* u, AST* v) {
	if(u->kind() == Kind::Addition || u->kind() == Kind::Subtraction) {
		AST* k = new AST(u->kind());

		for(unsigned int i=0; i<u->numberOfOperands(); i++) 
		{
			AST* z_ = div(u->operand(i)->copy(), v->copy());
			AST* z  = algebraicExpand(z_);

			delete z_;

			if(wasSimplified(z))
			{
				k->includeOperand(z);
			}
			else
			{
				delete z;
			}
		}

		if(k->numberOfOperands() == 0) {
			delete k;
			return integer(0);
		}
	
		return k;
	}

	AST* z_ = div(u->copy(), v->copy());

	AST* z = algebraicExpand(z_);

	delete z_;

	if(wasSimplified(z)) 
	{
		return z;
	}

	delete z;

	return integer(0);

	// printf("%i\n", wasSimplified(z));
	// printf("%i\n", k->numberOfOperands());

	// if(k->numberOfOperands() == 0) {
	// 	delete k;
	// 	return integer(0);
	// }

	// AST* r = algebraicExpand(k);

	// delete k;
	
	// return r;
}

AST* monomialPolyDiv(AST* u, AST* v, AST* L) {
	AST* q = integer(0);
	AST* r = u->copy();
	AST* vt = leadingMonomial(v, L);

	AST* f = G(r, vt);

	// printf("-> %s\n", r->toString().c_str());
	// printf("-> %s\n", vt->toString().c_str());
	// printf("-> %s\n", f->toString().c_str());
	// printf("-> %i\n", f->kind());

	while(f->kind() != Kind::Integer || f->value() != 0) {
		// printf("-> %s\n", f->toString().c_str());
		// printf("-> %i\n", f->kind());

		q = add({ q, f->copy() });

		AST* r_ = sub({ r, mul({ f, v->copy() }) });
		// printf("-> %s\n", r_->toString().c_str());
	
		r = algebraicExpand(r_);
		// printf("-> %s\n", r_->toString().c_str());

		delete r_;

		f = G(r, vt);
		// printf("\n");
	}
	
	AST* l = list({ reduceAST(q), reduceAST(r) });
	
	delete q;
	delete r;
	delete vt;
	delete f;

	return l;
}


// TODO

// monomialBasedPolyExpansion(a^2*b + 2*a*b^2 + b^3 + 2*a + 2*b + 3, a+b, [a, b], t) -> b*t^2 + 2*t + 3
AST* monomialBasedPolyExpansion(AST* u, AST* v, AST* L, AST* t) {
	if(u->kind() == Kind::Integer && u->value() == 0) {
		return integer(0);
	}

	AST* d = monomialPolyDiv(u, v, L);
	AST* q = d->operand(0)->copy();
	AST* r = d->operand(1)->copy();

	AST* k = add({
		mul({
			t->copy(),
			monomialBasedPolyExpansion(q, v, L, t),
		}),
		r
	});

	AST* x = algebraicExpand(k);

	delete k;
	delete q;
	delete d;

	return x;
}

// monomialPolyRem can be used for simplification, for example
// monomialPolyRem(a*i^3 + b*i^2 + c*i + d, i^2 + 1, [i]) -> -a*i - b + c*i + d
// simplification when i^2 + 1 = 0
// also
//monomialPolyRem(sin^4(x)+sin^3(x)+2*sin^2(x)cos^2(x)+cos^4(x), sin^2(x)+cos^2(x)-1, [cos(x), sin(x)]) -> 1 + sin^3(x)
AST* monomialPolyRem(AST* u, AST* v, AST* L) {
	AST* d = monomialPolyDiv(u,v,L);
	AST* r = d->operand(1)->copy();
	delete d;
	return r;
}

AST* monomialPolyQuo(AST* u, AST* v, AST* L)
{
	AST* d = monomialPolyDiv(u,v,L);
	AST* r = d->operand(0)->copy();
	delete d;
	return r;
}

AST* expandProduct(AST* r, AST* s) 
{
	if(r->kind() == Kind::Addition) 
	{
		AST* f = r->operand(0);

		AST* v = sub({
			r->copy(),
			f->copy(),
		});

		AST* k = reduceAST(v);
	
		delete v;

		AST* z = add({
			expandProduct(f, s),
			expandProduct(k, s),
		});

		delete k;

		AST* y = reduceAST(z);

		delete z;

		return y;
	}
	else if(s->kind() == Kind::Addition) 
	{
		return expandProduct(s, r);
	}

	AST* a = algebraicExpand(r);
	AST* b = algebraicExpand(s);

	if(a->kind() == Kind::Addition || b->kind() == Kind::Addition) 
	{
		AST* t = expandProduct(a, b);
		
		delete a;
		delete b;
		
		return t;
	}

	AST* t = mul({ a, b });

	AST* k = reduceAST(t);

	delete t;

	return k;
}

long fact(long i) 
{
	if(i==1 || i==0)
		return 1;
	return i * fact(i-1);
}

AST* expandPower(AST* u, AST* n)
{
	if(u->kind() == Kind::Addition) 
	{
		printf("u(x) = %s\n", u->toString().c_str());
		AST* f = u->operand(0);

		AST* o = sub({ u->copy(), f->copy() });

		AST* r = reduceAST(o);

		delete o;

		printf("1. %s\n", r->toString().c_str());

		
		// r = u->copy();

		// if(r->numberOfOperands() > 1)
		// {
		// 	r->removeOperand(0L);
		// }
		// else
		// {
		// 	delete r;
		// 	r = integer(0);
		// }
	
		printf("2. %s\n", r->toString().c_str());

		AST* s = integer(0);
	
		for(int k = 0; k <= n->value(); k++) {
			int a = n->value();
		
			int d = fact(a) / (fact(k) * fact(a - k));
	
			AST* z = mul({
				integer(d),
				power(f->copy(), integer(a - k))
			});

			AST* b = integer(k);
			AST* t = expandPower(r, b);
		
			s = add({ s, expandProduct(z, t) });

			delete b;
			delete z;
			delete t;
		}
		
		delete r;

		return s;
	}
	
	AST* v = power(u->copy(), n->copy());

	return v;
}

AST* algebraicExpand(AST* u) {
	if(u->isTerminal())
	{
		return reduceAST(u);
	}

	AST* u_ = reduceAST(u);
	
	if(u_->kind() == Kind::Addition) {
		AST* v = u_->operand(0)->copy();
		
		u_->deleteOperand(0);

		AST* t = add({
			algebraicExpand(v),
			algebraicExpand(u_)
		});

		delete u_;
		u_ = reduceAST(t);

		delete v;
	} 
	else if(u_->kind() == Kind::Multiplication) 
	{
		AST* v = u_->operand(0)->copy();
		u_->deleteOperand(0);

		AST* z = expandProduct(u_, v);

		delete u_;
		u_ = reduceAST(z);

		delete z;
	} 
	else if(u_->kind() == Kind::Power) 
	{
		AST* b = u_->operand(0);
		AST* e = u_->operand(1);

		if(e->kind() == Kind::Integer && e->value() >= 2)
		{
			AST* t = expandPower(b, e);

			delete u_;
			u_ = reduceAST(t);
			delete t;
		}
	
		if(e->kind() == Kind::Integer && e->value() <= -2) 
		{
			AST* p = power(b->copy(), integer(e->value() - 1));
			
			AST* b_ = p->operand(0);
			AST* e_ = p->operand(1);
	
			AST* t = expandPower(b_, e_);
			
			delete p;
			
			delete u_;

			u_ = power(t, integer(-1));
		}
	}

	AST* t = mapUnaryAST(u_, algebraicExpand);

	delete u_;

	u_ = reduceAST(t);

	delete t;

	return u_;
}



AST* expandProductRoot(AST* r, AST* s) {
	if(r->kind() == Kind::Addition) {
		AST* f = r->operand(0);

		AST* v = sub({
			r->copy(),
			f->copy(),
		});
	
		AST* k = reduceAST(v);
	
		AST* z = add({
			mul({f->copy(), s->copy()}),
			mul({k->copy(), s->copy()}),
		});

		AST* y = reduceAST(z);

		delete v;
		delete k;
		delete z;

		return y;
	}

	if(s->kind() == Kind::Addition) {
		return expandProductRoot(s, r);
	}

	// AST* a = algebraicExpand(r);
	// AST* b = algebraicExpand(s);

	// if(a->kind() == Kind::Addition || b->kind() == Kind::Addition) {
	// 	AST* t = expandProduct(a, b);
		
	// 	delete a;
	// 	delete b;
		
	// 	return t;
	// }

	AST* t = mul({ r->copy(), s->copy() });

	AST* k = reduceAST(t);

	delete t;

	return k;
}

AST* expandPowerRoot(AST* u, AST* n) {
	if(u->kind() == Kind::Addition) {
		AST* f = u->operand(0);

		AST* r_ = sub({ u->copy(), f->copy() });

		AST* r = reduceAST(r_);
		
		delete r_;

		AST* s = integer(0);
	
		for(int k_ = 0; k_ <= n->value(); k_++) {
			AST* k = integer(k_);

			AST* c_ = div(
				integer(fact(n->value())),
				integer(fact(k->value()) * fact(n->value() - k->value()))
			);
	
			AST* c = reduceAST(c_);
	
	
			AST* z_ = mul({
				c->copy(),
				power(
					f->copy(),
					integer(n->value() - k->value())
				)
			});

			AST* z = reduceAST(z_);


			AST* t = expandPowerRoot(r, k);
		
			s = add({ s, expandProductRoot(z, t) });

			delete c;
			delete k;
			delete z;
			delete t;
			delete z_;
			delete c_;
		}
		
		delete r;

		return s;
	}
	
	AST* v = power(
		u->copy(),
		n->copy()
	);

	AST* reduced = reduceAST(v);
	delete v;

	return reduced;
}

AST* algebraicExpandRoot(AST* u) {
	if(u->isTerminal())
		return reduceAST(u);
	
	AST* u_ = reduceAST(u);
	
	if(u_->kind() == Kind::Addition) {
		AST* v = u_->operand(0);
	
		AST* a = sub({
			u_->copy(),
			v->copy()
		});

		AST* k = reduceAST(a);
	
		AST* t = add({
			algebraicExpandRoot(v),
			algebraicExpandRoot(k)
		});

		delete u_;
		u_ = reduceAST(t);
		
		delete k;
		delete a;
		delete t;
	}

	if(u_->kind() == Kind::Multiplication) {

		AST* v = u_->operand(0);
		AST* e = div(
			u_->copy(),
			v->copy()
		);
		
		AST* t = reduceAST(e);

		AST* z = expandProductRoot(t, v);

		delete u_;
		u_ = reduceAST(z);

		delete t;
		delete e;
		delete z;

	}

	if(u_->kind() == Kind::Power) {

		AST* b = u_->operand(0)->copy();
		AST* e = u_->operand(1)->copy();

		if(e->kind() == Kind::Integer && e->value() >= 2) {
			AST* t = expandPowerRoot(b, e);

			delete u_;
			u_ = reduceAST(t);
			delete t;
		}
	
		if(e->kind() == Kind::Integer && e->value() <= -2) {
			AST* p_ = power(u_->copy(), integer(-1));
			AST* p = reduceAST(p_);

			delete p_;
			
			AST* b_ = p->operand(0);
			AST* e_ = p->operand(1);
	
			AST* t = expandPowerRoot(b_, e_);
			
			delete p;
			
			delete u_;
			u_ = power(t, integer(-1));
		}

		delete b;
		delete e;
	}
	
	AST* k = reduceAST(u_);

	delete u_;

	return k;
}



}
