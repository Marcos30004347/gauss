#include "Polynomial.hpp"
#include "Resultant.hpp"
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

AST* degree(AST* u, AST* v) {
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

AST* coeff(AST* u, AST* x, AST* j) 
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

AST* leadCoeff(AST* u, AST* x) {

	AST* d = degree(u, x);

	AST* lc = coeff(u, x, d);

	delete d;

	return lc;
}

AST* divideGPE(AST* u, AST* v, AST* x) {
	AST *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *t9, *t10;

	AST* q = integer(0);
	AST* r = u->copy();

	AST* m = degree(r, x);
	AST* n = degree(v, x);

	AST* lcv = leadCoeff(v, x);

	while(
		m->kind() != Kind::MinusInfinity &&
		(m->kind() == Kind::Integer && n->kind() == Kind::Integer &&
		m->value() >= n->value())
	) {
		AST* lcr = leadCoeff(r, x);

		AST* s = div(lcr->copy(), lcv->copy());

		t1 = power(x->copy(), sub({m->copy(), n->copy()}));
		t2 = mulPoly(s, t1);
		t3 = addPoly(q, t2);

		delete q;
	
		q = reduceAST(t3);
	
		delete t1;
		delete t2;
		delete t3;

		t1 = power(x->copy(), m->copy());
		t2 = mulPoly(lcr, t1);
		t3 = subPoly(r, t2);

		t4 = power(x->copy(), n->copy());
		t5 = mulPoly(lcv, t4);
		t6 = subPoly(v, t5);

		t7 = mulPoly(t6, s);
		t8 = power(x->copy(), sub({m->copy(), n->copy()}));
		t9 = mulPoly(t7, t8);
		t10 = subPoly(t3, t9);
	
		delete r;
	
		r = reduceAST(t10);

		delete t1;
		delete t2;
		delete t3;
		delete t4;
		delete t5;
		delete t6;
		delete t7;
		delete t8;
		delete t9;
		delete t10;
	
		delete m;
		delete lcr;
		delete s;
	
		m = degree(r, x);

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

	AST* e = mul({div(integer(1), leadCoeff(U,x)), U});
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

	AST* c = leadCoeff(U, x);

	AST* App_ = quotientGPE(App, c, x);
	// AST* App_ = algebraicExpand(App__);
	// delete App__;

	delete App;
	App = App_;

	AST* Bpp_ = quotientGPE(Bpp, c, x);
	// AST* Bpp_ = algebraicExpand(Bpp__);
	// delete Bpp__;

	delete Bpp;
	Bpp = Bpp_;

	AST* U_ = quotientGPE(U, c, x);
	// AST* U_ = algebraicExpand(U__);
	// delete U__;
	
	delete U;
	U = U_;

	delete c;
	delete Ap;
	delete Bp;
	delete V;

	return list({ U, App, Bpp });
}

AST* mulPolyRec(AST* p1, AST* p2)
{
	if(p1->kind() == Kind::Addition)
	{
		AST* res = add({});

		for(unsigned int i = 0; i < p1->numberOfOperands(); i++)
		{
			res->includeOperand(mulPoly(p1->operand(i), p2));
		}

		return res;
	}

	if(p2->kind() == Kind::Addition)
	{
		AST* res = add({});
		
		for(unsigned int i = 0; i < p2->numberOfOperands(); i++)
		{
			res->includeOperand(mulPoly(p2->operand(i), p1));
		}

		return res;
	}

	return mul({ p1->copy(), p2->copy() });
}

AST* mulPoly(AST* p1, AST* p2)
{
	AST* t1 = mulPolyRec(p1, p2);
	AST* t2 = reduceAST(t1);
	
	delete t1;

	return t2;
}

AST* subPoly(AST* p1, AST* p2)
{
	AST* s = sub({ p1->copy(), p2->copy() });
	AST* p = reduceAST(s);
	delete s;
	return p;
}

AST* addPoly(AST* p1, AST* p2)
{
	AST* s = add({ p1->copy(), p2->copy() });
	AST* p = reduceAST(s);
	delete s;
	return p;
}

AST* recPolyDiv(AST* u, AST* v, AST* L, AST* K) {
	assert(
		K->identifier() == "Z" || K->identifier() == "Q",
		"Field needs to be Z or Q"
	);

	if(L->numberOfOperands() == 0) 
	{
		AST* k = div(u->copy(), v->copy());
	
		AST* d = algebraicExpand(k);

		delete k;
		
		if(K->identifier() == "Z") 
		{
			if(d->kind() == Kind::Integer) 
			{
				return list({ d, integer(0) });
			}
			
			delete d;

			return list({ integer(0), u->copy() });
		}
	
		return list({ d, integer(0) });
	}

	AST* x = first(L);
	AST* r = u->copy();

	AST* m = degree(r, x);
	AST* n = degree(v, x);
	
	AST* q = integer(0);
	AST* lcv = leadCoeff(v, x);

	while(m->kind() != Kind::MinusInfinity && m->value() >= n->value()) 
	{
		AST* lcr = leadCoeff(r, x);
	
		AST* R = rest(L);

		AST* d = recPolyDiv(lcr, lcv, R, K);

		delete R;
		
		if(d->operand(1)->isNot(0)) 
		{
			AST* result = algebraicExpand(q);

			delete x;
			delete m;
			delete n;
			delete q;
			delete d;
			delete lcv;
			delete lcr;

			return list({ result, r });
		}

		AST* j = power(x->copy(), sub({ m->copy(), n->copy() }));

		q = add({q, mul({ d->operand(0)->copy(), j->copy()})});

		AST* t1 = mulPoly(v, d->operand(0));
		AST* t2 = mulPoly(t1, j);
		AST* t3 = subPoly(r, t2);

		delete r;
	
		r = reduceAST(t3);
	
		delete t1;
		delete t2;
		delete t3;

		delete m;
	
		m = degree(r, x);
	
		delete lcr;
		delete d;
		delete j;
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

AST* pdiv(AST* f, AST* g, AST* x)
{
	assert(g->isNot(0), "Division by zero!");

	AST *lg, *k, *q, *r, *t, *m, *n, *j;
	AST *t1, *t2, *t3, *t4, *t5, *t6;

	m = degree(f, x);
	n = degree(g, x);

	if(m->value() < n->value())
	{
		delete m;
		delete n;

		return list({ integer(0), f->copy() });
	}

	if(g->is(1))
	{
		delete m;
		delete n;

		return list({ f->copy(), integer(0)});
	}

	q = integer(0);
	r = f->copy();
	t = m->copy();

	k = add({ sub({ m, n }), integer(1) });

	lg = leadCoeff(g, x);

	// printf("df\n");
	// printf("%s\n", m->toString().c_str());

	// printf("dg\n");
	// printf("%s\n", n->toString().c_str());

	// printf("lc_g\n");
	// printf("%s\n", lg->toString().c_str());

	while(true)
	{
		t1 = leadCoeff(r, x);
		j = sub({ t, n->copy() });
		k = sub({ k, integer(1) });
		t3 = power(x->copy(), j);

		t2 = mulPoly(q, lg); //mul({ q->copy(), lg->copy() });
		t4 = mulPoly(t1, t3); 
		q  = addPoly(t2, t4);

		delete t2;
		delete t4;
		// q = add({ t2->copy(), mul({ t1->copy(), t3->copy()}) });
		
		t4 = mulPoly(r, lg);// mul({ r->copy(), lg->copy() });
		t5 = mulPoly(g, t1); //mul({ g->copy(), t1->copy(), t3->copy() });
		t6 = mulPoly(t5, t3);
		r  = subPoly(t4, t6);

		// t5 = sub({ t4, t5 });
		// r = algebraicExpand(t6);

		delete t3;
		delete t4;
		delete t5;
		delete t6;

		t = degree(r, x);

		if(t->kind() == Kind::MinusInfinity || t->value() < n->value())
		{
			break;
		}

		delete t1;
	}

	q = mul({q, power(lg->copy(), k->copy())});
	r = mul({r, power(lg->copy(), k->copy())});

	delete k;
	delete t;
	delete lg;

	t1 = algebraicExpand(q);
	t2 = algebraicExpand(r);
	
	delete q;
	delete r;
	
	return list({t1, t2});
}

AST* pseudoDivision(AST* u, AST* v, AST* x) 
{
	AST* p = integer(0);
	AST* s = u->copy();

	AST* m = degree(s, x);
	AST* n = degree(v, x);

	AST* delta = integer(max(m->value() - n->value() + 1, Int(0)));

	AST* lcv = leadCoeff(v, x);

	long tal = 0;

	while(m->kind() != Kind::MinusInfinity && m->value() >= n->value()) 
	{
		AST* lcs = leadCoeff(s, x);

		AST* j = power(x->copy(), sub({ m->copy(), n->copy() }));
		
		AST* t1 = mulPoly(lcv, p);
		AST* t2 = mulPoly(lcs, j);
		AST* t3 = addPoly(t1, t2);
		
		delete t1;
		delete t2;

		delete p;
		p = reduceAST(t3);
		
		delete t3;	

		AST* t4 = mulPoly(lcv, s);
		AST* t5 = mulPoly(lcs, v);
		AST* t6 = mulPoly(t5, j);

		delete s;
	
		s = subPoly(t4, t6);
	
		delete t4;
		delete t5;
		delete t6;

		tal = tal + 1;

		delete m;
	
		m = degree(s, x);
	
		delete lcs;
	
		delete j;
	}

	AST* k = power(lcv->copy(), integer(delta->value() - tal));

	AST* A = mulPoly(k, p);
	AST* B = mulPoly(k, s);

	AST* Q = reduceAST(A);
	AST* R = reduceAST(B);

	delete A;	
	delete B;
	
	delete p;
	delete k;
	delete s;
	delete m;
	delete n;

	delete delta;
	delete lcv;
	
	return list({ Q, R });
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
	assert(
		K->identifier() == "Z" || K->identifier() == "Q", 
		"field must be Z or Q"
	);

	if(u->is(0))
	{
		return integer(0);
	}

	if(isConstant(u)) 
	{
		if(isGreaterZero(u)) 
		{
			if(K->identifier() == "Z") 
			{
				return integer(1);
			}
	
			return power(u->copy(), integer(-1));
		}
		else 
		{
			if(K->identifier() == "Z") 
			{
				return integer(-1);
			}

			return mul({ integer(-1), power(u->copy(), integer(-1)) });
		}
	}

	if(L->numberOfOperands() == 0)
	{
		return undefined();
	}
	
	AST* lc = leadCoeff(u, L->operand(0));

	AST* rL = rest(L);

	AST* cf = getNormalizationFactor(lc, rL, K);

	delete rL;
	delete lc;

	return cf;
}

AST* normalizePoly(AST* u, AST* L, AST* K) 
{
	if(u->kind() == Kind::Integer && u->value() == 0)
	{
		return integer(0);
	}

	AST* u__ = mul({ getNormalizationFactor(u, L, K), u->copy() });

	AST* u_ = algebraicExpand(u__);

	delete u__;
	
	return u_;
}


AST* unitNormal(AST* v, AST* K)
{
	assert(
		K->identifier() == "Z" || K->identifier() == "Q", 
		"field must be Z or Q"
	);

	if(K->identifier() == "Z")
	{
		if(isLessZero(v))
		{
			return integer(-1);
		}
	
		return integer(1);
	}

	if(K->identifier() == "Q")
	{
		if(isLessZero(v))
		{
			return mul({ integer(-1), power(v->copy(), integer(-1)) });
		}

		return power(v->copy(), integer(-1));
	}

	return integer(1);
}

// Finds the content of u with respect to x using
// the auxiliary variables R with coeff domain K,
// with is Z or Q
AST* polynomialContent(AST* u, AST* x, AST* R, AST* K) 
{	
	if(u->is(0))
	{
		return integer(0);
	}

	AST* n = degree(u, x);

	AST* g = coeff(u, x, n);

	AST* k = sub({ u->copy(), mul({g->copy(), power(x->copy(), n->copy())}) });

	AST* v = algebraicExpand(k);

	delete k;

	if(v->is(0))
	{
		AST* un = unitNormal(g, K);
		AST* t = mul({un, g->copy()});
	
		delete g;
	
		g = reduceAST(t);
		
		delete t;
	}
	else
	{
		while(v->isNot(0))
		{
			AST* d = degree(v, x);
			AST* c = leadCoeff(v, x);
		
			AST* t = mvPolyGCD(g, c, R, K);

			delete g;

			g = t;

			k = sub({v, mul({c->copy(), power(x->copy(), d->copy())})});
			v = algebraicExpand(k);

			delete c;
			delete d;

			delete k;
		}
	}

	delete n;
	delete v;

	return g;
}


// Finds the content of u with respect to x using
// the auxiliary variables R with coeff domain K,
// with is Z or Q
AST* polynomialContentSubResultant(AST* u, AST* x, AST* R, AST* K) 
{
	if(u->is(0))
	{
		return integer(0);
	}

	AST* n = degree(u, x);

	AST* g = coeff(u, x, n);

	AST* k = sub({u->copy(), mul({g->copy(), power(x->copy(), n->copy())})});

	AST* v = algebraicExpand(k);

	delete k;

	if(v->is(0))
	{
		AST* un = unitNormal(g, K);
		AST* t = mul({un, g->copy()});
	
		delete g;
	
		g = reduceAST(t);
		
		delete t;
	}
	else
	{
		while(v->isNot(0))
		{
			AST* d = degree(v, x);
			AST* c = leadCoeff(v, x);

			AST* t = mvSubResultantGCD(g, c, R, K);

			delete g;

			g = t;

			k = sub({v, mul({c->copy(), power(x->copy(), d->copy())})});
			v = algebraicExpand(k);

			delete c;
			delete d;

			delete k;
		}
	}


	delete n;
	delete v;
	
	return g;
}

AST* subResultantGCDRec(AST* u, AST* v, AST* L, AST* K)
{
	if(L->numberOfOperands() == 0) 
	{
		if(K->identifier() == "Z") 
		{ 
			return integerGCD(u, v); 
		}
	
		if(K->identifier() == "Q") 
		{ 
			return integer(1); 
		}
	}

	AST* x = first(L);

	AST* du = degree(u, x);
	AST* dv = degree(v, x);
	
	AST* U = nullptr;
	AST* V = nullptr;

	if(du->value() >= dv->value())
	{
		U = u->copy();
		V = v->copy();
	}
	else
	{
		U = v->copy();
		V = u->copy();
	}

	delete du;
	delete dv;

	AST* R = rest(L);
	
	AST* contU = polynomialContentSubResultant(U, x, R, K);
	AST* contV = polynomialContentSubResultant(V, x, R, K);

	AST* d = subResultantGCDRec(contU, contV, R, K);

	AST* tmp1 = recQuotient(U, contU, L, K);
	AST* tmp2 = recQuotient(V, contV, L, K);

	delete U;
	U = tmp1;
	
	delete V;
	V = tmp2;

	AST* tmp3 = leadCoeff(U, x);
	AST* tmp4 = leadCoeff(V, x);

	AST* g = subResultantGCDRec(tmp3, tmp4, R, K);

	delete tmp3;
	delete tmp4;

	int i = 1;

	AST* delta = nullptr;
	AST* y = nullptr;
	AST* b = nullptr;
	AST* dp = nullptr;

	while (V->isNot(0))
	{
		AST* r = pseudoRemainder(U, V, x);
	
		if(r->isNot(0))
		{
			if(i == 1)
			{
	
				AST* tmp3 = add({
					degree(U, x),
					mul({integer(-1), degree(V, x) }),
					integer(1)
				});

				delta = algebraicExpand(tmp3);

				delete tmp3;

				y = integer(-1);
				
				AST* tmp4 = power(integer(-1), delta->copy());
				
				b = algebraicExpand(tmp4);
				
				delete tmp4;
			}
			else
			{
				dp = delta->copy();
			
				AST* tmp3 = add({
					degree(U, x),
					mul({integer(-1), degree(V, x)}),
					integer(1)
				});

				delete delta;
				delta = algebraicExpand(tmp3);

				delete tmp3;

				AST* f = leadCoeff(U, x);

				AST* tmp4 = power(mul({integer(-1), f->copy()}), sub({dp->copy(), integer(1)}));
				AST* tmp5 = power(y->copy(), sub({dp->copy(), integer(2)}));
				
				AST* tmp6 = algebraicExpand(tmp4);
				AST* tmp7 = algebraicExpand(tmp5);
				
				delete tmp4;
				delete tmp5;
				
				y = recQuotient(tmp6, tmp7, R, K);

				delete tmp6;
				delete tmp7;

				AST* tmp8 = mul({
					integer(-1),
					f->copy(),
					power(y->copy(), sub({ delta->copy(), integer(1) }))
				});
				
				b = algebraicExpand(tmp8);
				
				delete tmp8;		
			}
			
			delete U;
			U = V->copy();

			delete V;

			V = recQuotient(r, b, L, K);

			i = i + 1;
		}
		else
		{
			delete U;
			U = V->copy();

			delete V;
			V = r->copy();
		}

		delete r;
	}

	delete delta;
	delete y;
	delete b;

	AST* tmp5 = leadCoeff(U, x);

	AST* s = recQuotient(tmp5, g, R, K);

	delete tmp5;

	AST* W = recQuotient(U, s, L, K);

	delete s;

	AST* contW = polynomialContentSubResultant(W, x, R, K);
	AST* ppW = recQuotient(W, contW, L, K);
		
	AST* tmp6 = mul({d->copy(), ppW->copy()});
	AST* res = algebraicExpand(tmp6);
	
	delete tmp6;
	delete contW;
	delete ppW;

	delete W;
	delete U;
	delete V;

	delete contU;
	delete contV;

	delete g;
	delete d;

	delete R;
	delete x;

	return res;
}

AST* mvSubResultantGCD(AST* u, AST* v, AST* L, AST* K) {
	if(u->is(0)) {
		return normalizePoly(v, L, K);
	}

	if(v->is(0)) {
		return normalizePoly(u, L, K);
	}

	AST* gcd = subResultantGCDRec(u, v, L, K);

	AST* r = normalizePoly(gcd, L, K);

	delete gcd;
	
	return r;
}

AST* mvPolyGCDRec(AST* u, AST* v, AST* L, AST* K) 
{
	if(L->numberOfOperands() == 0) 
	{
		if(K->identifier() == "Z") 
		{ 
			return integerGCD(u, v); 
		}
	
		if(K->identifier() == "Q") 
		{ 
			return integer(1); 
		}
	}

	AST* x = first(L);
	AST* R = rest(L);

	AST* cont_u = polynomialContent(u, x, R, K);

	AST* cont_v = polynomialContent(v, x, R, K);

	AST* d = mvPolyGCDRec(cont_u, cont_v, R, K);

	AST* pp_u = recQuotient(u, cont_u, L, K);
	AST* pp_v = recQuotient(v, cont_v, L, K);

	while(pp_v->isNot(0)) 
	{

		AST* r = pseudoRemainder(pp_u, pp_v, x);
	
		AST* pp_r = nullptr;

		if(r->is(0)) {
			pp_r = integer(0);
		} 
		else 
		{
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
	if(u->is(0)) {
		return normalizePoly(v, L, K);
	}

	if(v->is(0)) {
		return normalizePoly(u, L, K);
	}

	AST* gcd = mvPolyGCDRec(u, v, L, K);
	AST* r = normalizePoly(gcd, L, K);

	delete gcd;
	
	return r;
}

AST* leadMonomial(AST* u, AST* L) {
	if(L->numberOfOperands() == 0) {
		return u->copy();
	}

	AST* x = first(L);
	AST* m = degree(u, x);

	AST* c = coeff(u, x, m);

	AST* restL = rest(L);

	AST* r_ = mul({
		power(x->copy(), m->copy()),
		leadMonomial(c, restL)
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
	AST* vt = leadMonomial(v, L);

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

AST* algebraicExpandRec(AST* u);

AST* expandProduct(AST* r, AST* s) 
{
	if(r->is(0)) return integer(0);
	if(s->is(0)) return integer(0);

	if(r->kind() == Kind::Addition && r->numberOfOperands() == 0) return integer(0);
	if(s->kind() == Kind::Addition && s->numberOfOperands() == 0) return integer(0);
	
	if(r->kind() == Kind::Addition) 
	{
		AST* f = r->operand(0);

		AST* k = r->copy();
	
		k->deleteOperand(0);

		AST* z = add({ expandProduct(f, s), expandProduct(k, s) });

		delete k;

		AST* y = reduceAST(z);

		delete z;

		return y;
	}
	else if(s->kind() == Kind::Addition) 
	{
		return expandProduct(s, r);
	}

	AST* a = algebraicExpandRec(r);
	AST* b = algebraicExpandRec(s);

	if(a->kind() == Kind::Addition || b->kind() == Kind::Addition) 
	{
		AST* t = expandProduct(a, b);

		delete a;
		delete b;
		
		return t;
	}

	return mul({ a, b });
}


AST* expandPower(AST* u, AST* n)
{
	if(u->kind() == Kind::Addition) 
	{
		AST* f = u->operand(0);
	
		AST* o = sub({ u->copy(), f->copy() });

		AST* s = integer(0);

		for(long k = 0; k <= n->value(); k++) 
		{
			Int a = n->value();
		
			Int d = fact(a) / (fact(k) * fact(a - k));
	
			AST* z = mul({
				integer(d),
				power(f->copy(), integer(a - k))
			});

			AST* b = integer(k);
			AST* t = expandPower(o, b);
		
			s = add({ s, expandProduct(z, t) });

			delete b;
			delete z;
			delete t;
		}
	
		delete o;

		return s;
	}

	AST* v = power(u->copy(), n->copy());

	return v;
}


AST* expandProductRoot(AST* r, AST* s) 
{
	if(r->kind() == Kind::Addition) 
	{
		AST* f = r->operand(0);

		AST* k = r->copy();
	
		k->deleteOperand(0);
	
		AST* z = add({
			mul({f->copy(), s->copy()}),
			mul({k->copy(), s->copy()}),
		});

		delete k;

		return z;
	}

	if(s->kind() == Kind::Addition) 
	{
		return expandProductRoot(s, r);
	}

	return mul({ r->copy(), s->copy() });
}

AST* expandPowerRoot(AST* u, AST* n) 
{
	if(u->kind() == Kind::Addition) {
		AST* f = u->operand(0);

		AST* r_ = sub({ u->copy(), f->copy() });

		AST* r = reduceAST(r_);
		
		delete r_;

		AST* s = integer(0);
	
		for(int k_ = 0; k_ <= n->value(); k_++) 
		{
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

AST* algebraicExpandRoot(AST* u) 
{
	if(u->isTerminal())
		return reduceAST(u);

	AST* u_ = reduceAST(u);

	if(u_->kind() == Kind::Addition)
	{
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

	if(u_->kind() == Kind::Multiplication)
	{
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


AST* algebraicExpandRec(AST* u) 
{
	if(u->isTerminal()) return u->copy();
	
	// Those kinds needs to be reduced to other
	// kinds before expansion, 
	// subtraction -> addition
	// division    -> multiplication by inverses
	// factorial   -> possible to an integer
	AST* z = nullptr;

	if(
		u->kind() == Kind::Subtraction ||
		u->kind() == Kind::Division    ||
		u->kind() == Kind::Factorial 
	)
	{
		z = reduceAST(u);
	}
	else
	{
		z = u->copy();
	}


	if(z->kind() == Kind::Power) 
	{
		// printf("pow -> : %s\n", z->toString().c_str());
		AST* b = z->operand(0);
		AST* e = z->operand(1);

		if(e->kind() == Kind::Integer)
		{
			AST* t = expandPower(b, e);

			delete z;

			z = t;
		}
		// printf("pow <- : %s\n", z->toString().c_str());
	}

	if(z->kind() == Kind::Multiplication) 
	{
		// printf("mul -> : %s\n", z->toString().c_str());

		AST* v = z->operand(0);

		z->removeOperand(0L);

		if(z->numberOfOperands() == 0)
		{
			delete z;
			z = algebraicExpandRec(v);

			delete v;
		}
		else
		{
			AST* q = expandProduct(z, v);

			delete z;
		
			z = q;

			delete v;
		}
		// printf("mul <- : %s\n", z->toString().c_str());

	}

	if(z->kind() == Kind::Addition) 
	{
		// printf("add -> : %s\n", z->toString().c_str());

		AST* v = z->operand(0);
	
		z->removeOperand(0L);

		if(z->numberOfOperands() == 0)
		{
			delete z;

			z = algebraicExpandRec(v);

			delete v;
		}
		else
		{
			AST* t = add({ algebraicExpandRec(v), algebraicExpandRec(z) });

			delete v;

			delete z;

			z = t;
		}
	
		// printf("add <- : %s\n", z->toString().c_str());
	} 

	AST* k = reduceAST(z);

	delete z;

	return k;
}

AST* algebraicExpand(AST* u)
{
	AST* t = algebraicExpandRec(u);
	
	AST* k = reduceAST(t);
	
	delete t;
	
	return k;
} 


AST* cont(AST* u, AST* x)
{
	AST* n, *c, *c1, *c2, *tmp;

	u = algebraicExpand(u);

	if(u->is(0))
	{
		delete u;

		return integer(0);
	}

	if(u->numberOfOperands() >= 2)
	{
		n = degree(u, x);
		
		c1 = coeff(u, x, n);

		tmp = sub({ u->copy(), mul({ c1->copy(), power(x->copy(), n->copy()) }) });
		
		delete u;
	
		u = algebraicExpand(tmp);
		
		delete tmp;

		delete n;

		n = degree(u, x);
		
		c2 = coeff(u, x, n);

		tmp = sub({ u->copy(), mul({ c2->copy(), power(x->copy(), n->copy()) }) });

		delete u;

		u = algebraicExpand(tmp);
		
		delete tmp;

		c = integerGCD(c1, c2);

		delete n;
		delete c1;
		delete c2;

		while(u->isNot(0))
		{
			n  = degree(u, x);

			c1 = coeff(u, x, n);

			tmp = sub({ u->copy(), mul({ c1->copy(), power(x->copy(), n->copy()) }) });

			delete u;
		
			u = algebraicExpand(tmp);

			delete tmp;

			c2 = integerGCD(c, c1);

			delete c;

			c = c2;

			delete n;
			delete c1;
		}

		delete u;

		return c;
	}

	if(u->numberOfOperands() == 1)
	{
		n = degree(u, x);
		c = coeff(u, x, n);

		delete n;

		delete u;
	
		return c;
	}

	delete u;

	return integer(0);
}

AST* cont(AST* u, AST* L, AST* K)
{
	AST* R = rest(L);
	AST* C = polynomialContentSubResultant(u, L->operand(0), R, K);

	delete R;

	return C;
}

AST* pp(AST* u, AST* L, AST* K)
{
	AST* c = cont(u, L, K);

	AST* p = recQuotient(u, c, L, K);
	
	delete c;

	return p;
}

AST* pp(AST* u, AST* c, AST* L, AST* K)
{
	AST* R = rest(L);

	AST* p = recQuotient(u, c, L, K);
	
	delete R;

	return p;
}


}
