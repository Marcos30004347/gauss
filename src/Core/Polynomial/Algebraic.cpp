#include "Algebraic.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Algebra.hpp"

using namespace ast;
using namespace algebra;

namespace polynomial
{

AST* algMulInverse(AST* v, AST* p, AST* a) {
	// TODO: assert following statements
	// 1. a is a symbol that represents an algebraic number
	// 2. p is a monoic irrecudctible polynomial in Q[a] with degree(p,a) >= 2
	// 3. v is a non zero polynomial in Q(a) with degree(v) < degree(p)

	AST* w = extendedEuclideanAlgGPE(v,p,a);
	AST* r = w->operand(1)->copy();
	delete w;
	return r;
}

AST* algDivide(AST* u, AST* v, AST* p, AST* a) {
	// TODO: assert following statements
	// a is a symbol that represents an algebraic number;
	// p is a monic, irrreducible polynomial in Q[α] with deg(p, α) ≥ 2;
	// u and v are both polynomials in Q(a) with degree < deg(p) and v != 0;
	AST* w = algMulInverse(v, p, a);

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


	AST* d = degree(u, x);

	if(d->value() == 0) {
		delete d;
		return remainderGPE(u, p, a);
	}

	AST* r = new AST(u->kind());

	for(int i=0; i <= d->value(); i++) {
		AST* d = integer(i);

		AST* coeff_ = coeff(u, x, d);
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

AST* algPolynomialDivision(AST* u, AST* v, AST* x, AST* p, AST* a) {
	// TODO: assert following statements
	// u, v : polynomials in Q(a)[x] with v != 0;
	// x : a symbol;
	// a : a symbol that represents an algebraic number;
	// p : a monic, irrreducible polynomial in Q[a] with degree ≥ 2;

	AST* q = integer(0);
	AST* r = u->copy();
	AST* m = degree(r, x);
	AST* n = degree(v, x);
	AST* lcv = leadCoeff(v, x);

	
	AST* p_ = deepReplace(p, x, a);


	while(
		m->kind() != Kind::MinusInfinity && (
			m->kind() == Kind::Integer && n->kind() == Kind::Integer &&
			m->value() >= n->value()
		)
	) {
	
		AST* lcr = leadCoeff(r, x);

		AST* s = algDivide(lcr, lcv, p_, a);
		
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

		m = degree(r, x);
		
		delete s;
		delete lcr;		
	}

	delete m;
	delete n;
	delete lcv;
	delete p_;

	return list({ q, r });
}

AST* algPolynomialRemainder(AST* u, AST* v, AST* x, AST* p, AST* a) {
	AST* res = algPolynomialDivision(u,v,x,p,a);
	AST* r = res->operand(1)->copy();
	delete res;
	return r;
}

AST* algPolynomialQuotient(AST* u, AST* v, AST* x, AST* p, AST* a) {
	AST* res = algPolynomialDivision(u,v,x,p,a);
	AST* r = res->operand(0)->copy();
	delete res;
	return r;
}

AST* algPolynomialGCD(AST* u, AST* v, AST* x, AST* p, AST* a) {
	AST* U = u->copy();
	AST* V = v->copy();

	while(
		V->kind() != Kind::Integer ||
		V->value() != 0
	) {
		AST* R = algPolynomialRemainder(U, V, x, p, a);

		delete U;

		U = V->copy();

		delete V;

		V = R->copy();

		delete R;
	}

	AST* r = algMonic(U, x, p, a);	
	
	delete U;
	delete V;
	
	return r;
}

AST* algMonic(AST* u,AST* x, AST* p,AST* a) {
	AST* lc = leadCoeff(u, x);
	AST* k_ = algPolynomialQuotient(u, lc, x, p, a);
	delete lc;
	return k_;
}
}
