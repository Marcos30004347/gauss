#include "Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"

#include <numeric>

using namespace ast;
using namespace expand;
using namespace simplification;
using namespace algebra;

namespace polynomial {

void includeVariable(std::vector<AST*>& vars, AST* u) {
	// TODO: optimize
	bool included = false;

	for(AST* k : vars) {
		if(k->match(u)){
			included = true;
			break;
		}
	}

	if(!included) {
		vars.push_back(u->deepCopy());
	}
}

std::vector<AST*> variables(AST* u) {
	std::vector<AST*> vars = std::vector<AST*>(0);

	if(
		u->kind() == Kind::Addition ||
		u->kind() == Kind::Subtraction ||
		u->kind() == Kind::Multiplication ||
		u->kind() == Kind::Division ||
		u->kind() == Kind::Factorial
	) {

		for(int i=0; i<u->numberOfOperands(); i++) {
			std::vector<AST*> vargs = variables(u->operand(i));

			for(AST* a : vargs)
				includeVariable(vars, a);

			for(AST* a : vargs)
				delete a;
		}

		return vars;
	}


	if(u->kind() == Kind::Power) {
		if(!isConstant(u->operand(0))) {
			includeVariable(vars, u->operand(0));
		}
	}

	if(u->kind() == Kind::Symbol) {
		includeVariable(vars, u);
	}

	if(u->kind() == Kind::FunctionCall) {
		std::vector<AST*> var_args = std::vector<AST*>(0);

		for(int j=0; j<u->numberOfOperands(); j++) {
			std::vector<AST*> vargs = variables(u->operand(j));
			for(AST* a : vargs)
				var_args.push_back(a);
		}

		if(var_args.size() > 0)
			includeVariable(vars, u);

		for(AST* a : var_args)
			delete a;
	}


	return vars;
}

bool isPolynomialGPE(AST* u, std::vector<AST*> vars) {
	if(u->kind() == Kind::Integer)
		return true;

	std::vector<AST*> vs = variables(u);
	bool inc = false;


	for(int j=0; j<vars.size(); j++) {
		for(int i=j; i<vs.size(); i++) {
			if(vars[j]->match(vs[i])) {
				inc = true;
				break;
			}
		}

		if(!inc) {
			for(AST* k : vs)
				delete k;
			return false;
		}
	}

	for(AST* k : vs)
		delete k;

	return true;
}


AST* degreeGPE(AST* u, AST* x) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return new AST(Kind::MinusInfinity);

	if(!isPolynomialGPE(u, { x })) {
		return inte(0);
	}

	if(
		u->kind() == Kind::Power &&
		// expoent is integer
		u->operand(1)->kind() == Kind::Integer &&
		// base is equal to x
		u->operand(0)->match(x)
	) {
		return u->operand(1)->deepCopy();
	}

	if(
		u->kind() == Kind::Multiplication ||
		u->kind() == Kind::Addition ||
		u->kind() == Kind::Subtraction
	) {

		AST* best = inte(0);
		AST* zero = inte(0);

		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* tmp = degreeGPE(u->operand(i), x);

			if(tmp->kind() != Kind::Integer || tmp->value() == 0) {
				delete tmp;
				continue;
			}

			if(best->match(zero)) {
				delete best;
				best = tmp;
				continue;
			}

			if(tmp->value() > best->value()) {
				delete best;
				best = tmp;
				continue;
			}

			delete tmp;
		}

		if(best->match(zero)) {
			delete zero;
			delete best;
			return mul({inte(-1), inf()});
		}

		delete zero;
		return best;
	}

	if(u->match(x)) {
		return inte(1);
	}

	return inte(0);
}

AST* coefficientXToZeroGPE(AST* u, AST* x) {

	if(
		u->kind() == Kind::Power &&
		u->operand(0)->match(x)
	) {
		return inte(0);
	} else if(u->match(x)) {
		return inte(0);
	}

	if(u->kind() == Kind::Addition || u->kind() == Kind::Subtraction) {
		
		AST* r = new AST(u->kind());
	
		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* c = coefficientXToZeroGPE(u->operand(i), x);
			
			if(c->kind() == Kind::Integer && c->value() == 0) {
				delete c;
				continue;
			}
			r->includeOperand(c);
		}

		if(r->numberOfOperands() == 0) {
			delete r;
			return inte(0);
		}

		if(r->numberOfOperands() == 1) {
			AST* z = r->operand(0)->deepCopy();
			delete r;
			return z;
		}

		return r;
	}

	if(u->kind() == Kind::Multiplication) {
		int x_index = -1;

		for(int i = 0; i<u->numberOfOperands(); i++) {
			if(
				u->operand(i)->kind() == Kind::Power &&
				u->operand(i)->operand(0)->match(x)
			) {
				x_index = i;
				break;
			}

			if(u->operand(i)->match(x)) {
				x_index = i;
				break;
			}
		}
	
		if(x_index == -1) {
			AST* r = new AST(Kind::Multiplication);
			for(int i=0; i<u->numberOfOperands(); i++) {
				r->includeOperand(u->operand(i)->deepCopy());
			}

			if(r->numberOfOperands() == 0) {
				delete r;
				return inte(0);
			}

			if(r->numberOfOperands() == 1) {
				AST* z = r->operand(0)->deepCopy();
				delete r;
				return z;
			}

			return r;
		} else {
			return inte(0);
		}
	}

	return u->deepCopy();
}

AST* coefficientXToOneGPE(AST* u, AST* x) { // x is just a symbol
	if(
		u->kind() == Kind::Power &&
		u->operand(1)->kind() == Kind::Integer &&
		u->operand(0)->match(x) &&
		u->operand(1)->value() == 1
	) {
		return inte(1);
	} else if(u->match(x)) {
		return inte(1);
	}

	if(u->kind() == Kind::Addition || u->kind() == Kind::Subtraction) {
		AST* r = new AST(u->kind());
	
		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* c = coefficientXToOneGPE(u->operand(i), x);
			if(c->kind() == Kind::Integer && c->value() == 0) {
				delete c;
				continue;
			}
			r->includeOperand(c);
		}
		if(r->numberOfOperands() == 0) {
			delete r;
			return inte(0);
		}

		if(r->numberOfOperands() == 1) {
			AST* z = r->operand(0)->deepCopy();
			delete r;
			return z;
		}
	
		return r;
	}

	if(u->kind() == Kind::Multiplication) {
		int x_index = -1;

		for(int i = 0; i<u->numberOfOperands(); i++) {

			if(
				u->operand(i)->kind() == Kind::Power &&
				u->operand(i)->operand(1)->kind() == Kind::Integer &&
				u->operand(i)->operand(0)->match(x) &&
				u->operand(i)->operand(1)->value() == 1
			) {
				x_index = i;
				break;
			} else 
			if(u->operand(i)->match(x)) {
				x_index = i;
				break;
			}
		}
	
		if(x_index != -1) {
			AST* r = new AST(Kind::Multiplication);

			for(int i=0; i < u->numberOfOperands(); i++) {
				if(i == x_index) continue;
				r->includeOperand(u->operand(i)->deepCopy());
			}
		
			if(r->numberOfOperands() == 0) {
				delete r;
				return inte(0);
			}

			if(r->numberOfOperands() == 1) {
				AST* z = r->operand(0)->deepCopy();
				delete r;
				return z;
			}
		
			return r;
		}
	}

	return inte(0);
}

AST* coefficientXToNGPE(AST* u, AST* x) {
	if(u->match(x)) {
		return inte(1);
	}

	if(u->kind() == Kind::Addition || u->kind() == Kind::Subtraction) {
		AST* r = new AST(u->kind());
	
		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* c = coefficientXToNGPE(u->operand(i), x);
			if(c->kind() == Kind::Integer && c->value() == 0) {
				delete c;
				continue;
			}
			r->includeOperand(c);
		}

		if(r->numberOfOperands() == 0) {
			delete r;
			return inte(0);
		}

		if(r->numberOfOperands() == 1) {
			AST* z = r->operand(0)->deepCopy();
			delete r;
			return z;
		}
		return r;
	}

	if(u->kind() == Kind::Multiplication) {
		int x_index = -1;

		for(int i = 0; i<u->numberOfOperands(); i++) {
			if(u->operand(i)->match(x)) {
				x_index = i;
				break;
			}
		}
	
		if(x_index != -1) {
			AST* r = new AST(Kind::Multiplication);

			for(int i=0; i<u->numberOfOperands(); i++) {
				if(i == x_index) continue;
				r->includeOperand(u->operand(i)->deepCopy());
			}

			if(r->numberOfOperands() == 0) {
				delete r;
				return inte(0);
			}

			if(r->numberOfOperands() == 1) {
				AST* z = r->operand(0)->deepCopy();
				delete r;
				return z;
			}

			return r;
		}

	}

	return inte(0);
}

AST* coefficientGPE(AST* u, AST* x) {
	if(x->kind() == Kind::Power && x->operand(1)->value() == 0) {
		return coefficientXToZeroGPE(u, x->operand(0));
	}

	if(x->kind() == Kind::Power && x->operand(1)->value() == 1) {
		return coefficientXToOneGPE(u, x->operand(0));
	}

	if(x->kind() == Kind::Power) {
		return coefficientXToOneGPE(u, x);
	}

	return coefficientXToOneGPE(u, x);

}

// AST* coefficientGPE(AST* u, AST* x) {
// 	// printf("coefficient\n");
// 	// printf("u %s\n", u->toString().c_str());
// 	// printf("x %s\n", x->toString().c_str());
// 	// assert(
// 	// 	x->kind() == Kind::Power,
// 	// 	"coefficientGPE: 'param(x)=%s' needs to be a power!",
// 	// 	x->toString().c_str()
// 	// );

// 	// assert(
// 	// 	!isConstant(x->operand(0)),
// 	// 	"coefficientGPE: base of 'param(x)=%s' "
// 	// 	"cant be constant!",
// 	// 	x->toString().c_str()
// 	// );
// 	// assert(
// 	// 	x->operand(1)->kind() == Kind::Integer &&
// 	// 	x->operand(1)->value() >= 0,
// 	// 	"coefficientGPE: expoent of 'param(x)=%s' "
// 	// 	"needs to be a non negative integer! ",
// 	// 	x->toString().c_str()
// 	// );



// 	if(u->kind() == Kind::Multiplication) {
// 		bool found = false;

// 		signed long count = 0;

// 		AST* res = new AST(Kind::Multiplication);

// 		for(int i=0; i<u->numberOfOperands(); i++) {
// 			AST* tmp = coefficientGPE(u->operand(i), x);

// 			if(tmp->kind() == Kind::Integer && tmp->value() == 1) {
// 				count++;

// 				assert(
// 					count < 2,
// 					"coefficientGPE: 'arg(u)=%s' have more than one 'arg(x)=%s'! "
// 					"Simplify 'arg(u)=%s' to solve the error!",
// 					u->toString().c_str(),
// 					x->toString().c_str(),
// 					u->toString().c_str()
// 				);

// 				delete tmp;
// 				found = true;
// 				continue;
// 			}

// 			delete tmp;
// 			res->includeOperand(u->operand(i)->deepCopy());
// 		}

// 		if(!found) {
// 			delete res;
// 			return inte(0);
// 		}

// 		if(res->numberOfOperands() == 0) {
// 			delete res;
// 			return inte(1);
// 		}

// 		if(res->numberOfOperands() == 1) {
// 			AST* r = res->operand(0);

// 			res->removeOperand(0L);
// 			delete res;

// 			return r;
// 		}

// 		return res;
// 	}

// 	if(u->kind() == Kind::Addition || u->kind() == Kind::Subtraction) {
// 		AST* res = new AST(u->kind());
		
// 		for(int i=0; i<u->numberOfOperands(); i++) {
// 			AST* coeff = coefficientGPE(u->operand(i), x);
// 			// printf("CO %s\n", u->operand(i)->toString().c_str());
// 			// printf("CO %s\n\n", coeff->operand(i)->toString().c_str());

// 			if(coeff->kind() == Kind::Integer && coeff->value() == 0) {
// 				delete coeff;
// 				continue;
// 			}

// 			res->includeOperand(coeff);
// 		}

// 		if(res->numberOfOperands() == 1) {
// 			AST* r = res->operand(0);
// 			res->removeOperand(0L);
// 			delete res;
// 			return r;
// 		}

// 		if(res->numberOfOperands() == 0) {
// 			delete res;
// 			return inte(0);
// 		}

// 		return res;
// 	}

// 	if(x->operand(1)->value() == 1) {
// 		if(u->match(x->operand(0))) {
// 			return inte(1);
// 		}
// 	}

// 	if(x->operand(1)->value() == 0) {
// 		// if(u->match(x->operand(0))) {
// 		// 	return inte(0);
// 		// }
	
// 		// printf("asdasd %s\n", u->toString().c_str());
// 		// printf("asdasa %s\n", x->toString().c_str());
// 		return u->deepCopy();
// 	}

// 	if(u->match(x)) {
// 		return inte(1);
// 	}

// 	return inte(0);
// }

AST* leadingCoefficientGPE(AST* u, AST* x) {
	// assert(
	// 	!isConstant(x),
	// 	"leadingCoefficientGPE: 'param(x)=%s' "
	// 	"cant be a constant expression",
	// 	x->toString().c_str()
	// );

	AST* po = pow(
		x->deepCopy(),
		degreeGPE(u, x)
	);

	// printf("deg = %s\n", degreeGPE(u, x)->toString().c_str());
	// printf("u = %s\n", u->toString().c_str());
	// printf("x = %s\n", po->toString().c_str());

	AST* lc = coefficientGPE(u, po);
	// printf("lc = %s\n", lc->toString().c_str());

	AST* r = expandAST(lc);

	delete po;
	delete lc;

	return r;
}

std::pair<AST*, AST*> divideGPE(AST* u, AST* v, AST* x) {
	assert(
		isPolynomialGPE(u, {x}),
		"'param(u)=%s' needs to be a "
		"GPE(General Polynomial Expression)! "
		"in 'param(x)=%s'",
		u->toString().c_str(),
		x->toString().c_str()
	);

	assert(
		isPolynomialGPE(v, {x}),
		"'param(v)=%s' needs to be a "
		"GPE(General Polynomial Expression)! "
		"in 'param(x)=%s'",
		v->toString().c_str(),
		x->toString().c_str()
	);

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

		AST* s = div(lcr->deepCopy(), lcv->deepCopy());

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
	
		q = expandAST(q_);
	
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

		r = expandAST(r_);
	
	
		delete r_;
		delete m;
		delete lcr;
		delete s;
	
	
		m = degreeGPE(r, x);

	}

	std::pair<AST*, AST*> res = std::make_pair(expandAST(q), expandAST(r));
	
	delete q;
	delete r;
	delete m;
	delete n;
	delete lcv;

	return res;
}

AST* quotientGPE(AST* u, AST* v, AST* x) {
	std::pair<AST*, AST*> res = divideGPE(u,v,x);
	delete res.second;
	return res.first;
}

AST* remainderGPE(AST* u, AST* v, AST* x) {
	std::pair<AST*, AST*> res = divideGPE(u,v,x);
	delete res.first;
	return res.second;
}

AST* expandGPE(AST* u, AST* v, AST* x, AST* t) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return inte(0);

	std::pair<AST*, AST*> d = divideGPE(u, v, x);

	AST* q = d.first;
	AST* r = d.second;

	AST* exp = add({
		mul({
			t->deepCopy(),
			expandGPE(q, v, x, t)
		}),
		r->deepCopy()
	});

	AST* res = expandAST(exp);

	delete exp;
	delete q;
	delete r;

	return res;
}


AST* gcdGPE(AST* u, AST* v, AST* x) {
	if(
		u->kind() == Kind::Integer && u->value() == 0 &&
		v->kind() == Kind::Integer && v->value() == 0
	) {
		return inte(0);
	}

	AST* U = u->deepCopy();
	AST* V = v->deepCopy();

	while (V->kind() != Kind::Integer ||(V->kind() == Kind::Integer && V->value() != 0)) {
		AST* R = remainderGPE(U, V, x);
		delete U;
		U = V->deepCopy();
		delete V;
		V = R->deepCopy();
		delete R;
	}

	AST* e = mul({div(inte(1), leadingCoefficientGPE(U,x)), U->deepCopy()});

	AST* res = expandAST(e);

	delete e;
	delete U;
	delete V;

	return res;
}

std::vector<AST*> extendedEuclideanAlgGPE(AST* u, AST* v, AST* x) {
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
		std::pair<AST*, AST*> d = divideGPE(U,V,x);
	
		AST* q = d.first;
		AST* r = d.second;
	
		AST* A = sub({ App->deepCopy(), mul({q->deepCopy(), Ap->deepCopy()}) });
		AST* B = sub({ Bpp->deepCopy(), mul({q->deepCopy(), Bp->deepCopy()}) });

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

		delete q;
		delete r;
	}

	AST* c = leadingCoefficientGPE(U, x);

	AST* App__ = div(App->deepCopy(), c->deepCopy());

	AST* App_ = expandAST(App__);
	delete App__;

	delete App;
	App = App_;

	AST* Bpp__ = div(Bpp->deepCopy(), c->deepCopy());

	AST* Bpp_ = expandAST(Bpp__);
	delete Bpp__;

	delete Bpp;
	Bpp = Bpp_;

	AST* U__ = div(U->deepCopy(), c->deepCopy());
	
	AST* U_ = expandAST(U__);
	delete U__;
	
	delete U;
	U = U_;

	delete c;
	delete Ap;
	delete Bp;
	delete V;

	return { U, App, Bpp };
}

AST* algMulInverseAST(AST* v, AST* p, AST* a) {
	// TODO: assert following statements
	// 1. a is a symbol that represents an algebraic number
	// 2. p is a monoic irrecudctible polynomial in Q[a] with degree(p,a) >= 2
	// 3. v is a non zero polynomial in Q(a) with degree(v) < degree(p)
	

	std::vector<AST*> w = extendedEuclideanAlgGPE(v,p,a);
	
	delete w[0];
	delete w[2];

	return w[1];
}

AST* algDivideAST(AST* u, AST* v, AST* p, AST* a) {
	// TODO: assert following statements
	// a is a symbol that represents an algebraic number;
	// p is a monic, irrreducible polynomial in Q[α] with deg(p, α) ≥ 2;
	// u and v are both polynomials in Q(a) with degree < deg(p) and v != 0;
	AST* w = algMulInverseAST(v, p, a);

	AST* e = mul({u->deepCopy(), w->deepCopy()});

	AST* k = expandAST(e);

	AST* r = remainderGPE(k, p, a);

	delete w;
	delete e;
	delete k;

	return r;
}

AST* algCoeffSimp(AST* u, AST* x, AST* p, AST* a) {
	assert(
		x->kind() == Kind::Symbol,
		"algCoeffSimp: 'param(x)=%s' needs to be a symbol",
		x->toString().c_str()
	);


	AST* d = degreeGPE(u, x);

	if(d->value() == 0) {
		delete d;
		return remainderGPE(u, p, a);
	}

	AST* r = new AST(u->kind());

	for(int i=0; i <= d->value(); i++) {
		AST* x_ = pow(x->deepCopy(), inte(i));
		AST* coeff_ = coefficientGPE(u, x_);
		AST* coeff = expandAST(coeff_);
		AST* k = remainderGPE(coeff, p, a);
		
		delete coeff_;
		delete coeff;
		
		r->includeOperand(mul({k, x_}));
	}
	
	AST* res = expandAST(r);
	
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

	AST* q = inte(0);
	AST* r = u->deepCopy();
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
			q->deepCopy(),
			mul({
				s->deepCopy(),
				pow(x->deepCopy(),
				sub({m->deepCopy(), n->deepCopy()}))
			})
		});

	
		delete q;
	
		q = expandAST(q_);
	
		delete q_;

		AST* e = sub({
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
					})
				}),
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
	
		AST* r_ = expandAST(e);
	
		// printf("r' = %s\n", r_->toString().c_str());
	
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
	AST* U = u->deepCopy();
	AST* V = v->deepCopy();

	while(
		V->kind() != Kind::Integer ||
		V->value() != 0
	) {
		AST* R = algPolynomialRemainderAST(U, V, x, p, a);

		delete U;

		U = V->deepCopy();

		delete V;

		V = R->deepCopy();

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
		AST* d_ = div(u->deepCopy(), v->deepCopy());
		AST* d = expandAST(d_);
		delete d_;
		
		if(K->identifier() == "Z") {
			if(d->kind() == Kind::Integer) {
				return list({ d, inte(0) });
			}
			
			delete d;
			
			return list({ inte(0), u->deepCopy() });
		}

		if(K->identifier() == "Q") {
			AST* L_ = list({ d, inte(0) });
			return L_;
		}

		// Undefined integral domain
		return new AST(Kind::Undefined);
	}

	AST* x = first(L);
	AST* r = u->deepCopy();

	AST* m = degreeGPE(r, x);
	AST* n = degreeGPE(v, x);

	AST* q = inte(0);
	
	AST* lcv = leadingCoefficientGPE(v, x);

	while(m->kind() != Kind::MinusInfinity && m->value() >= n->value()) {
		AST* lcr = leadingCoefficientGPE(r, x);
		
		AST* restL = rest(L);

		AST* d = recPolyDiv(lcr, lcv, restL, K);
		
		delete restL;
		
		if(d->operand(1)->kind() != Kind::Integer || d->operand(1)->value() != 0) {
			AST* result = expandAST(q);
			
			delete x;
			delete m;
			delete n;
			delete q;
			delete d;
			delete lcv;
			delete lcr;
			
			return list({ result, r });
		} else {
			AST* c = d->operand(0)->deepCopy();
			
			AST* q_ = add({
				q->deepCopy(),
				mul({
					c->deepCopy(),
					pow(
						x->deepCopy(),
						sub({ m->deepCopy(), n->deepCopy() })
					)
				})
			});

			delete q;
			q = expandAST(q_);
			delete q_;
			
			AST* r_ = sub({
				r->deepCopy(),
				mul({
					c->deepCopy(),
					v->deepCopy(),
					pow(
						x->deepCopy(),
						sub({ m->deepCopy(), n->deepCopy() })
					)
				})
			});

			delete r;
			r = expandAST(r_);
			delete r_;

			delete m;
			m = degreeGPE(r, x);
			delete c;
		}
	
		delete lcr;
		delete d;
	}

	AST* result = expandAST(q);

	delete x;
	delete m;
	delete n;
	delete q;
	delete lcv;

	return list({ result, r });
}

AST* recQuotient(AST* u, AST* v, AST* L, AST* K) {
	AST* r = recPolyDiv(u, v, L, K);
	AST* q = r->operand(0)->deepCopy();
	delete r;
	return q;
}

AST* recRemainder(AST* u, AST* v, AST* L, AST* K) {
	AST* r = recPolyDiv(u, v, L, K);
	AST* q = r->operand(1)->deepCopy();
	delete r;
	return q;
}


AST* pseudoDivision(AST* u, AST* v, AST* x) {
	AST* p = inte(0);
	AST* s = u->deepCopy();
	AST* m = degreeGPE(s, x);
	AST* n = degreeGPE(v, x);


	AST* e_ = add({
		sub({
			m->deepCopy(),
			n->deepCopy()
		}),
		inte(1)
	});
	AST* e = expandAST(e_);
	delete e_;

	AST* zero = inte(0);

	AST* delta = max(e, zero);

	delete e;
	delete zero;
	
	AST* ex = pow(x->deepCopy(), n->deepCopy());
	AST* lcv = coefficientGPE(v, ex);
	delete ex;

	int tal = 0;
	while(m->kind() != Kind::MinusInfinity && m->value() >= n->value()) {
		AST* ex_ = pow(x->deepCopy(), m->deepCopy());
		AST* lcs = coefficientGPE(s, ex_);
		
		delete ex_;

		p = add({
			mul({lcv->deepCopy(), p}),
			mul({lcs->deepCopy(), pow(x->deepCopy(), sub({m->deepCopy(), n->deepCopy()}))}),
		});

		AST* s_ = sub({
			mul({lcv->deepCopy(), s}),
			mul({
				lcs->deepCopy(),
				v->deepCopy(),
				pow(
					x->deepCopy(),
					sub({
						m->deepCopy(),
						n->deepCopy(),
					})
				)
			}),
		});
	
		s = expandAST(s_);
		delete s_;
	
		tal = tal + 1;

		delete m;
		m = degreeGPE(s, x);
		delete lcs;

	}

	AST* resQ = mul({
		pow(
			lcv->deepCopy(),
			sub({delta->deepCopy(), inte(tal)})
		),
		p->deepCopy()
	});

	AST* resR = mul({
		pow(
			lcv->deepCopy(),
			sub({ delta->deepCopy(), inte(tal) })
		),
		s->deepCopy()
	});


	AST* res = list({ expandAST(resQ), expandAST(resR) });
	
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
	AST* q = r->operand(0)->deepCopy();
	delete r;
	return q;
}

AST* pseudoRemainder(AST* u, AST* v, AST* x) {
	AST* r = pseudoDivision(u, v, x);
	AST* q = r->operand(1)->deepCopy();
	delete r;
	return q;
}


bool isGreatherThanZero(AST* k) {
	AST* k_ = expandAST(k);
	
	assert(isConstant(k_), "k is not constant!");
	bool r = false;
	
	if(k_->kind() == Kind::Integer) {
		r = k_->value() > 0;
	} 
	else if(k_->kind() == Kind::Fraction) {
		AST* n = k_->operand(0);
		AST* d = k_->operand(1);

		r = (n->value() > 0 && d->value() > 0) ||
				(n->value() < 0 && d->value() < 0);
	}

	delete k_;
	return r;
}

AST* getNormalizationFactor(AST* u, AST* L, AST* K) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return inte(0);

	if(isConstant(u)) {
		if(isGreatherThanZero(u)) {
			if(K->identifier() == "Z") {
				return inte(1);
				// return u->deepCopy();
			}
			if(K->identifier() == "Q") {
				return pow(u->deepCopy(), inte(-1));
			}

			// undefined integral domain
			return new AST(Kind::Undefined);
		} else {
			if(K->identifier() == "Z") {
				return inte(-1);
				// return mul({ inte(-1), u->deepCopy() });
			}
			if(K->identifier() == "Q") {
				return mul({ inte(-1), pow(u->deepCopy(), inte(-1)) });
			}

			// undefined integral domain
			return new AST(Kind::Undefined);
		}
	}

	if(L->numberOfOperands() == 0)
		return new AST(Kind::Undefined);
	
	AST* lc = leadingCoefficientGPE(u, L->operand(0));
	AST* restL = rest(L);

	AST* c = getNormalizationFactor(lc, restL, K);

	delete restL;
	delete lc;

	return c;
}

AST* normalizePoly(AST* u, AST* L, AST* K) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return inte(0);

	AST* u__ = mul({ getNormalizationFactor(u, L, K), u->deepCopy() });

	AST* u_ = expandAST(u__);

	delete u__;
	
	return u_;
}


AST* polynomialContent(AST* u, AST* x, AST* R, AST* K) {
	if(u->kind() == Kind::Integer && u->value() == 0) {
		return inte(0);
	}

	AST* deg = degreeGPE(u, x);

	AST* p = pow(x->deepCopy(), inte(deg->value()));

	AST* gcd = coefficientGPE(u, p);

	delete p;

	for(int i=deg->value() - 1; i>= 0; i--) {
		AST* p = pow(x->deepCopy(), inte(i));
		AST* coef = coefficientGPE(u, p);
		
		delete p;
		AST* gcd_ = mvPolyGCD(gcd, coef, R, K);

		delete gcd;
		delete coef;
		gcd = gcd_;
	}

	delete deg;

	return gcd;
}

AST* mvPolyGCDRec(AST* u, AST* v, AST* L, AST* K) {

	if(L->numberOfOperands() == 0) {
		if(K->identifier() == "Z") { return integerGCD(u, v); }
		if(K->identifier() == "Q") { return inte(1); }
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
			pp_r = inte(0);
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
	
	AST* k = mul({ d->deepCopy(), pp_u->deepCopy() });
	AST* result = expandAST(k);

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
	// printf("\nmvPolyGCD\n");
	// printf("u = %s\n", u->toString().c_str());
	// printf("v = %s\n", v->toString().c_str());

	AST* gcd = mvPolyGCDRec(u, v, L, K);

	AST* r = normalizePoly(gcd, L, K);

	delete gcd;
	
	return r;
}

AST* leadingMonomial(AST* u, AST* L) {
	if(L->numberOfOperands() == 0) {
		return u->deepCopy();
	}

	AST* x = first(L);
	AST* m = degreeGPE(u, x);
	AST* e = pow(x->deepCopy(), m->deepCopy());
	AST* c = coefficientGPE(u, e);
	

	AST* restL = rest(L);

	AST* r_ = mul({
		pow(x->deepCopy(), m->deepCopy()),
		leadingMonomial(c, restL)
	});

	delete c;
	delete e;
	delete m;
	delete x;
	delete restL;

	AST* r = expandAST(r_);

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
		for(int i=0; i<u->numberOfOperands(); i++) {
			if(u->operand(i)->kind() == Kind::Fraction)
				return false;

			if(
				u->operand(i)->kind() == Kind::Power &&
				u->operand(i)->operand(1)->kind() == Kind::Integer &&
				u->operand(i)->operand(1)->value() < 0
			)	return false;
		}

		return true;
	}

	return false;
}

AST* G(AST* u, AST* v) {
	// maybe wrong
	if(u->kind() == Kind::Addition || u->kind() == Kind::Subtraction) {
		AST* k = new AST(u->kind());

		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* z_ = div(u->operand(i)->deepCopy(), v->deepCopy());
			AST* z = expandAST(z_);
			
			// printf("===========================\n");
			// printf("u[i] = %s\n", u->operand(i)->toString().c_str());
			// printf("v = %s\n", v->toString().c_str());
			// printf("z = %s\n", z->toString().c_str());
			// printf("===========================\n");
			delete z_;
			if(wasSimplified(z)) {
				k->includeOperand(z->deepCopy());
			}

			delete z;
		}

		if(k->numberOfOperands() == 0) {
			delete k;
			return inte(0);
		}
	
		AST* r = expandAST(k);
	
		delete k;
		return r;
	}
	AST* k = new AST(u->kind());

	AST* z_ = div(u->deepCopy(), v->deepCopy());
	AST* z = expandAST(z_);

	if(wasSimplified(z)) {
		k->includeOperand(z->deepCopy());
	}

	if(k->numberOfOperands() == 0) {
		delete k;
		return inte(0);
	}

	AST* r = expandAST(k);
	delete k;
	return r;
}

AST* monomialPolyDiv(AST* u, AST* v, AST* L) {
	AST* q = inte(0);
	AST* r = u->deepCopy();
	AST* vt = leadingMonomial(v, L);
	AST* f = G(r, vt);

	while(f->kind() != Kind::Integer || f->value() != 0) {
		
		AST* q_ = add({ q->deepCopy(), f->deepCopy() });
		delete q;
		
		q = expandAST(q_);
	
		
		delete q_;

		AST* r_ = sub({
			r->deepCopy(),
			mul({ f->deepCopy(), v->deepCopy() })
		});
	
		delete r;
		r = expandAST(r_);
		delete r_;
	
		delete f;
		f = G(r, vt);

	}
	
	AST* l = list({ q->deepCopy(), r->deepCopy() });
	
	delete q;
	delete r;
	delete vt;
	delete f;

	return l;
}


// TODO
// monomialPolyExpansion can be used to rewrite the default operations
// like degreeGPE and coefficientGPE for multivariable polynomials
// reducing them to single variable polynomials

// and them using degreeGPE and coefficientGPE to make the query
// monomialPolyExpansion(a^2*b + 2*a*b^2 + b^3 + 2*a + 2*b + 3, a+b, [a, b], t) -> b*t^2 + 2*t + 3
AST* monomialPolyExpansion(AST* u, AST* v, AST* L, AST* t) {
	if(u->kind() == Kind::Integer && u->value() == 0) {
		return inte(0);
	}

	AST* d = monomialPolyDiv(u, v, L);
	
	AST* q = d->operand(0)->deepCopy();
	AST* r = d->operand(1)->deepCopy();

	AST* res_ = mul({
		t->deepCopy(),
		add({
			monomialPolyExpansion(q,v,L,t),
			r->deepCopy()
		})
	});

	AST* res = expandAST(res_);

	delete res_;
	delete q;
	delete r;
	delete d;

	return res;
}

// monomialPolyRem can be used for simplification, for example
// monomialPolyRem(a*i^3 + b*i^2 + c*i + d, i^2 + 1, [i]) -> -a*i - b + c*i + d
// simplification when i^2 + 1 = 0
// also
//monomialPolyRem(sin^4(x)+sin^3(x)+2*sin^2(x)cos^2(x)+cos^4(x), sin^2(x)+cos^2(x)-1, [cos(x), sin(x)]) -> 1 + sin^3(x)
AST* monomialPolyRem(AST* u, AST* v, AST* L) {
	AST* d = monomialPolyDiv(u,v,L);
	AST* r = d->operand(1)->deepCopy();
	delete d;
	return r;
}

AST* monomialPolyQuo(AST* u, AST* v, AST* L) {
	AST* d = monomialPolyDiv(u,v,L);
	AST* r = d->operand(0)->deepCopy();
	delete d;
	return r;
}




}
