#include "Trigonometry.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Polynomial/Polynomial.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace polynomial;
using namespace simplification;

namespace simplification {

AST* substituteTrig(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* g = mapUnaryAST(u, substituteTrig);

	if(g->kind() == Kind::FunctionCall) {
		
		if(g->funName() == "tan") {
			AST* k = div(
				funCall("sin", { g->operand(0)->deepCopy() }),
				funCall("cos", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "cot") {
			AST* k = div(
				funCall("cos", { g->operand(0)->deepCopy() }),
				funCall("sin", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "sec") {
			AST* k = div(
				integer(1),
				funCall("cos", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "csc") {
			AST* k = div(
				integer(1),
				funCall("sin", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}
	}

	return g;
}

bool isDivisionByZero(AST* k) {
	AST* d = denominator(k);
	if(d->kind() == Kind::Integer && d->value() == 0) {
		delete d;
		return true;
	}

	delete d;
	return false;
}

AST* expandExponentialRules(AST* A) {
	if(A->kind() == Kind::Addition) {
		AST* f = A->operand(0);
		
		AST* k = sub({
			A->deepCopy(),
			f->deepCopy()
		});

		AST* a_ = funCall(
			"exp", {
				f->deepCopy()
			}
		);

		AST* b_ = funCall(
			"exp", {
				reduceAST(k)
			}
		);

		AST* r_ =  mul({
			expandExponential(a_),
			expandExponential(b_),
		});
	
		AST* r = reduceAST(r_);
	
		delete k;
		delete a_;
		delete b_;
		delete r_;
	
		if(isDivisionByZero(r)) {
			delete r;
			return undefined();
		}

		return r;
	}

	if(A->kind() == Kind::Multiplication) {
		AST* f = A->operand(0);

		if(f->kind() == Kind::Integer) {
			AST* k = div(A->deepCopy(), f->deepCopy());
			
			AST* p_ = power(
				funCall(
					"exp", {
						reduceAST(k)
					}
				),
				f->deepCopy()
			);

			AST* p = reduceAST(p_);

			delete p_;
			delete k;

			if(isDivisionByZero(p)) {
				delete p;
				return undefined();
			}
	
			return p;
		}
	}

	return funCall("exp", { A->deepCopy() });
}

AST* expandExponential(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* u_ = algebraicExpand(u);
	AST* v = mapUnaryAST(u_, expandExponential);
	delete u_;

	if(
		v->kind() == Kind::FunctionCall &&
		v->funName() == "exp"
	) {
		AST* r = expandExponentialRules(v->operand(0));
		delete v;
		return r;
	}

	if(isDivisionByZero(v)) {
		delete v;
		return undefined();
	}

	return v->deepCopy();
}

signed long integer_fact(signed long i) {
	if(i == 1 || i == 0)
		return 1;

	return i * integer_fact(i - 1);
} 

AST* multipleAndlgeCos(AST* n, AST* theta) {
	AST* r = integer(0);

	for(int j=0; j <= n->value(); j++) {
		if(j%2 != 0)
			continue;
		
		signed long b = integer_fact(n->value())/(integer_fact(j) * integer_fact(n->value() - j));
		AST* c_ = funCall("cos", {
			theta->deepCopy()
		});

		AST* s_ = funCall("sin", {
			theta->deepCopy()
		});
		
		AST* c = c_;
		AST* s = s_;

		if(theta->kind() == Kind::Addition) {
			AST* e_s = expandTrigRules(s_);
			s = e_s->operand(0)->deepCopy();
			AST* e_c = expandTrigRules(c_);
			c = e_s->operand(1)->deepCopy();
			
			delete s_;
			delete e_s;
			delete c_;
			delete e_c;
		}


		AST* e = mul({
			power(integer(-1), div(integer(j), integer(2))),
			integer(b),
			power(c, sub({ n->deepCopy(), integer(j) })),
			power(s, integer(j))
		});
	
		r = add({r, e});
	}

	AST* expanded = reduceAST(r);

	delete r;

	return expanded;
}


AST* multipleAndlgeSin(AST* n, AST* theta) {
	AST* r = integer(0);

	for(int j=0; j <= n->value(); j++) {
		if(j%2 != 1)
			continue;
		
		signed long b = integer_fact(n->value())/(integer_fact(j) * integer_fact(n->value() - j));
		AST* c_ = funCall("cos", {
			theta->deepCopy()
		});

		AST* s_ = funCall("sin", {
			theta->deepCopy()
		});
		
		AST* c = c_;
		AST* s = s_;

		if(theta->kind() == Kind::Addition) {
			AST* e_s = expandTrigRules(s_);
			s = e_s->operand(0)->deepCopy();
			AST* e_c = expandTrigRules(c_);
			c = e_s->operand(1)->deepCopy();
			
			delete s_;
			delete e_s;
			delete c_;
			delete e_c;
		}


		AST* e = mul({
			power(integer(-1), div( sub({ integer(j), integer(1) }), integer(2))),
			integer(b),
			power(c, sub({ n->deepCopy(), integer(j) })),
			power(s, integer(j))
		});
	
		r = add({r, e});
	}

	AST* expanded = reduceAST(r);

	delete r;

	return expanded;
}


AST* expandTrigRules(AST* A) {
	if(A->kind() == Kind::Addition) {
		AST* f = expandTrigRules(A->operand(0));
		AST* A__ = sub({
			A->deepCopy(),
			A->operand(0)->deepCopy()
		});
		
		AST* A_ = reduceAST(A__);
		
		delete A__;

		AST* r = expandTrigRules(A_);
		delete A_;

		AST* s_ = add({
			mul({
				f->operand(0)->deepCopy(),
				r->operand(1)->deepCopy(),
			}),
			mul({
				f->operand(1)->deepCopy(),
				r->operand(0)->deepCopy(),
			}),
		});

		AST* s = reduceAST(s_);
		delete s_;

		AST* c_ = sub({
			mul({
				f->operand(1)->deepCopy(),
				r->operand(1)->deepCopy(),
			}),
			mul({
				f->operand(0)->deepCopy(),
				r->operand(0)->deepCopy(),
			}),
		});

		AST* c = reduceAST(s_);
		delete c_;

		delete f;
		delete r;

		return list({s,c});
	}

	if(A->kind() == Kind::Multiplication) {
		AST* f = A->operand(0);
		if(f->kind() == Kind::Integer) {
			AST* k_ = div(
				A->deepCopy(),
				f->deepCopy()
			);
			AST* k = reduceAST(k_);
			delete k_;

			AST* a = multipleAndlgeSin(f, k);
			AST* b = multipleAndlgeCos(f, k);
			delete k;

			return list({a, b});
		}
	}

	return list({
		funCall("sin", {A->deepCopy()}),
		funCall("cos", {A->deepCopy()}),
	});
}

AST* expandTrig(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* u_ = algebraicExpand(u);
	AST* v = mapUnaryAST(u_, expandExponential);
	delete u_;

	if(
		v->kind() == Kind::FunctionCall &&
		v->funName() == "sin"
	) {
		AST* a_ = expandTrigRules(v->operand(0));
		AST* a = reduceAST(a_);

		AST* r = a->operand(0)->deepCopy();
		
		delete a;
		delete a_;

		if(isDivisionByZero(r)) {
			delete r;
			return undefined();
		}

		return r;
	}

	if(
		v->kind() == Kind::FunctionCall &&
		v->funName() == "cos"
	) {
		AST* a_ = expandTrigRules(v->operand(0));
		AST* a = reduceAST(a_);

		AST* r = a->operand(1)->deepCopy();

		delete a;
		delete a_;

		if(isDivisionByZero(r)) {
			delete r;
			return undefined();
		}

		return r;
	}

	if(isDivisionByZero(v)) {
		delete v;
		return undefined();
	}

	return v->deepCopy();
}

}
