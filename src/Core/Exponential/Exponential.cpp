#include "Exponential.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Polynomial/Polynomial.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace polynomial;
using namespace simplification;

namespace exponential {

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

		if(f->kind() == Kind::Integer || f->kind() == Kind::Symbol) {
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

	return v;
}

AST* contractExponentialRules(AST* u) {
	AST* v = algebraicExpandRoot(u);

	if(v->kind() == Kind::Power) {
		AST* b = v->operand(0);
		AST* s = v->operand(1);

		if(b->kind() == Kind::FunctionCall && b->funName() == "exp") {
			AST* p = mul({
				b->operand(0)->deepCopy(),
				s->deepCopy()
			});
	
			if(
				p->kind() == Kind::Multiplication || 
				p->kind() == Kind::Power
			) {
				AST* p_ = contractExponentialRules(p);
				delete p;
				p = p_;
				delete v;

				AST* r = funCall("exp", { p });
				return r;
			}

			delete p;
			return v;
		}
	}

	if(v->kind() == Kind::Multiplication) {
		AST* p = integer(1);
		AST* s = integer(0);

		for(int i=0; i<v->numberOfOperands(); i++) {
			AST* y = v->operand(i);
			if(y->kind() == Kind::FunctionCall && y->funName() == "exp") {
				s = add({
					s, y->operand(0)->deepCopy()
				});
			} else {
				p = mul({
					p,
					y->deepCopy()
				});
			}
		}
	
		AST* r_;
	
		if(s->kind() == Kind::Integer && s->value() == 0) {
			delete s;
			r_ = p;
		} else {
			r_ = mul({ funCall("exp", {s}), p });
		}
	
		AST* r = reduceAST(r_);
		
		delete v;
		delete r_;
		
		return r;
	}

	if(v->kind() == Kind::Addition) {
		AST* s = integer(0);
		for(int i=0; i<v->numberOfOperands(); i++) {
			AST* y = v->operand(i);
			if(
				y->kind() == Kind::Multiplication ||
				y->kind() == Kind::Power 
			) {
				s = add({
					s, contractExponentialRules(y)
				});
			} else {
				s = add({
					s,
					y->deepCopy()
				});
			}
		}

		AST* r = reduceAST(s);

		delete v;
		delete s;
		
		return r;
	}

	return v;
}

AST* contractExponential(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* v_ = mapUnaryAST(u, contractExponential);
	AST* v = reduceAST(v_);
	delete v_;

	if(
		v->kind() == Multiplication ||
		v->kind() == Power
	) {
		AST* t = contractExponentialRules(v);
		delete v;
		return t;
	}

	return v;
}

}
