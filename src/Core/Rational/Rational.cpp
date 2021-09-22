#include "Rational.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Algebra/Set.hpp"


using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace rational {

bool isRationalExpression(AST* u, AST* S) {
	AST* n = numerator(u);
	AST* d = denominator(u);
	
	bool r = isGerenalPolynomial(n, S) &&
			 isGerenalPolynomial(d, S);

	delete n;
	delete d;

	return r;
}

AST* rationalVariables(AST* u) {
	AST* n = numerator(u);
	AST* d = denominator(u);

	AST* K = variables(n);
	AST* J = variables(d);

	AST* R = unification(K, J);

	delete n;
	delete d;
	delete K;
	delete J;

	return R;
}

AST* rationalizeSum(AST* u, AST* v) {
	AST* m = numerator(u);
	AST* r = denominator(u);
	AST* n = numerator(v);
	AST* s = denominator(v);

	if(
		r->kind() == Kind::Integer && r->value() == 1 &&
		s->kind() == Kind::Integer && s->value() == 1
	) {

		delete r;
		delete s;
		delete m;
		delete n;

		AST* t = add({ u->copy(), v->copy() });

		AST* k = reduceAST(t);

		delete t;

		return k;
	}

	AST* num_a 	= mul({ m->copy(), s->copy() });
	AST* num_b 	= mul({ n->copy(), r->copy() });
	AST* den 		= mul({ r->copy(), s->copy() });

	AST* num 	= rationalizeSum(num_a, num_b);
	
	AST* o 	= div(
		num,
		den
	);
	

	delete m;
	delete n;
	delete s;
	delete r;
	delete num_a;
	delete num_b;

	return o;
}

AST* rationalize(AST* u) {

	if(u->kind() == Kind::Power) {
		return power(
			rationalize(u->operand(0)),
			reduceAST(u->operand(1))
		);
	}

	if(u->kind() == Kind::Multiplication) {
		AST* f = u->operand(0);
		
		AST* k_ = div(u->copy(), f->copy());
		AST* k = reduceAST(k_);

		delete k_;
	
		AST* r = mul({
			rationalize(f),
			rationalize(k)
		});

		delete k;

		return r;
	}

	if(u->kind() == Kind::Addition) {
		AST* f = u->operand(0);
	
		AST* k_ = sub({ u->copy(), f->copy() });
		AST* k = reduceAST(k_);

		delete k_;
	
		AST* g = rationalize(f);
		AST* r = rationalize(k);

		delete k;
	
		AST* t = rationalizeSum(g, r);
	
		delete g;
		delete r;

		return t;
	}

	return u->copy();
}


AST* numerator(AST* u) {
	if(u->kind() == Kind::Fraction || u->kind() == Kind::Division)
		return u->operand(0)->copy();
	
	if(u->kind() == Kind::Power) {
		if(u->operand(1)->kind() == Kind::Integer && u->operand(1)->value() < 0) {
			return integer(1);
		}
		return u->copy();
	}

	if(u->kind() == Kind::Multiplication) {
		if(u->numberOfOperands() == 1) {
			return numerator(u->operand(0));
		}
		AST* v = u->operand(0);
		
		AST* h_ = div(
			u->copy(),
			v->copy()
		);
	
		AST* h = reduceAST(h_);
	
		AST* r_ = mul({
			numerator(v),
			numerator(h)
		});

		AST* r = reduceAST(r_);

		delete h;
		delete h_;
		delete r_;

		return r;
	}

	return u->copy();
}

AST* denominator(AST* u) {
	if(u->kind() == Kind::Fraction || u->kind() == Kind::Division)
		return u->operand(1)->copy();
	
	if(u->kind() == Kind::Power) {
		if(u->operand(1)->kind() == Kind::Integer && u->operand(1)->value() < 0) {
			AST* e = power(u->copy(), integer(-1));
		
			AST* r = reduceAST(e);
	
			delete e;
	
			return r;
		}
	
		return integer(1);
	}

	if(u->kind() == Kind::Multiplication) {
		if(u->numberOfOperands() == 1) {
			return denominator(u->operand(0));
		}

		AST* v = u->operand(0);
		
		AST* h_ = div(
			u->copy(),
			v->copy()
		);
	
		AST* h = reduceAST(h_);
	
		AST* r_ = mul({
			denominator(v),
			denominator(h)
		});

		AST* r = reduceAST(r_);

		delete h;
		delete h_;
		delete r_;

		return r;
	}

	return integer(1);
}

AST* expandRational(AST* u) {
	AST* n_ = numerator(u);
	AST* d_ = denominator(u);

	AST* n = algebraicExpand(n_);
	AST* d = algebraicExpand(d_);

	AST* k = div(n, d);
	
	delete n_;
	delete d_;

	return k;
}
 
}
