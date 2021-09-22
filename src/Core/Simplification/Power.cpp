#include "Power.hpp"
#include "Rationals.hpp"
#include "Multiplication.hpp"

#include <cstdio>

using namespace ast;
using namespace algebra;

namespace simplification {

AST* reduceIntegerPower(AST* v, AST* n) {
	//(a/b)^n
	if(v->kind() == Kind::Fraction || v->kind() == Kind::Integer) {
		AST* r = power(v->copy(), n->copy());
		AST* r_simp = reduceRNEAST(r);
		destroyASTs({ r });
		return r_simp;
	}

	// (v^0) = 1
	if(n->kind() == Kind::Integer && n->value() == 0)
		return integer(1);

	// (n^1) = n
	if(n->kind() == Kind::Integer && n->value() == 1)
		return v->copy();
	
	// (n^i)^k = n^(i*k)
	if(v->kind() == Kind::Power) {
		AST* r = v->operand(0);
		AST* s = v->operand(1);

		AST* p = mul({ n->copy(), s->copy() });

		AST* p_simp = reduceMultiplicationAST(p);

		delete p;

		if(p_simp->kind() == Kind::Integer) {
			AST* z = reduceIntegerPower(r, p_simp);
			
			delete p_simp;
			
			return z;
		}

		AST* v = power(r->copy(), p_simp);

		return v;
	}

	// (a*b*c)^n = a^n * b^n * c^n
	if(v->kind() == Kind::Multiplication) {
		AST* r = mapBinaryAST(v, n, reduceIntegerPower);
		AST* r_simp = reducePowerAST(r);
		destroyASTs({r});
		return r_simp;
	}

	return power(v->copy(), n->copy());
}

AST* reducePowerAST(AST* u) {
	
	AST* v = base(u);
	AST* n = expoent(u);

	if(v->kind() == Kind::Undefined) {
		destroyASTs({v, n});
		return undefined();
	}

	if(n->kind() == Kind::Undefined) {
		destroyASTs({v, n});
		return undefined();
	}

	if(n->kind() == Kind::Infinity) {
		destroyASTs({v, n});
		return new AST(Kind::Infinity);
	}

	if(n->kind() == Kind::MinusInfinity) {
		destroyASTs({v, n});
		return integer(0);
	}

	if(v->kind() == Kind::Integer && v->value() == 0) {

		if(n->kind() == Kind::Integer && n->value() > 0) {
			destroyASTs({v, n});
			return integer(0);
		}
		else {
			destroyASTs({v, n});
			return undefined();
		}
	}

	if(v->kind() == Kind::Integer && v->value() == 1) {
		destroyASTs({v, n});
		return integer(1);
	}

	if(n->kind() == Kind::Integer) {
		AST* res = reduceIntegerPower(v, n);
		destroyASTs({v, n});
		return res;
	}
	destroyASTs({v, n});
	return u->copy();
}


}
