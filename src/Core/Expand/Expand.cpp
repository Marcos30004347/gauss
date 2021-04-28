
#include <assert.h>

#include "Expand.hpp"
#include "Division.hpp"
#include "Polynomial.hpp"
#include "Multiplication.hpp"

using namespace ast;
using namespace algebra;

namespace expand {

/**
 * Tink aboud adding the following expansions:
 * 
 * sin(x+y) = sin(x)cos(y)+cos(x)sin(y)
 * cos(2x) = 2cos(x)^2 - 1
 * e^(a+ln(b)) = (e^a)b
 */
AST* expandAST(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction || 
		u->kind() == Kind::FunctionCall ||
		u->kind() == Kind::Undefined ||
		u->kind() == Kind::Symbol ||
		u->kind() == Kind::Infinity
	) {
		return u->deepCopy();
	}

	AST* k = mapUnaryAST(u, expandAST);

	if(k->kind() == Kind::Power) {
		AST* k_ = expandMultinomialAST(k);
		delete k;
		return k_;
	}

	if(k->kind() == Kind::Multiplication) {
		AST* k_ = expandMultiplicationAST(k);
		delete k;
		return k_;
	}

	if(k->kind() == Kind::Division) {
		AST* k_ = expandDivisionAST(k);
		delete k;
		return k_;
	}

	// TODO: expand polynomial functions 

	return k;
}
}
