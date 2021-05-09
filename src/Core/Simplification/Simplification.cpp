#include "Simplification.hpp"
#include "Addition.hpp"
#include "Subtraction.hpp"
#include "Division.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"
#include "Rationals.hpp"
#include "Factorial.hpp"

#include <cstdio>

using namespace ast;
using namespace algebra;

namespace simplification {

AST* reduceAST(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Symbol ||
		u->kind() == Kind::Infinity ||
		u->kind() == Kind::MinusInfinity ||
		u->kind() == Kind::Undefined ||
		u->kind() == Kind::List ||
		u->kind() == Kind::Set
	) return u->deepCopy();


	if(u->kind() == Kind::Fraction)
		return reduceRNEAST(u);


	AST* v = mapUnaryAST(u, reduceAST);

	if(v->kind() == Kind::Addition) {
		AST* res1 = reduceAdditionAST(v);
		delete v;
		return res1;
	}
	if(v->kind() == Kind::Subtraction) {
		
		AST* res2 = reduceSubtractionAST(v);
		delete v;
		return res2;
	}
	if(v->kind() == Kind::Division) {

		AST* res3 = reduceDivisionAST(v);
		delete v;
		return res3;
	}
	if(v->kind() == Kind::Multiplication) {

		AST* res4 = reduceMultiplicationAST(v);
		delete v;
		return res4;
	}
	if(v->kind() == Kind::Power) {

		AST* res5 = reducePowerAST(v);
		delete v;
		return res5;
	}
	if(v->kind() == Kind::Factorial) {

		AST* res6 = reduceFactorialAST(v);
		delete v;
		return res6;
	}

	return v;
}


}
