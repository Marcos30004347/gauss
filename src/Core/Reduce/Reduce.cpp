#include "Reduce.hpp"
#include "Addition.hpp"
#include "Subtraction.hpp"
#include "Division.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"
#include "Rationals.hpp"

#include <cstdio>

using namespace ast;
using namespace algebra;

namespace reduce {

AST* reduceAST(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* v = mapUnaryST(u, reduceAST);

	switch (v->kind()) {
	case Kind::Fraction:
		AST* res = reduceRNEAST(v);
		delete v;
		return res;
	case Kind::Addition:
		AST* res = reduceAdditionAST(v);
		delete v;
		return res;
	case Kind::Subtraction:
		AST* res = reduceSubtractionAST(v);
		delete v;
		return res;
	case Kind::Division:
		AST* res = reduceDivisionAST(v);
		delete v;
		return res;
	case Kind::Multiplication:
		AST* res = reduceMultiplicationAST(v);
		delete v;
		return res;
	case Kind::Power:
		AST* res = reducePowerAST(v);
		delete v;
		return res;
	}

	return v;
}

// AST* reduce(AST* u) {
//     if (u == nullptr)
// 		return nullptr;

// 	if(
// 		u->kind() == Kind::Integer || 
// 		u->kind() == Kind::Symbol
// 	) return reduceAdditionAST(u);

// 	std::list<AST*> tems;

// 	for(int i=0; i<u->numberOfOperands(); i++) {
// 		tems.insert(tems.end(), reduceAST(u->operand(i)));
// 	}

// 	AST* reduced = new AST(u->kind());

//     cout << node->data << " ";
// }

}
