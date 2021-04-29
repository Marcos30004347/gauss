
#include <assert.h>

#include "Expand.hpp"
#include "Core/Reduce/Reduce.hpp"
#include "Division.hpp"
#include "Polynomial.hpp"
#include "Multiplication.hpp"

using namespace ast;
using namespace reduce;
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

	// if(!k->match(u)) {
	// 	printf("* %s → %s\n", u->toString().c_str(), k->toString().c_str());
	// 	printf("\n");
	// }

	AST* k_ = reduceAST(k);

	// if(!k->match(k_)) {
	// 	printf("* %s → %s\n", k->toString().c_str(), k_->toString().c_str());
	// 	printf("\n");
	// }

	delete k;
	k = k_;


	if(k->kind() == Kind::Power) {
		AST* k_ = expandMultinomialAST(k);
		// if(!k->match(k_)) {
		// 	printf("* %s → %s\n", k->toString().c_str(), k_->toString().c_str());
		// 	printf("\n");
		// }
		delete k;
		k = k_;
	} else if(k->kind() == Kind::Multiplication) {
		AST* k_ = expandMultiplicationAST(k); 
		// if(!k->match(k_)) {
		// 	printf("* %s → %s\n", k->toString().c_str(), k_->toString().c_str());
		// 	printf("\n");
		// }
		delete k;
		k = k_;
	} else if(k->kind() == Kind::Division) {
		AST* k_ = expandDivisionAST(k);
		// if(!k->match(k_)) {
		// 	printf("* %s → %s\n", k->toString().c_str(), k_->toString().c_str());
		// 	printf("\n");
		// }
		delete k;
		k = k_;
	}

	AST* res = reduceAST(k);

	// if(!k->match(res)) {
	// 	printf("* %s → %s\n", k->toString().c_str(), res->toString().c_str());
	// 	printf("\n");
	// }
	
	delete k;

	return res;
}


}
