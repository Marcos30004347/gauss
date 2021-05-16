#include "Trigonometry.hpp"
#include "Core/Simplification/Division.hpp"
#include "Core/Trigonometry/Trigonometry.hpp"
#include "Core/Rational/Rational.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace trigonometry;

namespace simplification {

AST* reduceTrigonometricAST(AST* u) {
	AST* v 	= substituteTrig(u);
	AST* w 	= rationalize(v);
	AST* n_ = numerator(w);
	AST* d_ = denominator(w);
	AST* n 	= expandTrig(n_);
	AST* d 	= expandTrig(d_);
	AST* k 	= contractTrig(d);
	
	delete v;
	delete w;
	delete d;
	delete n_;
	delete d_;
	
	if(k->kind() == Kind::Integer && k->value() == 0) {
		delete k;
		delete n;
		return undefined();
	}

	AST* r_ = div(n, k);
	AST* r = reduceDivisionAST(r);

	delete r_;

	return r;
}

}
