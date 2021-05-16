#include "Exponential.hpp"
#include "Core/Simplification/Division.hpp"
#include "Core/Trigonometry/Trigonometry.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace trigonometry;

namespace simplification {

AST* reduceExponentialAST(AST* u) {
	AST* u_ = rationalize(u);
	AST* n_ = numerator(u_);
	AST* d_ = denominator(u_);
	AST* n = contractExponential(n_);
	AST* d = contractExponential(d_);
	
	delete u_;
	delete n_;
	delete d_;

	if(d->kind() == Kind::Integer && d->value() == 0) {
		delete n;
		delete d;
		return undefined();
	}

	return div(n, d);
}

}
