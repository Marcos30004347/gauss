#include "Trigonometry.hpp"
#include "gauss/Simplification/Division.hpp"
#include "gauss/Trigonometry/Trigonometry.hpp"
#include "gauss/Rational/Rational.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace trigonometry;

namespace simplification {

Expr reduceTrigonometricAST(Expr&& u) {
	Expr v 	= substituteTrig(u);
	Expr w 	= rationalize(v);
	Expr n_ = numerator(w);
	Expr d_ = denominator(w);
	Expr n 	= expandTrig(n_);
	Expr d 	= expandTrig(d_);
	Expr k 	= contractTrig(d);

	if(k.kind() == Kind::Integer && k.value() == 0) {
		return undefined();
	}

	return reduceDivisionAST(n / k);
}

Expr reduceTrigonometricAST(Expr& u) {
	return reduceTrigonometricAST(std::forward<Expr>(u));
}

}
