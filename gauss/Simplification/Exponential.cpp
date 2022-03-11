#include "Exponential.hpp"
#include "gauss/AST/AST.hpp"
#include "gauss/Simplification/Division.hpp"
#include "gauss/Exponential/Exponential.hpp"
#include "gauss/Rational/Rational.hpp"
#include "gauss/Simplification/Simplification.hpp"
#include <utility>

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace exponential;

namespace simplification {

Expr reduceExponentialAST(Expr&& u) {
	Expr u_ = rationalize(u);
	Expr n_ = numerator(u_);
	Expr d_ = denominator(u_);
	Expr n = contractExponential(n_);
	Expr d = contractExponential(d_);

	if(d.kind() == Kind::Integer && d.value() == 0) {
		return undefined();
	}

	return div(n, d);
}

Expr reduceExponentialAST(Expr& u) {
	return reduceExponentialAST(std::forward<Expr>(u));
}

}
