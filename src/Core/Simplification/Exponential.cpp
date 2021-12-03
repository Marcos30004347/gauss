#include "Exponential.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Simplification/Division.hpp"
#include "Core/Exponential/Exponential.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Simplification/Simplification.hpp"
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
