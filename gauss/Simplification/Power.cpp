#include "Power.hpp"
#include "gauss/AST/AST.hpp"
#include "Rationals.hpp"
#include "Multiplication.hpp"

#include <cstdio>
#include <tuple>

using namespace ast;
using namespace algebra;

namespace simplification {

Expr reduceIntegerPower(Expr& v, Expr& n) {
	//(a/b)^n
	if(v.kind() == Kind::Fraction || v.kind() == Kind::Integer) {
		return reduceRNEAST(power(v, n));
	}

	// (v^0) = 1
	if(n == 0)
		return integer(1);

	// (n^1) = n
	if(n == 1)
		return v;

	// (n^i)^k = n^(i*k)
	if(v.kind() == Kind::Power) {
		Expr r = v[0];
		Expr s = v[1];

		Expr p = mul({ n, s });

		Expr p_simp = reduceMultiplicationAST(p);

		if(p_simp.kind() == Kind::Integer) {
			return reduceIntegerPower(r, p_simp);
		}

		return power(r, p_simp);
	}

	// (a*b*c)^n = a^n * b^n * c^n
	if(v.kind() == Kind::Multiplication) {
		Expr r = mapBinaryAST(v, n, reduceIntegerPower);
		Expr r_simp = reducePowerExpr(r);
		return r_simp;
	}

	return power(v, n);
}

Expr reducePowerExpr(Expr&& u) {
	Expr v = u[0];
	Expr n = u[1];

	if(v.kind() == Kind::Undefined) {
		return undefined();
	}

	if(n == undefined()) {
		return n;
	}

	if(n == inf()) {
		return n;
	}

	if(n == -inf()) {
		return 0;
	}

	if(v == 0)
	{
		if(n > 0) {
			return integer(0);
		}
		else {
			return undefined();
		}
	}

	if(v == 1) {
		return 1;
	}

	if(n.kind() == Kind::Integer) {
		Expr res = reduceIntegerPower(v, n);
		return res;
	}

	return Expr(u);
}

Expr reducePowerExpr(Expr& u) {
	return reducePowerExpr(std::forward<Expr>(u));
}



}
