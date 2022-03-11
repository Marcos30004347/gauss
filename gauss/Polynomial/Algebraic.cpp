#include "Algebraic.hpp"
#include "gauss/Algebra/List.hpp"
#include "gauss/Algebra/Algebra.hpp"

using namespace ast;
using namespace algebra;

namespace polynomial
{

Expr algMulInverse(Expr v, Expr p, Expr a) {
	// TODO: assert following statements
	// 1. a is a symbol that represents an algebraic number
	// 2. p is a monoic irrecudctible polynomial in Q[a] with degree(p,a) >= 2
	// 3. v is a non zero polynomial in Q(a) with degree(v) < degree(p)

	Expr w = extendedEuclideanAlgGPE(v,p,a);
	Expr r = w[1];

	return r;
}

Expr algDivide(Expr u, Expr v, Expr p, Expr a) {
	// TODO: assert following statements
	// a is a symbol that represents an algebraic number;
	// p is a monic, irrreducible polynomial in Q[α] with deg(p, α) ≥ 2;
	// u and v are both polynomials in Q(a) with degree < deg(p) and v != 0;
	Expr w = algMulInverse(v, p, a);
	Expr e = mul({u, w});
	Expr k = algebraicExpand(e);
	Expr r = remainderGPE(k, p, a);

	return r;
}

Expr algCoeffSimp(Expr u, Expr x, Expr p, Expr a) {
	// assert(
	// 	x.kind() == Kind::Symbol,
	// 	"algCoeffSimp: 'param(x)=%s' needs to be a symbol",
	// 	x->toString().c_str()
	// );

	Expr d = degree(u, x);

	if(d.value() == 0) {
		return remainderGPE(u, p, a);
	}

	Expr r = Expr(u.kind());

	for(int i=0; i <= d.value(); i++) {
		Expr d = integer(i);

		Expr coeff_ = coeff(u, x, d);
		Expr coeff = algebraicExpand(coeff_);
		Expr k = remainderGPE(coeff, p, a);

		r.insert(mul({k, power(x, d)}));
	}

	Expr res = algebraicExpand(r);




	return res;
}

Expr algPolynomialDivision(Expr u, Expr v, Expr x, Expr p, Expr a) {
	// TODO: assert following statements
	// u, v : polynomials in Q(a)[x] with v != 0;
	// x : a symbol;
	// a : a symbol that represents an algebraic number;
	// p : a monic, irrreducible polynomial in Q[a] with degree ≥ 2;

	Expr q = integer(0);
	Expr r = u;
	Expr m = degree(r, x);
	Expr n = degree(v, x);
	Expr lcv = leadCoeff(v, x);

	Expr p_ = deepReplace(p, x, a);

	while(
		m.kind() != Kind::MinusInfinity && (
			m.kind() == Kind::Integer && n.kind() == Kind::Integer &&
			m.value() >= n.value()
		)
	) {

		Expr lcr = leadCoeff(r, x);

		Expr s = algDivide(lcr, lcv, p_, a);

		Expr q_ = add({
			q,
			mul({
				s,
				power(x,
				sub({m, n}))
			})
		});




		q = algebraicExpand(q_);



		Expr e = sub({
			sub({
				r,
				mul({
					lcr,
					power(x, m)
				})
			}),
			mul({
				sub({
					v,
					mul({
						lcv,
						power(x, n)
					})
				}),
				s,
				power(
					x,
					sub({
						m,
						n
					})
				)
			})
		});

		Expr r_ = algebraicExpand(e);




		r = algCoeffSimp(r_, x, p_, a);





		m = degree(r, x);



	}






	return list({ q, r });
}

Expr algPolynomialRemainder(Expr u, Expr v, Expr x, Expr p, Expr a) {
	Expr res = algPolynomialDivision(u,v,x,p,a);
	Expr r = res[1];

	return r;
}

Expr algPolynomialQuotient(Expr u, Expr v, Expr x, Expr p, Expr a) {
	Expr res = algPolynomialDivision(u,v,x,p,a);
	Expr r = res[0];

	return r;
}

Expr algPolynomialGCD(Expr u, Expr v, Expr x, Expr p, Expr a) {
	Expr U = u;
	Expr V = v;

	while(
		V.kind() != Kind::Integer ||
		V.value() != 0
	) {
		Expr R = algPolynomialRemainder(U, V, x, p, a);



		U = V;



		V = R;


	}

	Expr r = algMonic(U, x, p, a);




	return r;
}

Expr algMonic(Expr u,Expr x, Expr p,Expr a) {
	Expr lc = leadCoeff(u, x);
	Expr k_ = algPolynomialQuotient(u, lc, x, p, a);

	return k_;
}
}
