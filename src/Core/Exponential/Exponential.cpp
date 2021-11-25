#include "Exponential.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Polynomial/Polynomial.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace polynomial;
using namespace simplification;

namespace exponential {

Expr expandExponentialRules(Expr A) {
	if(A.kind() == Kind::Addition) {
		Expr f = A[0];
	
		Expr k = sub({
			A,
			f
		});

		Expr a_ = funCall(
			"exp", {
				f
			}
		);

		Expr b_ = funCall(
			"exp", {
				reduceAST(k)
			}
		);

		Expr r_ =  mul({
			expandExponential(a_),
			expandExponential(b_),
		});
	
		Expr r = reduceAST(r_);
	
		
		
		
		
	
		if(isDivisionByZero(r)) {
			
			return undefined();
		}

		return r;
	}

	if(A.kind() == Kind::Multiplication) {
		Expr f = A[0];

		if(f.kind() == Kind::Integer || f.kind() == Kind::Symbol) {
			Expr k = div(A, f);
			
			Expr p_ = power(
				funCall(
					"exp", {
						reduceAST(k)
					}
				),
				f
			);

			Expr p = reduceAST(p_);

			
			

			if(isDivisionByZero(p)) {
				
				return undefined();
			}
	
			return p;
		}
	}

	return funCall("exp", { A });
}

Expr expandExponential(Expr u) {
	if(
		u.kind() == Kind::Integer ||
		u.kind() == Kind::Fraction ||
		u.kind() == Kind::Symbol
	) return u;

	Expr u_ = algebraicExpand(u);
	
	Expr v = mapUnaryAST(u_, expandExponential);
	
	

	if(
		v.kind() == Kind::FunctionCall &&
		v.funName() == "exp"
	) {
		Expr r = expandExponentialRules(v[0]);
		
		return r;
	}

	if(isDivisionByZero(v)) {
		
		return undefined();
	}

	return v;
}

Expr contractExponentialRules(Expr u) {
	Expr v = algebraicExpandRoot(u);

	if(v.kind() == Kind::Power) {
		Expr b = v[0];
		Expr s = v[1];

		if(b.kind() == Kind::FunctionCall && b.funName() == "exp") {
			Expr p = mul({
				b[0],
				s
			});
	
			if(
				p.kind() == Kind::Multiplication || 
				p.kind() == Kind::Power
			) {
				Expr p_ = contractExponentialRules(p);
				
				p = p_;
				

				Expr r = funCall("exp", { p });
				return r;
			}

			
			return v;
		}
	}

	if(v.kind() == Kind::Multiplication) {
		Expr p = integer(1);
		Expr s = integer(0);

		for(unsigned int i=0; i<v.size(); i++) {
			Expr y = v[i];
			if(y.kind() == Kind::FunctionCall && y.funName() == "exp") {
				s = add({
					s, y[0]
				});
			} else {
				p = mul({
					p,
					y
				});
			}
		}
	
		Expr r_;
	
		if(s.kind() == Kind::Integer && s.value() == 0) {
			
			r_ = p;
		} else {
			r_ = mul({ funCall("exp", {s}), p });
		}
	
		Expr r = reduceAST(r_);
		
		
		
		
		return r;
	}

	if(v.kind() == Kind::Addition) {
		Expr s = integer(0);
		for(unsigned int i=0; i<v.size(); i++) {
			Expr y = v[i];
			if(
				y.kind() == Kind::Multiplication ||
				y.kind() == Kind::Power 
			) {
				s = add({
					s, contractExponentialRules(y)
				});
			} else {
				s = add({
					s,
					y
				});
			}
		}

		Expr r = reduceAST(s);

		
		
		
		return r;
	}

	return v;
}

Expr contractExponential(Expr u) {
	if(
		u.kind() == Kind::Integer ||
		u.kind() == Kind::Fraction ||
		u.kind() == Kind::Symbol
	) return u;

	Expr v_ = mapUnaryAST(u, contractExponential);
	Expr v = reduceAST(v_);
	

	if(
		v.kind() == Multiplication ||
		v.kind() == Power
	) {
		Expr t = contractExponentialRules(v);
		
		return t;
	}

	return v;
}

}
