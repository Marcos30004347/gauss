#include "Rationals.hpp"
#include "Core/Rational/Rational.hpp"

#include <cstdio>

using namespace ast;
using namespace algebra;
using namespace rational;

namespace simplification {

Expr reduceRationalNumber(Expr u);

Expr evaluateProduct(Expr u, Expr v) {
	if(isConstant(u) && isConstant(v)) {
		if(u.kind() == Kind::Integer && v.kind() == Kind::Integer)
			return integer(u.value() * v.value());

		Expr num_u = numerator(u);	
		Expr den_u = denominator(u);	
		Expr num_v = numerator(v);	
		Expr den_v = denominator(v);	

		Expr f = fraction(
			evaluateProduct(num_u, num_v),
			evaluateProduct(den_u, den_v)  
		);

		Expr simp_f = reduceRNEAST(f);
		
		
		
		return simp_f;
	}

	return mul({u, v});
}

Expr evaluateAddition(Expr u, Expr v) {
	if(u.kind() == Kind::Integer && v.kind() == Kind::Integer)
		return integer(u.value() + v.value());

	// x/y + r = (x+ry)/y
	if(u.kind() == Kind::Fraction && v.kind() == Kind::Integer) {
		Expr num_u 	= numerator(u);
		Expr den_u 	= denominator(u);
		Expr k 			= evaluateProduct(v, den_u);

		Expr f = fraction(
			evaluateAddition(num_u, k),
			denominator(u)
		);
		
		Expr simp_f = reduceRNEAST(f);
		
		
	
		return simp_f;
	}

	// r + x/y= (ry + x)/y
	if(v.kind() == Kind::Fraction && u.kind() == Kind::Integer) {
		Expr num_v 	= numerator(v);
		Expr den_v 	= denominator(v);
		Expr k 			= evaluateProduct(u, den_v);

		Expr f = fraction(
			evaluateAddition(k, num_v),
			denominator(v)
		);

		Expr simp_f = reduceRNEAST(f);
		
		
	
		return simp_f;
	}
	// x/y + a/b = (xb + ay)/ab
	if(v.kind() == Kind::Fraction && u.kind() == Kind::Fraction) {
		Expr num_u 					= numerator(u);
		Expr den_u 					= denominator(u);
		Expr num_v 					= numerator(v);
		Expr den_v 					= denominator(v);
		Expr num_u_t_dun_v 	= evaluateProduct(num_u, den_v);
		Expr num_v_t_dun_u 	= evaluateProduct(num_v, den_u);

		Expr f = fraction(
			evaluateAddition(num_u_t_dun_v, num_v_t_dun_u),
			evaluateProduct(den_u, den_v)
		);
	
		Expr simp_f = reduceRNEAST(f);
		
		
		
		return simp_f;
	}

	return add({ u, v });
}


Expr evaluateSubtraction(Expr u, Expr v) {
	if(u.kind() == Kind::Integer && v.kind() == Kind::Integer)
		return integer(u.value() - v.value());

	// x/y - r = (x-ry)/y
	if(u.kind() == Kind::Fraction && v.kind() == Kind::Integer) {
		Expr num_u 	= numerator(u);
		Expr den_u 	= denominator(u);
		Expr k 			= evaluateProduct(v, den_u);

		Expr f = fraction(
			evaluateSubtraction(num_u, k),
			denominator(u)
		);
	
		Expr simp_f = reduceRNEAST(f);
		
		
	
		return simp_f;
	}

	// r - x/y= (ry - x)/y
	if(v.kind() == Kind::Fraction && u.kind() == Kind::Integer) {
		Expr num_v 	= numerator(v);
		Expr den_v 	= denominator(v);
		Expr k 			= evaluateProduct(u, den_v);

		Expr f = fraction(
			evaluateSubtraction(k, num_v),
			denominator(v)
		);
		
		Expr simp_f = reduceRNEAST(f);
		
		
	
		return simp_f;
	}
	// x/y - a/b = (xb - ay)/ab
	if(v.kind() == Kind::Fraction && u.kind() == Kind::Fraction) {
		Expr num_u 					= numerator(u);
		Expr den_u 					= denominator(u);
		Expr num_v 					= numerator(v);
		Expr den_v 					= denominator(v);
		Expr num_u_t_dun_v 	= evaluateProduct(num_u, den_v);
		Expr num_v_t_dun_u 	= evaluateProduct(num_v, den_u);

		Expr f = fraction(
			evaluateSubtraction(num_u_t_dun_v, num_v_t_dun_u),
			evaluateProduct(den_u, den_v)
		);
	
		Expr simp_f = reduceRNEAST(f);
		
		
		
		return simp_f;
	}

	return sub({ u, v });
}

Expr evaluateDivision(Expr v, Expr w) {
	if(!isConstant(v) || !isConstant(w))
		return div(v, w);

	Expr num_w = numerator(w);
	Expr den_w = denominator(w);

	Expr num_v = numerator(v);
	Expr den_v = denominator(v);

	if(num_w.kind() == Kind::Integer && num_w.value() == 0) {
		
		return undefined();
	}

	// if(num_w.kind() == Kind::Integer && num_w.value() == 1) {
	// 	
	// 	return v;
	// }

	Expr f = fraction(
		evaluateProduct(num_v, den_w),
		evaluateProduct(num_w, den_v)
	);	

	Expr simp_f = reduceRNEAST(f);

	

	return simp_f;
}

Expr evaluatePower(Expr v, Expr w) {
	if(w.kind() != Kind::Integer)
			return power(v ,w);
	
	Int n = w.value();

	if(!isConstant(v))
			return power(v, integer(n));

	Expr num_v = numerator(v);
	if(num_v.kind() == Kind::Integer && num_v.value() != 0) {
		
		if(n > 0) {
			Expr e 		= integer(n - Int(1));
			Expr s 		= evaluatePower(v, e);
			Expr res 	= evaluateProduct(s, v);
			
			
			return res;
		} else if(n == 0) {
				return integer(1);
		} else if(n == -1) {
				Expr f 				= fraction(denominator(v), numerator(v));
				Expr f_simp 	= reduceRNEAST(f);
				
				
				
				return f_simp;
		} else {
				Expr f 			= fraction(denominator(v), numerator(v));
				Expr f_simp = reduceRNEAST(f);
				Expr e 			= integer(-1 * n);
				Expr res 		= evaluatePower(f_simp, e);

				

				return res;
		}
	}

	
		
	if(n >= 1) 
		return integer(0);

	return undefined();
}



Expr simplyfyQuotient(Expr u, Expr v) {
    return integer(u.value() / v.value());
}

Expr reduceRationalNumber(Expr u) {
    if(u.kind() == Kind::Integer) {
		return u;
	}

    if(u.kind() == Kind::Fraction) {

		Expr n = u[0];
		Expr d = u[1];
		if(d.value() == 1) {

				return n;
		}

		if(n.value() % d.value() == 0) {
				return simplyfyQuotient(n, d);
		}
		else {
			Expr g = integerGCD(n, d);

			if(d.value() > 0) {

				if(g.value() > 0) {
					Expr a = simplyfyQuotient(n, g);
					Expr b = simplyfyQuotient(d, g);

					Int nume = a.value();
					Int deno = b.value();

					

					return fraction(nume,deno);
				}
				else {
					Expr min_one = integer(-1);
					Expr a = evaluateProduct(min_one, g);

					Expr f = fraction(simplyfyQuotient(n, a), simplyfyQuotient(d, a));
					
					return f;
				}
			}
			else if(d.value() < 0) {

				Expr min_one = integer(-1);
				Expr prodn = evaluateProduct(n, min_one);
				Expr prodd = evaluateProduct(d, min_one);
			
				Expr a = simplyfyQuotient(prodn, g);
				Expr b = simplyfyQuotient(prodd, g);
				
				Expr f = fraction(simplyfyQuotient(prodn, g), simplyfyQuotient(prodd, g));

				

				return f;
			}
		}
    }

    return u;    
}

Expr reduceRationalNumberASTRec(Expr u) {
	if(u.kind() == Kind::Integer)
		return u;

	if(u.kind() == Kind::Fraction) {
		Expr deno = denominator(u);
		
		if(deno == 0) {
			return undefined();
		}

		return u;
	}

	if(u.size() == 1) {
		Expr v = reduceRationalNumberASTRec(u[0]);

		if(v.kind() == Kind::Undefined) {
			
			return undefined();
		} else if(u.kind() == Kind::Addition) {
				return v;
		} else if(u.kind() == Kind::Subtraction) {
			Expr mi_one = integer(-1);
			Expr res = evaluateProduct(mi_one, v);
			
			
			
			return res;
		}
		else
				return u;
	}

	// if(u.size() < 2)
	// 		return u[0];

	if(
		u.kind() == Kind::Addition ||
		u.kind() == Kind::Multiplication ||
		u.kind() == Kind::Subtraction ||
		u.kind() == Kind::Division
	) {
		Expr v = reduceRationalNumberASTRec(u[0]);
		Expr w = reduceRationalNumberASTRec(u[1]);

		if(
			v.kind() == Kind::Undefined ||
			w.kind() == Kind::Undefined
		) {
			
			return undefined();
		}

		Expr res = nullptr;

		if(u.kind() == Kind::Addition) {

			res = evaluateAddition(v, w);
		}

		if(u.kind() == Kind::Subtraction)
			res = evaluateSubtraction(v, w);

		if(u.kind() == Kind::Multiplication)
			res = evaluateProduct(v, w);

		if(u.kind() == Kind::Division)
			res = evaluateDivision(v, w);
	
		

		return res;

	} else if(u.kind() == Kind::Power) {
		Expr v = reduceRationalNumberASTRec(u[0]);

		if(v.kind() == Kind::Undefined) {
			
			return undefined();
		}
		else if(u[1].kind() == Kind::Integer) {
			Expr res = evaluatePower(v, u[1]);
			
			return res;
		}

		
		return u;
	}

	return u;
}

Expr reduceRNEAST(Expr u) {

	if(!isRNE(u))
		return u;

	Expr v = reduceRationalNumberASTRec(u);

	if(v.kind() == Kind::Undefined) {
		
		return undefined();
	}

	Expr res = reduceRationalNumber(v);
	


	return res;
}

}
