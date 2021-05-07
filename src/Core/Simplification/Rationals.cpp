#include "Rationals.hpp"

#include <cstdio>

using namespace ast;
using namespace algebra;

namespace simplification {

AST* reduceRationalNumber(AST* u);

AST* evaluateProduct(AST* u, AST* v) {
	if(isConstant(u) && isConstant(v)) {
		if(u->kind() == Kind::Integer && v->kind() == Kind::Integer)
			return integer(u->value() * v->value());

		AST* num_u = numerator(u);	
		AST* den_u = denominator(u);	
		AST* num_v = numerator(v);	
		AST* den_v = denominator(v);	

		AST* f = fraction(
			evaluateProduct(num_u, num_v),
			evaluateProduct(den_u, den_v)  
		);

		AST* simp_f = reduceRNEAST(f);
		
		destroyASTs({ num_u, num_v, den_u, den_v, f });
		
		return simp_f;
	}

	return mul({u->deepCopy(), v->deepCopy()});
}

AST* evaluateAddition(AST* u, AST* v) {
	if(u->kind() == Kind::Integer && v->kind() == Kind::Integer)
		return integer(u->value() + v->value());

	// x/y + r = (x+ry)/y
	if(u->kind() == Kind::Fraction && v->kind() == Kind::Integer) {
		AST* num_u 	= numerator(u);
		AST* den_u 	= denominator(u);
		AST* k 			= evaluateProduct(v, den_u);

		AST* f = fraction(
			evaluateAddition(num_u, k),
			denominator(u)
		);
		
		AST* simp_f = reduceRNEAST(f);
		
		destroyASTs({num_u, den_u, k, f});
	
		return simp_f;
	}

	// r + x/y= (ry + x)/y
	if(v->kind() == Kind::Fraction && u->kind() == Kind::Integer) {
		AST* num_v 	= numerator(v);
		AST* den_v 	= denominator(v);
		AST* k 			= evaluateProduct(u, den_v);

		AST* f = fraction(
			evaluateAddition(k, num_v),
			denominator(v)
		);

		AST* simp_f = reduceRNEAST(f);
		
		destroyASTs({ num_v, den_v, k, f });
	
		return simp_f;
	}
	// x/y + a/b = (xb + ay)/ab
	if(v->kind() == Kind::Fraction && u->kind() == Kind::Fraction) {
		AST* num_u 					= numerator(u);
		AST* den_u 					= denominator(u);
		AST* num_v 					= numerator(v);
		AST* den_v 					= denominator(v);
		AST* num_u_t_dun_v 	= evaluateProduct(num_u, den_v);
		AST* num_v_t_dun_u 	= evaluateProduct(num_v, den_u);

		AST* f = fraction(
			evaluateAddition(num_u_t_dun_v, num_v_t_dun_u),
			evaluateProduct(den_u, den_v)
		);
	
		AST* simp_f = reduceRNEAST(f);
		
		destroyASTs({num_u, den_u, num_v, den_v, num_u_t_dun_v, num_v_t_dun_u, f });
		
		return simp_f;
	}

	return add({ u->deepCopy(), v->deepCopy() });
}


AST* evaluateSubtraction(AST* u, AST* v) {
	if(u->kind() == Kind::Integer && v->kind() == Kind::Integer)
		return integer(u->value() - v->value());

	// x/y - r = (x-ry)/y
	if(u->kind() == Kind::Fraction && v->kind() == Kind::Integer) {
		AST* num_u 	= numerator(u);
		AST* den_u 	= denominator(u);
		AST* k 			= evaluateProduct(v, den_u);

		AST* f = fraction(
			evaluateSubtraction(num_u, k),
			denominator(u)
		);
	
		AST* simp_f = reduceRNEAST(f);
		
		destroyASTs({ num_u, den_u, k, f});
	
		return simp_f;
	}

	// r - x/y= (ry - x)/y
	if(v->kind() == Kind::Fraction && u->kind() == Kind::Integer) {
		AST* num_v 	= numerator(v);
		AST* den_v 	= denominator(v);
		AST* k 			= evaluateProduct(u, den_v);

		AST* f = fraction(
			evaluateSubtraction(k, num_v),
			denominator(v)
		);
		
		AST* simp_f = reduceRNEAST(f);
		
		destroyASTs({num_v, den_v, k, f});
	
		return simp_f;
	}
	// x/y - a/b = (xb - ay)/ab
	if(v->kind() == Kind::Fraction && u->kind() == Kind::Fraction) {
		AST* num_u 					= numerator(u);
		AST* den_u 					= denominator(u);
		AST* num_v 					= numerator(v);
		AST* den_v 					= denominator(v);
		AST* num_u_t_dun_v 	= evaluateProduct(num_u, den_v);
		AST* num_v_t_dun_u 	= evaluateProduct(num_v, den_u);

		AST* f = fraction(
			evaluateSubtraction(num_u_t_dun_v, num_v_t_dun_u),
			evaluateProduct(den_u, den_v)
		);
	
		AST* simp_f = reduceRNEAST(f);
		
		destroyASTs({num_u, den_u, num_v, den_v, num_u_t_dun_v, num_v_t_dun_u, f });
		
		return simp_f;
	}

	return sub({ u->deepCopy(), v->deepCopy() });
}

AST* evaluateDivision(AST* v, AST* w) {
	if(!isConstant(v) || !isConstant(w))
		return div(v, w);

	AST* num_w = numerator(w);
	AST* den_w = denominator(w);

	AST* num_v = numerator(v);
	AST* den_v = denominator(v);

	if(num_w->kind() == Kind::Integer && num_w->value() == 0) {
		destroyASTs({ num_w, den_w, num_v, den_v });
		return new AST(Kind::Undefined);
	}

	// if(num_w->kind() == Kind::Integer && num_w->value() == 1) {
	// 	destroyASTs({ num_w, den_w, num_v, den_v });
	// 	return v->deepCopy();
	// }

	AST* f = fraction(
		evaluateProduct(num_v, den_w),
		evaluateProduct(num_w, den_v)
	);	

	AST* simp_f = reduceRNEAST(f);

	destroyASTs({ num_w, den_w, num_v, den_v, f });

	return simp_f;
}

AST* evaluatePower(AST* v, AST* w) {
	if(w->kind() != Kind::Integer)
			return power(v->deepCopy() ,w->deepCopy());
	
	long long n = w->value();

	if(!isConstant(v))
			return power(v->deepCopy(), integer(n));

	AST* num_v = numerator(v);
	if(num_v->kind() == Kind::Integer && num_v->value() != 0) {
		destroyASTs({ num_v });

		if(n > 0) {
			AST* e 		= integer(n - 1);
			AST* s 		= evaluatePower(v, e);
			AST* res 	= evaluateProduct(s, v);
			
			destroyASTs({ e, s });
			return res;
		} else if(n == 0) {
				return integer(1);
		} else if(n == -1) {
				AST* f 				= fraction(denominator(v), numerator(v));
				AST* f_simp 	= reduceRNEAST(f);
				
				destroyASTs({f});
				
				return f_simp;
		} else {
				AST* f 			= fraction(denominator(v), numerator(v));
				AST* f_simp = reduceRNEAST(f);
				AST* e 			= integer(-1 * n);
				AST* res 		= evaluatePower(f_simp, e);

				destroyASTs({ f, f_simp, e });

				return res;
		}
	}

	destroyASTs({ num_v });
		
	if(n >= 1) 
		return integer(0);

	return new AST(Kind::Undefined);
}



AST* simplyfyQuotient(AST* u, AST* v) {
    return integer(u->value() / v->value());
}

AST* reduceRationalNumber(AST* u) {

    if(u->kind() == Kind::Integer) {
			return u->deepCopy();
		}

    if(u->kind() == Kind::Fraction) {

			AST* n = u->operand(0);
			AST* d = u->operand(1);
			if(d->value() == 1) {

					return n->deepCopy();
			}

			if(n->value() % d->value() == 0) {
					return simplyfyQuotient(n, d);
			}
			else {
					AST* g = gcd(n, d);

					if(d->value() > 0) {

							if(g->value() > 0) {
									AST* a = simplyfyQuotient(n, g);
									AST* b = simplyfyQuotient(d, g);

									signed long nume = a->value();
									signed long deno	= b->value();

									destroyASTs({a, b, g});

									return fraction(nume,deno);
							}
							else {
								AST* min_one = integer(-1);
								AST* a = evaluateProduct(min_one, g);

								AST* f = fraction(simplyfyQuotient(n, a), simplyfyQuotient(d, a));
								destroyASTs({a, g, min_one});
								return f;
							}
					}
					else if(d->value() < 0) {

						AST* min_one = integer(-1);
						AST* prodn = evaluateProduct(n, min_one);
						AST* prodd = evaluateProduct(d, min_one);
					
						AST* a = simplyfyQuotient(prodn, g);
						AST* b = simplyfyQuotient(prodd, g);
						
						AST* f = fraction(simplyfyQuotient(prodn, g), simplyfyQuotient(prodd, g));

						destroyASTs({a, b, prodn, prodd, min_one, g});

						return f;
					}
			}
    }

    return u->deepCopy();    
}

AST* reduceRationalNumberASTRec(AST* u) {
	if(u->kind() == Kind::Integer)
		return u->deepCopy();

	if(u->kind() == Kind::Fraction) {
		AST* deno = denominator(u);
		
		bool is_den_zero = deno->kind() == Kind::Integer && deno->value() == 0;
		delete deno;

		if(is_den_zero) {
			return new AST(Kind::Undefined);
		}

		return u->deepCopy();
	}

	if(u->numberOfOperands() == 1) {
		AST* v = reduceRationalNumberASTRec(u->operand(0));

		if(v->kind() == Kind::Undefined) {
			delete v;
			return new AST(Kind::Undefined);
		} else if(u->kind() == Kind::Addition) {
				return v;
		} else if(u->kind() == Kind::Subtraction) {
			AST* mi_one = integer(-1);
			AST* res = evaluateProduct(mi_one, v);
			
			destroyASTs({ v, mi_one });
			
			return res;
		}
		else
				return u->deepCopy();
	}

	// if(u->numberOfOperands() < 2)
	// 		return u->operand(0)->deepCopy();

	if(
			u->kind() == Kind::Addition ||
			u->kind() == Kind::Multiplication ||
			u->kind() == Kind::Subtraction ||
			u->kind() == Kind::Division
	) {
			AST* v = reduceRationalNumberASTRec(u->operand(0));
			AST* w = reduceRationalNumberASTRec(u->operand(1));

			if(
				v->kind() == Kind::Undefined ||
				w->kind() == Kind::Undefined
			) {
				destroyASTs({v,w});
				return new AST(Kind::Undefined);
			}

			AST* res = nullptr;

			if(u->kind() == Kind::Addition) {

				res = evaluateAddition(v, w);
			}

			if(u->kind() == Kind::Subtraction)
				res = evaluateSubtraction(v, w);

			if(u->kind() == Kind::Multiplication)
				res = evaluateProduct(v, w);

			if(u->kind() == Kind::Division)
				res = evaluateDivision(v, w);
		
			destroyASTs({ v, w });

			return res;

	} else if(u->kind() == Kind::Power) {
			AST* v = reduceRationalNumberASTRec(u->operand(0));

			if(v->kind() == Kind::Undefined) {
				destroyASTs({v});
				return new AST(Kind::Undefined);
			}
			else if(u->operand(1)->kind() == Kind::Integer) {
				AST* res = evaluatePower(v, u->operand(1));
				destroyASTs({v});
				return res;
			}

			destroyASTs({v});
			return u->deepCopy();
	}

	return u->deepCopy();
}

AST* reduceRNEAST(AST* u) {
	if(!isRNE(u))
		return u->deepCopy();

	AST* v = reduceRationalNumberASTRec(u);

	if(v->kind() == Kind::Undefined) {
		delete v;
		return new AST(Kind::Undefined);
	}

	AST* res = reduceRationalNumber(v);
	delete v;
	return res;
}

}
