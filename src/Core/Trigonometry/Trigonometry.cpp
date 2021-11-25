#include "Trigonometry.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Polynomial/Polynomial.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace polynomial;
using namespace simplification;

namespace trigonometry {

Expr expandTrigRules(Expr A);
Expr contractTrigRules(Expr u);

Expr substituteTrig(Expr u) {
	if(
		u.kind() == Kind::Integer ||
		u.kind() == Kind::Fraction ||
		u.kind() == Kind::Symbol
	) return u;

	Expr g = mapUnaryAST(u, substituteTrig);

	if(g.kind() == Kind::FunctionCall) {
		if(g.funName() == "tan") {
			Expr k = div(
				funCall("sin", { g[0] }),
				funCall("cos", { g[0] })
			);
			
			return k;
		}

		if(g.funName() == "cot") {
			Expr k = div(
				funCall("cos", { g[0] }),
				funCall("sin", { g[0] })
			);
			
			return k;
		}

		if(g.funName() == "sec") {
			Expr k = div(
				integer(1),
				funCall("cos", { g[0] })
			);
			
			return k;
		}

		if(g.funName() == "csc") {
			Expr k = div(
				integer(1),
				funCall("sin", { g[0] })
			);
			
			return k;
		}
	}

	return g;
}

Expr multipleAndlgeCos(Expr n, Expr theta) {
	Expr r = integer(0);

	for(Int j=0; j <= n.value(); j++) {
		if(j%2 != 0)
			continue;
		
		Int b = 
			fact(n.value()) / 
			(fact(j) * fact(n.value() - j));
		
		Expr c_ = funCall("cos", {
			theta
		});

		Expr s_ = funCall("sin", {
			theta
		});
		
		Expr c = c_;
		Expr s = s_;

		if(theta.kind() == Kind::Addition) {
			Expr e_s = expandTrigRules(s_);
			s = e_s[0];
			Expr e_c = expandTrigRules(c_);
			c = e_s[1];
			
			
			
			
			
		}


		Expr e = mul({
			power(integer(-1), div(integer(j), integer(2))),
			integer(b),
			power(c, sub({ n, integer(j) })),
			power(s, integer(j))
		});

		r = add({r, e});

	}
	Expr d = reduceAST(r);
	
	return d;
}


Expr multipleAndlgeSin(Expr n, Expr theta) {
	Expr r = integer(0);

	for(int j=0; j <= n.value(); j++) {
		if(j%2 != 1)
			continue;
		
		Int b = fact(n.value())/(fact(j) * fact(n.value() - j));
		Expr c_ = funCall("cos", {
			theta
		});

		Expr s_ = funCall("sin", {
			theta
		});
		
		Expr c = c_;
		Expr s = s_;

		if(theta.kind() == Kind::Addition) {
			Expr e_s = expandTrigRules(s_);
			s = e_s[0];
			Expr e_c = expandTrigRules(c_);
			c = e_s[1];
			
			
			
			
			
		}


		Expr e = mul({
			power(integer(-1), div( sub({ integer(j), integer(1) }), integer(2))),
			integer(b),
			power(c, sub({ n, integer(j) })),
			power(s, integer(j))
		});
	
		r = add({r, e});
	}

	Expr expanded = reduceAST(r);

	

	return expanded;
}


Expr expandTrigRules(Expr A) {
	if(A.kind() == Kind::Addition) {
		Expr f = expandTrigRules(A[0]);
		
		Expr A__ = sub({
			A,
			A[0]
		});
		
		Expr A_ = reduceAST(A__);

		Expr r = expandTrigRules(A_);

		Expr s_ = add({
			mul({
				f[0],
				r[1],
			}),
			mul({
				f[1],
				r[0],
			}),
		});

		Expr s = reduceAST(s_);

		Expr c_ = sub({
			mul({
				f[1],
				r[1],
			}),
			mul({
				f[0],
				r[0],
			}),
		});

		Expr c = reduceAST(s_);

		
		
		
		

		
		

		return list({s,c});
	}

	if(A.kind() == Kind::Multiplication) {
		Expr f = A[0];
		if(f.kind() == Kind::Integer) {
			Expr k_ = div(
				A,
				f
			);

			Expr k = reduceAST(k_);
	
			Expr a = multipleAndlgeSin(f, k);
			Expr b = multipleAndlgeCos(f, k);
	

			
			

			return list({a, b});
		}
	}

	return list({
		funCall("sin", { A }),
		funCall("cos", { A }),
	});
}

Expr expandTrig(Expr u) {
	if(
		u.kind() == Kind::Integer ||
		u.kind() == Kind::Fraction ||
		u.kind() == Kind::Symbol
	) return u;

	Expr u_ = algebraicExpand(u);
	Expr v = mapUnaryAST(u_, expandTrig);
	

	if(
		v.kind() == Kind::FunctionCall &&
		v.funName() == "sin"
	) {
		Expr a_ = expandTrigRules(v[0]);
		Expr a = reduceAST(a_);

		Expr r = a[0];
		
		
		
		

		if(isDivisionByZero(r)) {
			
			return undefined();
		}

		return r;
	}

	if(
		v.kind() == Kind::FunctionCall &&
		v.funName() == "cos"
	) {
		Expr a_ = expandTrigRules(v[0]);
		Expr a = reduceAST(a_);

		Expr r = a[1];

		
		
		

		if(isDivisionByZero(r)) {
			
			return undefined();
		}

		return r;
	}

	if(isDivisionByZero(v)) {
		
		return undefined();
	}

	return v;
}

Int floor(double x) {
	Int xi = (Int)x;
	return x < xi ? xi - 1 : xi;
}

Expr contractTrigPower(Expr u) {

	Expr b = u[0];
	Expr n = u[1];

	if(n.kind() == Kind::Integer && n.value() > 0) {
		if(
			b.kind() == Kind::FunctionCall && b.funName() == "cos"
		) {
			Expr theta = b[0];
	
			if(n.value() % 2 == 0) {
				// n even
				Expr n_ = integer(n.value()/2);

				Expr p0 = div(
					div(
						integer(fact(n.value())),
						integer(fact(n_.value()) * fact(n.value() - n_.value()))
					),
					power(integer(2), n)
				);

				Expr p1 = div(
					integer(1),
					power(integer(2), sub({n, integer(1)}))
				);

				Expr p2 = integer(0);

				for(int j=0; j <= n_.value() - 1; j++) {
					Expr b = div(
						integer(fact(n.value())),
						integer(fact(j) * fact(n.value() - j))
					);

					Expr c = funCall("cos", {
						mul({
							sub({
								n,
								mul({integer(2), integer(j)})
							}),
							theta
						})
					});

					p2 = add({ p2, mul({ b, c })});
				}

				

				Expr r_ = add({
					p0,
					mul({p1, p2})
				});
				Expr r = reduceAST(r_);
	
				

				return r;
			}


			if(n.value() % 2 == 1) {
				// n ood
				Expr n_ = integer(floor(n.value().longValue()/2.0f));

				Expr p1 = div(
					integer(1),
					power(integer(2), sub({n, integer(1)}))
				);

				Expr p2 = integer(0);
				for(int j=0; j<=n_.value() - 1; j++) {
					Expr b = div(
						integer(fact(n.value())),
						integer(fact(j) * fact(n.value() - j))
					);
					Expr c = funCall("cos", {
						mul({
							sub({
								n,
								mul({integer(2), integer(j)})
							}),
							theta
						})
					});
					p2 = add({
						p2, mul({
							b, c
						})
					});
				}
				

				Expr r_ = mul({p1, p2});

				Expr r = reduceAST(r_);
	
				

				return r;
			}
		}
		if(
			b.kind() == Kind::FunctionCall && b.funName() == "sin" 
		) {
			
			Expr theta = b[0];
	
			if(n.value() % 2 == 0) {
				// n even
				Expr n_ = integer(n.value()/2);

				Expr p0 = div(
					mul({
						power(integer(-1), n),
						div(
							integer(fact(n.value())),
							integer(fact(n.value()/2) * fact(n.value() - (n.value()/2)))
						)
					}),
					power(integer(2), n)
				);
	
				Expr p1 = div(
					power(integer(-1), integer(n.value()/2)),
					power(integer(2), sub({ n, integer(1) }))
				);

				Expr p2 = integer(0);
				for(int j=0; j<=n_.value() - 1; j++) {
					Expr a = power(
						integer(-1),
						integer(j)
					);
				
					Expr b = div(
						integer(fact(n.value())),
						integer(fact(j) * fact(n.value() - j))
					);
				
					Expr c = funCall("cos", {
						mul({
							sub({
								n,
								mul({integer(2), integer(j)})
							}),
							theta
						})
					});
	
					p2 = add({
						p2, mul({
							a, b, c
						})
					});
				}
				

				Expr r_ = add({
					p0,
					mul({p1, p2})
				});

				Expr r = reduceAST(r_);
	
				

				return r;
			}

			if(n.value() % 2 == 1) {
				// n odd
				Expr n_ = integer(floor(n.value().longValue() / 2.f));

				Expr p1 = div(
					power(integer(-1), integer((n.value() - 1)/2)),
					power(integer(2), sub({ n, integer(1) }))
				);
	
				Expr p2 = integer(0);
				for(int j=0; j<=n_.value() - 1; j++) {
					Expr a = power(
						integer(-1),
						integer(j)
					);
				
					Expr b = div(
						integer(fact(n.value())),
						integer(fact(j) * fact(n.value() - j))
					);
				
					Expr c = funCall("sin", {
						mul({
							sub({
								n,
								mul({integer(2), integer(j)})
							}),
							theta
						})
					});
	
					p2 = add({
						p2, mul({
							b, a, c
						})
					});
				}

				

				Expr r_ = mul({ p1, p2 });

				Expr r = reduceAST(r_);
	
				

				return r;
			}
		}
	}
	return u;
}

Expr separateSinCos(Expr u) {
	if(u.kind() == Kind::Multiplication) {
		Expr s = integer(1);
		Expr r = integer(1);

		for(unsigned int i=0; i<u.size(); i++) {
			Expr y = u[i];
			
			if(
				(y.kind() == Kind::FunctionCall && y.funName() == "sin") ||
				(y.kind() == Kind::FunctionCall && y.funName() == "cos")
			) {
				s = mul({
					s,
					y
				});
			} else
			if(
				y.kind() == Kind::Power && 
				y[0].kind() == Kind::FunctionCall &&
				(
					y[0].funName() == "sin" ||
					y[0].funName() == "cos"
				) &&
				y[1].kind() == Kind::Integer &&
				y[1].value() > 0
			) {
				s = mul({
					s,
					y
				});
			} else {
				r = mul({
					r,
					y
				});
			}
		}
	
		Expr L = list({reduceAST(r), reduceAST(s)});
		
		
		
		
		return L;
	}

	if(
		(u.kind() == Kind::FunctionCall && u.funName() == "sin") ||
		(u.kind() == Kind::FunctionCall && u.funName() == "cos")
	) {
		return list({integer(1), u});
	}

	if(
		u.kind() == Kind::Power && 
		u[0].kind() == Kind::FunctionCall &&
		(
			u[0].funName() == "sin" ||
			u[0].funName() == "cos"
		) &&
		u[1].kind() == Kind::Integer &&
		u[1].value() > 0
	) {
		return list({integer(1), u});
	}

	return list({u, integer(1)});
}

Expr contractTrigProduct(Expr u) {

	if(u.kind() == Kind::Integer) {
		return u;
	}

	if(u.size() == 1) {
		return u[0];
	}

	if(u.size() == 2) {
		Expr A = u[0];
		Expr B = u[1];

		if(A.kind() == Kind::Power) {
			A = contractTrigPower(A);
			
			Expr C = mul({
				A,
				B,
			});
			Expr r = contractTrigRules(C);
		
			
			
		
			return r;
		}

		if(B.kind() == Kind::Power) {
			B = contractTrigPower(B);
			Expr C = mul({
				A,
				B,
			});
			Expr r = contractTrigRules(C);
		
			
			
		
			return r;
		}

		Expr theta 	= A[0];
		Expr phi 		= B[0];

		if(
			A.kind() == Kind::FunctionCall &&
			A.funName() == "sin" &&
			B.kind() == Kind::FunctionCall &&
			B.funName() == "sin"
		) {
			Expr t = sub({
				div(
					funCall("cos", { sub({theta, phi})}),
					integer(2)
				),
				div(
					funCall("cos", { add({theta, phi})}),
					integer(2)
				)
			});

			Expr r = reduceAST(t);

			
			return r;
		}

		if(
			A.kind() == Kind::FunctionCall &&
			A.funName() == "cos" &&
			B.kind() == Kind::FunctionCall &&
			B.funName() == "cos"
		) {
			Expr t = add({
				div(
					funCall("cos", { add({theta, phi})}),
					integer(2)
				),
				div(
					funCall("cos", { sub({theta, phi})}),
					integer(2)
				)
			});

			Expr r = reduceAST(t);

			
			return r;
		}

		if(
			A.kind() == Kind::FunctionCall &&
			A.funName() == "sin" &&
			B.kind() == Kind::FunctionCall &&
			B.funName() == "cos"
		) {
			Expr t = add({
				div(
					funCall("sin", { add({theta, phi})}),
					integer(2)
				),
				div(
					funCall("sin", { sub({theta, phi})}),
					integer(2)
				)
			});

			Expr r = reduceAST(t);

			
			return r;
		}

		if(
			A.kind() == Kind::FunctionCall &&
			A.funName() == "cos" &&
			B.kind() == Kind::FunctionCall &&
			B.funName() == "sin"
		) {
			Expr t = add({
				div(
					funCall("sin", { add({theta, phi})}),
					integer(2)
				),
				div(
					funCall("sin", { sub({phi, theta})}),
					integer(2)
				)
			});

			Expr r = reduceAST(t);

			
			return r;
		}
	}

	Expr A = u[0];

	Expr k_ = div(
		u,
		A
	);

	Expr k = reduceAST(k_);

	

	Expr B = contractTrigProduct(k);

	

	Expr r = mul({
		A,
		B,
	});

	Expr d = contractTrigRules(r);

	

	return d;
}

Expr contractTrigRules(Expr u) 
{
	Expr v = algebraicExpandRoot(u);

	if(v.kind() == Kind::Power) {
		Expr t = contractTrigPower(v);

		
		return t;
	}

	if(v.kind() == Kind::Multiplication) {
		Expr s = separateSinCos(v);
	

		Expr c = s[0];
		Expr d = s[1];
		
		if(d.kind() == Kind::Integer && d.value() == 1) {
			
			return v;
		}

		if(
			(d.kind() == Kind::FunctionCall && d.funName() == "sin") ||
			(d.kind() == Kind::FunctionCall && d.funName() == "cos")
		) {
			
			return v;
		}

		if(d.kind() == Kind::Power) {
			Expr k_ = mul({
				c,
				contractTrigPower(d)
			});

			Expr k = reduceAST(k_);
			Expr r = algebraicExpandRoot(k);
			
			
			
			
			

			return r;
		}

		Expr k_ = mul({
			c,
			contractTrigProduct(d)
		});

		Expr k = reduceAST(k_);
		Expr r = algebraicExpandRoot(k);

		
		
		
		

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
					s, contractTrigRules(y)
				});

			} else {
				s = add({
					s, y
				});
			}
		}

		Expr r = reduceAST(s);

		
		
		return r;
	}

	return v;
}

Expr contractTrig(Expr u) {
	if(
		u.kind() == Kind::Integer ||
		u.kind() == Kind::Fraction ||
		u.kind() == Kind::Symbol
	) return u;

	Expr v_ = mapUnaryAST(u, contractTrig);
	Expr v = reduceAST(v_);
	

	if(
		v.kind() == Kind::Multiplication ||
		v.kind() == Kind::Power
	) 
	{
		Expr t = contractTrigRules(v);
		
		return t;
	}

	return v;
}

}
