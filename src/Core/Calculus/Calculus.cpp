
#include "Calculus.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Simplification/Simplification.hpp"
using namespace ast;
using namespace algebra;
using namespace simplification;

namespace calculus {

// AST* integralTable(AST* u, AST* x) {

// 	if(u->freeOf(x))
// 		return integer(0);
	
// 	if(u->kind() == Kind::Power) {
// 		if(u->operand(1)->freeOf(x)) {
// 			if(u->operand(0)->match(x)) {
// 				if(
// 					u->operand(1)->kind() == Kind::Integer &&
// 					u->operand(1)->value() == -1
// 				) return funCall("ln", { u->operand(0)->deepCopy() });
			
// 				return mul({
// 					div(
// 						integer(1),
// 						add({ u->operand(1)->deepCopy(), integer(1) })
// 					),
// 					power(
// 						u->operand(0)->deepCopy(),
// 						add({ u->operand(1)->deepCopy(), integer(1) })
// 					),
// 				});
// 			}
// 		}

// 		if(u->operand(0)->freeOf(x)) {
// 			if(u->operand(1)->match(x)) {
// 				return div(
// 					power(
// 						u->operand(0)->deepCopy(),
// 						u->operand(1)->deepCopy()
// 					),
// 					funCall(
// 						"ln",
// 						{ 
// 							u->operand(0)->deepCopy() 
// 						}
// 					)
// 				);
// 			}
// 		}
// 	}

// 	if(u->kind() == Kind::FunctionCall) {
// 		if(u->funName() == "exp") {
// 			return funCall(
// 				"exp", {
// 					u->operand(0)->deepCopy()
// 				}
// 			);
// 		}

// 		if(u->funName() == "sin") {
// 			return mul({
// 				integer(-1),
// 				funCall(
// 					"cos", {
// 						u->operand(0)->deepCopy()
// 					}
// 				)
// 			});
// 		}
	
// 		if(u->funName() == "cos") {
// 			return funCall(
// 				"sin", {
// 					u->operand(0)->deepCopy()
// 				}
// 			);
// 		}

// 		if(u->funName() == "tan") {
// 			return mul({
// 				integer(-1),
// 				funCall(
// 					"ln",  {
// 						funCall(
// 							"cos", {
// 								u->operand(0)->deepCopy()
// 							}
// 						)
// 					}
// 				)
// 			});
// 		}
	
// 		if(u->funName() == "sinh") {
// 			return mul({
// 				integer(-1),
// 				funCall(
// 					"cosh", {
// 						u->operand(0)->deepCopy()
// 					}
// 				)
// 			});
// 		}
	
// 		if(u->funName() == "cosh") {
// 			return funCall(
// 				"sinh", {
// 					u->operand(0)->deepCopy()
// 				}
// 			);
// 		}

// 		if(u->funName() == "tanh") {
// 			return mul({
// 				integer(-1),
// 				funCall(
// 					"ln",  {
// 						funCall(
// 							"cosh", {
// 								u->operand(0)->deepCopy()
// 							}
// 						)
// 					}
// 				)
// 			});
// 		}
// 	}

// 	return new AST(Kind::Fail);
// }


// AST* linearProperties(AST* u, AST* x) {
// 	if(u->kind() == Kind::Integral) {
// 		return integrate(u, x);
// 	}

// 	if(u->kind() == Kind::Multiplication) {
// 		AST* a_ = integer(1);	 
// 		AST* b_ = integer(1);

// 		for(int i=0; i<u->numberOfOperands(); i++) {
// 			if(u->operand(i)->freeOf(x)) {
// 				a_ = mul({ a_, u->operand(i)->deepCopy()});
// 			} else {
// 				b_ = mul({ b_, u->operand(i)->deepCopy()});
// 			}
// 		}

// 		AST* a = reduceAST(a_);
// 		delete a_;
// 		AST* b = reduceAST(b_);
// 		delete b_;
		
// 		AST* u_ = mul({
// 			a->deepCopy(),
// 			integral(b->deepCopy(), x->deepCopy())
// 		});

// 		delete a;
// 		delete b;

// 		return u_;
// 	}

// 	if(u->kind() == Kind::Addition || u->kind() == Kind::Subtraction) {
// 		AST* u_ = new AST(u->kind());

// 		for(int i=0; i<u->numberOfOperands(); i++) {
// 			u_->includeOperand(linearProperties(u->operand(i), x));
// 		}

// 		return u_;
// 	}

// 	if(u->kind() == Kind::FunctionCall)
// 		return integral(u->deepCopy(), x->deepCopy());
	
// 	if(u->kind() == Kind::Power) {
// 		if(
// 			!u->operand(0)->freeOf(x) &&
// 			 u->operand(1)->freeOf(x) &&
// 			 u->operand(0)->kind() == Kind::Multiplication
// 		) {
// 			AST* a_ = integer(1);	 
// 			AST* b_ = integer(1);

// 			for(int i=0; i<u->numberOfOperands(); i++) {
// 				if(u->operand(i)->freeOf(x)) {
// 					a_ = mul({ a_, u->operand(i)->deepCopy()});
// 				} else {
// 					b_ = mul({ b_, u->operand(i)->deepCopy()});
// 				}
// 			}

// 			AST* a = reduceAST(a_);
// 			delete a_;

// 			AST* b = reduceAST(b_);
// 			delete b_;
			
// 			AST* a = power(a, u->operand(1)->deepCopy());
// 			AST* b = power(b, u->operand(1)->deepCopy());
			
// 			return mul({a, integral(b, x->deepCopy())});
// 		}
// 	}

// 	return u->deepCopy();
// }

AST* derivative(AST* u, AST* x) {
	return new AST(
		Kind::Derivative,
		{ u, x }
	);
}

// AST* integral(AST* u, AST* x) {
// 	return new AST(
// 		Kind::Integral,
// 		{ u, x }
// 	);
// }

// AST* integrate(AST* u, AST* x) {
// 	throw 'Not implemented';
// }

AST* derivateInverseTrig(AST* u, AST* x)
{
if(
		u->kind() == Kind::Power && 
		u->operand(1)->kind() == Integer && 
		u->operand(1)->value() == -1
	) {
		if(u->operand(0)->kind() == Kind::FunctionCall)
		{
			if(u->operand(0)->funName() == "sin")
			{
				AST* dx_ = mul({
					div(
						integer(1),
						power( 
							sub({
								integer(1),
								power(u->operand(0)->operand(0)->deepCopy(), integer(2))
							}),
							fraction(integer(1),integer(2))
						)
					),
					derivate(u->operand(0)->deepCopy(), x)
				});
				
				AST* dx = reduceAST(dx_);
				delete dx_;
				return dx;
			}


			if(u->operand(0)->funName() == "cos")
			{
				AST* dx_ = mul({
					integer(-1),
					div(
						integer(1),
						power( 
							sub({
								integer(1),
								power(u->operand(0)->operand(0)->deepCopy(), integer(2))
							}),
							fraction(integer(1),integer(2))
						)
					),
					derivate(u->operand(0)->deepCopy(), x)
				});
				
				AST* dx = reduceAST(dx_);
				delete dx_;
				return dx;
			}


			if(u->operand(0)->funName() == "tan")
			{
				AST* dx_ = mul({
					div(
						integer(1),
						add({
							integer(1),
							power(u->operand(0)->operand(0)->deepCopy(), integer(2))
						})
					),
					derivate(u->operand(0)->deepCopy(), x)
				});
				
				AST* dx = reduceAST(dx_);
				delete dx_;
				return dx;
			}
	
			if(u->operand(0)->funName() == "cot")
			{
				AST* dx_ = mul({
					integer(-1),
					div(
						integer(1),
						add({
							integer(1),
							power(u->operand(0)->operand(0)->deepCopy(), integer(2))
						})
					),
					derivate(u->operand(0)->deepCopy(), x)
				});
				
				AST* dx = reduceAST(dx_);
				delete dx_;
				return dx;
			}


			if(u->operand(0)->funName() == "sec")
			{
				AST* dx_ = mul({
					div(
						integer(1),
						mul({
							abs(u->operand(0)->operand(0)->deepCopy()),
							power( 
								sub({
									power(u->operand(0)->operand(0)->deepCopy(), integer(2)),
									integer(1)
								}),
								fraction(integer(1),integer(2))
							)
						})
					),
					derivate(u->operand(0)->deepCopy(), x)
				});
				
				AST* dx = reduceAST(dx_);
				delete dx_;
				return dx;
			}


			if(u->operand(0)->funName() == "csc")
			{
				AST* dx_ = mul({
					integer(-1),
					div(
						integer(1),
						mul({
							algebra::abs(u->operand(0)->operand(0)->deepCopy()),
							power( 
								sub({
									power(u->operand(0)->operand(0)->deepCopy(), integer(2)),
									integer(1)
								}),
								fraction(integer(1),integer(2))
							)
						})
					),
					derivate(u->operand(0)->deepCopy(), x)
				});
				
				AST* dx = reduceAST(dx_);
				delete dx_;
				return dx;
			}
		}
	}

	return nullptr;
}



AST* derivatePower(AST* u, AST* x)
{
	if(u->kind() == Kind::Power) 
	{
		AST* v = base(u);
		AST* w = expoent(u);

		AST* d_ = add({
			mul({
				expoent(u),
				power(
					base(u),
					sub({
						expoent(u),
						integer(1)
					})
				),
				derivate(v, x)
			}),

			mul({
				derivate(w, x),
				power(
					base(u),
					expoent(u)
				),
				funCall("ln", {
					base(u)
				})
			})
		});

		AST* d = reduceAST(d_);

		delete v;
		delete w;

		return d;
	}

	return nullptr;
}

AST* derivateSumsAndSubs(AST* u, AST* x)
{
	if(
		u->kind() == Kind::Addition || u->kind() == Kind::Subtraction
	) {
		AST* dx = new AST(u->kind());
		
		for(unsigned int i=0; i<u->numberOfOperands(); i++) 
		{
			dx->includeOperand(derivate(u->operand(i), x));
		}
		
		return dx;
	}

	return nullptr;
}	

AST* derivateMul(AST* u, AST* x)
{
	if(u->kind() == Kind::Multiplication) {
		AST* v = u->operand(0)->deepCopy();
		AST* w_ = div(u->deepCopy(), v->deepCopy());
		AST* w = reduceAST(w_);

		AST* d_ = add({
			mul({
				derivate(v, x),
				w->deepCopy()
			}),
			mul({
				v->deepCopy(),
				derivate(w, x)
			})
		});

		AST* d = reduceAST(d_);
	
		delete w_;
		delete v;
		delete w;
		delete d_;

		return d;
	}

	return nullptr;
}

AST* derivateFuncs(AST* u, AST* x)
{
	if(u->kind() == Kind::FunctionCall) {

		if(u->funName() == "abs") {
			return algebra::abs(derivate(u->operand(0)->deepCopy(), x));
		}

		if(u->funName() == "ln") {
			return reduceAST(div(
				integer(1),
				u->operand(0)->deepCopy()
			));
		}
	
		if(u->funName() == "log") {
			if(u->numberOfOperands() == 2) {
				AST* dx_ = div(
					integer(1),
					mul({
						// x
						u->operand(0)->deepCopy(),
						funCall(
							"ln", {
							// base
							u->operand(1)->deepCopy()
						})
					})
				);
				AST* dx = reduceAST(dx_);
				delete dx_;
				return dx;
			} else {
				AST* dx_ = div(
					integer(1),
					mul({
						// x
						u->operand(0)->deepCopy(),
						funCall(
							"ln", {
							integer(2)
						})
					})
				);
				AST* dx = reduceAST(dx_);
				delete dx_;
				return dx;
			}
		}

		if(u->funName() == "exp") {
			AST* dx_ = mul({
				funCall("exp", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}
	
		if(u->funName() == "tan") {
			AST* dx_ = mul({
				power(
					funCall("sec", { u->operand(0)->deepCopy() }),
					integer(2)
				),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "sinh") {
			AST* dx_ = mul({
				funCall("cosh", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "cosh") {
			AST* dx_ = mul({
				funCall("sinh", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "tanh") {
			AST* dx_ = mul({
				power(
					funCall("sech", { u->operand(0)->deepCopy() }),
					integer(2)
				),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "sec") {
			AST* dx_ = mul({
				funCall("sec", { u->operand(0)->deepCopy() }),
				funCall("tan", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}
	
		if(u->funName() == "csc") {
			AST* dx_ = mul({
				integer(-1),
				funCall("cot", { u->operand(0)->deepCopy() }),
				funCall("csc", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "cot") {
			AST* dx_ = mul({
				integer(-1),
				power(
					funCall("csc", { u->operand(0)->deepCopy() }),
					integer(2)
				),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "coth") {
			AST* dx_ = mul({
				integer(-1),
				power(
					funCall("csch", { u->operand(0)->deepCopy() }),
					integer(2)
				),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "sech") {
			AST* dx_ = mul({
				integer(-1),
				funCall("tanh", { u->operand(0)->deepCopy() }),
				funCall("sech", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "csch") {
			AST* dx_ = mul({
				integer(-1),
				funCall("coth", { u->operand(0)->deepCopy() }),
				funCall("csch", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "sin") {
			AST* d_ = mul({
				funCall("cos", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});

			AST* d = reduceAST(d_);
	
			delete d_;
	
			return d;
		}
	
		if(u->funName() == "cos") {
			AST* d_ = mul({
				integer(-1),
				funCall("sin", { u->operand(0)->deepCopy() }),
				derivate(u->operand(0)->deepCopy(), x)
			});

			AST* d = reduceAST(d_);
	
			delete d_;
	
			return d;
		}

		if(u->funName() == "arcsin")
		{
			AST* d_ = mul({
				div(
					integer(1),
					power(
						sub({
							integer(1),
							power(symbol("x"), integer(2))
						}),
						fraction(integer(1), integer(2))
					)
				),
				derivate(u->operand(0)->deepCopy(), x)
			});


			AST* d = reduceAST(d_);
			delete d_;
			return d;
		}

		if(u->funName() == "arccos")
		{
			AST* d_ = mul({
				integer(-1),
				div(
					integer(1),
					power(
						sub({
							integer(1),
							power(symbol("x"), integer(2))
						}),
						fraction(integer(1), integer(2))
					)
				),
				derivate(u->operand(0)->deepCopy(), x)
			});

			AST* d = reduceAST(d_);
			delete d_;
			return d;
		}


		if(u->funName() == "arctan")
		{
			AST* d_ = mul({
				div(
					integer(1),
					add({
						integer(1),
						power(symbol("x"), integer(2))
					})
				),
				derivate(u->operand(0)->deepCopy(), x)
			});


			AST* d = reduceAST(d_);
			delete d_;
			return d;
		}

		if(u->funName() == "arccot")
		{
			AST* d_ = mul({
				integer(-1),
				div(
					integer(1),
					add({
						integer(1),
						power(symbol("x"), integer(2))
					})
				),
				derivate(u->operand(0)->deepCopy(), x)
			});

			AST* d = reduceAST(d_);
			delete d_;
			return d;
		}


		if(u->funName() == "arcsec")
		{
			AST* d_ = mul({
				div(
					integer(1),
					mul({
						algebra::abs(u->operand(0)->deepCopy()),
						power(
							sub({
								power(symbol("x"), integer(2)),
								integer(1)
							}),
							fraction(integer(1), integer(2))
						)
					})
				),
				derivate(u->operand(0)->deepCopy(), x),
			});

			AST* d = reduceAST(d_);
			delete d_;
			return d;
		}


		if(u->funName() == "arccsc")
		{
			AST* d_ = mul({
				integer(-1),
				div(
					integer(1),
					mul({
						algebra::abs(u->operand(0)->deepCopy()),
						power(
							sub({
								power(symbol("x"), integer(2)),
								integer(1)
							}),
							fraction(integer(1), integer(2))
						)
					})
				),
				derivate(u->operand(0)->deepCopy(), x),
			});

			AST* d = reduceAST(d_);
			delete d_;
			return d;
		}

		if(u->funName() == "arccosh")
		{
			AST* d_ = mul({
				div(
					integer(1),
					power(
						sub({
							power(symbol("x"), integer(2)),
							integer(1)
						}),
						fraction(integer(1), integer(2))
					)
				),
				derivate(u->operand(0)->deepCopy(), x),
			});

			AST* d = reduceAST(d_);
			delete d_;
			return d;
		}
	
		if(u->funName() == "arctanh")
		{
			AST* d_ = mul({
				div(
					integer(1),
					sub({
						integer(1),
						power(symbol("x"), integer(2))
					})
				),
				derivate(u->operand(0)->deepCopy(), x),
			});

			AST* d = reduceAST(d_);
			delete d_;
			return d;
		}
	}

	return nullptr;
}

AST* derivate(AST* u, AST* x) {
	AST* dx;

	if(u->match(x))
	{
		return integer(1);
	}

	dx = derivateInverseTrig(u, x);
	if(dx) 
	{
		return dx;
	}

	dx = derivatePower(u, x);
	if(dx) 
	{
		return dx;
	}

	dx = derivateSumsAndSubs(u, x);
	if(dx) 
	{
		return dx;
	}

	dx = derivateMul(u, x);
	if(dx) 
	{
		return dx;
	}

	dx = derivateFuncs(u, x);
	if(dx) 
	{
		return dx;
	}

	if(u->freeOf(x)) 
	{
		return integer(0);
	}

	return derivative(
		u->deepCopy(), 
		x->deepCopy()
	);
}

}
