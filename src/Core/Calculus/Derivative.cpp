
#include "Calculus.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Simplification/Simplification.hpp"
using namespace ast;
using namespace algebra;
using namespace simplification;

namespace calculus {

AST* derivative(AST* u, AST* x) {
	return new AST(
		Kind::Derivative,
		{ u, x }
	);
}

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
								power(u->operand(0)->operand(0)->copy(), integer(2))
							}),
							fraction(integer(1),integer(2))
						)
					),
					derivate(u->operand(0)->copy(), x)
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
								power(u->operand(0)->operand(0)->copy(), integer(2))
							}),
							fraction(integer(1),integer(2))
						)
					),
					derivate(u->operand(0)->copy(), x)
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
							power(u->operand(0)->operand(0)->copy(), integer(2))
						})
					),
					derivate(u->operand(0)->copy(), x)
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
							power(u->operand(0)->operand(0)->copy(), integer(2))
						})
					),
					derivate(u->operand(0)->copy(), x)
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
							abs(u->operand(0)->operand(0)->copy()),
							power( 
								sub({
									power(u->operand(0)->operand(0)->copy(), integer(2)),
									integer(1)
								}),
								fraction(integer(1),integer(2))
							)
						})
					),
					derivate(u->operand(0)->copy(), x)
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
							algebra::abs(u->operand(0)->operand(0)->copy()),
							power( 
								sub({
									power(u->operand(0)->operand(0)->copy(), integer(2)),
									integer(1)
								}),
								fraction(integer(1),integer(2))
							)
						})
					),
					derivate(u->operand(0)->copy(), x)
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

		delete d_;
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
		AST* v = u->operand(0)->copy();
		AST* w_ = div(u->copy(), v->copy());
		AST* w = reduceAST(w_);

		AST* d_ = add({
			mul({
				derivate(v, x),
				w->copy()
			}),
			mul({
				v->copy(),
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
			return algebra::abs(derivate(u->operand(0)->copy(), x));
		}

		if(u->funName() == "ln") {
			return reduceAST(div(
				integer(1),
				u->operand(0)->copy()
			));
		}
	
		if(u->funName() == "log") {
			if(u->numberOfOperands() == 2) {
				AST* dx_ = div(
					integer(1),
					mul({
						// x
						u->operand(0)->copy(),
						funCall(
							"ln", {
							// base
							u->operand(1)->copy()
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
						u->operand(0)->copy(),
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
				funCall("exp", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}
	
		if(u->funName() == "tan") {
			AST* dx_ = mul({
				power(
					funCall("sec", { u->operand(0)->copy() }),
					integer(2)
				),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "sinh") {
			AST* dx_ = mul({
				funCall("cosh", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "cosh") {
			AST* dx_ = mul({
				funCall("sinh", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "tanh") {
			AST* dx_ = mul({
				power(
					funCall("sech", { u->operand(0)->copy() }),
					integer(2)
				),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "sec") {
			AST* dx_ = mul({
				funCall("sec", { u->operand(0)->copy() }),
				funCall("tan", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}
	
		if(u->funName() == "csc") {
			AST* dx_ = mul({
				integer(-1),
				funCall("cot", { u->operand(0)->copy() }),
				funCall("csc", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "cot") {
			AST* dx_ = mul({
				integer(-1),
				power(
					funCall("csc", { u->operand(0)->copy() }),
					integer(2)
				),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "coth") {
			AST* dx_ = mul({
				integer(-1),
				power(
					funCall("csch", { u->operand(0)->copy() }),
					integer(2)
				),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "sech") {
			AST* dx_ = mul({
				integer(-1),
				funCall("tanh", { u->operand(0)->copy() }),
				funCall("sech", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "csch") {
			AST* dx_ = mul({
				integer(-1),
				funCall("coth", { u->operand(0)->copy() }),
				funCall("csch", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
			});
			AST* dx = reduceAST(dx_);
			delete dx_;
			return dx;
		}

		if(u->funName() == "sin") {
			AST* d_ = mul({
				funCall("cos", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
			});

			AST* d = reduceAST(d_);
	
			delete d_;
	
			return d;
		}
	
		if(u->funName() == "cos") {
			AST* d_ = mul({
				integer(-1),
				funCall("sin", { u->operand(0)->copy() }),
				derivate(u->operand(0)->copy(), x)
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
				derivate(u->operand(0)->copy(), x)
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
				derivate(u->operand(0)->copy(), x)
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
				derivate(u->operand(0)->copy(), x)
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
				derivate(u->operand(0)->copy(), x)
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
						algebra::abs(u->operand(0)->copy()),
						power(
							sub({
								power(symbol("x"), integer(2)),
								integer(1)
							}),
							fraction(integer(1), integer(2))
						)
					})
				),
				derivate(u->operand(0)->copy(), x),
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
						algebra::abs(u->operand(0)->copy()),
						power(
							sub({
								power(symbol("x"), integer(2)),
								integer(1)
							}),
							fraction(integer(1), integer(2))
						)
					})
				),
				derivate(u->operand(0)->copy(), x),
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
				derivate(u->operand(0)->copy(), x),
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
				derivate(u->operand(0)->copy(), x),
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
		u->copy(), 
		x->copy()
	);
}

}
