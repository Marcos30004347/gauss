
#include "Calculus.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Simplification/Simplification.hpp"
using namespace ast;
using namespace algebra;
using namespace simplification;

namespace calculus {

Expr derivative(Expr u, Expr x) {
	return Expr(
		Kind::Derivative,
		{ u, x }
	);
}

Expr derivateInverseTrig(Expr u, Expr x)
{
if(
		u.kind() == Kind::Power && 
		u[1].kind() == Integer && 
		u[1].value() == -1
	) {
		if(u[0].kind() == Kind::FunctionCall)
		{
			if(u[0].funName() == "sin")
			{
				Expr dx_ = mul({
					div(
						integer(1),
						power( 
							sub({
								integer(1),
								power(u[0][0], integer(2))
							}),
							fraction(integer(1),integer(2))
						)
					),
					derivate(u[0], x)
				});
				
				Expr dx = reduceAST(dx_);
				
				return dx;
			}


			if(u[0].funName() == "cos")
			{
				Expr dx_ = mul({
					integer(-1),
					div(
						integer(1),
						power( 
							sub({
								integer(1),
								power(u[0][0], integer(2))
							}),
							fraction(integer(1),integer(2))
						)
					),
					derivate(u[0], x)
				});
				
				Expr dx = reduceAST(dx_);
				
				return dx;
			}


			if(u[0].funName() == "tan")
			{
				Expr dx_ = mul({
					div(
						integer(1),
						add({
							integer(1),
							power(u[0][0], integer(2))
						})
					),
					derivate(u[0], x)
				});
				
				Expr dx = reduceAST(dx_);
				
				return dx;
			}
	
			if(u[0].funName() == "cot")
			{
				Expr dx_ = mul({
					integer(-1),
					div(
						integer(1),
						add({
							integer(1),
							power(u[0][0], integer(2))
						})
					),
					derivate(u[0], x)
				});
				
				Expr dx = reduceAST(dx_);
				
				return dx;
			}


			if(u[0].funName() == "sec")
			{
				Expr dx_ = mul({
					div(
						integer(1),
						mul({
							abs(u[0][0]),
							power( 
								sub({
									power(u[0][0], integer(2)),
									integer(1)
								}),
								fraction(integer(1),integer(2))
							)
						})
					),
					derivate(u[0], x)
				});
				
				Expr dx = reduceAST(dx_);
				
				return dx;
			}


			if(u[0].funName() == "csc")
			{
				Expr dx_ = mul({
					integer(-1),
					div(
						integer(1),
						mul({
							algebra::abs(u[0][0]),
							power( 
								sub({
									power(u[0][0], integer(2)),
									integer(1)
								}),
								fraction(integer(1),integer(2))
							)
						})
					),
					derivate(u[0], x)
				});
				
				Expr dx = reduceAST(dx_);
				
				return dx;
			}
		}
	}

	return undefined();
}



Expr derivatePower(Expr u, Expr x)
{
	if(u.kind() == Kind::Power) 
	{
		Expr v = base(u);
		Expr w = expoent(u);

		Expr d_ = add({
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

		Expr d = reduceAST(d_);

		return d;
	}

	return undefined();
}

Expr derivateSumsAndSubs(Expr u, Expr x)
{
	if(
		u.kind() == Kind::Addition || u.kind() == Kind::Subtraction
	) {
		Expr dx = Expr(u.kind());
		
		for(unsigned int i=0; i<u.size(); i++) 
		{
			dx.insert(derivate(u[i], x));
		}
		
		return dx;
	}

	return undefined();
}	

Expr derivateMul(Expr u, Expr x)
{
	if(u.kind() == Kind::Multiplication) {
		Expr v = u[0];
		Expr w_ = div(u, v);
		Expr w = reduceAST(w_);

		Expr d_ = add({
			mul({
				derivate(v, x),
				w
			}),
			mul({
				v,
				derivate(w, x)
			})
		});

		Expr d = reduceAST(d_);
	
		
		
		
		

		return d;
	}

	return undefined();
}

Expr derivateFuncs(Expr u, Expr x)
{
	if(u.kind() == Kind::FunctionCall) {

		if(u.funName() == "abs") {
			return algebra::abs(derivate(u[0], x));
		}

		if(u.funName() == "ln") {
			return reduceAST(div(
				integer(1),
				u[0]
			));
		}
	
		if(u.funName() == "log") {
			if(u.size() == 2) {
				Expr dx_ = div(
					integer(1),
					mul({
						// x
						u[0],
						funCall(
							"ln", {
							// base
							u[1]
						})
					})
				);
				Expr dx = reduceAST(dx_);
				
				return dx;
			} else {
				Expr dx_ = div(
					integer(1),
					mul({
						// x
						u[0],
						funCall(
							"ln", {
							integer(2)
						})
					})
				);
				Expr dx = reduceAST(dx_);
				
				return dx;
			}
		}

		if(u.funName() == "exp") {
			Expr dx_ = mul({
				funCall("exp", { u[0] }),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}
	
		if(u.funName() == "tan") {
			Expr dx_ = mul({
				power(
					funCall("sec", { u[0] }),
					integer(2)
				),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "sinh") {
			Expr dx_ = mul({
				funCall("cosh", { u[0] }),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "cosh") {
			Expr dx_ = mul({
				funCall("sinh", { u[0] }),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "tanh") {
			Expr dx_ = mul({
				power(
					funCall("sech", { u[0] }),
					integer(2)
				),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "sec") {
			Expr dx_ = mul({
				funCall("sec", { u[0] }),
				funCall("tan", { u[0] }),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}
	
		if(u.funName() == "csc") {
			Expr dx_ = mul({
				integer(-1),
				funCall("cot", { u[0] }),
				funCall("csc", { u[0] }),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "cot") {
			Expr dx_ = mul({
				integer(-1),
				power(
					funCall("csc", { u[0] }),
					integer(2)
				),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "coth") {
			Expr dx_ = mul({
				integer(-1),
				power(
					funCall("csch", { u[0] }),
					integer(2)
				),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "sech") {
			Expr dx_ = mul({
				integer(-1),
				funCall("tanh", { u[0] }),
				funCall("sech", { u[0] }),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "csch") {
			Expr dx_ = mul({
				integer(-1),
				funCall("coth", { u[0] }),
				funCall("csch", { u[0] }),
				derivate(u[0], x)
			});
			Expr dx = reduceAST(dx_);
			
			return dx;
		}

		if(u.funName() == "sin") {
			Expr d_ = mul({
				funCall("cos", { u[0] }),
				derivate(u[0], x)
			});

			Expr d = reduceAST(d_);
	
			
	
			return d;
		}
	
		if(u.funName() == "cos") {
			Expr d_ = mul({
				integer(-1),
				funCall("sin", { u[0] }),
				derivate(u[0], x)
			});

			Expr d = reduceAST(d_);
	
			
	
			return d;
		}

		if(u.funName() == "arcsin")
		{
			Expr d_ = mul({
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
				derivate(u[0], x)
			});


			Expr d = reduceAST(d_);
			
			return d;
		}

		if(u.funName() == "arccos")
		{
			Expr d_ = mul({
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
				derivate(u[0], x)
			});

			Expr d = reduceAST(d_);
			
			return d;
		}


		if(u.funName() == "arctan")
		{
			Expr d_ = mul({
				div(
					integer(1),
					add({
						integer(1),
						power(symbol("x"), integer(2))
					})
				),
				derivate(u[0], x)
			});


			Expr d = reduceAST(d_);
			
			return d;
		}

		if(u.funName() == "arccot")
		{
			Expr d_ = mul({
				integer(-1),
				div(
					integer(1),
					add({
						integer(1),
						power(symbol("x"), integer(2))
					})
				),
				derivate(u[0], x)
			});

			Expr d = reduceAST(d_);
			
			return d;
		}


		if(u.funName() == "arcsec")
		{
			Expr d_ = mul({
				div(
					integer(1),
					mul({
						algebra::abs(u[0]),
						power(
							sub({
								power(symbol("x"), integer(2)),
								integer(1)
							}),
							fraction(integer(1), integer(2))
						)
					})
				),
				derivate(u[0], x),
			});

			Expr d = reduceAST(d_);
			
			return d;
		}


		if(u.funName() == "arccsc")
		{
			Expr d_ = mul({
				integer(-1),
				div(
					integer(1),
					mul({
						algebra::abs(u[0]),
						power(
							sub({
								power(symbol("x"), integer(2)),
								integer(1)
							}),
							fraction(integer(1), integer(2))
						)
					})
				),
				derivate(u[0], x),
			});

			Expr d = reduceAST(d_);
			
			return d;
		}

		if(u.funName() == "arccosh")
		{
			Expr d_ = mul({
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
				derivate(u[0], x),
			});

			Expr d = reduceAST(d_);
			
			return d;
		}
	
		if(u.funName() == "arctanh")
		{
			Expr d_ = mul({
				div(
					integer(1),
					sub({
						integer(1),
						power(symbol("x"), integer(2))
					})
				),
				derivate(u[0], x),
			});

			Expr d = reduceAST(d_);
			
			return d;
		}
	}

	return undefined();
}

Expr derivate(Expr u, Expr x) {
	Expr dx;

	if(u==(x))
	{
		return integer(1);
	}

	dx = derivateInverseTrig(u, x);
	if(dx != undefined()) 
	{
		return dx;
	}

	dx = derivatePower(u, x);
	if(dx != undefined()) 
	{
		return dx;
	}

	dx = derivateSumsAndSubs(u, x);
	if(dx != undefined()) 
	{
		return dx;
	}

	dx = derivateMul(u, x);
	if(dx != undefined()) 
	{
		return dx;
	}

	dx = derivateFuncs(u, x);
	if(dx != undefined()) 
	{
		return dx;
	}

	if(u.freeOf(x)) 
	{
		return integer(0);
	}

	return derivative(
		u, 
		x
	);
}

}
