
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

AST* derivate(AST* u, AST* x) {
	if(u->match(x))
		return integer(1);

	if(u->kind() == Kind::Power) {
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

	if(
		u->kind() == Kind::Addition ||
		u->kind() == Kind::Subtraction
	) {
		AST* d = new AST(u->kind());
		for(int i=0; i<u->numberOfOperands(); i++) {
			d->includeOperand(derivate(u->operand(i), x));
		}
		return d;
	}

	if(u->kind() == Kind::Multiplication) {
		AST* v = u->operand(0)->deepCopy();
		
		AST* w_ = div(u->deepCopy(), v->deepCopy());
		AST* w = reduceAST(w_);
		delete w_;

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

		delete v;
		delete w;
		delete d_;

		return d;
	}

	if(u->kind() == Kind::FunctionCall) {
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
	}

	if(u->freeOf(x)) {
		return integer(0);
	}

	return derivative(u->deepCopy(), x->deepCopy());
}

}
