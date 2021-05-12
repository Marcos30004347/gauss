
#include "Calculus.hpp"
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

AST* integral(AST* u, AST* x) {
	return new AST(
		Kind::Integral,
		{ u, x }
	);
}

AST* integrate(AST* u, AST* x) {
	throw 'Not implemented';
}

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
