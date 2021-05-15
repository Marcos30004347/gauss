#include <assert.h>

#include "Core/Algebra/Algebra.hpp"
#include "Core/Trigonometry/Trigonometry.hpp"

using namespace ast;
using namespace algebra;
using namespace trigonometry;

void should_expand_exp() {
	AST* u0 = funCall("exp", {
		add({
			symbol("u"),		
			symbol("v"),		
		})
	});

	AST* r0 = expandExponential(u0);

	AST* k0 = mul({
		funCall("exp", { symbol("u") }),
		funCall("exp", { symbol("v") }),
	});

	assert(r0->match(k0));

	AST* u1 = funCall("exp", {
		mul({symbol("w"), symbol("u")})
	});

	AST* r1 = expandExponential(u1);
	AST* k1 = power(
		funCall("exp", {symbol("w")}),
		symbol("u")
	);

	assert(r1->match(k1));

	AST* u2 = funCall("exp", {
		mul({
			integer(2),
			add({
				symbol("x"),
				symbol("y")
			})
		})
	});

	AST* r2 = expandExponential(u2);
	AST* k2 = mul({
		power(
			funCall("exp", {symbol("x")}),
			integer(2)
		),
		power(
			funCall("exp", {symbol("y")}),
			integer(2)
		),
	});

	assert(r2->match(k2));

	delete k0;
	delete k1;
	delete k2;
	delete r0;
	delete r1;
	delete r2;
	delete u0;
	delete u1;
	delete u2;
}

void should_expand_trig() {
	AST* u0 = funCall("cos", {
		mul({integer(5), symbol("x")}),
	});
	AST* r0 = expandTrig(u0);
	AST* k0 = add({
		power(funCall("cos", {symbol("x")}), integer(5)),
		mul({
			integer(-10),
			power(funCall("cos", {symbol("x")}), integer(3)),
			power(funCall("sin", {symbol("x")}), integer(2)),
		}),
		mul({
			integer(5),
			funCall("cos", {symbol("x")}),
			power(funCall("sin", {symbol("x")}), integer(4)),
		})
	});

	assert(r0->match(k0));

	AST* u1 = funCall("sin", {
		add({
			mul({integer(2), symbol("x")}),
			mul({integer(3), symbol("y")}),
		})
	});

	AST* r1 = expandTrig(u1);
	AST* k1 = add({
		mul({
			integer(2),
			funCall("cos", {symbol("x")}),
			funCall("sin", {symbol("x")}),
			add({
				power(
					funCall("cos", {symbol("y")}),
					integer(3)
				),
				mul({
					integer(-3),
					funCall("cos", {symbol("y")}),
					power(
						funCall("sin", {symbol("y")}),
						integer(2)
					),
				})
			})
		}),
		mul({
			add({
				power(
					funCall("cos", {symbol("x")}),
					integer(2)
				),
				mul({
					integer(-1),
					power(
						funCall("sin", {symbol("x")}),
						integer(2)
					),
				})
			}),
			add({
				mul({
					integer(3),
					power(
						funCall("cos", {symbol("y")}),
						integer(2)
					),
					funCall("sin", {symbol("y")}),
				}),
				mul({
					integer(-1),
					power(
						funCall("sin", {symbol("y")}),
						integer(3)
					),
				})
			})
		})
	});

	assert(r1->match(k1));
	
	delete u0;
	delete u1;
	delete r0;
	delete r1;
	delete k0;
	delete k1;
}
void should_substitute_trig() {}

int main() {
	should_expand_exp();
	should_expand_trig();
	return 0;
}
