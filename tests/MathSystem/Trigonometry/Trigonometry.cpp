#include <assert.h>

#include "MathSystem/Algebra/Algebra.hpp"
#include "MathSystem/Trigonometry/Trigonometry.hpp"

using namespace ast;
using namespace algebra;
using namespace trigonometry;

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

void should_contract_trig() {

	AST* u0 = mul({
		add({
			funCall("sin", {symbol("x")}),
			funCall("cos", {symbol("y")}),
		}),
		funCall("cos", {symbol("y")}),
	});

	AST* r0 = contractTrig(u0);
	AST* k0 = add({
		fraction(1,2),
		mul({
			fraction(1,2),
			funCall("cos", {
				mul({ integer(2), symbol("y") })
			})
		}),
		mul({
			fraction(1,2),
			funCall("sin", {
				add({ symbol("x"), symbol("y") })
			})
		}),
		mul({
			fraction(1,2),
			funCall("sin", {
				add({ symbol("x"), mul({ integer(-1), symbol("y")}) })
			})
		}),
	});

	assert(r0->match(k0));

	delete u0;
	delete r0;
	delete k0;
}

void should_substitute_trig() {
	AST* u0 = funCall("tan", { symbol("x") });
	AST* u1 = funCall("cot", { symbol("x") });
	AST* u2 = funCall("sec", { symbol("x") });
	AST* u3 = funCall("csc", { symbol("x") });
	AST* r0 = substituteTrig(u0);
	AST* r1 = substituteTrig(u1);
	AST* r2 = substituteTrig(u2);
	AST* r3 = substituteTrig(u3);
	AST* k0 = div(
		funCall("sin", { symbol("x") }),
		funCall("cos", { symbol("x") })
	);
	AST* k1 = div(
		funCall("cos", { symbol("x") }),
		funCall("sin", { symbol("x") })
	);
	AST* k2 = div(
		integer(1),
		funCall("cos", { symbol("x") })
	);
	AST* k3 = div(
		integer(1),
		funCall("sin", { symbol("x") })
	);

	assert(r0->match(k0));
	assert(r1->match(k1));
	assert(r2->match(k2));
	assert(r3->match(k3));
	
	delete u0;
	delete u1;
	delete u2;
	delete u3;
	delete r0;
	delete r1;
	delete r2;
	delete r3;
	delete k0;
	delete k1;
	delete k2;
	delete k3;
}

int main() {

	should_substitute_trig();
	should_expand_trig();
	should_contract_trig();

	return 0;
}
