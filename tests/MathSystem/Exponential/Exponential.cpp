#include <assert.h>

#include "MathSystem/Algebra/Algebra.hpp"
#include "MathSystem/Exponential/Exponential.hpp"

using namespace ast;
using namespace algebra;
using namespace exponential;

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

void should_contract_exponential() {
	
	AST* u0 = mul({
		funCall("exp", { symbol("u") }),
		funCall("exp", { symbol("v") }),
	});

	AST* r0 = contractExponential(u0);
	AST* k0 = funCall("exp", { add({ symbol("u"), symbol("v") })});
	
	assert(r0->match(k0));

	AST* u1 = power(
		funCall("exp", { symbol("u") }),
		symbol("w")
	);

	AST* r1 = contractExponential(u1);

	
	AST* k1 = funCall("exp", {
		mul({symbol("u"), symbol("w")})
	});

	assert(r1->match(k1));


	AST* u2 = power(
		funCall("exp", {
			funCall("exp", {symbol("x")})
		}),
		funCall("exp", {symbol("y")})
	);

	AST* r2 = contractExponential(u2);

	AST* k2 = funCall("exp", {
		funCall("exp", {
			add({
				symbol("x"),
				symbol("y"),
			})
		})
	});

	assert(r2->match(k2));

	delete u0;
	delete u1;
	delete u2;
	delete r0;
	delete r1;
	delete r2;
	delete k0;
	delete k1;
	delete k2;
}


int main() {

	should_expand_exp();
	should_contract_exponential();

	return 0;
}
