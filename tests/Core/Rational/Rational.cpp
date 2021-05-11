#include <assert.h>
#include "Core/Rational/Rational.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Algebra/Set.hpp"

using namespace ast;
using namespace algebra;
using namespace rational;
using namespace polynomial;
using namespace simplification;

void should_get_numerator() {
	AST* u0 = add({
		symbol("a"),
		symbol("b"),
		symbol("c"),
	});

	AST* u0_ = numerator(u0);

	AST* u1 = div(
		add({
			symbol("a"),
			symbol("b"),
			symbol("c"),
		}),
		add({
			symbol("d"),
			symbol("e"),
			symbol("f"),
		})
	);

	AST* u1_ = numerator(u1);

	AST* abc = add({
		symbol("a"),
		symbol("b"),
		symbol("c"),
	});

	assert(u0_->match(abc));
	assert(u1_->match(abc));
	
	delete u0;
	delete u1;
	delete u0_;
	delete u1_;
	delete abc;
}

void should_get_denominators() {
	AST* u0 = add({
		symbol("a"),
		symbol("b"),
		symbol("c"),
	});

	AST* u0_ = denominator(u0);

	AST* u1 = div(
		add({
			symbol("a"),
			symbol("b"),
			symbol("c"),
		}),
		add({
			symbol("d"),
			symbol("e"),
			symbol("f"),
		})
	);

	AST* u1_ = denominator(u1);

	AST* def = add({
		symbol("d"),
		symbol("e"),
		symbol("f"),
	});

	AST* one = integer(1);

	assert(u0_->match(one));
	assert(u1_->match(def));
	
	delete u0;
	delete u1;
	delete u0_;
	delete u1_;
	delete def;
	delete one;
}
void should_expand_rational_expressions() {
	AST* u = add({
		div(symbol("a"), symbol("b")),
		div(symbol("c"), symbol("d")),
		div(symbol("e"), symbol("f")),
	});

	AST* u_ = rationalize(u);
	AST* v = expandRational(u_);

	AST* k = div(
		add({
			mul({symbol("b"), symbol("d"), symbol("e")}),
			mul({symbol("b"), symbol("c"), symbol("f")}),
			mul({symbol("a"), symbol("d"), symbol("f")}),
		}),
		mul({symbol("b"), symbol("d"), symbol("f") })
	);

	assert(v->match(k));

	delete u;
	delete u_;
	delete v;
	delete k;
}

int main() {
	should_get_numerator();
	should_expand_rational_expressions();
	should_get_denominators();
}
