#include <assert.h>

#include "Core/Factorization/Utils.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Zp.hpp"
#include "Core/Algebra/List.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_repeat_square_integer()
{
	AST* i = integer(8);
	AST* x = symbol("x");

	AST* p = repeatSquaring(i, x, 13, 17);

	assert(p->kind() == Kind::Integer);
	assert(p->value() == -8);
	
	delete i;
	delete x;
	delete p;

		AST* ax = add({
		power(
			symbol("x"),
			integer(6)
		),
		mul({
			integer(-3),
			power(
				symbol("x"),
				integer(5)
			)
		}),
		power(
			symbol("x"),
			integer(4)
		),
		mul({
			integer(-3),
			power(
				symbol("x"),
				integer(3)
			)
		}),
		mul({
			integer(-1),
			power(
				symbol("x"),
				integer(2)
			)
		}),
		mul({
			integer(-3),
			symbol("x"),
		}),
		integer(1)
	});
}

int main()
{
	should_repeat_square_integer();
}
