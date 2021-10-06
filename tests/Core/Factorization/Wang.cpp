#include <assert.h>

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Factorization.hpp"
#include "Core/Factorization/Wang.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_get_nondivisors()
{
	AST* F = list({
		integer(-14),
		integer(3),
		integer(-11),
		integer(-17),
	});

	AST* L = list({
		symbol("y"),
		symbol("z"),
	});

	AST* K = symbol("Z");

	long d[4];

	int success = nondivisors(4, F, 1, L, K, d);

	assert(success = 1);

	printf("%li\n", d[0]);
	printf("%li\n", d[1]);
	printf("%li\n", d[2]);
	printf("%li\n", d[3]);

	delete F;

}

int main()
{
	should_get_nondivisors();
	return 0;
}
