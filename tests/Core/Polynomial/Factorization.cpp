#include "Core/Algebra/Set.hpp"	
#include "Core/Algebra/Algebra.hpp"	
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Factorization.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;
using namespace polynomial;

void should_get_r_matrix() {

	AST* u = add({
		power(symbol("x"), integer(6)),
		power(symbol("x"), integer(5)),
		symbol("x"),
		integer(4)
	});

	AST* x = symbol("x");

	AST* n = integer(6);

	RMatrix(u, x, n, 5);

	// 0 0 1 0 1 0
	// 0 4 3 3 1 1
	// 0 0 1 1 2 0
	// 0 0 3 4 2 0
	// 0 0 2 4 1 0
	// 0 1 3 3 1 4

	assert(getRMatrixValue(0,0) == 0);
	assert(getRMatrixValue(0,1) == 0);
	assert(getRMatrixValue(0,2) == 1);
	assert(getRMatrixValue(0,3) == 0);
	assert(getRMatrixValue(0,4) == 1);
	assert(getRMatrixValue(0,5) == 0);

	assert(getRMatrixValue(1,0) == 0);
	assert(getRMatrixValue(1,1) == 4);
	assert(getRMatrixValue(1,2) == 3);
	assert(getRMatrixValue(1,3) == 3);
	assert(getRMatrixValue(1,4) == 1);
	assert(getRMatrixValue(1,5) == 1);

	assert(getRMatrixValue(2,0) == 0);
	assert(getRMatrixValue(2,1) == 0);
	assert(getRMatrixValue(2,2) == 1);
	assert(getRMatrixValue(2,3) == 1);
	assert(getRMatrixValue(2,4) == 2);
	assert(getRMatrixValue(2,5) == 0);

	assert(getRMatrixValue(3,0) == 0);
	assert(getRMatrixValue(3,1) == 0);
	assert(getRMatrixValue(3,2) == 3);
	assert(getRMatrixValue(3,3) == 4);
	assert(getRMatrixValue(3,4) == 2);
	assert(getRMatrixValue(3,5) == 0);

	assert(getRMatrixValue(4,0) == 0);
	assert(getRMatrixValue(4,1) == 0);
	assert(getRMatrixValue(4,2) == 2);
	assert(getRMatrixValue(4,3) == 4);
	assert(getRMatrixValue(4,4) == 1);
	assert(getRMatrixValue(4,5) == 0);

	assert(getRMatrixValue(5,0) == 0);
	assert(getRMatrixValue(5,1) == 1);
	assert(getRMatrixValue(5,2) == 3);
	assert(getRMatrixValue(5,3) == 3);
	assert(getRMatrixValue(5,4) == 1);
	assert(getRMatrixValue(5,5) == 4);

	destroyRMatrix(6);

	delete u;
	delete x;
	delete n;

}

void should_get_berlekamp_factors() {
	AST* u = add({
		power(symbol("x"), integer(6)),
		power(symbol("x"), integer(5)),
		symbol("x"),
		integer(4)
	});

	AST* x = symbol("x");
	AST* factors = berlekampFactor(u, x, 5);

	AST* F = set({
		add({
			integer(3), 
			mul({
				integer(2),
				symbol("x"),
			}),
			power(
				symbol("x"),
				integer(2)
			)
		}),
		add({
			integer(4), 
			mul({
				integer(3),
				symbol("x"),
			}),
			power(
				symbol("x"),
				integer(2)
			)
		}),
		add({
			integer(2), 
			symbol("x"),
			power(
				symbol("x"),
				integer(2)
			)
		}),
	});

	assert(factors->match(F));

	delete u;
	delete x;
	delete F;
	delete factors;
}

int main() {

	should_get_r_matrix();
	should_get_berlekamp_factors();
	return 0;
}
