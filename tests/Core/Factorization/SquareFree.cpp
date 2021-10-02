#include <assert.h>

#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Factorization/SquareFree.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_get_square_free_factorization()
{
	AST* ax = add({
		power(symbol("x"), integer(8)),
		mul({
			integer(-2),
			power(symbol("x"), integer(6))
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(2))
		}),
		integer(-1)
	});
	AST* x = symbol("x");

	AST* sf_ax0 = squareFreeFactorization(ax, x);
	AST* sf_ax3 = squareFreeFactorization2(ax, x);

	AST* sf_ax_r = mul({
		power(
			add({
				integer(-1),
				power(symbol("x"), integer(2))
			}),
			integer(3)
		),
		add({
			integer(1),
			power(symbol("x"), integer(2))
		})
	});

	assert(sf_ax0->match(sf_ax_r));
	assert(sf_ax3->match(sf_ax_r));

	AST* sf_ax1 = squareFreeFactorization2(ax, x);

	assert(sf_ax1->match(sf_ax_r));

	AST* bx = add({
		power(symbol("x"), integer(11)),
		mul({integer(2), power(symbol("x"), integer(9))}),
		mul({integer(2), power(symbol("x"), integer(8))}),
		power(symbol("x"), integer(6)),
		power(symbol("x"), integer(5)),
		mul({integer(2), power(symbol("x"), integer(3))}),
		mul({integer(2), power(symbol("x"), integer(2))}),
		integer(1)
	});

	AST* q = power(integer(3), integer(1));

	AST* sf_bx0 = squareFreeFactorizationFiniteField(bx, x, q);

	AST* sf_bx_r = mul({
		add({symbol("x"), integer(1)}),
		power(
			add({
				power(symbol("x"), integer(2)),
				integer(1)
			}), 
			integer(3)
		),
		power(add({symbol("x"), integer(2)}), integer(4)),
	});

	assert(sf_bx0->match(sf_bx_r));

	delete x;
	delete q;
	
	delete ax;
	delete bx;

	delete sf_ax3;
	delete sf_bx0;
	delete sf_ax0;
	delete sf_ax1;
	delete sf_ax_r;
	delete sf_bx_r;
}

int main()
{
	should_get_square_free_factorization();
}
