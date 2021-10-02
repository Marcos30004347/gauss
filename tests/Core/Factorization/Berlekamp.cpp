#include <assert.h>

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Polynomial/Factorization.hpp"
#include "Core/Polynomial/Zp.hpp"
#include "Core/Factorization/Berlekamp.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_get_berlekamp_factors() 
{
	AST* u = add({
		power(symbol("x"), integer(6)),
		power(symbol("x"), integer(5)),
		symbol("x"),
		integer(4)
	});

	AST* x = symbol("x");
	AST* p = integer(5);

	AST* factors = berlekampFactors(u, x, 5);

  factorization::berlekamp(u, x, p);

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
	delete p;
	delete x;
	delete F;
	delete factors;
}

void should_factorize_with_berlekamp()
{

	factorization::buildBerkelampBasisMatrix(add({power(symbol("x"), integer(4)), power(symbol("x"), integer(2)), symbol("x"), integer(1) }), symbol("x"), integer(2));
	factorization::buildBerkelampBasisMatrix(add({power(symbol("x"), integer(8)), power(symbol("x"), integer(7)), power(symbol("x"), integer(4)), power(symbol("x"), integer(3)), symbol("x"), integer(1) }), symbol("x"), integer(3));

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

	AST* x = symbol("x");
	AST* q = integer(11);

	AST* f = berlekampFactors(ax, x, 11);
  factorization::berlekamp(ax, x, q);

	AST* F = set({
		add({
			symbol("x"),
			integer(1),
		}),
		add({
			power(symbol("x"), integer(2)),
			mul({
				integer(5),
				symbol("x")
			}),
			integer(3)
		}),
		add({
			power(symbol("x"), integer(3)),
			mul({
				integer(2),
				power(symbol("x"), integer(2))
			}),
			mul({
				integer(3),
				symbol("x")
			}),
			integer(4)
		}),
	});

	assert(f->match(F));

	delete ax;
	delete f;
	delete x;
	delete q;
	delete F;
}

void should_form_berkelamp_basis_matrix()
{
	AST* p = add({
		power(symbol("x"), integer(4)),
		power(symbol("x"), integer(2)),
		symbol("x"),
		integer(1)
	});

	AST* t = add({
		power(symbol("x"), integer(8)),
		power(symbol("x"), integer(7)),
		power(symbol("x"), integer(4)),
		power(symbol("x"), integer(3)),
		symbol("x"),
		integer(1)
	});

	AST* x = symbol("x");
	AST* r = integer(2);
	AST* q = integer(3);
	
	AST* A = factorization::buildBerkelampBasisMatrix(p, x, r);

	assert(A->operand(0)->kind() == Kind::List);
	assert(A->numberOfOperands() == 4);

	assert(A->operand(0)->numberOfOperands() == 4);
	assert(A->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(A->operand(0)->operand(0)->value() == 1);
	assert(A->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(A->operand(0)->operand(1)->value() == 0);
	assert(A->operand(0)->operand(2)->kind() == Kind::Integer);
	assert(A->operand(0)->operand(2)->value() == 0);
	assert(A->operand(0)->operand(3)->kind() == Kind::Integer);
	assert(A->operand(0)->operand(3)->value() == 0);

	assert(A->operand(1)->numberOfOperands() == 4);
	assert(A->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(A->operand(1)->operand(0)->value() == 0);
	assert(A->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(A->operand(1)->operand(1)->value() == 0);
	assert(A->operand(1)->operand(2)->kind() == Kind::Integer);
	assert(A->operand(1)->operand(2)->value() == 1);
	assert(A->operand(1)->operand(3)->kind() == Kind::Integer);
	assert(A->operand(1)->operand(3)->value() == 0);

	assert(A->operand(2)->numberOfOperands() == 4);
	assert(A->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(A->operand(2)->operand(0)->value() == 1);
	assert(A->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(A->operand(2)->operand(1)->value() == 1);
	assert(A->operand(2)->operand(2)->kind() == Kind::Integer);
	assert(A->operand(2)->operand(2)->value() == 1);
	assert(A->operand(2)->operand(3)->kind() == Kind::Integer);
	assert(A->operand(2)->operand(3)->value() == 0);

	assert(A->operand(3)->numberOfOperands() == 4);
	assert(A->operand(3)->operand(0)->kind() == Kind::Integer);
	assert(A->operand(3)->operand(0)->value() == 1);
	assert(A->operand(3)->operand(1)->kind() == Kind::Integer);
	assert(A->operand(3)->operand(1)->value() == 1);
	assert(A->operand(3)->operand(2)->kind() == Kind::Integer);
	assert(A->operand(3)->operand(2)->value() == 0);
	assert(A->operand(3)->operand(3)->kind() == Kind::Integer);
	assert(A->operand(3)->operand(3)->value() == 1);
	R = new int*[8];
	R[0] = new int[8];
	R[1] = new int[8];
	R[2] = new int[8];
	R[3] = new int[8];
	R[4] = new int[8];
	R[5] = new int[8];
	R[6] = new int[8];
	R[7] = new int[8];

	R[0][0] = 0;
	R[0][1] = 0;
	R[0][2] = 1;
	R[0][3] = 1;

	R[1][0] = 0;
	R[1][1] = 1;
	R[1][2] = 1;
	R[1][3] = 1;
	
	R[2][0] = 0;
	R[2][1] = 1;
	R[2][2] = 0;
	R[2][3] = 0;

	R[3][0] = 0;
	R[3][1] = 0;
	R[3][2] = 0;
	R[3][3] = 0;

	printf("%s\n", auxiliaryBasis(x, integer(4), 2)->toString().c_str());

	std::vector<std::vector<int>> ASD = {
		{0, 0, 0, 0, 0, 0, 0, 0},
		{0, 2, 0, 1, 0, 0, 0, 0},
		{0, 0, 2, 0, 0, 0, 1, 0},
		{1, 0, 2, 0, 0, 2, 0, 1},
		{0, 1, 0, 0, 0, 2, 0, 0},
		{1, 1, 0, 1, 2, 2, 0, 2},
		{1, 0, 0, 0, 1, 0, 1, 0},
		{2, 0, 1, 0, 0, 1, 0, 2},
	};
	printf("**asdasdasdasd*\n");
	
	for(int i = 0; i < 8; i++)
		for(int j = 0; j < 8; j++)
			R[i][j] = ASD[j][i];
	
	printf("**asdasdasdasd*\n");
	printf("%s\n", auxiliaryBasis(x, integer(8), 3)->toString().c_str());
	
	printf("%s\n", A->toString().c_str());
	// printf("%s\n", nullSpace_sZp(A, r->value())->toString().c_str());
	printf("%s\n", nullSpace_Zp(A, r->value())->toString().c_str());
	printf("***********\n");

	AST* B = factorization::buildBerkelampBasisMatrix(t, x, q);

	assert(B->operand(0)->kind() == Kind::List);
	assert(B->numberOfOperands() == 8);

	assert(B->operand(0)->numberOfOperands() == 8);
	assert(B->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(B->operand(0)->operand(0)->value() == 1);
	assert(B->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(0)->operand(1)->value() == 0);
	assert(B->operand(0)->operand(2)->kind() == Kind::Integer);
	assert(B->operand(0)->operand(2)->value() == 0);
	assert(B->operand(0)->operand(3)->kind() == Kind::Integer);
	assert(B->operand(0)->operand(3)->value() == 0);
	assert(B->operand(0)->operand(4)->kind() == Kind::Integer);
	assert(B->operand(0)->operand(4)->value() == 0);
	assert(B->operand(0)->operand(5)->kind() == Kind::Integer);
	assert(B->operand(0)->operand(5)->value() == 0);
	assert(B->operand(0)->operand(6)->kind() == Kind::Integer);
	assert(B->operand(0)->operand(6)->value() == 0);
	assert(B->operand(0)->operand(7)->kind() == Kind::Integer);
	assert(B->operand(0)->operand(7)->value() == 0);

	assert(B->operand(1)->numberOfOperands() == 8);
	assert(B->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(0)->value() == 0);
	assert(B->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(1)->value() == 0);
	assert(B->operand(1)->operand(2)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(2)->value() == 0);
	assert(B->operand(1)->operand(3)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(3)->value() == 1);
	assert(B->operand(1)->operand(4)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(4)->value() == 0);
	assert(B->operand(1)->operand(5)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(5)->value() == 0);
	assert(B->operand(1)->operand(6)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(6)->value() == 0);
	assert(B->operand(1)->operand(7)->kind() == Kind::Integer);
	assert(B->operand(1)->operand(7)->value() == 0);

	assert(B->operand(2)->numberOfOperands() == 8);
	assert(B->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(B->operand(2)->operand(0)->value() == 0);
	assert(B->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(2)->operand(1)->value() == 0);
	assert(B->operand(2)->operand(2)->kind() == Kind::Integer);
	assert(B->operand(2)->operand(2)->value() == 0);
	assert(B->operand(2)->operand(3)->kind() == Kind::Integer);
	assert(B->operand(2)->operand(3)->value() == 0);
	assert(B->operand(2)->operand(4)->kind() == Kind::Integer);
	assert(B->operand(2)->operand(4)->value() == 0);
	assert(B->operand(2)->operand(5)->kind() == Kind::Integer);
	assert(B->operand(2)->operand(5)->value() == 0);
	assert(B->operand(2)->operand(6)->kind() == Kind::Integer);
	assert(B->operand(2)->operand(6)->value() == 1);
	assert(B->operand(2)->operand(7)->kind() == Kind::Integer);
	assert(B->operand(2)->operand(7)->value() == 0);

	assert(B->operand(3)->numberOfOperands() == 8);
	assert(B->operand(3)->operand(0)->kind() == Kind::Integer);
	assert(B->operand(3)->operand(0)->value() == 1);
	assert(B->operand(3)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(3)->operand(1)->value() == 0);
	assert(B->operand(3)->operand(2)->kind() == Kind::Integer);
	assert(B->operand(3)->operand(2)->value() == 2);
	assert(B->operand(3)->operand(3)->kind() == Kind::Integer);
	assert(B->operand(3)->operand(3)->value() == 1);
	assert(B->operand(3)->operand(4)->kind() == Kind::Integer);
	assert(B->operand(3)->operand(4)->value() == 0);
	assert(B->operand(3)->operand(5)->kind() == Kind::Integer);
	assert(B->operand(3)->operand(5)->value() == 2);
	assert(B->operand(3)->operand(6)->kind() == Kind::Integer);
	assert(B->operand(3)->operand(6)->value() == 0);
	assert(B->operand(3)->operand(7)->kind() == Kind::Integer);
	assert(B->operand(3)->operand(7)->value() == 1);

	assert(B->operand(4)->numberOfOperands() == 8);
	assert(B->operand(4)->operand(0)->kind() == Kind::Integer);
	assert(B->operand(4)->operand(0)->value() == 0);
	assert(B->operand(4)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(4)->operand(1)->value() == 1);
	assert(B->operand(4)->operand(2)->kind() == Kind::Integer);
	assert(B->operand(4)->operand(2)->value() == 0);
	assert(B->operand(4)->operand(3)->kind() == Kind::Integer);
	assert(B->operand(4)->operand(3)->value() == 0);
	assert(B->operand(4)->operand(4)->kind() == Kind::Integer);
	assert(B->operand(4)->operand(4)->value() == 1);
	assert(B->operand(4)->operand(5)->kind() == Kind::Integer);
	assert(B->operand(4)->operand(5)->value() == 2);
	assert(B->operand(4)->operand(6)->kind() == Kind::Integer);
	assert(B->operand(4)->operand(6)->value() == 0);
	assert(B->operand(4)->operand(7)->kind() == Kind::Integer);
	assert(B->operand(4)->operand(7)->value() == 0);

	assert(B->operand(5)->numberOfOperands() == 8);
	assert(B->operand(5)->operand(0)->kind() == Kind::Integer);
	assert(B->operand(5)->operand(0)->value() == 1);
	assert(B->operand(5)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(5)->operand(1)->value() == 1);
	assert(B->operand(5)->operand(2)->kind() == Kind::Integer);
	assert(B->operand(5)->operand(2)->value() == 0);
	assert(B->operand(5)->operand(3)->kind() == Kind::Integer);
	assert(B->operand(5)->operand(3)->value() == 1);
	assert(B->operand(5)->operand(4)->kind() == Kind::Integer);
	assert(B->operand(5)->operand(4)->value() == 2);
	assert(B->operand(5)->operand(5)->kind() == Kind::Integer);
	assert(B->operand(5)->operand(5)->value() == 0);
	assert(B->operand(5)->operand(6)->kind() == Kind::Integer);
	assert(B->operand(5)->operand(6)->value() == 0);
	assert(B->operand(5)->operand(7)->kind() == Kind::Integer);
	assert(B->operand(5)->operand(7)->value() == 2);

	assert(B->operand(6)->numberOfOperands() == 8);
	assert(B->operand(6)->operand(0)->kind() == Kind::Integer);
	assert(B->operand(6)->operand(0)->value() == 1);
	assert(B->operand(6)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(6)->operand(1)->value() == 0);
	assert(B->operand(6)->operand(2)->kind() == Kind::Integer);
	assert(B->operand(6)->operand(2)->value() == 0);
	assert(B->operand(6)->operand(3)->kind() == Kind::Integer);
	assert(B->operand(6)->operand(3)->value() == 0);
	assert(B->operand(6)->operand(4)->kind() == Kind::Integer);
	assert(B->operand(6)->operand(4)->value() == 1);
	assert(B->operand(6)->operand(5)->kind() == Kind::Integer);
	assert(B->operand(6)->operand(5)->value() == 0);
	assert(B->operand(6)->operand(6)->kind() == Kind::Integer);
	assert(B->operand(6)->operand(6)->value() == 2);
	assert(B->operand(6)->operand(7)->kind() == Kind::Integer);
	assert(B->operand(6)->operand(7)->value() == 0);

	assert(B->operand(7)->numberOfOperands() == 8);
	assert(B->operand(7)->operand(0)->kind() == Kind::Integer);
	assert(B->operand(7)->operand(0)->value() == 2);
	assert(B->operand(7)->operand(1)->kind() == Kind::Integer);
	assert(B->operand(7)->operand(1)->value() == 0);
	assert(B->operand(7)->operand(2)->kind() == Kind::Integer);
	assert(B->operand(7)->operand(2)->value() == 1);
	assert(B->operand(7)->operand(3)->kind() == Kind::Integer);
	assert(B->operand(7)->operand(3)->value() == 0);
	assert(B->operand(7)->operand(4)->kind() == Kind::Integer);
	assert(B->operand(7)->operand(4)->value() == 0);
	assert(B->operand(7)->operand(5)->kind() == Kind::Integer);
	assert(B->operand(7)->operand(5)->value() == 1);
	assert(B->operand(7)->operand(6)->kind() == Kind::Integer);
	assert(B->operand(7)->operand(6)->value() == 0);
	assert(B->operand(7)->operand(7)->kind() == Kind::Integer);
	assert(B->operand(7)->operand(7)->value() == 0);

	delete p;
	delete t;
	delete x;
	delete r;
	delete q;
	delete A;
	delete B;
}

int main()
{
	should_form_berkelamp_basis_matrix();
}
