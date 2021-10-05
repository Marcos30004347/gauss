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

// void should_get_berlekamp_factors() 
// {
// 	AST* u = add({
// 		power(symbol("x"), integer(6)),
// 		power(symbol("x"), integer(5)),
// 		symbol("x"),
// 		integer(4)
// 	});

// 	AST* x = symbol("x");
// 	AST* p = integer(5);

// 	AST* factors = berlekampFactors(u, x, 5);

//   factorization::berlekamp(u, x, p);

// 	AST* F = set({
// 		add({
// 			integer(3),
// 			mul({
// 				integer(2),
// 				symbol("x"),
// 			}),
// 			power(
// 				symbol("x"),
// 				integer(2)
// 			)
// 		}),
// 		add({
// 			integer(4),
// 			mul({
// 				integer(3),
// 				symbol("x"),
// 			}),
// 			power(
// 				symbol("x"),
// 				integer(2)
// 			)
// 		}),
// 		add({
// 			integer(2),
// 			symbol("x"),
// 			power(
// 				symbol("x"),
// 				integer(2)
// 			)
// 		}),
// 	});

// 	assert(factors->match(F));

// 	delete u;
// 	delete p;
// 	delete x;
// 	delete F;
// 	delete factors;
// }

// void should_factorize_with_berlekamp()
// {

// 	factorization::buildBerkelampMatrix(add({power(symbol("x"), integer(4)), power(symbol("x"), integer(2)), symbol("x"), integer(1) }), symbol("x"), integer(2));
// 	factorization::buildBerkelampMatrix(add({power(symbol("x"), integer(8)), power(symbol("x"), integer(7)), power(symbol("x"), integer(4)), power(symbol("x"), integer(3)), symbol("x"), integer(1) }), symbol("x"), integer(3));

// 	AST* ax = add({
// 		power(
// 			symbol("x"),
// 			integer(6)
// 		),
// 		mul({
// 			integer(-3),
// 			power(
// 				symbol("x"),
// 				integer(5)
// 			)
// 		}),
// 		power(
// 			symbol("x"),
// 			integer(4)
// 		),
// 		mul({
// 			integer(-3),
// 			power(
// 				symbol("x"),
// 				integer(3)
// 			)
// 		}),
// 		mul({
// 			integer(-1),
// 			power(
// 				symbol("x"),
// 				integer(2)
// 			)
// 		}),
// 		mul({
// 			integer(-3),
// 			symbol("x"),
// 		}),
// 		integer(1)
// 	});

// 	AST* x = symbol("x");
// 	AST* q = integer(11);

// 	AST* f = berlekampFactors(ax, x, 11);
//   factorization::berlekamp(ax, x, q);

// 	AST* F = set({
// 		add({
// 			symbol("x"),
// 			integer(1),
// 		}),
// 		add({
// 			power(symbol("x"), integer(2)),
// 			mul({
// 				integer(5),
// 				symbol("x")
// 			}),
// 			integer(3)
// 		}),
// 		add({
// 			power(symbol("x"), integer(3)),
// 			mul({
// 				integer(2),
// 				power(symbol("x"), integer(2))
// 			}),
// 			mul({
// 				integer(3),
// 				symbol("x")
// 			}),
// 			integer(4)
// 		}),
// 	});

// 	assert(f->match(F));

// 	delete ax;
// 	delete f;
// 	delete x;
// 	delete q;
// 	delete F;
// }

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

	AST* k = add({
		power(symbol("x"), integer(5)),
		mul({integer(6), power(symbol("x"), integer(4))}),
		mul({integer(4), power(symbol("x"), integer(3))}),
		mul({integer(4), power(symbol("x"), integer(2))}),
		mul({integer(4), symbol("x")}),
		integer(6)
	});

	AST* x = symbol("x");
	AST* r = integer(2);
	AST* q = integer(3);
	AST* z = integer(7);
	
	AST* A = factorization::buildBerkelampMatrix(p, x, r);

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
	
	AST* Ab = buildBerlekampBasis(A, r);

	assert(Ab->kind() == Kind::List);
	assert(Ab->numberOfOperands() == 2);

	assert(Ab->operand(0)->kind() == Kind::List);
	assert(Ab->operand(0)->numberOfOperands() == 4);
	assert(Ab->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(Ab->operand(0)->operand(0)->value() == 1);
	assert(Ab->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(Ab->operand(0)->operand(1)->value() == 0);
	assert(Ab->operand(0)->operand(2)->kind() == Kind::Integer);
	assert(Ab->operand(0)->operand(2)->value() == 0);
	assert(Ab->operand(0)->operand(3)->kind() == Kind::Integer);
	assert(Ab->operand(0)->operand(3)->value() == 0);

	assert(Ab->operand(1)->kind() == Kind::List);
	assert(Ab->operand(1)->numberOfOperands() == 4);
	assert(Ab->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(Ab->operand(1)->operand(0)->value() == 0);
	assert(Ab->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(Ab->operand(1)->operand(1)->value() == 0);
	assert(Ab->operand(1)->operand(2)->kind() == Kind::Integer);
	assert(Ab->operand(1)->operand(2)->value() == 1);
	assert(Ab->operand(1)->operand(3)->kind() == Kind::Integer);
	assert(Ab->operand(1)->operand(3)->value() == 1);


	delete A;
	delete Ab;


	AST* B = factorization::buildBerkelampMatrix(t, x, q);

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


	AST* Bb = buildBerlekampBasis(B, q);

	assert(Bb->kind() == Kind::List);
	assert(Bb->numberOfOperands() == 3);

	assert(Bb->operand(0)->kind() == Kind::List);
	assert(Bb->operand(0)->numberOfOperands() == 8);
	assert(Bb->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(Bb->operand(0)->operand(0)->value() == 1);
	assert(Bb->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(Bb->operand(0)->operand(1)->value() == 0);
	assert(Bb->operand(0)->operand(2)->kind() == Kind::Integer);
	assert(Bb->operand(0)->operand(2)->value() == 0);
	assert(Bb->operand(0)->operand(3)->kind() == Kind::Integer);
	assert(Bb->operand(0)->operand(3)->value() == 0);
	assert(Bb->operand(0)->operand(4)->kind() == Kind::Integer);
	assert(Bb->operand(0)->operand(4)->value() == 0);
	assert(Bb->operand(0)->operand(5)->kind() == Kind::Integer);
	assert(Bb->operand(0)->operand(5)->value() == 0);
	assert(Bb->operand(0)->operand(6)->kind() == Kind::Integer);
	assert(Bb->operand(0)->operand(6)->value() == 0);
	assert(Bb->operand(0)->operand(7)->kind() == Kind::Integer);
	assert(Bb->operand(0)->operand(7)->value() == 0);

	assert(Bb->operand(2)->kind() == Kind::List);
	assert(Bb->operand(2)->numberOfOperands() == 8);
	assert(Bb->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(Bb->operand(2)->operand(0)->value() == 0);
	assert(Bb->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(Bb->operand(2)->operand(1)->value() == 0);
	assert(Bb->operand(2)->operand(2)->kind() == Kind::Integer);
	assert(Bb->operand(2)->operand(2)->value() == 0);
	assert(Bb->operand(2)->operand(3)->kind() == Kind::Integer);
	assert(Bb->operand(2)->operand(3)->value() == 1);
	assert(Bb->operand(2)->operand(4)->kind() == Kind::Integer);
	assert(Bb->operand(2)->operand(4)->value() == 0);
	assert(Bb->operand(2)->operand(5)->kind() == Kind::Integer);
	assert(Bb->operand(2)->operand(5)->value() == 0);
	assert(Bb->operand(2)->operand(6)->kind() == Kind::Integer);
	assert(Bb->operand(2)->operand(6)->value() == 0);
	assert(Bb->operand(2)->operand(7)->kind() == Kind::Integer);
	assert(Bb->operand(2)->operand(7)->value() == 1);

	assert(Bb->operand(1)->kind() == Kind::List);
	assert(Bb->operand(1)->numberOfOperands() == 8);
	assert(Bb->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(Bb->operand(1)->operand(0)->value() == 0);
	assert(Bb->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(Bb->operand(1)->operand(1)->value() == 2);
	assert(Bb->operand(1)->operand(2)->kind() == Kind::Integer);
	assert(Bb->operand(1)->operand(2)->value() == 2);
	assert(Bb->operand(1)->operand(3)->kind() == Kind::Integer);
	assert(Bb->operand(1)->operand(3)->value() == 1);
	assert(Bb->operand(1)->operand(4)->kind() == Kind::Integer);
	assert(Bb->operand(1)->operand(4)->value() == 1);
	assert(Bb->operand(1)->operand(5)->kind() == Kind::Integer);
	assert(Bb->operand(1)->operand(5)->value() == 1);
	assert(Bb->operand(1)->operand(6)->kind() == Kind::Integer);
	assert(Bb->operand(1)->operand(6)->value() == 1);
	assert(Bb->operand(1)->operand(7)->kind() == Kind::Integer);
	assert(Bb->operand(1)->operand(7)->value() == 0);

	delete B;
	delete Bb;


	AST* C = factorization::buildBerkelampMatrix(k, x, z);


	assert(C->operand(0)->kind() == Kind::List);
	assert(C->numberOfOperands() == 5);

	assert(C->operand(0)->numberOfOperands() == 5);
	assert(C->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(C->operand(0)->operand(0)->value() == 1);
	assert(C->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(C->operand(0)->operand(1)->value() == 0);
	assert(C->operand(0)->operand(2)->kind() == Kind::Integer);
	assert(C->operand(0)->operand(2)->value() == 0);
	assert(C->operand(0)->operand(3)->kind() == Kind::Integer);
	assert(C->operand(0)->operand(3)->value() == 0);
	assert(C->operand(0)->operand(4)->kind() == Kind::Integer);
	assert(C->operand(0)->operand(4)->value() == 0);

	assert(C->operand(1)->numberOfOperands() == 5);
	assert(C->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(C->operand(1)->operand(0)->value() == 4);
	assert(C->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(C->operand(1)->operand(1)->value() == 6);
	assert(C->operand(1)->operand(2)->kind() == Kind::Integer);
	assert(C->operand(1)->operand(2)->value() == 2);
	assert(C->operand(1)->operand(3)->kind() == Kind::Integer);
	assert(C->operand(1)->operand(3)->value() == 4);
	assert(C->operand(1)->operand(4)->kind() == Kind::Integer);
	assert(C->operand(1)->operand(4)->value() == 3);

	assert(C->operand(2)->numberOfOperands() == 5);
	assert(C->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(C->operand(2)->operand(0)->value() == 2);
	assert(C->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(C->operand(2)->operand(1)->value() == 3);
	assert(C->operand(2)->operand(2)->kind() == Kind::Integer);
	assert(C->operand(2)->operand(2)->value() == 6);
	assert(C->operand(2)->operand(3)->kind() == Kind::Integer);
	assert(C->operand(2)->operand(3)->value() == 1);
	assert(C->operand(2)->operand(4)->kind() == Kind::Integer);
	assert(C->operand(2)->operand(4)->value() == 4);

	assert(C->operand(3)->numberOfOperands() == 5);
	assert(C->operand(3)->operand(0)->kind() == Kind::Integer);
	assert(C->operand(3)->operand(0)->value() == 6);
	assert(C->operand(3)->operand(1)->kind() == Kind::Integer);
	assert(C->operand(3)->operand(1)->value() == 3);
	assert(C->operand(3)->operand(2)->kind() == Kind::Integer);
	assert(C->operand(3)->operand(2)->value() == 5);
	assert(C->operand(3)->operand(3)->kind() == Kind::Integer);
	assert(C->operand(3)->operand(3)->value() == 3);
	assert(C->operand(3)->operand(4)->kind() == Kind::Integer);
	assert(C->operand(3)->operand(4)->value() == 1);

	assert(C->operand(4)->numberOfOperands() == 5);
	assert(C->operand(4)->operand(0)->kind() == Kind::Integer);
	assert(C->operand(4)->operand(0)->value() == 1);
	assert(C->operand(4)->operand(1)->kind() == Kind::Integer);
	assert(C->operand(4)->operand(1)->value() == 5);
	assert(C->operand(4)->operand(2)->kind() == Kind::Integer);
	assert(C->operand(4)->operand(2)->value() == 5);
	assert(C->operand(4)->operand(3)->kind() == Kind::Integer);
	assert(C->operand(4)->operand(3)->value() == 6);
	assert(C->operand(4)->operand(4)->kind() == Kind::Integer);
	assert(C->operand(4)->operand(4)->value() == 6);

	AST* Cb = buildBerlekampBasis(C, z);

	assert(Cb->kind() == Kind::List);
	assert(Cb->numberOfOperands() == 2);

	assert(Cb->operand(0)->kind() == Kind::List);
	assert(Cb->operand(0)->numberOfOperands() == 5);
	assert(Cb->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(Cb->operand(0)->operand(0)->value() == 1);
	assert(Cb->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(Cb->operand(0)->operand(1)->value() == 0);
	assert(Cb->operand(0)->operand(2)->kind() == Kind::Integer);
	assert(Cb->operand(0)->operand(2)->value() == 0);
	assert(Cb->operand(0)->operand(3)->kind() == Kind::Integer);
	assert(Cb->operand(0)->operand(3)->value() == 0);
	assert(Cb->operand(0)->operand(4)->kind() == Kind::Integer);
	assert(Cb->operand(0)->operand(4)->value() == 0);

	assert(Cb->operand(1)->kind() == Kind::List);
	assert(Cb->operand(1)->numberOfOperands() == 5);
	assert(Cb->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(Cb->operand(1)->operand(0)->value() == 0);
	assert(Cb->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(Cb->operand(1)->operand(1)->value() == 5);
	assert(Cb->operand(1)->operand(2)->kind() == Kind::Integer);
	assert(Cb->operand(1)->operand(2)->value() == 6);
	assert(Cb->operand(1)->operand(3)->kind() == Kind::Integer);
	assert(Cb->operand(1)->operand(3)->value() == 5);
	assert(Cb->operand(1)->operand(4)->kind() == Kind::Integer);
	assert(Cb->operand(1)->operand(4)->value() == 1);

	delete C;
	delete Cb;

	delete p;
	delete t;
	delete k;
	delete x;
	delete r;
	delete q;
	delete z;
}

void should_factorize_square_free_poly_with_berlekamp()
{

	AST* ux = add({
		power(symbol("x"), integer(6)),
		power(symbol("x"), integer(5)),
		power(symbol("x"), integer(4)),
		power(symbol("x"), integer(3)),
		integer(1)
	});

	AST* x = symbol("x");
	
	AST* p0 = integer(2);

	berlekampFactors(ux, x, p0);

	delete ux;
	delete p0;
	delete x;
}

int main()
{

	// AST* fx = add({
	// 	power(symbol("x"), integer(6)),
	// 	power(symbol("x"), integer(5)),
	// 	power(symbol("x"), integer(4)),
	// 	power(symbol("x"), integer(3)),
	// 	integer(1)
	// });
	// printf("%s\n", remainderGPE_Zp(power(symbol("x"), integer(0)), fx, symbol("x"), 2)->toString().c_str());
	// printf("%s\n", remainderGPE_Zp(power(symbol("x"), integer(2)), fx, symbol("x"), 2)->toString().c_str());
	// printf("%s\n", remainderGPE_Zp(power(symbol("x"), integer(4)), fx, symbol("x"), 2)->toString().c_str());
	// printf("%s\n", remainderGPE_Zp(power(symbol("x"), integer(6)), fx, symbol("x"), 2)->toString().c_str());
	// printf("%s\n", remainderGPE_Zp(power(symbol("x"), integer(8)), fx, symbol("x"), 2)->toString().c_str());
	// printf("%s\n", remainderGPE_Zp(power(symbol("x"), integer(10)), fx, symbol("x"), 2)->toString().c_str());
	// printf("%s\n", remainderGPE_Zp(power(symbol("x"), integer(12)), fx, symbol("x"), 2)->toString().c_str());

	should_form_berkelamp_basis_matrix();
	should_factorize_square_free_poly_with_berlekamp();

	return 0;
}
