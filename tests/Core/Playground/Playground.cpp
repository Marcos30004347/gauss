#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Expand/Expand.hpp"

#include <assert.h>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace polynomial;


int main() {

// -2x + ax - (-2 + a)(x^1) - (ax + a - a(x^1))(1 + -1a)(x^(1 - 1))
	// expandAST(
	// 	sub({
	// 		add({
	// 			mul({
	// 				inte(-2),
	// 				symb("x")
	// 			}),
	// 			mul({
	// 				symb("a"),
	// 				symb("x")
	// 			})
	// 		}),
	// 		mul({
	// 			add({
	// 				inte(-2),
	// 				symb("a")
	// 			}),
	// 			pow(symb("x"), inte(1))
	// 		})
	// 	})
	// );
	// AST* lc0 = leadingCoefficientGPE(inte(-2), symb("x"));
	// printf("lc = %s\n\n", lc0->toString().c_str());

	// AST* lc1 = leadingCoefficientGPE(add({ pow(symb("x"), inte(2)) }), symb("x"));
	// printf("lc = %s\n\n", lc1->toString().c_str());

	// AST* lc2 = leadingCoefficientGPE(add({ pow(symb("x"), inte(2)), mul({inte(7), pow(symb("x"), inte(3))}) }), symb("x"));
	// printf("lc = %s\n\n", lc2->toString().c_str());
	AST* u = add({
			mul({inte(2), pow(symb("𝑥"), inte(2))}),
			mul({symb("𝛼"), symb("𝑥")})
		});
	
	AST* v = add({
			mul({symb("𝛼"), symb("𝑥")}),
			symb("𝛼")
		});

	AST* x = symb("𝑥");

	AST* p = add({
			pow(symb("𝛼"), inte(2)),
			inte(-2)
		});

	AST* a = symb("𝛼");

	std::vector<AST*> res = algPolynomialDivisionAST(u, v, x, p, a);

	printf("%s\n", res[0]->toString().c_str());
	printf("%s\n", res[1]->toString().c_str());

	delete u;
	delete v;
	delete x;
	delete p;
	delete a;

	for(AST* l : res)
		delete l;

}
