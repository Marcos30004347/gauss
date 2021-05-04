#include "Core/Polynomial/Zp.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Simplification/Addition.hpp"

#include <assert.h>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace simplification;
using namespace polynomial;



int main() {

	
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
			pow(symb("𝑥"), inte(2)),
			inte(-2)
		});

	AST* a = symb("𝛼");

	std::vector<AST*> res = algPolynomialDivisionAST(u, v, x, p, a);


	printf("%s\n", res[0]->toString().c_str());
	printf("%s\n", res[1]->toString().c_str());
	printf("\n\n\n");

	AST* u_ = add({
		pow(symb("x"), inte(2)),
		mul({
			sub({
				inte(-1),
				symb("a")
			}),
			symb("x")
		})
	});

	AST* v_ = add({
		pow(symb("x"), inte(2)),
		mul({
			sub({
				inte(-2),
				mul({inte(2), symb("a")})
			}),
			symb("x")
		}),
		inte(3),
		mul({inte(2), symb("a")})
	});

	AST* p_ = sub({
		pow(symb("x"), inte(2)),
		inte(2),
	});

	AST* x_ = symb("x");
	AST* a_ = symb("a");
	AST* k = algPolynomialGCDAST(u_, v_, x_, p_, a_);

	printf("RES = %s\n", k->toString().c_str());

	// delete u;
	// delete v;
	// delete x;
	// delete p;
	// delete a;

	// for(AST* l : res)
	// 	delete l;

}
