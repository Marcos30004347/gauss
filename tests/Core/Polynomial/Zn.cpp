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
	AST* x = symb("x");
	AST* u = add({
		pow(symb("x"), inte(5)),
		mul({inte(2), pow(symb("x"), inte(3))}),
		mul({inte(2), pow(symb("x"), inte(2))}),
		symb("x"),
		inte(2)
	});

	AST* v = add({
		pow(symb("x"), inte(4)),
		mul({inte(2), pow(symb("x"), inte(3))}),
		mul({inte(2), pow(symb("x"), inte(2))})
	});

	std::vector<AST*> rr = extendedEuclideanAlgGPE_Zp(u, v, symb("x"), 3);

	printf("gcd = %s\n", rr[0]->toString().c_str());
	printf("A = %s\n", rr[1]->toString().c_str());
	printf("B = %s\n", rr[2]->toString().c_str());

	AST* k_ = add({
		mul({rr[1]->deepCopy(), u->deepCopy()}),
		mul({rr[2]->deepCopy(), v->deepCopy()})
	});


	AST* k = Tnn(k_, x, 3);
	printf("A*U + B*V in Z_%i = %s\n", 3, k->toString().c_str());

	AST* u_ = add({
		pow(symb("x"), inte(2)),
		mul({inte(5), symb("x")}),
		inte(8)
	});

	AST* u_Znn5 = Tnn(u_, x, 5);
	AST* u_Zs5 = Ts(u_, x, 5);

	printf("%s\n", u_Znn5->toString().c_str());
	printf("%s\n", u_Zs5->toString().c_str());
	return 0;
}

