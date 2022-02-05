#include "MathSystem/Algebra/Algebra.hpp"
#include "MathSystem/Expand/Expand.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;
using namespace expand;

void should_expand_products() {
	Expr exp0 = (Expr("x") + 4) * (Expr("x") + 3);
	Expr exp1 = (Expr("x") + 4) * (Expr("x") + 3) * (Expr("x") + 5);

	assert(expandAST(exp0) == 12 + 7*Expr("x") + power(Expr("x"), 2));
	assert(expandAST(exp1) == 60 + 47*Expr("x") + 12*power(Expr("x"), 2) + power(Expr("x"), 3));
}

int main() {
	should_expand_products();

	return 0;
}
