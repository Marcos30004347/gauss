#include <assert.h>

#include "Factorial.hpp"
#include "Core/Expand/Expand.hpp"

using namespace ast;
using namespace expand;
using namespace algebra;

namespace simplification {

// signed long factorial(signed long n) {
// 	if(n == 0 || n == 1)
// 		return 1;
	
// 	return n * factorial(n - 1);
// }

Expr reduceFactorialAST(Expr u) {
	assert(isConstant(u[0]));

	Expr u_ = expandAST(u[0]);
	
	assert(u_.kind() == Kind::Integer);
	assert(u_.value() >= 0);

	return Expr(Kind::Integer, fact(u_.value()));
}

}
