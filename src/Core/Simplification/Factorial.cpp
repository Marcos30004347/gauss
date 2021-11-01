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

AST* reduceFactorialAST(AST* u) {
	assert(isConstant(u->operand(0)));

	AST* u_ = expandAST(u->operand(0));
	
	assert(u_->kind() == Kind::Integer);
	assert(u_->value() >= 0);

	return new AST(Kind::Integer, fact(u_->value()));
}

}
