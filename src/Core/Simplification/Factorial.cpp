#include <assert.h>

#include "Factorial.hpp"

using namespace ast;
using namespace algebra;

namespace simplification {

signed long factorial(signed long n) {
	if(n == 0 || n == 1)
		return 1;
	
	return n * factorial(n - 1);
}

AST* reduceFactorialAST(AST* u) {
	assert(u->operand(0)->kind() == Kind::Integer);
	assert(u->operand(0)->value() >= 0);
	return new AST(Kind::Integer, factorial(u->operand(0)->value()));
}

}
