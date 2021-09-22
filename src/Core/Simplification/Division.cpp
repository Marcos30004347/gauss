#include "Division.hpp"
#include "Multiplication.hpp"
#include "Power.hpp"

using namespace ast;
using namespace algebra;

namespace simplification {

AST* reduceDivisionAST(AST* u) {
	AST* p = power(u->operand(1)->copy(), integer(-1));
	AST* m = mul({u->operand(0)->copy(), reducePowerAST(p)});
	
	AST* res = reduceMultiplicationAST(m);
	
	destroyASTs({ p, m });
	
	return res;
}

}

