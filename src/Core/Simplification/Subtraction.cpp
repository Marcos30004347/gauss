#include "Subtraction.hpp"
#include "Addition.hpp"
#include "Multiplication.hpp"

#include "Core/Expand/Expand.hpp"

using namespace ast;
using namespace algebra;

namespace simplification {
AST* reduceSubtractionAST(AST* u) {
    AST* sum = new AST(Kind::Addition);
		
		sum->includeOperand(u->operand(0)->deepCopy());

		for(int i=1; i<u->numberOfOperands(); i++) {
			AST* t = mul({integer(-1), u->operand(i)->deepCopy()});
			sum->includeOperand(reduceMultiplicationAST(t));
			delete t;
    }

		AST* res = reduceAdditionAST(sum);
    delete sum;
	
		return res;
}

}
