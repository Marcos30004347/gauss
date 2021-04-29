#include "Subtraction.hpp"
#include "Addition.hpp"
#include "Multiplication.hpp"

using namespace ast;
using namespace algebra;

namespace reduce {
AST* reduceSubtractionAST(AST* u) {
    AST* sum = new AST(Kind::Addition);
		
		sum->includeOperand(u->operand(0)->deepCopy());
    // printf("reducing sub: \n");
		// printf("%s\n", u->toString().c_str());
		// printf("%s\n", u->operand(0)->toString().c_str());
	
		for(int i=1; i<u->numberOfOperands(); i++) {
			AST* pro = mul({inte(-1), u->operand(i)->deepCopy()});
			sum->includeOperand(reduceMultiplicationAST(pro));
			delete pro;
    }

		AST* res = reduceAdditionAST(sum);
	
    delete sum;
	
		return res;
}

}
