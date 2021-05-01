#include <assert.h>

#include "Multiplication.hpp"
// #include "Core/Simplification/Simplification.hpp"
// #include "Core/Simplification/Subtraction.hpp"
// #include "Core/Simplification/Multiplication.hpp"
// #include "Core/Simplification/Addition.hpp"

using namespace ast;
using namespace algebra;
// using namespace simplification;

namespace expand {

AST* expandProductOfSums(AST* a, AST* b) {
	assert(a->kind() == Kind::Addition || b->kind() == Kind::Addition);

	if(a->kind() == Kind::Addition && b->kind() == Kind::Addition) {
		AST* u = new AST(Kind::Addition);
	
		for(int i=0; i<a->numberOfOperands(); i++) {
			for(int j=0; j<b->numberOfOperands(); j++) {
				u->includeOperand(mul({
					a->operand(i)->deepCopy(),
					b->operand(j)->deepCopy()
				}));
			}
		}

		return { u };
	}

	if(a->kind() != Kind::Addition && b->kind() == Kind::Addition) {
		AST* u = new AST(Kind::Addition);
		
		for(int j=0; j<b->numberOfOperands(); j++) {
			AST* prod = mul({ a->deepCopy(), b->operand(j)->deepCopy() });				
			u->includeOperand(mul({
				a->deepCopy(),
				b->operand(j)->deepCopy()
			}));			
		}

		return u;
	}

	AST* u = new AST(Kind::Addition);

	for(int j=0; j<a->numberOfOperands(); j++) {
		u->includeOperand(mul({
			a->operand(j)->deepCopy(),
			b->deepCopy()
		}));			
	}

	return u;
}

AST* expandMultiplicationAST(AST* u){

	if(u->numberOfOperands() == 1) {
		AST* res = u->operand(0)->deepCopy();
		return res;
	}

	AST* expanded = u->deepCopy();

	for(int i=0; i<expanded->numberOfOperands(); i++) {
		if(expanded->operand(i)->kind() == Kind::Addition) {
			signed long no = expanded->numberOfOperands();
			
			int k = (i+1) % no;
	
			if(i==k) continue;

			AST* a = expanded->operand(i);
			AST* b = expanded->operand(k);

			if(k>i) {
				expanded->removeOperand(k);
				expanded->removeOperand(i);
			} else {
				expanded->removeOperand(i);
				expanded->removeOperand(k);
			}

			i = -1;

			expanded->includeOperand(expandProductOfSums(a, b), 0);
			
			delete a;
			delete b;
		} 
	}

	return expanded;
}

}
