#include <assert.h>

#include "Multiplication.hpp"
#include "Core/Reduce/Reduce.hpp"
#include "Core/Reduce/Subtraction.hpp"
// #include "Core/Reduce/Multiplication.hpp"
// #include "Core/Reduce/Addition.hpp"

using namespace ast;
using namespace algebra;
using namespace reduce;

namespace expand {

AST* expandProductOfSums(AST* a, AST* b) {
	// assert(a->kind() == Kind::Addition || b->kind() == Kind::Addition);

	if(a->kind() == Kind::Addition && b->kind() == Kind::Addition) {
		AST* u = new AST(Kind::Addition);
	
		for(int i=0; i<a->numberOfOperands(); i++) {
			for(int j=0; j<b->numberOfOperands(); j++) {

				AST* prod = mul({ a->operand(i)->deepCopy(), b->operand(j)->deepCopy() });				
				// AST* reduced_prod = reduceMultiplicationAST(prod);

				u->includeOperand(prod);
				
				// delete prod;
			}
		}
	
		AST* r = reduceAST(u);
		delete u;
		return { u };
	}

	if(a->kind() != Kind::Addition && b->kind() == Kind::Addition) {
		AST* u = new AST(Kind::Addition);
		
		for(int j=0; j<b->numberOfOperands(); j++) {
			AST* prod = mul({ a->deepCopy(), b->operand(j)->deepCopy() });				
			// AST* reduced_prod = reduceMultiplicationAST(prod);
			u->includeOperand(prod);			
			// delete prod;
		}

		// AST* r = reduceAST(u);
		// delete u;
		// return { r };
		return u;
	}

	AST* u = new AST(Kind::Addition);

	for(int j=0; j<a->numberOfOperands(); j++) {
		AST* prod = mul({  a->operand(j)->deepCopy(), b->deepCopy() });				
		// AST* reduced_prod = reduceMultiplicationAST(prod);
		u->includeOperand(prod);			
		// delete prod;
	}

		// AST* r = reduceAST(u);
		// delete u;
		// return { r };
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
			// printf("ADDITION\n");
			signed long no = expanded->numberOfOperands();
			signed long k = (i+1) % no;

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

			expanded->includeOperand(expandProductOfSums(a, b), 0);
			
			delete a;
			delete b;
		} else if(expanded->operand(i)->kind() == Kind::Subtraction) {
			// printf("SUBTRACTION\n");
			signed long no = expanded->numberOfOperands();
			signed long k = (i+1) % no;
			
			if(i == k) continue;

			AST* a = expanded->operand(i);
			AST* b = expanded->operand(k);

			if(k>i) {
				expanded->removeOperand(k);
				expanded->removeOperand(i);
			} else {
				expanded->removeOperand(i);
				expanded->removeOperand(k);
			}

			AST* a_ = reduceSubtractionAST(a);

			// printf("a' = %s\n", a_->toString().c_str());
			// printf("b = %s\n", b->toString().c_str());

			expanded->includeOperand(expandProductOfSums(a_, b), 0);
			
			delete a_;
			delete a;
			delete b;
		}
	}

	// printf("reducing: %s\n", expanded->toString().c_str());
	// AST* res = reduceAST(expanded);
	// // printf("reduced: %s\n", res->toString().c_str());

	// delete expanded;

	return expanded;
}

}
