#include "List.hpp"
#include "Core/Debug/Assert.hpp"

using namespace ast;

namespace algebra {

AST* list(std::vector<AST*> e) {
	return new AST(Kind::List, e);
}

AST* listDifference(AST* L, AST* M) {
	assert(L->kind() == Kind::List,"L is not a List!\n");
	assert(M->kind() == Kind::List,"M is not a List!\n");

	AST* l = new AST(Kind::List);

	for(int i=0; i<L->numberOfOperands(); i++) {
		bool inc = false;
		for(int j=i; j<M->numberOfOperands(); j++) {
			if(L->operand(i)->match(M->operand(j))) {
				inc = true;
				break;
			}
		}
		
		if(inc) continue;
		l->includeOperand(L->operand(i)->deepCopy());
	}

	return l;
}

AST* listJoin(AST* L, AST* M) {
	assert(L->kind() == Kind::List,"L is not a list!\n");
	assert(M->kind() == Kind::List,"M is not a list!\n");

	AST* l = new AST(Kind::List);

	for(int i=0; i<L->numberOfOperands(); i++) {
		l->includeOperand(L->operand(i)->deepCopy());
	}

	for(int j=0; j<M->numberOfOperands(); j++) {
		bool inc = false;
		
		for(int i=0; i<L->numberOfOperands(); i++) {
			if(L->operand(i)->match(M->operand(j))) {
				inc = true;
				break;
			}
		}
		
		if(inc) continue;
		
		l->includeOperand(M->operand(j)->deepCopy());
	}

	return l;
}

AST* listAdjoin(AST* x, AST* L) {
	assert(L->kind() == Kind::List,"L is not a list!\n");

	AST* l = new AST(Kind::List);

	l->includeOperand(x->deepCopy());

	for(int i=0; i<L->numberOfOperands(); i++) {
		l->includeOperand(L->operand(i)->deepCopy());
	}

	return l;
}

AST* listFirst(AST* L) {
	return L->operand(0)->deepCopy();
}

AST* listRest(AST* L, int i) {
	AST* l = new AST(Kind::List);
	for(int j=i; j<L->numberOfOperands(); j++) {
		l->includeOperand(L->operand(j)->deepCopy());
	}
	return l;
}

}
