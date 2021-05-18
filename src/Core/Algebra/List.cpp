#include "List.hpp"
#include "Core/Debug/Assert.hpp"

using namespace ast;

namespace algebra {

AST* list(std::vector<AST*> e) {
	return new AST(Kind::List, e);
}

AST* remove(AST* L, AST* M) {
	assert(L->kind() == Kind::List,"L is not a List!\n");
	assert(M->kind() == Kind::List,"M is not a List!\n");

	AST* l = new AST(Kind::List);

	int p = 0;

	for(unsigned int i=0; i<L->numberOfOperands(); i++) {
		bool inc = false;
		
		for(unsigned int j=p; j<M->numberOfOperands(); j++) {
			if(L->operand(i)->match(M->operand(j))) {
				p = j+1;
				inc = true;
				break;
			}
		}
		
		if(inc)
			continue;

		l->includeOperand(L->operand(i)->deepCopy());
	}

	return l;
}

AST* join(AST* L, AST* M) {
	assert(L->kind() == Kind::List,"L is not a list!\n");
	assert(M->kind() == Kind::List,"M is not a list!\n");

	AST* l = new AST(Kind::List);

	for(unsigned int i=0; i<L->numberOfOperands(); i++) {
		l->includeOperand(L->operand(i)->deepCopy());
	}

	for(unsigned int j=0; j<M->numberOfOperands(); j++) {
		bool inc = false;
		
		for(unsigned int i=0; i<L->numberOfOperands(); i++) {
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

AST* adjoin(AST* x, AST* L, AST* (*f)(AST* const)) {
	assert(L->kind() == Kind::List,"L is not a list!\n");
	if(f) {
		if(L->numberOfOperands() > 0) {
			AST* K = list({ x->deepCopy(), L->operand(0)->deepCopy() });
			AST* r = f(K);
			delete K;

			AST* l = new AST(Kind::List);

			if(r->kind() == Kind::List) {
				for(unsigned int i=0; i<r->numberOfOperands(); i++) {
					l->includeOperand(r->operand(i)->deepCopy());
				}
				for(unsigned int i=1; i<L->numberOfOperands(); i++) {
					l->includeOperand(L->operand(i)->deepCopy());
				}
		
				delete r;

				return l;
			}

			l->includeOperand(r->deepCopy());
			delete r;
			for(unsigned int i=1; i<L->numberOfOperands(); i++) {
				l->includeOperand(L->operand(i)->deepCopy());
			}

			return l;
		}
	
		AST* l = new AST(Kind::List);
		l->includeOperand(x->deepCopy());

		return l;
	}

	AST* l = new AST(Kind::List);

	l->includeOperand(x->deepCopy());

	for(unsigned int i=0; i<L->numberOfOperands(); i++) {
		l->includeOperand(L->operand(i)->deepCopy());
	}

	return l;
}

AST* first(AST* L) {
	return L->operand(0)->deepCopy();
}

AST* rest(AST* L, int i) {
	AST* l = new AST(Kind::List);
	for(unsigned int j=i; j<L->numberOfOperands(); j++) {
		l->includeOperand(L->operand(j)->deepCopy());
	}
	return l;
}

}
