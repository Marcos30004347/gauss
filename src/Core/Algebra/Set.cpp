#include "Set.hpp"
#include "Core/Debug/Assert.hpp"

using namespace ast;

namespace algebra {

AST* set(std::vector<AST*> e) {
	AST* S = new AST(Kind::Set);

	for(long unsigned int i=0; i<e.size(); i++) {
		S->includeOperand(e[i]);
	}

	return S;
}

void combUtil(AST* ans, AST* tmp, AST* n, int left, int k) {
	if (k == 0) {
		ans->includeOperand(tmp->deepCopy());
		return;
	}

	for (unsigned int i = left; i < n->numberOfOperands(); ++i) {
		tmp->includeOperand(n->operand(i)->deepCopy());
		combUtil(ans, tmp, n, i + 1, k - 1);
		
		AST* b = tmp->operand(tmp->numberOfOperands() - 1);
		tmp->removeOperand(tmp->numberOfOperands() - 1);
		delete b;
	}
}


AST* combination(AST* n, AST* k) {
	assert(k->kind() == Kind::Integer, "k needs to be an integer!");

	AST* ans = new AST(Kind::Set);
	AST* tmp = new AST(Kind::Set);

	combUtil(ans, tmp, n, 0, k->value());

	delete tmp;
	return ans;
}

AST* difference(AST* L, AST* M) {
	assert(L->kind() == Kind::Set, "L is not a Set!\n");
	assert(M->kind() == Kind::Set, "M is not a Set!\n");

	AST* S = new AST(Kind::Set);

	for(unsigned int i=0; i<L->numberOfOperands(); i++) {
		bool inc = false;

		for(unsigned int j=0; j<M->numberOfOperands(); j++) {
			if(L->operand(i)->match(M->operand(j))) {
				inc = true;
				break;
			}
		}
		
		if(inc) continue;
		S->includeOperand(L->operand(i)->deepCopy());
	}

	return S;
}

AST* unification(AST* L, AST* M) {
	AST* D = difference(M, L);

	AST* S = new AST(Kind::Set);

	for(unsigned int i=0; i<L->numberOfOperands(); i++) {
		S->includeOperand(L->operand(i)->deepCopy());
	}

	for(unsigned int i=0; i<D->numberOfOperands(); i++) {
		S->includeOperand(D->operand(i)->deepCopy());
	}

	delete D;
	return S;
}

AST* intersection(AST* L, AST* M) {
	AST* D_ = difference(L, M);
	AST* D = difference(L, D_);
	delete D_;
	return D;
}

bool exists(AST* L, AST* e) {
	for(unsigned int i=0; i<L->numberOfOperands(); i++) {
		if(L->operand(i)->match(e)) 
			return true;
	}

	return false;
}


AST* cleanUp(AST* C, AST* t) {
	assert(C->kind() == Kind::Set, "C is not a set!");
	assert(t->kind() == Kind::Set, "t is not a set!");
	
	AST* C_ = new AST(Kind::Set);
	
	for(unsigned int i=0; i<C->numberOfOperands(); i++) {
		bool inc = false;
	
		for(unsigned int j=0; j<C->operand(i)->numberOfOperands(); j++) {
	
			for(unsigned int k=0; k<t->numberOfOperands(); k++) {
				if(t->operand(i)->match(C->operand(i)->operand(j))) {
					inc = true;
					break;
				}
			}
	
			if(inc) break;
		}
	
		if(inc) continue;
	
		C_->includeOperand(C->operand(i)->deepCopy());
	}

	return C_;
}

}
