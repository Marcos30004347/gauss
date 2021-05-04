#include "Set.hpp"
#include "Core/Debug/Assert.hpp"
using namespace ast;

namespace algebra {

AST* set(std::vector<AST*> e) {
	return new AST(Kind::Set);
}

void combUtil(AST* ans, AST* tmp, AST* n, int left, int k) {
	if (k == 0) {
		ans->includeOperand(tmp->deepCopy());
		return;
	}

	for (int i = left; i < n->numberOfOperands(); ++i) {
		tmp->includeOperand(n->operand(i)->deepCopy());
		combUtil(ans, tmp, n, i + 1, k - 1);
		
		AST* b = tmp->operand(tmp->numberOfOperands() - 1);
		tmp->removeOperand(tmp->numberOfOperands() - 1);
		delete b;
	}
}


AST* combination(AST* n, int k) {
	AST* ans = new AST(Kind::Set);
	AST* tmp = new AST(Kind::Set);

	combUtil(ans, tmp, n, 0, k);

	delete tmp;
	return ans;
}

AST* setDifference(AST* L, AST* M) {
	assert(L->kind() == Kind::Set,"L is not a Set!\n");
	assert(M->kind() == Kind::Set,"M is not a Set!\n");

	AST* S = new AST(Kind::Set);

	for(int i=0; i<L->numberOfOperands(); i++) {
		bool inc = false;
		for(int j=0; j<M->numberOfOperands(); j++) {
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

AST* setUnion(AST* L, AST* M) {
	AST* D = setDifference(L, M);
	AST* s = new AST(Kind::Set);

	for(int i=0; i<L->numberOfOperands(); i++) {
		s->includeOperand(L->operand(i)->deepCopy());
	}

	for(int i=0; i<D->numberOfOperands(); i++) {
		s->includeOperand(D->operand(i)->deepCopy());
	}

	delete D;
	return s;
}

AST* setIntersection(AST* L, AST* M) {
	AST* D_ = setDifference(L, M);
	AST* D = setDifference(L, D_);
	delete D_;
	return D;
}

bool setExists(AST* L, AST* e) {
	for(int i=0; i<L->numberOfOperands(); i++) {
		if(L->operand(i)->match(e)) 
			return true;
	}

	return false;
}

}
