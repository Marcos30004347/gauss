#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;

void should_create_lists() {
	AST* L = list({ inte(0), inte(1), inte(2) });
	
	assert(L->kind() == Kind::List);

	assert(L->numberOfOperands() == 3);

	for(int i=0; i<3; i++) {
		assert(L->operand(i)->kind() == Kind::Integer);
		assert(L->operand(i)->value() == i);
	}
	
	delete L;
}

void should_remove_elements_from_list() {
	AST* L = list({ inte(0), inte(1), inte(2) });
	AST* M = list({ inte(0), inte(2) });

	AST* K = remove(L, M);

	assert(K->kind() == Kind::List);


	assert(K->numberOfOperands() == 1);
	assert(K->operand(0)->kind() == Kind::Integer);
	assert(K->operand(0)->value() == 1);

	AST* L_ = list({ inte(0), inte(0), inte(1), inte(2) });
	AST* M_ = list({ inte(0), inte(2) });
	AST* K_ = remove(L_, M_);

	assert(K_->kind() == Kind::List);

	assert(K_->numberOfOperands() == 2);
	assert(K_->operand(0)->kind() == Kind::Integer);
	assert(K_->operand(0)->value() == 0);
	assert(K_->operand(1)->kind() == Kind::Integer);
	assert(K_->operand(1)->value() == 1);

	delete L;
	delete L_;
	delete M;
	delete M_;
	delete K;
	delete K_;
}

void should_join_elements_from_list() {
	AST* L = list({ inte(0), inte(1), inte(2) });
	AST* M = list({ inte(3), inte(4) });

	AST* K = join(L, M);

	assert(K->kind() == Kind::List);
	assert(K->numberOfOperands() == L->numberOfOperands() + M->numberOfOperands());

	for(int i=0; i<K->numberOfOperands(); i++) {
		assert(K->operand(i)->kind() == Kind::Integer);
		assert(K->operand(i)->value() == i);
	}

	delete L;
	delete M;
	delete K;
}

void should_adjoin_elements_from_list() {
	AST* L = list({inte(1), inte(2) });
	AST* x = inte(0);

	AST* K = adjoin(x, L);

	assert(K->kind() == Kind::List);
	assert(K->numberOfOperands() == L->numberOfOperands() + 1);

	for(int i=0; i<K->numberOfOperands(); i++) {
		assert(K->operand(i)->kind() == Kind::Integer);
		assert(K->operand(i)->value() == i);
	}

	delete x;
	delete L;
	delete K;
}

void should_get_rest_elements_from_list() {
	AST* L = list({ inte(0), inte(1), inte(2) });
	AST* L_ = rest(L);

	assert(L_->kind() == Kind::List);
	assert(L_->numberOfOperands() == L->numberOfOperands() - 1);
	assert(L_->operand(0)->kind() == Kind::Integer);
	assert(L_->operand(0)->value() == 1);
	assert(L_->operand(1)->kind() == Kind::Integer);
	assert(L_->operand(1)->value() == 2);

	delete L;
	delete L_;
}

void should_get_first_elements_from_list() {
	AST* L = list({ inte(0), inte(1), inte(2) });
	AST* f = first(L);

	assert(f->kind() == Kind::Integer);
	assert(f->value() == 0);
	
	delete L;
	delete f;
}

void should_compare_lists() {
	AST* L = list({ inte(0), inte(1), inte(2) });
	AST* L_ = list({ inte(0), inte(1), inte(2) });
	AST* L__ = list({ inte(2), inte(1), inte(0) });

	assert(L->match(L_));
	assert(!L->match(L__));

	delete L;
	delete L_;
	delete L__;
}

int main() {
	should_create_lists();
	should_remove_elements_from_list();
	should_join_elements_from_list();
	should_adjoin_elements_from_list();
	should_get_rest_elements_from_list();
	should_get_first_elements_from_list();
	should_compare_lists();
	return 0;
}
