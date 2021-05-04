#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/Algebra.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;

void should_create_sets() {
	AST* S = set({ inte(0), inte(1), inte(2) });

	assert(S->kind() == Kind::Set);

	assert(S->operand(0)->kind() == Kind::Integer);
	assert(S->operand(0)->value() == 0);
	assert(S->operand(1)->kind() == Kind::Integer);
	assert(S->operand(1)->value() == 1);
	assert(S->operand(2)->kind() == Kind::Integer);
	assert(S->operand(2)->value() == 2);
	
	delete S;
}

void should_not_include_duplicates_on_set() {
	AST* S = set({inte(0), inte(1), inte(2)});
	AST* zero = inte(0);
	
	S->includeOperand(zero);
	assert(S->numberOfOperands() == 3);
	
	delete S;
	delete zero;
}

void should_get_set_difference() {
	AST* S = set({inte(0), inte(1), inte(2)});
	AST* S_ = set({inte(0), inte(2)});
	AST* M = difference(S, S_);

	assert(M->kind() == Kind::Set);
	assert(M->numberOfOperands() == 1);
	assert(M->operand(0)->kind() == Kind::Integer);
	assert(M->operand(0)->value() == 1);

	delete S;
	delete S_;
	delete M;
}

void should_get_set_unnification() {
	AST* S = set({inte(0), inte(1), inte(2)});
	AST* S_ = set({inte(0), inte(3)});
	AST* M = unification(S, S_);

	assert(M->numberOfOperands() == 4);
	assert(M->kind() == Kind::Set);
	assert(M->operand(0)->kind() == Kind::Integer);
	assert(M->operand(0)->value() == 0);
	assert(M->operand(1)->kind() == Kind::Integer);
	assert(M->operand(1)->value() == 1);
	assert(M->operand(2)->kind() == Kind::Integer);
	assert(M->operand(2)->value() == 2);
	assert(M->operand(3)->kind() == Kind::Integer);
	assert(M->operand(3)->value() == 3);

	delete S;
	delete S_;
	delete M;
}

void should_get_set_intersection() {
	AST* S = set({inte(0), inte(1), inte(2)});
	AST* S_ = set({inte(0), inte(3)});
	AST* M = intersection(S, S_);

	assert(M->numberOfOperands() == 1);
	assert(M->kind() == Kind::Set);
	assert(M->operand(0)->kind() == Kind::Integer);
	assert(M->operand(0)->value() == 0);
	
	delete S;
	delete S_;
	delete M;
}

void should_get_if_element_exists_inside_set() {
	AST* S = set({inte(0), inte(1), inte(2)});
	AST* one = inte(1);
	
	assert(exists(S, one));
	
	delete S;
	delete one;
}


void should_get_combinations_of_elements() {
	AST* S = set({inte(0), inte(1), inte(3)});
	AST* two = inte(2);
	AST* S_ = combination(S, two);

	assert(S_->kind() == Kind::Set);
	assert(S_->numberOfOperands() == 3); // {0,1}, {0,3}, {1, 3}

	assert(S_->operand(0)->kind() == Kind::Set);
	assert(S_->operand(0)->numberOfOperands() == two->value()); // {0,1}
	assert(S_->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(S_->operand(0)->operand(0)->value() == 0);
	assert(S_->operand(0)->operand(1)->kind() == Kind::Integer);
	assert(S_->operand(0)->operand(1)->value() == 1);

	assert(S_->operand(1)->kind() == Kind::Set);
	assert(S_->operand(1)->numberOfOperands() == two->value()); // {0,3}
	assert(S_->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(S_->operand(1)->operand(0)->value() == 0);
	assert(S_->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(S_->operand(1)->operand(1)->value() == 3);

	assert(S_->operand(2)->kind() == Kind::Set);
	assert(S_->operand(2)->numberOfOperands() == two->value()); // {1,3}
	assert(S_->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(S_->operand(2)->operand(0)->value() == 1);
	assert(S_->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(S_->operand(2)->operand(1)->value() == 3);

	delete two;
	delete S;
	delete S_;
}

int main() {
	should_create_sets();
	should_not_include_duplicates_on_set();
	should_get_set_difference();
	should_get_set_unnification();
	should_get_set_intersection();
	should_get_if_element_exists_inside_set();
	should_get_combinations_of_elements();

	return 0;
}
