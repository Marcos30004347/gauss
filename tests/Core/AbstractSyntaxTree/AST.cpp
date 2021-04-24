#include "Core/AbstractSyntaxTree/AST.hpp"

#include <assert.h>

using namespace ast;

void should_create_ast_nodes() {
	// Undefined
	AST* ast0 = new AST(Kind::Undefined);

	// 3 + 4 + 5
	AST* ast1 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Integer, 4),
		new AST(Kind::Integer, 5),
	});

	// 3 + 4 + 5/6
	AST* ast2 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Integer, 4),
		new AST(Kind::Fraction, {
			new AST(Kind::Integer, 5),
			new AST(Kind::Integer, 6),
		}),
	});

	assert(ast0->kind() == Kind::Undefined);

	assert(ast1->kind() == Kind::Addition);
	assert(ast1->operand(0)->kind() == Kind::Integer);
	assert(ast1->operand(0)->value() == 3);
	assert(ast1->operand(1)->kind() == Kind::Integer);
	assert(ast1->operand(1)->value() == 4);
	assert(ast1->operand(2)->kind() == Kind::Integer);
	assert(ast1->operand(2)->value() == 5);

	assert(ast2->kind() == Kind::Addition);
	assert(ast2->operand(0)->kind() == Kind::Integer);
	assert(ast2->operand(0)->value() == 3);
	assert(ast2->operand(1)->kind() == Kind::Integer);
	assert(ast2->operand(1)->value() == 4);
	assert(ast2->operand(2)->kind() == Kind::Fraction);
	assert(ast2->operand(2)->operand(0)->kind() == Kind::Integer);
	assert(ast2->operand(2)->operand(0)->value() == 5);
	assert(ast2->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(ast2->operand(2)->operand(1)->value() == 6);

	delete ast0;
	delete ast1;
	delete ast2;
}

void should_match_ast_nodes() {
	AST* ast0 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Integer, 4),
		new AST(Kind::Integer, 5),
	});

	AST* ast1 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Integer, 4),
		new AST(Kind::Integer, 5),
	});

	AST* ast2 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Integer, 4),
	});

	assert(ast0->match(ast1));
	assert(!ast0->match(ast2));
	
	delete ast0;
	delete ast1;
	delete ast2;
}


void should_deep_copy_ast_nodes() {
	AST* ast0 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Integer, 4),
		new AST(Kind::Integer, 5),
	});
	
	AST* ast1 = ast0->deepCopy();
	
	assert(ast0 != ast1);
	assert(ast0->match(ast1));
	
	delete ast0;
	delete ast1;
}

int main() {
	should_create_ast_nodes();
	should_match_ast_nodes();
	should_deep_copy_ast_nodes();
	return 0;
}
