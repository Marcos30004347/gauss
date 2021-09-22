#include "Core/AST/AST.hpp"

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
	
	AST* ast1 = ast0->copy();
	
	assert(ast0 != ast1);
	assert(ast0->match(ast1));
	
	delete ast0;
	delete ast1;
}

// void should_get_if_asts_are_analogous() {
// 	AST* a0 = new AST(Kind::Addition, {
// 		new AST(Kind::Integer, 3),
// 		new AST(Kind::Integer, 4),
// 		new AST(Kind::Integer, 5),
// 	});

// 	AST* b0 = new AST(Kind::Addition, {
// 		new AST(Kind::Integer, 1),
// 		new AST(Kind::Integer, 2),
// 		new AST(Kind::Integer, 3),
// 	});

// 	assert(a0->analogous(b0));

// 	delete a0;
// 	delete b0;


// 	AST* a1 = new AST(Kind::Addition, {
// 		new AST(Kind::Multiplication, {
// 			new AST(Kind::Integer, 2),
// 			new AST(Kind::Symbol, "x")
// 		}),
// 		new AST(Kind::Integer, 1),
// 		new AST(Kind::Integer, 5),
// 	});

// 	AST* b1 = new AST(Kind::Addition, {
// 		new AST(Kind::Multiplication, {
// 			new AST(Kind::Integer, 5),
// 			new AST(Kind::Symbol, "y")
// 		}),
// 		new AST(Kind::Integer, 7),
// 		new AST(Kind::Integer, 8),
// 	});

// 	assert(a1->analogous(b1));

// 	delete a1;
// 	delete b1;	

// 	AST* a2 = new AST(Kind::Division, {
// 		new AST(Kind::Multiplication, {
// 			new AST(Kind::Integer, 2),
// 			new AST(Kind::Symbol, "x")
// 		}),
// 		new AST(Kind::Integer, 1),
// 	});

// 	AST* b2 = new AST(Kind::Division, {
// 		new AST(Kind::Multiplication, {
// 			new AST(Kind::Integer, 5),
// 			new AST(Kind::Symbol, "y")
// 		}),
// 		new AST(Kind::Integer, 7),
// 	});

// 	assert(a2->analogous(b2));

// 	delete a2;
// 	delete b2;

// }

void should_return_if_ast_is_of_form()
{
	AST* a0 = new AST(
		Kind::Addition, {
			new AST(Kind::Integer, 1),
			new AST(Kind::Integer, 2),
			new AST(Kind::Integer, 3),
			new AST(Kind::Symbol, "x"),
	});

	AST* b0 = new AST(
		Kind::Addition, {
			new AST(Kind::Integer, 5),
			new AST(Kind::Symbol, "y"),
			new AST(Kind::Integer, 7),
			new AST(Kind::Integer, 6),
	});

	assert(a0->isOfForm(b0));

	delete a0;
	delete b0;

	AST* a1 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Integer, 4),
		new AST(Kind::Division, {
			new AST(Kind::Symbol, "x"),
			new AST(Kind::Integer, 6),
		}),
	});

	AST* b1 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Integer, 4),
		new AST(Kind::Division, {
			new AST(Kind::Symbol, "y"),
			new AST(Kind::Integer, 1),
		}),
	});

	assert(a1->isOfForm(b1));

	delete a1;
	delete b1;


	AST* a2 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::FunctionCall, {
			new AST(Kind::Symbol, "f"),
			new AST(Kind::Integer, 3),
			new AST(Kind::Symbol, "x"),
		}),
		new AST(Kind::Division, {
			new AST(Kind::Symbol, "x"),
			new AST(Kind::Integer, 6),
		}),
	});

	AST* b2 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 1),
		new AST(Kind::FunctionCall, {
			new AST(Kind::Symbol, "f"),
			new AST(Kind::Integer, 6),
			new AST(Kind::Symbol, "z"),
		}),
		new AST(Kind::Division, {
			new AST(Kind::Symbol, "y"),
			new AST(Kind::Integer, 1),
		}),
	});

	assert(a2->isOfForm(b2));

	delete a2;
	delete b2;

	AST* a3 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::FunctionCall, {
			new AST(Kind::Symbol, "f"),
			new AST(Kind::Integer, 3),
			new AST(Kind::Symbol, "x"),
		}),
		new AST(Kind::Division, {
			new AST(Kind::Symbol, "x"),
			new AST(Kind::Integer, 6),
		}),
	});

	AST* b3 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 1),
		new AST(Kind::FunctionCall, {
			new AST(Kind::Symbol, "g"),
			new AST(Kind::Integer, 6),
			new AST(Kind::Symbol, "z"),
		}),
		new AST(Kind::Division, {
			new AST(Kind::Symbol, "y"),
			new AST(Kind::Integer, 1),
		}),
	});

	// function call name dont match
	assert(!a3->isOfForm(b3));

	delete a3;
	delete b3;


	AST* a4 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::FunctionCall, {
			new AST(Kind::Symbol, "f"),
			new AST(Kind::Integer, 3),
			new AST(Kind::Symbol, "x"),
		}),
		new AST(Kind::Division, {
			new AST(Kind::Symbol, "x"),
			new AST(Kind::Integer, 6),
		}),
	});

	AST* b4 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 1),
		new AST(Kind::FunctionCall, {
			new AST(Kind::Symbol, "f"),
			new AST(Kind::Symbol, "z"),
			new AST(Kind::Integer, 6),
		}),
		new AST(Kind::Division, {
			new AST(Kind::Symbol, "y"),
			new AST(Kind::Integer, 1),
		}),
	});

	// function call argument order dont match
	assert(!a4->isOfForm(b4));

	delete a4;
	delete b4;

	AST* a5 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 3),
		new AST(Kind::Derivative, {
			new AST(Kind::Symbol, "f"),
			new AST(Kind::Symbol, "x"),
		}),
		new AST(Kind::Power, {
			new AST(Kind::Symbol, "x"),
			new AST(Kind::Integer, 6),
		}),
	});

	AST* b5 = new AST(Kind::Addition, {
		new AST(Kind::Integer, 1),
		new AST(Kind::Derivative, {
			new AST(Kind::Symbol, "z"),
			new AST(Kind::Symbol, "y"),
		}),
		new AST(Kind::Power, {
			new AST(Kind::Symbol, "y"),
			new AST(Kind::Integer, 9),
		}),
	});

	assert(a5->isOfForm(b5));
	
	delete a5;
	delete b5;

}

int main() {
	should_create_ast_nodes();
	should_match_ast_nodes();
	should_deep_copy_ast_nodes();
	should_return_if_ast_is_of_form();
	// should_get_if_asts_are_analogous();
	return 0;
}
