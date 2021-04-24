#ifndef MATH_ABSTRACT_SYNTAX_TREE_AST_H
#define MATH_ABSTRACT_SYNTAX_TREE_AST_H

#include <list>
#include <string>

namespace ast {

enum Kind {
	// Types
	Undefined = 0,
	Integer,
	Fraction,
	Function,

	// Operators
	Symbol,
	Addition,
	Subtraction,
	Multiplication,
	Division,
	Power,
	Factorial
};

class AST {
public:
	AST(Kind kind);
	AST(Kind kind, signed long value);
	AST(Kind kind, const char* identifier);
	AST(Kind kind, std::list<AST*> operands);

	~AST();
	
	Kind kind();
	void print();

	AST* deepCopy();

	bool match(AST* const other);

	AST* operand(signed long i);

	void includeOperand(AST* expr);
	void includeOperand(AST* expr, signed long i);

	void removeOperand(AST* expr);
	void removeOperand(signed long i);

	unsigned numberOfOperands();

	const signed long value();
	const std::string identifier();

	AST(AST&) 	= delete;
	AST(AST&&) 	= delete;

private:
	const Kind 										_kind;
	const signed long 						_value;
	const std::string 						_identifier;
	
	std::list<AST*> 	_operands;

	AST(Kind kind, const signed long value, const std::string identifier);
};

void destroyASTs(std::list<AST*>);
AST* mapBinaryAST(AST* a, AST* n, AST*(*)(AST*, AST*));


}// ast



#endif
