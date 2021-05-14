#ifndef MATH_ABSTRACT_SYNTAX_TREE_AST_H
#define MATH_ABSTRACT_SYNTAX_TREE_AST_H

#include <vector>
#include <string>

namespace ast {

enum Kind {
	// Terminals
	Undefined = 0,
	Integer,
	Symbol,
	Infinity,
	MinusInfinity,
	Fraction,
	Fail,

	// Operators
	Addition,
	Subtraction,
	Multiplication,
	Division,
	Power,
	Factorial,

	// Functions
	FunctionCall,
	Derivative,

	// Data Types
	List,
	Set,

};

class AST {
public:
	AST(Kind kind);
	AST(Kind kind, signed long value);
	AST(Kind kind, const char* identifier);
	AST(Kind kind, std::vector<AST*> operands);

	~AST();
	
	Kind kind();
	std::string toString();

	AST* deepCopy();

	bool match(AST* const other);
	bool analogous(AST* const other);
	bool freeOf(AST* const other);
	bool freeOfElementsInSet(AST* const set);
	bool isTerminal();

	AST* operand(signed long i);

	bool includeOperand(AST* expr);
	bool includeOperand(AST* expr, signed long i);

	bool removeOperand(AST* expr);
	bool removeOperand(signed long i);

	unsigned numberOfOperands();

	const signed long value();
	const std::string identifier();
	const std::string funName();

	AST(AST&) 	= delete;
	AST(AST&&) 	= delete;

private:
	const Kind 										_kind;
	const signed long 						_value;
	const std::string 						_identifier;
	
	std::vector<AST*> 	_operands;

	AST(Kind kind, const signed long value, const std::string identifier);
};

void destroyASTs(std::vector<AST*>);
AST* mapBinaryAST(AST* a, AST* n, AST*(*)(AST*, AST*));
AST* mapUnaryAST(AST* u, AST*(*f)(AST*));
AST* deepReplace(AST* tree, AST* subtree, AST* v);
}// ast



#endif
