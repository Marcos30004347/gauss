#ifndef MATH_ABSTRACT_SYNTAX_TREE_AST_H
#define MATH_ABSTRACT_SYNTAX_TREE_AST_H

#include <vector>
#include <string>

namespace ast {

enum Kind {
	Undefined = 0,
	Integer,
	Symbol,
	Infinity,
	MinusInfinity,
	Fraction,
	Tensor,
	Fail,

	Addition,
	Subtraction,
	Multiplication,
	Division,
	Power,
	Factorial,

	FunctionCall,

	Derivative,

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

	AST* operand(unsigned long i);
	AST* operandList();

	bool includeOperand(AST* expr);
	bool includeOperand(AST* expr, signed long i);

	bool removeOperand(AST* expr);
	bool removeOperand(signed long i);

	unsigned numberOfOperands();

	signed long value();
	const std::string identifier();
	const std::string funName();

	AST(AST&) 	= delete;
	AST(AST&&) 	= delete;

private:
	std::vector<AST*> _operands;
	Kind 							_kind;
	signed long 			_value;
	std::string 			_identifier;
	

	AST(Kind kind, const signed long value, const std::string identifier);
};

void destroyASTs(std::vector<AST*>);
AST* mapBinaryAST(AST* a, AST* n, AST*(*)(AST*, AST*));
AST* mapUnaryAST(AST* u, AST*(*f)(AST*));
AST* deepReplace(AST* tree, AST* subtree, AST* v);
AST* construct(Kind kind, AST* L);

template<typename... types>
AST* mapAST(AST*(*f)(AST*, types ... args), AST* u, types ... params) {
	if(
			u->kind() == Kind::Integer ||
			u->kind() == Kind::Fraction ||
			u->kind() == Kind::Symbol ||
			u->kind() == Kind::Infinity ||
			u->kind() == Kind::MinusInfinity
	){
		return f(u, params...);
	}

	if(u->numberOfOperands() == 0) {
		return f(u, params...);
	}

	AST* t = new AST(u->kind());

	if(u->kind() == Kind::FunctionCall) {
		t->includeOperand(new AST(Kind::Symbol, u->funName().c_str()));
	}

	for(unsigned int i=0; i< u->numberOfOperands(); i++) {
		t->includeOperand(f(u->operand(i), params...));
	}

	return t;
}

}// ast

#endif
