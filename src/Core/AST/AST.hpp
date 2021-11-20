#ifndef MATH_ABSTRACT_SYNTAX_TREE_AST_H
#define MATH_ABSTRACT_SYNTAX_TREE_AST_H

#include <vector>
#include <string>

#include "Integer.hpp"

namespace ast {

enum Kind {
	Undefined = 0,
	Integer,
	Symbol,
	Infinity,
	MinusInfinity,
	Fraction,
	Tensor,
	Matrix,
	Fail,

	Addition,
	Subtraction,
	Multiplication,
	Division,
	Power,
	Factorial,

	FunctionCall,

	Integral,
	Derivative,

	List,
	Set,

};

class AST {
public:
	AST(Kind kind);
	AST(Kind kind, Int value);
	AST(Kind kind, const char* identifier);
	AST(Kind kind, std::vector<AST*> operands);

	~AST();
	
	Kind kind() const;
	std::string toString();

	AST* copy();

	bool match(AST* const other);
	// bool analogous(AST* const other);
	bool freeOf(AST* const other);
	bool freeOfElementsInSet(AST* const set);
	bool isTerminal();

	bool is(signed long i);

	bool isNot(signed long i);

	bool isOfForm(AST* const form, AST* const syms = nullptr);
	
	AST* operand(unsigned long i);
	AST* operand(Int i);
	AST* operandList();
	AST* symbols();
	
	bool includeOperand(AST* expr);
	bool includeOperand(AST* expr, signed long i);
	bool includeOperand(AST* expr, Int i);

	bool removeOperand(AST* expr);
	bool removeOperand(signed long i);
	bool removeOperand(Int i);
	bool deleteOperand(signed long i);
	bool deleteOperand(Int i);

	unsigned numberOfOperands() const;

	Int value() const;

	const std::string identifier();
	const std::string funName();

	AST(AST&) 	= delete;
	AST(AST&&) 	= delete;

private:
	std::vector<AST*> _operands;
	Kind 							_kind;
	Int								_value;
	std::string 			_identifier;

	AST(Kind kind, const Int value, const std::string identifier);
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
