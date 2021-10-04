#include "AST.hpp"
#include <assert.h>

namespace ast {

AST::AST(Kind kind) {
	this->_kind = kind;
	this->_operands = std::vector<AST*>(0);
	this->_identifier = "";
	this->_value = 0;
}

AST::AST(Kind kind, signed long value) {
	this->_kind = kind;
	this->_operands = std::vector<AST*>(0);
	this->_value = value;
	this->_identifier = "";
}

AST::AST(Kind kind, const char* identifier) {
	this->_kind = kind;
	this->_operands = std::vector<AST*>(0);
	this->_identifier = identifier;
	this->_value = 0;
}

AST::AST(Kind kind, std::vector<AST*> operands) {
	this->_kind = kind;
	this->_operands = operands;
	this->_identifier = "";
	this->_value = 0;
}

AST::AST(Kind kind, const signed long value, const std::string identifier) {
	this->_kind = kind;
	this->_operands = std::vector<AST*>(0);
	this->_identifier = identifier;
	this->_value = value;
}


AST::~AST() {
	for(AST* e : this->_operands) {
		delete e;
	}
}

Kind AST::kind() const {
	return this->_kind;
}

bool AST::is(signed long i)
{
	return (this->kind() == Kind::Integer && this->value() == i);
}

bool AST::isNot(signed long i)
{
	return this->kind() != Kind::Integer ||
				 (this->kind() == Kind::Integer && this->value() != i);
}

AST* AST::operand(unsigned long i) {
	if(
		this->kind() == Kind::Integer ||
		this->kind() == Kind::Symbol ||
		this->kind() == Kind::Infinity ||
		this->kind() == Kind::MinusInfinity
	) {
		return this;
	}
	
	if(this->kind() == Kind::FunctionCall) {
		i = i+1;
	}

	return this->_operands[i];
}

bool AST::includeOperand(AST* expr) {
	if(this->kind() == Kind::Set) {
		for(unsigned int i=0; i<this->numberOfOperands(); i++) {
			if(this->operand(i)->match(expr))
				return false;
		}
	}

	this->_operands.push_back(expr);

	return true;
}

bool AST::includeOperand(AST* expr,signed long i) {
	if(expr->kind() == Kind::Set) {
		for(unsigned int i=0; i<this->numberOfOperands(); i++) {
			if(this->operand(i)->match(expr))
				return false;
		}
	}

	std::vector<AST*>::iterator it = this->_operands.begin();
	std::advance(it, i);
	this->_operands.insert(it, expr);
	return true;
}

bool AST::removeOperand(AST* u) {
	for(unsigned int i=0; i<this->numberOfOperands(); i++) {
		if(this->operand(i)->match(u)) {
			this->_operands.erase(this->_operands.begin() + i);
			return true;
		}
	}
	return false;
}

bool AST::removeOperand(signed long i) {
	this->_operands.erase(this->_operands.begin() + i);
	return true;
}

bool AST::deleteOperand(signed long i) {
	AST* k = this->_operands[i];
	
	this->_operands.erase(this->_operands.begin() + i);

	delete k;

	return true;
}

unsigned AST::numberOfOperands() const {
	switch (this->kind()) {
	case Kind::Integer:
	case Kind::Fraction:
	case Kind::Symbol:
	case Kind::Infinity:
	case Kind::MinusInfinity:
		return 1;

	case Kind::FunctionCall:
		return this->_operands.size() - 1;

	default:
		return this->_operands.size();
	}
}

signed long AST::value() const {
	return this->_value;
}

const std::string AST::identifier() {
	return this->_identifier;
}

AST* AST::copy() {
	AST* u = new AST(this->kind(), this->value(), this->identifier());

	switch (this->kind()) {
	case Kind::Integer:
	case Kind::Symbol:
	case Kind::Infinity:
	case Kind::MinusInfinity:
		break;
	case Kind::FunctionCall:
		u->includeOperand(new AST(Kind::Symbol, this->funName().c_str()));
		for(unsigned int i=0; i<this->numberOfOperands(); i++)
			u->includeOperand(this->operand(i)->copy());
		break;
	case Kind::Fraction:
		for(unsigned int i=0; i<2; i++)
			u->includeOperand(this->operand(i)->copy());
		break;
	default:
		for(unsigned int i=0; i<this->numberOfOperands(); i++)
			u->includeOperand(this->operand(i)->copy());
		break;
	}

	return u;
}

bool AST::match(AST* const other) {
	if(this->kind() != other->kind())
	{
		return false;
	}
    
	if(this->numberOfOperands() != other->numberOfOperands())
	{
		return false;
	}

	if(this->kind() == Kind::FunctionCall) 
	{
		if(this->funName() != other->funName())
		{
			return false;
		}
		// TODO: match arguments
	}

	if(this->kind() == Kind::Fraction)
	{
		return this->operand(0)->match(other->operand(0)) &&
					 this->operand(1)->match(other->operand(1));
	}

	if(this->kind() == Kind::Integer)	
	{
		return this->value() == other->value();
	}
	
	if(this->kind() == Kind::Symbol)
	{
		return this->identifier() == other->identifier();
	}
	
	if(this->kind() == Kind::Undefined)
	{
		return this->value() == other->value();
	}
	
	if(this->kind() == Kind::Factorial)
	{
		return this->operand(0)->match(other->operand(0));
	}
	
	if(this->kind() == Kind::Division)
	{
		return this->operand(0)->match(other->operand(0)) &&
					 this->operand(1)->match(other->operand(1));
	}

	if(
		this->kind() == Kind::Infinity 			||
		this->kind() == Kind::MinusInfinity
	)
	{
	 return this->kind() == other->kind();
	}

	if(this->kind() == Kind::Subtraction) 
	{
		unsigned int matches = 0;

		if(!this->operand(0)->match(other->operand(0))) 
		{
			return false;
		}

		matches++;

		for(unsigned int i=1; i < this->numberOfOperands(); i++)
		{
			for(unsigned int j=1; j < other->numberOfOperands(); j++)
			{
				if(this->operand(i)->match(other->operand(j))) 
				{
					matches++;
					break;
				}
			}
		}

		return matches == this->numberOfOperands();    
	}


	if(
		this->kind() == Kind::Addition 				||
		this->kind() == Kind::Multiplication 	||
		this->kind() == Kind::Set
	) {

		unsigned int matches = 0;

		for(unsigned int i=0; i < this->numberOfOperands(); i++) 
		{
			for(unsigned int j=0; j < other->numberOfOperands(); j++) 
			{

				if(this->operand(i)->match(other->operand(j))) 
				{
					matches++;
					break;
				}
			}
		}
	
		return matches == this->numberOfOperands();    
	}

	// order of the operators does matter
	for(unsigned int i=0; i < this->numberOfOperands(); i++) 
	{
		if(!this->operand(i)->match(other->operand(i)))
		{
			return false;
		}
	}

	return true;
}

bool AST::isTerminal() {
	if(
		this->kind() == Kind::Integer 			||
		this->kind() == Kind::Fraction 			||
		this->kind() == Kind::Infinity 			||
		this->kind() == Kind::MinusInfinity ||
		this->kind() == Kind::Symbol 				||
		this->kind() == Kind::Tensor
	) return true;

	return false;
}


bool AST::freeOf(AST* const other) {
	if(this->match(other))
		return false;

	if(
		this->kind() == Kind::Integer ||
		this->kind() == Kind::Fraction ||
		this->kind() == Kind::Symbol
	) return true;

	for(unsigned int i=0; i<this->numberOfOperands(); i++) {
		if(!this->operand(i)->freeOf(other))
			return false;
	}

	return true;
}

bool AST::freeOfElementsInSet(AST* const set) {
	for(unsigned int i=0; i<set->numberOfOperands(); i++) {
		if(!this->freeOf(set->operand(i))) {
			return false;
		}
	}

	return true;
}

// bool AST::analogous(AST* const other) {
// 	if(this->kind() != other->kind())
// 		return false;
    
// 	if(this->numberOfOperands() != other->numberOfOperands())
// 		return false;

// 	// compare expressions that have meaningfull data
// 	if(this->kind() == Kind::Fraction)
// 		return this->operand(0)->analogous(other->operand(0)) &&
// 					 this->operand(1)->analogous(other->operand(1));

// 	if(
// 		this->kind() == Kind::Integer ||
// 		this->kind() == Kind::Symbol ||
// 		this->kind() == Kind::Undefined ||
// 		this->kind() == Kind::Infinity ||
// 		this->kind() == Kind::MinusInfinity
// 	)	return true;
	
// 	if(this->kind() == Kind::Factorial)
// 		return this->operand(0)->analogous(other->operand(0));
	
// 	if(this->kind() == Kind::Division)
// 		return this->operand(0)->analogous(other->operand(0)) &&
// 					 this->operand(1)->analogous(other->operand(1));

// 	if(this->kind() == Kind::Subtraction) {
// 		long matches = 0;
// 		long match = false;

// 		if(!this->operand(0)->analogous(other->operand(0))) {
// 			return false;
// 		}

// 		matches++;

// 		for(unsigned int i=1; i < this->numberOfOperands(); i++) {
// 			for(unsigned int j=1; j < other->numberOfOperands(); j++) {
// 				if(this->operand(i)->analogous(other->operand(j))) {
// 					matches++;
// 					match = true;
// 					break;
// 				}
// 			}

// 			if(match) {
// 				match = false;
// 				continue;
// 			}
// 		}

// 		return matches == this->numberOfOperands();    
// 	}

// 	if(
// 		this->kind() == Kind::Addition ||
// 		this->kind() == Kind::Multiplication
// 	) {

// 		if(other->kind() != this->kind())
// 			return false;

// 		long matches = 0;
// 		long match = false;

// 		for(unsigned int i=0; i < this->numberOfOperands(); i++) {
// 			for(unsigned int j=0; j < other->numberOfOperands(); j++) {
// 				if(this->operand(i)->analogous(other->operand(j))) {
// 					matches++;
// 					match = true;
// 					break;
// 				}
// 			}

// 			if(match) {
// 				match = false;
// 				continue;
// 			}
// 		}
// 			return matches == this->numberOfOperands();    
// 	}


// 	// order of the operators does matter
// 	for(unsigned int i=0; i < this->numberOfOperands(); i++) 
// 		if(!this->operand(i)->analogous(other->operand(i)))
// 			return false;

// 	return true;
// }

std::string AST::toString() {
	std::string res = "";
	
	switch(this->kind()) {
		case Kind::Fail:
			res += "Fail";
			break;

		case Kind::Addition:
			// res += "(";
			for(unsigned int i=0; i<this->numberOfOperands(); i++) {
				res += this->operand(i)->toString();
				if(i != this->numberOfOperands() -1)
					res += " + ";
			}
			// res += ")";
			break;
		
		case Kind::Subtraction:
			// res += "(";
			for(unsigned int i=0; i<this->numberOfOperands(); i++) {
					res += this->operand(i)->toString();
					if(i != this->numberOfOperands() -1)
						res += " - ";
			}
			// res += ")";
			break;
		
		case Kind::Power:
			// res += "(";
			if(this->operand(0)->numberOfOperands() > 1) {
				res += "(";
			}
			res += this->operand(0)->toString();
			if(this->operand(0)->numberOfOperands() > 1) {
				res += ")";
			}
			if(this->operand(1)->numberOfOperands() > 1 || this->operand(1)->kind() == Kind::Fraction) {
				res += "^(";
			} else {
				res += "^";
			}
			res += this->operand(1)->toString();
			if(this->operand(1)->numberOfOperands() > 1 || this->operand(1)->kind() == Kind::Fraction) {
				res += ")";
			}
			// res += ")";
			break;
		
		case Kind::Multiplication:
			// res += "(";
			for(unsigned int i=0; i<this->numberOfOperands(); i++) {
				// if(i == 0 && this->operand(i)->kind() == Kind::Integer && this->operand(i)->value() == -1) {
				// 	res += "-";
				// 	continue;
				// } 
				if(
					this->operand(i)->kind() == Kind::Addition || 
					this->operand(i)->kind() == Kind::Subtraction ||
					this->operand(i)->kind() == Kind::Power ||
					this->operand(i)->kind() == Kind::Division ||
					this->operand(i)->kind() == Kind::Fraction
				) res += "(";

				res += this->operand(i)->toString();

				if(
					this->operand(i)->kind() == Kind::Addition || 
					this->operand(i)->kind() == Kind::Subtraction ||
					this->operand(i)->kind() == Kind::Power ||
					this->operand(i)->kind() == Kind::Division ||
					this->operand(i)->kind() == Kind::Fraction
				) res += ")";
				
				if(
					i != this->numberOfOperands() -1
				) res += "*";
			}
			// res += ")";
			break;
		
		case Kind::Division:
			res += "(";
			res += this->operand(0)->toString();
			res += ")";
			res += "/";
			res += "(";
			res += this->operand(1)->toString();
			res += ")";
			// res += ")";
			break;
		case Kind::Fraction:
			// res += "(";
			res += this->operand(0)->toString();
			res += "/";
			res += this->operand(1)->toString();
			// res += ")";
			break;
		
		case Kind::Factorial:
			res += "!";
			// res += "(";
			res += this->operand(0)->toString();
			// res += ")";
			break;
			
		case Kind::Integer:
			res += std::to_string(this->value());
			break;
		case Kind::Infinity:
			res += "∞";
			break;
		case Kind::MinusInfinity:
			res += "-∞";
			break;
		case Kind::Symbol:
			res += this->identifier();
			break;

		case Kind::FunctionCall:
			res += this->funName();
			res += "(";
			for(unsigned int i=0; i<this->numberOfOperands(); i++) {
				res += this->operand(i)->toString();
				if(i != this->numberOfOperands() - 1)
					res += ", ";	
			}
			res += ")";
			break;
		case Kind::List:
			res += "[";
			for(unsigned int i=0; i<this->numberOfOperands(); i++) {
				res += this->operand(i)->toString();
				if(i!= this->numberOfOperands() - 1)
					res += ", ";
			}
			res += "]";
			break;
		case Kind::Set:
			res += "{";
			for(unsigned int i=0; i<this->numberOfOperands(); i++) {
				res += this->operand(i)->toString();
				if(i!= this->numberOfOperands() - 1)
					res += ", ";
			}
			res += "}";
			break;
		case Kind::Undefined:
			res += "Undefined";
			break;

		case Kind::Derivative:
			res += "diff(";
			res += this->operand(0)->toString();
			res += ", ";
			res += this->operand(1)->toString();
			res += ")";
			break;

		case Kind::Integral:
			res += "inte(";
			res += this->operand(0)->toString();
			res += ", ";
			res += this->operand(1)->toString();
			res += ")";
			break;


		case Kind::Matrix:
			res += "[";
			for(unsigned int i = 0; i < this->numberOfOperands(); i++)
			{
				res += "[";
					for(unsigned int j = 0; j < this->numberOfOperands(); j++)
					{
						res += this->operand(i)->operand(j)->toString();
						if(j < this->operand(i)->numberOfOperands() - 1)
						{
							res += ",";
						}
					}
				res += "]";
			
				if(i < this->numberOfOperands() - 1)
				{
					res += ",";
				}
			}
			res += "]";
			break;

		default:
		  res += "Not implemented(" + std::to_string(this->kind()) + ")" ;
			break;
	}

	return res;
}



AST* AST::operandList() {
	AST* L = new AST(Kind::List);

	for(unsigned int i=0;i<this->numberOfOperands(); i++) {
		L->includeOperand(this->operand(i)->copy());
	}

	return L;
}


void destroyASTs(std::vector<AST*> l) {
	for(AST* a : l)
		delete a;
}

const std::string AST::funName() {
	std::vector<AST*>::iterator it = this->_operands.begin();
	return (*it)->identifier();
}


AST* mapUnaryAST(AST* u, AST*(*f)(AST*)) {
	if(
			u->kind() == Kind::Integer ||
			u->kind() == Kind::Fraction ||
			u->kind() == Kind::Symbol ||
			u->kind() == Kind::Infinity ||
			u->kind() == Kind::MinusInfinity
	)
	{
		return f(u);
	}

	if(u->numberOfOperands() == 0) {
		return f(u);
	}

	AST* t = new AST(u->kind());

	if(u->kind() == Kind::FunctionCall) 
	{
		t->includeOperand(new AST(Kind::Symbol, u->funName().c_str()));
	}

	for(unsigned int i=0; i < u->numberOfOperands(); i++) {
		t->includeOperand(f(u->operand(i)));
	}

	return t;
}

// AST* AST::symbols()
// {
// 	if(this->kind() == Kind::Symbol)
// 	{
// 		return new AST(Kind::Set, { this->copy() });
// 	}

// 	AST* syms = new AST(Kind::Set);

// 	if(
// 		this->kind() == Kind::Addition       ||
// 		this->kind() == Kind::Subtraction    ||
// 		this->kind() == Kind::Power          ||
// 		this->kind() == Kind::Division		   ||
// 		this->kind() == Kind::Multiplication ||
// 		this->kind() == Kind::Matrix			   ||
// 		this->kind() == Kind::Set			   		 ||
// 		this->kind() == Kind::List
// 	)
// 	{
// 		for(unsigned int i = 0; i < this->numberOfOperands(); i++)
// 		{
// 			AST* s = this->operand(i)->symbols();

// 			if(s->numberOfOperands() > 0)
// 			{
// 				for(unsigned int k = 0; k < s->numberOfOperands(); k++)
// 				{
// 					syms->includeOperand(s->operand(k)->copy());
// 				}
// 			}

// 			delete s;
// 		}
// 	}

// 	if(
// 		this->kind() == Kind::Derivative ||
// 		this->kind() == Kind::Integral   ||
// 		this->kind() == Kind::Factorial
// 	)
// 	{
// 		AST* s = this->operand(0)->symbols();
// 		if(s->numberOfOperands() > 0)
// 		{
// 			for(unsigned int k = 0; k < s->numberOfOperands(); k++)
// 			{
// 				syms->includeOperand(s->operand(k)->copy());
// 			}
// 		}

// 		delete s;
// 	}

// 	if(
// 		this->kind() == Kind::Derivative ||
// 		this->kind() == Kind::Integral   ||
// 		this->kind() == Kind::Factorial
// 	)
// 	{
// 		AST* s = this->operand(0)->symbols();
// 		if(s->numberOfOperands() > 0)
// 		{
// 			for(unsigned int k = 0; k < s->numberOfOperands(); k++)
// 			{
// 				syms->includeOperand(s->operand(k)->copy());
// 			}
// 		}

// 		delete s;
// 	}


// 	return syms;
// }

AST* mapBinaryAST(AST* u, AST* v, AST*(*f)(AST*, AST*)) {
	if(
			u->kind() == Kind::Integer ||
			u->kind() == Kind::Symbol ||
			u->kind() == Kind::Infinity || 
			u->kind() == Kind::MinusInfinity
	) return f(u, v);

	if(u->numberOfOperands() == 0)
		return f(u, v);

	AST* t = new AST(u->kind());

	if(u->kind() == Kind::FunctionCall) {
		t->includeOperand(new AST(Kind::Symbol, u->funName().c_str()));
	}

	for(unsigned int i=0; i< u->numberOfOperands(); i++)
			t->includeOperand(f(u->operand(i), v));

	return t;
}

AST* deepReplace(AST* tree, AST* subtree, AST* v) {
	if(tree->kind() == subtree->kind()) {
		if(tree->match(subtree)) {
			return v->copy();		
		}
	}
	
	if(tree->numberOfOperands() > 1) {
		AST* t = new AST(tree->kind());

		for(unsigned int i=0; i<tree->numberOfOperands(); i++) {
			t->includeOperand(deepReplace(tree->operand(i), subtree, v));
		}

		return t;
	}

	return tree->copy();
}

AST* construct(Kind kind, AST* L) {
	AST* u = new AST(kind);

	for(unsigned int i=0; i<L->numberOfOperands(); i++) {
		u->includeOperand(L->operand(i)->copy());
	}

	return u;
}

bool AST::isOfForm(AST* const v, AST* const syms) {
	// assert syms is a list of symbols
	if(syms != nullptr)
	{
		assert(syms->kind() == Kind::List);

		for(unsigned int i=0; i < syms->numberOfOperands(); i++)
		{
			assert(syms->operand(i)->kind() == Kind::Symbol);
		}

	}


	AST* u = this;

	if(u->kind() != v->kind())
	{
		return false;
	}


	if(u->isTerminal() && v->isTerminal())
	{

		if(syms != nullptr && u->kind() == Kind::Symbol)
		{
			for(unsigned int i=0; i < syms->numberOfOperands(); i++)
			{
				if(u->identifier() == syms->operand(i)->identifier())
				{
					if(u->identifier() != v->identifier())
					{
						return false;
					}
				}
			}
		}
	
		return true;
	}

	if(u->numberOfOperands() != v->numberOfOperands())
	{
		return false;
	}

	if(u->kind() == Kind::FunctionCall)
	{
		if(u->funName() != v->funName())
		{
			return false;
		}

		if(u->numberOfOperands() != v->numberOfOperands())
		{
			return false;
		}

		for(unsigned int i = 0; i < u->numberOfOperands(); i++)
		{
			if(!u->operand(i)->isOfForm(v->operand(i), syms)) 
			{
				return false;
			}
		}
	}

	if(u->kind() == Kind::Factorial )
	{
		return u->operand(0)->isOfForm(v->operand(0), syms);
	}

	if(
		u->kind() == Kind::Power   ||
		u->kind() == Kind::Division||
		u->kind() == Kind::Derivative 
	)
	{
		return u->operand(0)->isOfForm(v->operand(0), syms) && 
	 		u->operand(1)->isOfForm(v->operand(1), syms);
	}

	unsigned int m = 0;
	unsigned int n = u->numberOfOperands();

	for(unsigned int i = 0; i < n; i++)
	{
		bool match = false;

		for(unsigned int j = 0; j < n; j++)
		{
			match = u->operand(i)->isOfForm(v->operand(j), syms);
		
			if(match) 
			{
				break;
			}
		}
	
		if(match)
		{
			m++;
		}
	}

	return n == m;
}

}
