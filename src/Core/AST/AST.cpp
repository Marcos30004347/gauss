#include "AST.hpp"

namespace ast {

AST::AST(Kind kind)
: _operands{}, _kind{kind}, _identifier{""}, _value{0} {}

AST::AST(Kind kind, signed long value)
: _operands{}, _kind{kind}, _identifier{""}, _value{value} {}

AST::AST(Kind kind, const char* identifier)
: _operands{}, _kind{kind}, _identifier{identifier}, _value{0} {}

AST::AST(Kind kind, std::vector<AST*> operands)
: _operands{operands}, _kind{kind}, _identifier{""}, _value{0} {}

AST::AST(Kind kind, const signed long value, const std::string identifier)
: _operands{}, _kind{kind}, _identifier{identifier}, _value{value} {}


AST::~AST() {
	for(AST* e : this->_operands) {
		delete e;
	}
}

Kind AST::kind() {
	return this->_kind;
}

AST* AST::operand(signed long i) {
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
		for(int i=0; i<this->numberOfOperands(); i++) {
			if(this->operand(i)->match(expr))
				return false;
		}
	}

	this->_operands.push_back(expr);
	return true;
}

bool AST::includeOperand(AST* expr,signed long i) {
	if(expr->kind() == Kind::Set) {
		for(int i=0; i<this->numberOfOperands(); i++) {
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
	for(int i=0; i<this->numberOfOperands(); i++) {
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

unsigned AST::numberOfOperands() {
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

const signed long AST::value() {
	return this->_value;
}

const std::string AST::identifier() {
	return this->_identifier;
}

AST* AST::deepCopy() {
	AST* u = new AST(this->kind(), this->value(), this->identifier());

	switch (this->kind()) {
	case Kind::Integer:
	case Kind::Symbol:
	case Kind::Infinity:
	case Kind::MinusInfinity:
		break;
	case Kind::FunctionCall:
		u->includeOperand(new AST(Kind::Symbol, this->funName().c_str()));
		for(int i=0; i<this->numberOfOperands(); i++)
			u->includeOperand(this->operand(i)->deepCopy());
		break;
	case Kind::Fraction:
		for(int i=0; i<2; i++)
			u->includeOperand(this->operand(i)->deepCopy());
		break;
	default:
		for(int i=0; i<this->numberOfOperands(); i++)
			u->includeOperand(this->operand(i)->deepCopy());
		break;
	}

	return u;
}

bool AST::match(AST* const other) {
	if(this->kind() != other->kind())
        return false;
    
	if(this->numberOfOperands() != other->numberOfOperands())
			return false;

	// compare expressions that have meaningfull data
	if(this->kind() == Kind::FunctionCall) {
		if(this->funName() != other->funName())
			return false;
	}

	if(this->kind() == Kind::Fraction)
		return this->operand(0)->match(other->operand(0)) &&
					 this->operand(1)->match(other->operand(1));

	if(this->kind() == Kind::Integer)	
		return this->value() == other->value();
	
	if(this->kind() == Kind::Symbol)
		return this->identifier() == other->identifier();
	
	if(this->kind() == Kind::Undefined)
		return this->value() == other->value();
	
	if(this->kind() == Kind::Factorial)
		return this->operand(0)->match(other->operand(0));
	
	if(this->kind() == Kind::Division)
		return this->operand(0)->match(other->operand(0)) &&
					 this->operand(1)->match(other->operand(1));

	if(
		this->kind() == Kind::Infinity ||
		this->kind() == Kind::MinusInfinity
	) return this->kind() == other->kind();

	if(this->kind() == Kind::Subtraction) {

		long matches = 0;
		long match = false;

		if(!this->operand(0)->match(other->operand(0))) {
			return false;
		}

		matches++;

		for(int i=1; i < this->numberOfOperands(); i++) {
			for(int j=1; j < other->numberOfOperands(); j++) {
				if(this->operand(i)->match(other->operand(j))) {
					matches++;
					match = true;
					break;
				}
			}

			if(match) {
				match = false;
				continue;
			}
		}

		return matches == this->numberOfOperands();    
	}


	if(
		this->kind() == Kind::Addition ||
		this->kind() == Kind::Multiplication ||
		this->kind() == Kind::Set
	) {

		long matches = 0;
		long match = false;

		for(int i=0; i < this->numberOfOperands(); i++) {
			for(int j=0; j < other->numberOfOperands(); j++) {
				if(this->operand(i)->match(other->operand(j))) {
					matches++;
					match = true;
					break;
				}
			}

			if(match) {
				match = false;
				continue;
			}
		}

		return matches == this->numberOfOperands();    
	}


	// order of the operators does matter
	for(int i=0; i < this->numberOfOperands(); i++) 
		if(!this->operand(i)->match(other->operand(i)))
			return false;

	return true;
}

bool AST::isTerminal() {
	if(
		this->kind() == Kind::Integer ||
		this->kind() == Kind::Fraction ||
		this->kind() == Kind::Infinity ||
		this->kind() == Kind::MinusInfinity ||
		this->kind() == Kind::Symbol
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

	for(int i=0; i<this->numberOfOperands(); i++) {
		if(!this->operand(i)->freeOf(other))
			return false;
	}

	return true;
}

bool AST::freeOfElementsInSet(AST* const set) {
	for(int i=0; i<set->numberOfOperands(); i++) {
		if(!this->freeOf(set->operand(i))) {
			return false;
		}
	}

	return true;
}

bool AST::analogous(AST* const other) {
	if(this->kind() != other->kind())
		return false;
    
	if(this->numberOfOperands() != other->numberOfOperands())
		return false;

	// compare expressions that have meaningfull data
	if(this->kind() == Kind::Fraction)
		return this->operand(0)->analogous(other->operand(0)) &&
					 this->operand(1)->analogous(other->operand(1));

	if(
		this->kind() == Kind::Integer ||
		this->kind() == Kind::Symbol ||
		this->kind() == Kind::Undefined ||
		this->kind() == Kind::Infinity ||
		this->kind() == Kind::MinusInfinity
	)	return true;
	
	if(this->kind() == Kind::Factorial)
		return this->operand(0)->analogous(other->operand(0));
	
	if(this->kind() == Kind::Division)
		return this->operand(0)->analogous(other->operand(0)) &&
					 this->operand(1)->analogous(other->operand(1));

	if(this->kind() == Kind::Subtraction) {
		long matches = 0;
		long match = false;

		if(!this->operand(0)->analogous(other->operand(0))) {
			return false;
		}

		matches++;

		for(int i=1; i < this->numberOfOperands(); i++) {
			for(int j=1; j < other->numberOfOperands(); j++) {
				if(this->operand(i)->analogous(other->operand(j))) {
					matches++;
					match = true;
					break;
				}
			}

			if(match) {
				match = false;
				continue;
			}
		}

		return matches == this->numberOfOperands();    
	}

	if(
		this->kind() == Kind::Addition ||
		this->kind() == Kind::Multiplication
	) {

		if(other->kind() != this->kind())
			return false;

		long matches = 0;
		long match = false;

		for(int i=0; i < this->numberOfOperands(); i++) {
			for(int j=0; j < other->numberOfOperands(); j++) {
				if(this->operand(i)->analogous(other->operand(j))) {
					matches++;
					match = true;
					break;
				}
			}

			if(match) {
				match = false;
				continue;
			}
		}
			return matches == this->numberOfOperands();    
	}


	// order of the operators does matter
	for(int i=0; i < this->numberOfOperands(); i++) 
		if(!this->operand(i)->analogous(other->operand(i)))
			return false;

	return true;
}

std::string AST::toString() {
	std::string res = "";
	
	switch(this->kind()) {
		case Kind::Addition:
			// res += "(";
			for(int i=0; i<this->numberOfOperands(); i++) {
				res += this->operand(i)->toString();
				if(i != this->numberOfOperands() -1)
					res += " + ";
			}
			// res += ")";
			break;
		
		case Kind::Subtraction:
			// res += "(";
			for(int i=0; i<this->numberOfOperands(); i++) {
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
			for(int i=0; i<this->numberOfOperands(); i++) {
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
			for(int i=0; i<this->numberOfOperands(); i++) {
				res += this->operand(i)->toString();
				if(i != this->numberOfOperands() - 1)
					res += ", ";	
			}
			res += ")";
			break;
		case Kind::List:
			res += "[";
			for(int i=0; i<this->numberOfOperands(); i++) {
				res += this->operand(i)->toString();
				if(i!= this->numberOfOperands() - 1)
					res += ", ";
			}
			res += "]";
			break;
		case Kind::Set:
			res += "{";
			for(int i=0; i<this->numberOfOperands(); i++) {
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
			res += "derivative(";
			res += this->operand(0)->toString();
			res += ", ";
			res += this->operand(1)->toString();
			res += ")";

		default:
		  res += "Not implemented(" + std::to_string(this->kind()) + ")" ;
			break;
	}

	return res;
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
	){
		return f(u);
	}

	if(u->numberOfOperands() == 0) {
		return f(u);
	}

	AST* t = new AST(u->kind());

	if(u->kind() == Kind::FunctionCall) {
		t->includeOperand(new AST(Kind::Symbol, u->funName().c_str()));
	}

	for(int i=0; i< u->numberOfOperands(); i++) {
		t->includeOperand(f(u->operand(i)));
	}

	return t;
}


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

	for(int i=0; i< u->numberOfOperands(); i++)
			t->includeOperand(f(u->operand(i), v));

	return t;
}

AST* deepReplace(AST* tree, AST* subtree, AST* v) {
	if(tree->kind() == subtree->kind()) {
		if(tree->match(subtree)) {
			return v->deepCopy();		
		}
	}
	
	if(tree->numberOfOperands() > 1) {
		AST* t = new AST(tree->kind());

		for(int i=0; i<tree->numberOfOperands(); i++) {
			t->includeOperand(deepReplace(tree->operand(i), subtree, v));
		}

		return t;
	}

	return tree->deepCopy();
}


}
