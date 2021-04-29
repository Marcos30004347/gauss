#include "AST.hpp"

namespace ast {

AST::AST(Kind kind)
: _operands{}, _kind{kind}, _identifier{""}, _value{0} {}

AST::AST(Kind kind, signed long value)
: _operands{}, _kind{kind}, _identifier{""}, _value{value} {}

AST::AST(Kind kind, const char* identifier)
: _operands{}, _kind{kind}, _identifier{identifier}, _value{0} {}

AST::AST(Kind kind, std::list<AST*> operands)
: _operands{operands}, _kind{kind}, _identifier{""}, _value{0} {}

AST::AST(Kind kind, const signed long value, const std::string identifier)
: _operands{}, _kind{kind}, _identifier{identifier}, _value{value} {}


AST::~AST() {
	for(AST* e : this->_operands)
		delete e;
}

Kind AST::kind() {
	return this->_kind;
}

AST* AST::operand(signed long i) {
	if(
		this->kind() == Kind::Integer ||
		this->kind() == Kind::Symbol ||
		this->kind() == Kind::Infinity
	) return this;
	
	if(this->kind() == Kind::FunctionCall) {
		i = i+1;
	}
	
	std::list<AST*>::iterator it = this->_operands.begin();
	std::advance(it, i);
	return *it;
}

void AST::includeOperand(AST* expr) {
	this->_operands.insert(this->_operands.end(), expr);
}

void AST::includeOperand(AST* expr,signed long i) {
	std::list<AST*>::iterator it = this->_operands.begin();
	std::advance(it, i);
	this->_operands.insert(it, expr);
}

void AST::removeOperand(AST* u) {
	this->_operands.remove(u);
}

void AST::removeOperand(signed long i) {
	std::list<AST*>::iterator it = this->_operands.begin();
	std::advance(it, i);
	this->_operands.remove(*it);
}

unsigned AST::numberOfOperands() {
	switch (this->kind()) {
	case Kind::Integer:
	case Kind::Fraction:
	case Kind::Symbol:
	case Kind::Infinity:
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
    if(this->kind() == Kind::Integer)
        return this->value() == other->value();
    if(this->kind() == Kind::Symbol)
        return this->identifier() == other->identifier();
    if(this->kind() == Kind::Undefined)
        return this->value() == other->value();
    if(this->kind() == Kind::Infinity)
        return this->kind() == other->kind();
    // order of the operators dont matter
    if(
        this->kind() == Kind::Addition ||
        this->kind() == Kind::Subtraction ||
        this->kind() == Kind::Multiplication
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

bool AST::freeOf(AST* const other) {
	for(int i=0; i<this->numberOfOperands(); i++) {
		if(this->operand(i)->match(other))
			return false;
	}

	return true;
}


std::string AST::toString() {
	std::string res = "";

	switch(this->kind()) {
		case Kind::Addition:
			for(int i=0; i<this->numberOfOperands(); i++) {
				res += this->operand(i)->toString();
				if(i != this->numberOfOperands() -1)
					res += " + ";
			}
			break;
		
		case Kind::Subtraction:
			for(int i=0; i<this->numberOfOperands(); i++) {
					res += this->operand(i)->toString();
					if(i != this->numberOfOperands() -1)
						res += " - ";
			}
			break;
		
		case Kind::Power:
			res += this->operand(0)->toString();
			if(this->operand(1)->numberOfOperands() > 1) {
				res += "^(";
			} else {
				res += "^";
			}
			res += this->operand(1)->toString();
			if(this->operand(1)->numberOfOperands() > 1) {
				res += ")";
			}
			break;
		
		case Kind::Multiplication:
			for(int i=0; i<this->numberOfOperands(); i++) {
				if(
					this->operand(i)->kind() == Kind::Addition || 
					this->operand(i)->kind() == Kind::Subtraction 
				) res += "(";

				res += this->operand(i)->toString();

				if(
					this->operand(i)->kind() == Kind::Addition || 
					this->operand(i)->kind() == Kind::Subtraction 
				) res += ")";
				
				if(
					i != this->numberOfOperands() -1 &&
					this->operand(i+1)->kind() != this->operand(i)->kind()
				) res += "⋅";
			}
			break;
		
		case Kind::Division:
		case Kind::Fraction:
			res += this->operand(0)->toString();
			res += "/";
			res += this->operand(1)->toString();
			break;
		
		case Kind::Factorial:
			res += "!";
			res += this->operand(0)->toString();
			break;
			
		case Kind::Integer:
			res += std::to_string(this->value());
			break;
		case Kind::Infinity:
			res += "∞";
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
		case Kind::Undefined:
			res += "Undefined";
			break;

		default:
		  res += "Not implemented";
			break;
	}

	return res;
}

void destroyASTs(std::list<AST*> l) {
	for(AST* a : l)
		delete a;
}

const std::string AST::funName() {
	std::list<AST*>::iterator it = this->_operands.begin();
	return (*it)->identifier();
}

AST* mapUnaryAST(AST* u, AST*(*f)(AST*)) {
	if(
			u->kind() == Kind::Integer ||
			u->kind() == Kind::Symbol ||
			u->kind() == Kind::Infinity
	) return f(u);

	if(u->numberOfOperands() == 0)
		return f(u);

	AST* t = new AST(u->kind());

	for(int i=0; i< u->numberOfOperands(); i++) {
		// printf("map %s\n", u->operand(i)->toString().c_str());
			t->includeOperand(f(u->operand(i)));
	}

	return t;
}


AST* mapBinaryAST(AST* u, AST* v, AST*(*f)(AST*, AST*)) {
	if(
			u->kind() == Kind::Integer ||
			u->kind() == Kind::Symbol ||
			u->kind() == Kind::Infinity
	) return f(u, v);

	if(u->numberOfOperands() == 0)
		return f(u, v);

	AST* t = new AST(u->kind());

	for(int i=0; i< u->numberOfOperands(); i++)
			t->includeOperand(f(u->operand(i), v));

	return t;
}

	

}
