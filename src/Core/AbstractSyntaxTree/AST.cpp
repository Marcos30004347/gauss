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
	switch (this->kind())
	{
	case Kind::Integer:
	case Kind::Fraction:
	case Kind::Symbol:
	case Kind::Function:
		return 1;
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

	switch (this->kind())
	{
	case Kind::Integer:
	case Kind::Symbol:
	case Kind::Function:
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

void AST::print() {
	switch(this->kind()) {
		case Kind::Addition:
			for(int i=0; i<this->numberOfOperands(); i++) {
				this->operand(i)->print();
				if(i != this->numberOfOperands() -1)
					printf(" + ");
			}
			break;
		
		case Kind::Subtraction:
			for(int i=0; i<this->numberOfOperands(); i++) {
					this->operand(i)->print();
					if(i != this->numberOfOperands() -1)
							printf(" - ");
			}
			break;
		
		case Kind::Power:
			this->operand(0)->print();
			printf("^(");
			this->operand(1)->print();
			printf(")");
			break;
		
		case Kind::Multiplication:
			for(int i=0; i<this->numberOfOperands(); i++) {
				if(
					this->operand(i)->kind() == Kind::Addition || 
					this->operand(i)->kind() == Kind::Subtraction 
				) printf("(");

				this->operand(i)->print();

				if(
					this->operand(i)->kind() == Kind::Addition || 
					this->operand(i)->kind() == Kind::Subtraction 
				) printf(")");

				if(i != this->numberOfOperands() -1)
					printf("⋅");
			}
			break;
		
		case Kind::Division:
		case Kind::Fraction:
			printf("(");
			this->operand(0)->print();
			printf(" / ");
			this->operand(1)->print();
			printf(")");
			break;
		
		case Kind::Factorial:
			printf("!");
			this->operand(0)->print();
			break;
			
		case Kind::Integer:
			printf("%ld", this->value());
			break;
		
		case Kind::Symbol:
			printf("%s", this->identifier().c_str());
			break;

		default:
			printf("Undefined");
			break;
	}
}

void destroyASTs(std::list<AST*> l) {
	for(AST* a : l)
		delete a;
}

AST* mapBinaryAST(AST* u, AST* v, AST*(*f)(AST*, AST*)) {
	if(
			u->kind() == Kind::Integer ||
			u->kind() == Kind::Symbol ||
			u->kind() == Kind::Function
	) return f(u, v);

	if(u->numberOfOperands() == 0)
		return f(u, v);

	AST* t = new AST(u->kind());

	for(int i=0; i< u->numberOfOperands(); i++)
			t->includeOperand(f(u->operand(i), v));
}

	

}
