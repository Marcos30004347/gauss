#include <assert.h>
#include <algorithm>
#include <string.h>
#include <cstdio>
#include "Algebra.hpp"

#include "Core/Reduce/Rationals.hpp"

using namespace ast;
using namespace reduce;

namespace algebra {

AST* inte(signed long val) {
	return new AST(Kind::Integer, val);
}

AST* symb(const char* identifier) {
	return new AST(Kind::Symbol, identifier);
}

AST* frac(signed long n, signed long d) {
	return new AST(Kind::Fraction, {
		new AST(Kind::Integer, n),
		new AST(Kind::Integer, d)
	});
}

ast::AST* frac(AST* n, AST* d) {
	assert(isConstant(n));
	assert(isConstant(d));
	return new AST(Kind::Fraction, { n, d });
}

AST* add(std::list<AST*> terms) {
	return new AST(Kind::Addition, terms);
}

AST* sub(std::list<AST*> terms) {
	return new AST(Kind::Subtraction, terms);
}

AST* mul(std::list<AST*> terms) {
	return new AST(Kind::Multiplication, terms);
}

AST* div(AST* num, AST* den) {
	return new AST(Kind::Division, { num, den });
}

AST* pow(ast::AST* bas, ast::AST* exp) {
	return new AST(Kind::Power, { bas, exp });
}

AST* fact(AST* u) {
	return new AST(Kind::Factorial, {
		u,
	});
}

bool isConstant(AST* u) {
	if(
		u->kind() == Kind::Symbol ||
		u->kind() == Kind::Undefined
	) return false;
	
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction
	) return true;
	
	for(int i=0; i<u->numberOfOperands(); i++) {
		if(isConstant(u->operand(i)))
			return false;
	}

	return false; 
}

int gcd_rec(int a, int b) {
	if (b == 0)
		return a;

	return gcd_rec(b, a % b);
}

AST* base(AST* u) {
	if(u->kind() == Kind::Power) 
		return u->operand(0)->deepCopy();
	
	return u->deepCopy();
}

AST* exp(AST* u) {
	if(u->kind() == Kind::Power) 
		return u->operand(1)->deepCopy();

	return inte(1);
}

AST* gcd(AST* a, AST* b) {
    return inte(gcd_rec(a->value(), b->value()));
}

AST* num(AST* u) {
    if(u->kind() != Kind::Fraction && u->kind() != Kind::Division)
        return u->deepCopy();
    
    return u->operand(0)->deepCopy();
}

AST* den(AST* u) {
    if(u->kind() != Kind::Fraction && u->kind() != Kind::Division)
        return inte(1);
    
    return u->operand(1)->deepCopy();
}

bool isRNE(AST* u) {
	if(u->kind() == Kind::Integer)
		return true;

	if(u->kind() == Kind::Fraction)
		return isConstant(u->operand(0)) && isConstant(u->operand(1));

	if(u->kind() == Kind::Addition && u->numberOfOperands() <= 2) {
		for(int i=0; i < u->numberOfOperands(); i++) {
	
			// printf("KKKKK\n");
			// u->operand(i)->print();
			// printf("\n");
			// printf("is rne: %i\n", isRNE(u->operand(i)));
		
			if(!isRNE(u->operand(i)))
				return false;
		}
		return true;
	}

	if(u->kind() == Kind::Subtraction && u->numberOfOperands() <= 2) {
		for(int i=0; i < u->numberOfOperands(); i++)
			if(!isRNE(u->operand(i)))
				return false;
		return true;
	}

	if(u->kind() == Kind::Multiplication && u->numberOfOperands() == 2) {
		for(int i=0; i < u->numberOfOperands(); i++)
			if(!isRNE(u->operand(i)))
				return false;
		return true;
	}

	if(u->kind() == Kind::Division) {
		for(int i=0; i < u->numberOfOperands(); i++)
			if(!isRNE(u->operand(i)))
				return false;
		return true;
	}

	if(u->kind() == Kind::Power) {
		AST* b = base(u);
		AST* e = exp(u);
	
		bool is_rne = isRNE(b) && e->kind() == Kind::Integer;

		destroyASTs({b, e});
	
		return is_rne;
	}

	return false;
}

bool compareSymbols(std::string a, std::string b) {
    return std::lexicographical_compare(a.c_str(), a.c_str() + a.length(), b.c_str(), b.c_str() + b.length());
}

bool compareConstants(AST* u, AST* v) {
	if(u->kind() == Kind::Integer && v->kind() == Kind::Integer)
			return u->value() < v->value();

	AST* d = gcd(u, v);
	AST* num_u = num(u);
	AST* num_v = num(v);

	if(
			d->kind() == Kind::Integer &&
			num_u->kind() == Kind::Integer &&
			num_v->kind() == Kind::Integer
	) {
		AST* o_e = mul({d, num_u});
		AST* o_f = mul({d, num_v});
	
		AST* e = reduceRNEAST(o_e);
		AST* f = reduceRNEAST(o_f);
		
		bool res = e->value() < f->value();
		
		destroyASTs({ d, num_u, num_v, o_e, o_f, e, f });
		
		return  res;
	}

	destroyASTs({ d, num_u, num_v });

	return false;
}

bool compareProductsAndSummations(AST* u, AST* v) {
	long m = u->numberOfOperands() - 1;
	long n = v->numberOfOperands() - 1;

	for(int k=0; k <= std::min(m, n); k++) {
		if(!u->operand(m - k)->match(v->operand(n - k))) {
			return orderRelation(u->operand(m - k), v->operand(n - k));
		}
	}
	
	return m < n;
}

bool comparePowers(AST* u, AST* v) {
	AST* b_u = base(u);
	AST* b_v = base(v);

	if(!b_u->match(b_v)) {
		bool res = orderRelation(b_u, b_v);
		destroyASTs({b_u, b_v});
		return res; 
	}
	
	AST* e_u = exp(u);
	AST* e_v = exp(v);

	bool res = orderRelation(e_u, e_v);
	destroyASTs({e_u, e_v, b_u, b_v});
	return res;
}

bool compareFactorials(AST* u, AST* v) {
    return orderRelation(u->operand(0), v->operand(0));
}

bool compareFunctions(AST* u, AST* v) {
	if(u->operand(0)->identifier() == v->operand(0)->identifier())
		return orderRelation(u->operand(0), v->operand(0));

	AST* argsu = u->operand(1);
	AST* argsv = v->operand(1);

	if(argsu->numberOfOperands() >= 1 && argsv->numberOfOperands() >= 1) {
		long m = argsu->numberOfOperands() - 1;
		long n = argsv->numberOfOperands() - 1;

		for(int k=0; k <= std::min(m, n); k++) {
			if(!argsu->operand(m - k)->match(argsv->operand(n - k))) {
				bool res = orderRelation(argsu->operand(m - k), argsv->operand(n - k));
				destroyASTs({ argsu, argsv });
				return res;
			}
		}
		
		destroyASTs({ argsu, argsv });
		return m < n;
	}

	destroyASTs({ argsu, argsv });
	return true;
}

bool orderRelation(AST* u, AST* v) {
 if(isConstant(u) && isConstant(v))
	return compareConstants(u, v);

	if(u->kind() == Kind::Symbol && v->kind() == Kind::Symbol)
			return compareSymbols(u->identifier(), v->identifier());

	if(u->kind() == Kind::Addition && v->kind() == Kind::Addition)
			return compareProductsAndSummations(u, v);

	if(u->kind() == Kind::Multiplication && v->kind() == Kind::Multiplication)
			return compareProductsAndSummations(u, v);
	
	if(u->kind() == Kind::Power && v->kind() == Kind::Power)
			return comparePowers(u, v);

	if(u->kind() == Kind::Factorial && v->kind() == Kind::Factorial)
			return compareFactorials(u, v);

	if(u->kind() == Kind::FunctionCall && v->kind() == Kind::FunctionCall)
			return compareFunctions(u, v);

	if(isConstant(u))
			return true;

	// if(IsConstant(v))
	//     return false;

	if(
		u->kind() == Kind::Multiplication && (
		v->kind() == Kind::Power ||
		v->kind() == Kind::Addition ||
		v->kind() == Kind::Factorial ||
		v->kind() == Kind::FunctionCall ||
		v->kind() == Kind::Symbol
	)) {
		AST* m = mul({v->deepCopy()});
		bool res = orderRelation(u, m);
		destroyASTs({ m });
		return res;
	} 

	if(
		u->kind() == Kind::Power && (
		v->kind() == Kind::Addition ||
		v->kind() == Kind::Factorial ||
		v->kind() == Kind::FunctionCall ||
		v->kind() == Kind::Symbol
	)) {

		AST* m = pow(v->deepCopy(), inte(1));
		bool res = orderRelation(u, m);
		delete m;
		return res;
	} 

	if(
		u->kind() == Kind::Addition && (
		v->kind() == Kind::Factorial ||
		v->kind() == Kind::FunctionCall ||
		v->kind() == Kind::Symbol
	)) {
		AST* m = add({v->deepCopy()});
		bool res = orderRelation(u, m);
		destroyASTs({ m });
		return res;
	}

	if(
		u->kind() == Kind::Factorial && (
		v->kind() == Kind::FunctionCall ||
		v->kind() == Kind::Symbol
	)) {
		if(u->operand(0)->match(v)) {
			return false;
		} else {
			AST* m = fact(v->deepCopy());
			bool res = orderRelation(u, m);
			destroyASTs({m});
			return res;
		}
	}

	if(u->kind() == Kind::FunctionCall && v->kind() == Kind::Symbol) {
		if(u->operand(0)->identifier() == v->identifier()) {
			return false;
		} else {
			return orderRelation(u->operand(0), v);
		}
	}

	return !orderRelation(v, u);
}


AST* binomial(signed long n, std::vector<signed long> ks) {
	AST* p = new AST(Kind::Multiplication);
	for(signed long k : ks)
		p->includeOperand(fact(inte(k)));
	return div(fact(inte(n)), p);
}


AST* funCall(const char* id, std::list<AST*> args) {
	AST* f = new AST(Kind::FunctionCall);
	f->includeOperand(symb(id));
	for(AST* a : args)
		f->includeOperand(a);
	return f;
}

AST* inf() {
	return new AST(Kind::Infinity);
}

} // algebra
