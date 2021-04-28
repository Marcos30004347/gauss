#include "Polynomial.hpp"
#include "Core/Debug/Assert.hpp"

using namespace ast;
using namespace algebra;

namespace polynomial {

void includeVariable(std::vector<AST*>& vars, AST* u) {
	// TODO: optimize
	bool included = false;
	
	for(AST* k : vars) {
		if(k->match(u)){
			included = true;
			break;
		}
	}
	
	if(!included) {
		vars.push_back(u->deepCopy()); 
	}	
}

std::vector<AST*> variables(AST* u) {
	std::vector<AST*> vars = std::vector<AST*>(0);
	
	if(
		u->numberOfOperands() > 1 &&
		u->kind() != Kind::Power &&
		u->kind() != Kind::Division
	) {

		for(int i=0; i<u->numberOfOperands(); i++) {
			std::vector<AST*> vargs = variables(u->operand(i));
			
			for(AST* a : vargs)
				includeVariable(vars, a);
			
			for(AST* a : vargs)
				delete a;
		}

		return vars;
	}

	if(u->kind() == Kind::Power) {
		if(!isConstant(u->operand(0))) {
			includeVariable(vars, u->operand(0));
		}
	}

	if(u->kind() == Kind::Symbol) {
		includeVariable(vars, u);
	}	

	if(u->kind() == Kind::FunctionCall) {
		std::vector<AST*> var_args = std::vector<AST*>(0);

		for(int j=0; j<u->numberOfOperands(); j++) {
			std::vector<AST*> vargs = variables(u->operand(j));
			for(AST* a : vargs)
				var_args.push_back(a);				
		}

		if(var_args.size() > 0)
			includeVariable(vars, u);

		for(AST* a : var_args)
			delete a;
	}


	return vars;
}

bool isPolynomialGPE(ast::AST* u, std::vector<AST*> vars) {
	std::vector<AST*> vs = variables(u);
	
	bool inc = false;

	for(int j=0; j<vars.size(); j++) {
		for(int i=j; i<vs.size(); i++) {
			if(vars[j]->match(vs[i])) {
				inc = true;
				break;
			}
		}
		if(!inc) {
			for(AST* k : vs)
				delete k;
			return false;
		}
	}	

	for(AST* k : vs)
		delete k;

	return true;
}


AST* degreeGPE(AST* u, AST* x) {

	if(!isPolynomialGPE(u, { x }))
		return mul({inte(-1), inf()}); 

	if(
		u->kind() == Kind::Power &&
		// expoent is integer
		u->operand(1)->kind() == Kind::Integer &&
		// base is equal to x
		u->operand(0)->match(x)
	) {
		return u->operand(1)->deepCopy();
	}

	if(
		u->kind() == Kind::Multiplication ||
		u->kind() == Kind::Addition ||
		u->kind() == Kind::Subtraction
	) {
		AST* best = mul({inte(-1), inf()});
		AST* minus_inf = mul({inte(-1), inf()});

		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* tmp = degreeGPE(u->operand(i), x);

			if(tmp->kind() != Kind::Integer) {
				delete tmp;
				continue;
			}

			if(best->match(minus_inf)) {
				delete best;
				best = tmp;
				continue;
			}

			if(tmp->value() > best->value()) {
				delete best;
				best = tmp;
				continue;
			}

			delete tmp;
		}	

		if(best->match(minus_inf)) {
			delete minus_inf;
			delete best;
			return mul({inte(-1), inf()});
		}

		delete minus_inf;
		return best;
	}

	if(u->match(x)) {
		return inte(1);
	}

	return mul({inte(-1), inf()});
}

AST* coefficientGPE(AST* u, AST* x) {
	assert(
		x->kind() == Kind::Power,
		"coefficientGPE: 'param(x)=%s' needs to be a power!",
		x->toString().c_str()
	);

	assert(
		!isConstant(x->operand(0)),
		"coefficientGPE: base of 'param(x)=%s' "
		"cant be constant!",
		x->toString().c_str()
	);

	assert(
		x->operand(1)->kind() == Kind::Integer &&
		x->operand(1)->value() >= 0,
		"coefficientGPE: expoent of 'param(x)=%s' "
		"needs to be a non negative integer! ",
		x->toString().c_str()
	);

	if(u->match(x)) {
		return inte(1);
	}

	if(u->kind() == Kind::Multiplication) {
		bool found = false;

		signed long count = 0;

		AST* res = new AST(Kind::Multiplication);
		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* tmp = coefficientGPE(u->operand(i), x);

			if(tmp->kind() == Kind::Integer && tmp->value() == 1) {
				count++;

				assert(
					count < 2,
					"coefficientGPE: 'arg(u)=%s' have more than one 'arg(x)=%s'! "
					"Simplify 'arg(u)=%s' to solve the error!",
					u->toString().c_str(),
					x->toString().c_str(),
					u->toString().c_str()
				);

				delete tmp;
				found = true;
				continue;
			}

			delete tmp;
			res->includeOperand(u->operand(i)->deepCopy());
		}	

		if(!found) {
			delete res;
			return inte(0);
		}
	
		if(res->numberOfOperands() == 1) {
			AST* r = res->operand(0);
			res->removeOperand(0L);
			delete res;
			return r;
		}
	
		return res;
	}

	if(u->kind() == Kind::Addition || u->kind() == Kind::Subtraction) {
		AST* res = new AST(u->kind());
		
		for(int i=0; i<u->numberOfOperands(); i++) {
			AST* coeff = coefficientGPE(u->operand(i), x);
			
			if(coeff->kind() == Kind::Integer && coeff->value() == 0) {
				delete coeff;
				continue;
			}
			
			res->includeOperand(coeff);
		}

		if(res->numberOfOperands() == 1) {
			AST* r = res->operand(0);
			res->removeOperand(0L);
			delete res;
			return r;
		}
		return res;
	}

	return inte(0);
}

AST* leadingCoefficientGPE(AST* u, AST* x) {
	assert(
		!isConstant(x),
		"leadingCoefficientGPE: 'param(x)=%s' "
		"cant be a constant expression",
		x->toString().c_str()
	);


	AST* po = pow(
		x->deepCopy(),
		degreeGPE(u, x)
	);

	AST* lc = coefficientGPE(u, po);

	delete po;

	return lc;
}



}
