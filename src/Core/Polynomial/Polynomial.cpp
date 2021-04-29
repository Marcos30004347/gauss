#include "Polynomial.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Reduce/Reduce.hpp"
#include "Core/Expand/Expand.hpp"

using namespace ast;
using namespace expand;
using namespace reduce;
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
	if(u->kind() == Kind::Integer)
		return true;

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

	if(x->operand(1)->value() == 1) {
		if(u->match(x->operand(0))) {
			return inte(1);
		}
	}

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

		if(res->numberOfOperands() == 0) {
			delete res;
			return inte(1);
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
	
		if(res->numberOfOperands() == 0) {
			delete res;
			return inte(0);
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

std::pair<ast::AST*, ast::AST*> divideGPE(AST* u, AST* v, AST* x) {
	assert(
		isPolynomialGPE(u, {x}),
		"'param(u)=%s' needs to be a "
		"GPE(General Polynomial Expression)! "
		"in 'param(x)=%s'", 
		u->toString().c_str(),
		x->toString().c_str()
	);

	assert(
		isPolynomialGPE(v, {x}),
		"'param(v)=%s' needs to be a "
		"GPE(General Polynomial Expression)! "
		"in 'param(x)=%s'", 
		v->toString().c_str(),
		x->toString().c_str()
	);

	AST* q = inte(0);
	AST* r = u->deepCopy();

	AST* m = degreeGPE(r, x);
	AST* n = degreeGPE(v, x);

	AST* lcv = leadingCoefficientGPE(v, x);

	while(m->value() >= n->value()) {
		AST* lcr = leadingCoefficientGPE(r, x);

		AST* q_ = add({
			q->deepCopy(),
			mul({
				div(lcr->deepCopy(), lcv->deepCopy()),
				pow(
					x->deepCopy(),
					sub({
						m->deepCopy(),
						n->deepCopy()
					})
				)
			})
		});

		delete q;
		q = expandAST(q_);

		delete q_;

		AST* r_ = sub({
			sub({
				r->deepCopy(),
				mul({
					lcr->deepCopy(),
					pow(x->deepCopy(), m->deepCopy())
				})
			}),
			mul({
				sub({
					v->deepCopy(),
					mul({
						lcv->deepCopy(),
						pow(x->deepCopy(), n->deepCopy())
					}),
				}),
				div(lcr->deepCopy(), lcv->deepCopy()),
				pow(
					x->deepCopy(),
					sub({m->deepCopy(), n->deepCopy()})
				)
			})
		});

		delete r;

		r = expandAST(r_);

		delete r_;
		delete m;
		delete lcr;

		m = degreeGPE(r, x);
	}

	std::pair<AST*, AST*> res = std::make_pair(expandAST(q), expandAST(r));
	
	delete r;
	delete q;
	delete m;
	delete n;
	delete lcv;
	
	return res;
}

AST* quotientGPE(AST* u, AST* v, AST* x) {
	std::pair<ast::AST*, ast::AST*> res = divideGPE(u,v,x);
	delete res.second;
	return res.first;
}

AST* remainderGPE(AST* u, AST* v, AST* x) {
	std::pair<ast::AST*, ast::AST*> res = divideGPE(u,v,x);
	delete res.first;
	return res.second;
}

AST* expandGPE(AST* u, AST* v, AST* x, AST* t) {
	if(u->kind() == Kind::Integer && u->value() == 0)
		return inte(0);

	std::pair<AST*, AST*> d = divideGPE(u, v, x);
	
	AST* q = d.first;
	AST* r = d.second;

	AST* exp = add({
		mul({
			t->deepCopy(),
			expandGPE(q, v, x, t)
		}),
		r->deepCopy()
	});

	AST* res = expandAST(exp);
	
	delete exp;
	delete q;
	delete r;
	
	return res;
}	

}
