#include "Addition.hpp"
#include "Multiplication.hpp"
#include "Rationals.hpp"
#include "Power.hpp"
#include "Core/Expand/Expand.hpp"

#include <vector>

using namespace ast;
using namespace expand;
using namespace algebra;

namespace simplification {

std::vector<AST*> simplifyAdditionRec(std::vector<AST*> L);

std::vector<AST*> restAddition(std::vector<AST*> p, int from = 1) {
    std::vector<AST*> a;

    for(int i=from; i <= p.size() - 1; i++) 
        a.push_back(p[i]->deepCopy());

    return a;
}

std::vector<AST*> adjoinAdditions(AST* p, std::vector<AST*> q) {
	std::vector<AST*> tmp = std::vector<AST*>(0);
	
	if(q.size() == 0) {
		return {p->deepCopy()};
	}

	std::vector<AST*> u = simplifyAdditionRec({ p, q[0] });
	
	for(AST* k : u) {
		tmp.push_back(k);
	}

	for(int i=1; i<q.size(); i++) {
		tmp.push_back(q[i]->deepCopy());
	}

	return tmp;
}


std::vector<AST*> mergeAdditions(std::vector<AST*> p, std::vector<AST*> q) {
	// return a copy of q
	if(p.size() == 0) {
		std::vector<AST*> r;
		for(AST* k:q) r.push_back(k->deepCopy());
		return r;
	}

	// return a copy of p
	if(q.size() == 0) {
		std::vector<AST*> r;
		for(AST* k:p) r.push_back(k->deepCopy());
		return r;
	}
	
	std::vector<AST*> H = simplifyAdditionRec({ p[0], q[0] });
	
	if(H.size() == 0) {
		std::vector<AST*> a = restAddition(p);
		std::vector<AST*> b = restAddition(q);

		std::vector<AST*> R = mergeAdditions(a, b);
		
		for(AST* k: a) delete k;
		for(AST* k: b) delete k;
		
		return R;
	}
	
	if(H.size() == 1) {
	
		std::vector<AST*> a = restAddition(p);
		std::vector<AST*> b = restAddition(q);

		std::vector<AST*> R = mergeAdditions(a, b);
		
		for(AST* k: a) delete k;
		for(AST* k: b) delete k;
		
		std::vector<AST*> _R = R;

		R = adjoinAdditions(H[0], R);
	
		for(AST* k : _R)
			delete k;

		for(AST* k : H)
			delete k;

		return R;
	}

	if(H[0]->match(p[0])) {
		std::vector<AST*> restP = restAddition(p);
		std::vector<AST*> mer = mergeAdditions(restP, q);
		std::vector<AST*> res = adjoinAdditions(p[0], mer);
	
		for(AST* k : mer)
			delete k;
		
		for(AST* k : restP)
			delete k;

		for(AST* k : H)
			delete k;

		return res;
	}

	std::vector<AST*> restQ = restAddition(q);
	std::vector<AST*> mer = mergeAdditions(p, restQ);
	std::vector<AST*> res = adjoinAdditions(q[0], mer);

	for(AST* k : mer)
		delete k;
	
	for(AST* k : restQ)
		delete k;

	for(AST* k : H)
			delete k;
	return res;

}


AST* nonConstantCoefficient(AST* a) {
	if(a->kind() == Kind::Power) {
		if(!isConstant(a->operand(0)) || !isConstant(a->operand(1)))
			return a->deepCopy();
	}

	AST* res = new AST(Kind::Multiplication);

	for(int i=0; i<a->numberOfOperands(); i++) {
		if(!isConstant(a->operand(i)))
			res->includeOperand(a->operand(i)->deepCopy());
	}

	if(res->numberOfOperands() == 0) {
		delete res;
		return new AST(Kind::Undefined);
	}

	if(res->numberOfOperands() == 1) {
		AST* r = res->operand(0)->deepCopy();
		delete res;
		return r;
	}

	return res;
}

AST* constantCoefficient(AST* a) {
	if(a->kind() == Kind::Power) {
		if(!isConstant(a->operand(0)) || !isConstant(a->operand(1)))
			return inte(1);

		return a->deepCopy();
	}

	AST* res = new AST(Kind::Multiplication);
	
	for(int i=0; i<a->numberOfOperands(); i++)
		if(isConstant(a->operand(i)))
			res->includeOperand(a->operand(i)->deepCopy());

	if(res->numberOfOperands() == 0) {
		delete res;
		return inte(1);
	}

	if(res->numberOfOperands() == 1) {
		AST* r = res->operand(0)->deepCopy();
		// printf("ASDSADASDASDAS %s\n", r->toString().c_str());
		delete res;
		return r;
	}

	if(res->numberOfOperands() > 1) {
		AST* old = res;
		res = reduceRNEAST(res);
		delete old;
	}

	return res;
}

std::vector<AST*> simplifyAdditionRec(std::vector<AST*> L) {
	if(
		L.size() == 2 &&
		L[0]->kind() != Kind::Addition &&
		L[1]->kind() != Kind::Addition
	) {
		AST* u1 = L[0];
		AST* u2 = L[1];

		if(isConstant(u1) && isConstant(u2)) {
			AST* P_ = add({u1->deepCopy(), u2->deepCopy()});
			AST* P = reduceRNEAST(P_);

			delete P_;
		
			if(P->kind() == Kind::Integer && P->value() == 0) {
				delete P;
				return {};
			}
			
			return {P};
		}

		if(u2->kind() == Kind::Infinity) {
			if(u1->kind() == Kind::MinusInfinity)
				return {new AST(Kind::Undefined)};
			return {new AST(Kind::Infinity)};
		} 

		if(u2->kind() == Kind::MinusInfinity) {
			if(u1->kind() == Kind::Infinity)
				return {new AST(Kind::Undefined)};
			return {new AST(Kind::MinusInfinity)};
		} 

		if(u1->kind() == Kind::Integer && u1->value() == 0) {
			return {u2->deepCopy()};
		}
	
		if(u2->kind() == Kind::Integer && u2->value() == 0) {
			return {u1->deepCopy()};
		}

		AST* nc_u1 = nonConstantCoefficient(u1);
		AST* nc_u2 = nonConstantCoefficient(u2);


		if(nc_u1->match(nc_u2)) {

	
			AST* S_ = add({
				constantCoefficient(u1),
				constantCoefficient(u2)
			});
	
			AST* P_ = mul({reduceAdditionAST(S_), nonConstantCoefficient(u1)});
			AST* P = reduceMultiplicationAST(P_);
			
			delete S_;
			delete P_;
			delete nc_u1;
			delete nc_u2;

			if(P->kind() == Kind::Integer && P->value() == 0) {
				delete P;
				return {};
			}
			
			return { P };
		}
	
		delete nc_u1;
		delete nc_u2;

		if(orderRelation(u2, u1))
			return simplifyAdditionRec({u2, u1});

		std::vector<AST*> L_;

		for(AST* k : L)
			L_.push_back(k->deepCopy());
		
		return L_;
	}

	if(
		L.size() == 2 &&
		(
			L[0]->kind() == Kind::Addition ||
			L[1]->kind() == Kind::Addition
		)
	) {
		AST* u1 = L[0];
		AST* u2 = L[1];

		if(
			u1->kind() == Kind::Addition &&
			u2->kind() == Kind::Addition
		) {
			std::vector<AST*> U1;
			std::vector<AST*> U2;
			
			for(int i=0; i<u1->numberOfOperands(); i++)
				U1.push_back(u1->operand(i)->deepCopy());
	
			for(int i=0; i<u2->numberOfOperands(); i++)
				U2.push_back(u2->operand(i)->deepCopy());
			

			std::vector<AST*> L_ = mergeAdditions(U1, U2);
			
			for(AST* k : U1)
				delete k;

			for(AST* k : U2)
				delete k;

			return L_;
		}
	
		if(u1->kind() == Kind::Addition) {
			std::vector<AST*> U1;
			
			for(int i=0; i<u1->numberOfOperands(); i++)
				U1.push_back(u1->operand(i)->deepCopy());
	

			std::vector<AST*> L_ = mergeAdditions(U1, { u2 });
			
			for(AST* k : U1)
				delete k;

			return L_;
		}
	
		if(u2->kind() == Kind::Addition) {
			std::vector<AST*> U2;
			
			for(int i=0; i<u2->numberOfOperands(); i++)
				U2.push_back(u2->operand(i)->deepCopy());
	

			std::vector<AST*> L_ = mergeAdditions({u1}, U2);
			
			for(AST* k : U2)
				delete k;

			return L_;
		}
	}

	AST* u1 = L[0];

	std::vector<AST*> restL = restAddition(L);
	
	std::vector<AST*> w = simplifyAdditionRec(restL);
	
	for(AST* k : restL)
		delete k;

	if(u1->kind() == Kind::Addition) {
		std::vector<AST*> U1;
		
		for(int i=0; i<u1->numberOfOperands(); i++)
			U1.push_back(u1->deepCopy());


		std::vector<AST*> L_ = mergeAdditions(U1, w);

		for(AST* k : U1)
			delete k;
		for(AST* k : w)
			delete k;
		
		return L_;
	}

	std::vector<AST*> L_ = mergeAdditions({u1}, w);
	
	for(AST* k : w)
		delete k;
	
	return L_;
}

AST* reduceAdditionAST(AST* u) {
	// printf("u = %s\n", u->toString().c_str());

	if(u->kind() == Kind::Undefined)
		return new AST(Kind::Undefined);
	
	if(u->numberOfOperands() == 1)
		return u->operand(0)->deepCopy();

	std::vector<AST*> L;
	for(int i=0; i<u->numberOfOperands(); i++)
		L.push_back(u->operand(i)->deepCopy());

	std::vector<AST*> R = simplifyAdditionRec(L);
	
	for(AST* t : L)
		delete t;

	if(R.size() == 0)
			return inte(0);

	if(R.size() == 1)
			return R[0];

	AST* res = new AST(Kind::Addition);
	
	for(AST* t : R)
		res->includeOperand(t);
	
	return res;
}

}

