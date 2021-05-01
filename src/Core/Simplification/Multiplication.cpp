#include "Multiplication.hpp"
#include "Addition.hpp"
#include "Rationals.hpp"
#include "Power.hpp"
#include "Core/Expand/Expand.hpp"

#include <vector>

using namespace ast;
using namespace expand;
using namespace algebra;

namespace simplification {

std::vector<AST*> simplifyProductRec(std::vector<AST*> L);

std::vector<AST*> restMultiplication(std::vector<AST*> p, int from = 1) {
    std::vector<AST*> a;

    for(int i=from; i <= p.size() - 1; i++) 
        a.push_back(p[i]->deepCopy());

    return a;
}

std::vector<AST*> adjoinProducts(AST* p, std::vector<AST*> q) {
	std::vector<AST*> tmp = std::vector<AST*>(0);
	
	if(q.size() == 0) {
		return {p->deepCopy()};
	}

	std::vector<AST*> u = simplifyProductRec({ p, q[0] });
	
	for(AST* k : u) {
		tmp.push_back(k);
	}

	for(int i=1; i<q.size(); i++) {
		tmp.push_back(q[i]->deepCopy());
	}
	
	return tmp;
}

std::vector<AST*> mergeProducts(std::vector<AST*> p, std::vector<AST*> q) {
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
	
	std::vector<AST*> H = simplifyProductRec({ p[0], q[0] });
	
	if(H.size() == 0) {
		std::vector<AST*> a = restMultiplication(p);
		std::vector<AST*> b = restMultiplication(q);

		std::vector<AST*> R = mergeProducts(a, b);
		
		for(AST* k: a) delete k;
		for(AST* k: b) delete k;
		
		return R;
	}
	
	if(H.size() == 1) {
	
		std::vector<AST*> a = restMultiplication(p);
		std::vector<AST*> b = restMultiplication(q);

		std::vector<AST*> R = mergeProducts(a, b);
		
		for(AST* k: a) delete k;
		for(AST* k: b) delete k;
		
		std::vector<AST*> _R = R;

		R = adjoinProducts(H[0], R);
	
		for(AST* k : _R)
			delete k;

		for(AST* k : H)
			delete k;

		return R;
	}

	if(H[0]->match(p[0])) {
		std::vector<AST*> restP = restMultiplication(p);
		std::vector<AST*> mer = mergeProducts(restP, q);
		std::vector<AST*> res = adjoinProducts(p[0], mer);
	
		for(AST* k : mer)
			delete k;
		
		for(AST* k : restP)
			delete k;

		for(AST* k : H)
			delete k;

		return res;
	}

	std::vector<AST*> restQ = restMultiplication(q);
	std::vector<AST*> mer = mergeProducts(p, restQ);
	std::vector<AST*> res = adjoinProducts(q[0], mer);

	for(AST* k : mer)
		delete k;
	
	for(AST* k : restQ)
		delete k;

	for(AST* k : H)
			delete k;
	return res;

}


std::vector<AST*> simplifyProductRec(std::vector<AST*> L) {
	// printf("L[0]=%s \n", L[0]->toString().c_str());
	// printf("L[1]=%s \n", L[1]->toString().c_str());

	if(
		L.size() == 2 &&
		L[0]->kind() != Kind::Multiplication &&
		L[1]->kind() != Kind::Multiplication
	) {
	
		AST* u1 = L[0];
		AST* u2 = L[1];

		if(isConstant(u1) && isConstant(u2)) {
			AST* P_ = mul({u1->deepCopy(), u2->deepCopy()});
			AST* P = reduceRNEAST(P_);
			
			delete P_;
		
			if(P->kind() == Kind::Integer && P->value() == 1) {
				delete P;
				return {};
			}
			
			return { P };
		}

		if(u1->kind() == Kind::Infinity) {
			if(u2->kind() == Kind::Integer && u2->value() == 0) {
				return {new AST(Kind::Undefined)};
			}
			else if(u2->kind() == Kind::Integer && u2->value() == -1) {
				return {new AST(Kind::MinusInfinity)};
			} else {
				return {new AST(Kind::Infinity)};
			}
		} 

		if(u1->kind() == Kind::MinusInfinity) {
			if(u2->kind() == Kind::Integer && u2->value() == 0) {
				return {new AST(Kind::Undefined)};
			}
			else if(u2->kind() == Kind::Integer && u2->value() == -1) {
				return {new AST(Kind::Infinity)};
			} else {
				return {new AST(Kind::MinusInfinity)};
			}
		} 

		if(u2->kind() == Kind::Infinity) {
			if(u1->kind() == Kind::Integer && u1->value() == 0) {
				return {new AST(Kind::Undefined)};
			}
			else if(u1->kind() == Kind::Integer && u1->value() == -1) {
				return {new AST(Kind::MinusInfinity)};
			} else {
				return {new AST(Kind::Infinity)};
			}
		} 

		if(u2->kind() == Kind::MinusInfinity) {
			if(u1->kind() == Kind::Integer && u1->value() == 0) {
				return {new AST(Kind::Undefined)};
			}
			else if(u1->kind() == Kind::Integer && u1->value() == -1) {
				return {new AST(Kind::Infinity)};
			} else {
				return {new AST(Kind::MinusInfinity)};
			}
		} 

		if(u1->kind() == Kind::Integer && u1->value() == 1) {
			return {u2->deepCopy()};
		}

		if(u1->kind() == Kind::Integer && u1->value() == 0) {
			return {inte(0)};
		}

		if(u2->kind() == Kind::Integer && u2->value() == 1) {
			return {u1->deepCopy()};
		}

		if(u2->kind() == Kind::Integer && u2->value() == 0) {
			return {inte(0)};
		}

		AST* base_u1 = base(u1);
		AST* base_u2 = base(u2);
	
		if(base_u1->match(base_u2)) {

			AST* S_ = add({ exp(u1), exp(u2) });

			AST* P_ = pow(base(u1), reduceAdditionAST(S_));
			AST* P = reducePowerAST(P_);
			
			delete S_;
			delete P_;
			delete base_u1;
			delete base_u2;

			if(P->kind() == Kind::Integer && P->value() == 1) {
				delete P;
				return {};
			}
			
			return {P};
		}
	
		delete base_u1;
		delete base_u2;

		if(orderRelation(u2, u1))
			return simplifyProductRec({u2, u1});

		std::vector<AST*> L_;
		for(AST* k : L)
			L_.push_back(k->deepCopy());
		
		return L_;
	}

	if(
		L.size() == 2 &&
		(
			L[0]->kind() == Kind::Multiplication ||
			L[1]->kind() == Kind::Multiplication
		)
	) {
		AST* u1 = L[0];
		AST* u2 = L[1];

		if(
			u1->kind() == Kind::Multiplication &&
			u2->kind() == Kind::Multiplication
		) {
			std::vector<AST*> U1;
			std::vector<AST*> U2;
			
			for(int i=0; i<u1->numberOfOperands(); i++)
				U1.push_back(u1->operand(i)->deepCopy());
	
			for(int i=0; i<u2->numberOfOperands(); i++)
				U2.push_back(u2->operand(i)->deepCopy());
			

			std::vector<AST*> L_ = mergeProducts(U1, U2);
			
			for(AST* k : U1)
				delete k;

			for(AST* k : U2)
				delete k;

			return L_;
		}
	
		if(u1->kind() == Kind::Multiplication) {
			std::vector<AST*> U1;
			
			for(int i=0; i<u1->numberOfOperands(); i++)
				U1.push_back(u1->operand(i)->deepCopy());
	

			std::vector<AST*> L_ = mergeProducts(U1, { u2 });
			
			for(AST* k : U1)
				delete k;

			return L_;
		}
	
		if(u2->kind() == Kind::Multiplication) {
			std::vector<AST*> U2;
			
			for(int i=0; i<u2->numberOfOperands(); i++)
				U2.push_back(u2->operand(i)->deepCopy());
	

			std::vector<AST*> L_ = mergeProducts({u1}, U2);
			
			for(AST* k : U2)
				delete k;

			return L_;
		}
	}

	AST* u1 = L[0];

	std::vector<AST*> restL = restMultiplication(L);
	
	std::vector<AST*> w = simplifyProductRec(restL);
	
	for(AST* k : restL)
		delete k;

	if(u1->kind() == Kind::Multiplication) {
		std::vector<AST*> U1;
		
		for(int i=0; i<u1->numberOfOperands(); i++)
			U1.push_back(u1->deepCopy());


		std::vector<AST*> L_ = mergeProducts(U1, w);

		for(AST* k : U1)
			delete k;
		for(AST* k : w)
			delete k;
		
		return L_;
	}

	std::vector<AST*> L_ = mergeProducts({u1}, w);
	
	for(AST* k : w)
		delete k;
	
	return L_;
}

AST* reduceMultiplicationAST(AST* u) {
	// printf("mul %s\n", u->toString().c_str());
	if(u->kind() == Kind::Undefined)
		return new AST(Kind::Undefined);
	
	for(int i=0; i<u->numberOfOperands(); i++) {
		AST* o = u->operand(i);
		if(o->kind() == Kind::Integer && o->value() == 0)
			return inte(0);
	}

	for(int i=0; i<u->numberOfOperands(); i++) {
		AST* o = u->operand(i);
	}

	if(u->numberOfOperands() == 1)
		return u->operand(0)->deepCopy();

	std::vector<AST*> L;

	for(int i=0; i<u->numberOfOperands(); i++)
		L.push_back(u->operand(i)->deepCopy());

	std::vector<AST*> R = simplifyProductRec(L);
	
	for(AST* t : L)
		delete t;


	if(R.size() == 1)
			return R[0];

	if(R.size() == 0)
		return inte(1);
	
	AST* res = new AST(Kind::Multiplication);
	
	for(AST* t : R)
		res->includeOperand(t);
	
	return res;
}

}

