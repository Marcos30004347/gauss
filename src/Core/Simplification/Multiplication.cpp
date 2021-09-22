#include "Multiplication.hpp"
#include "Addition.hpp"
#include "Rationals.hpp"
#include "Power.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Algebra/List.hpp"

#include <vector>

using namespace ast;
using namespace expand;
using namespace algebra;

namespace simplification {

AST* simplifyProductRec(AST* L);


AST* mergeProducts(AST* p, AST* q) {
	if(p->numberOfOperands() == 0) 
		return q->copy();

	if(q->numberOfOperands() == 0) 
		return p->copy();
	
	AST* L = list({ p->operand(0)->copy(), q->operand(0)->copy() });

	AST* H = simplifyProductRec(L);
	
	delete L;
	
	if(H->numberOfOperands() == 0) {
		AST* a = rest(p);
		AST* b = rest(q);

		AST* R = mergeProducts(a, b);
		
		delete H;
		delete a;
		delete b;
		
		return R;
	}
	
	if(H->numberOfOperands() == 1) {
	
		AST* a = rest(p);
		AST* b = rest(q);

		AST* R = mergeProducts(a, b);
		
		delete a;
		delete b;
		
		AST* R_ = R;

		R = adjoin(H->operand(0), R, simplifyProductRec);
	
		delete R_;
		delete H;

		return R;
	}

	if(H->operand(0)->match(p->operand(0))) {
		AST* restP 	= rest(p);
		AST* mer 		= mergeProducts(restP, q);
		AST* res 		= adjoin(p->operand(0), mer, simplifyProductRec);
	
		delete mer;
		delete restP;
		delete H;

		return res;
	}

	AST* restQ 	= rest(q);
	AST* mer 		= mergeProducts(p, restQ);
	AST* res 		= adjoin(q->operand(0), mer, simplifyProductRec);

	delete mer;
	delete restQ;
	delete H;

	return res;
}


AST* simplifyProductRec(AST* L) {

	if(
		L->numberOfOperands() == 2 &&
		L->operand(0)->kind() != Kind::Multiplication &&
		L->operand(1)->kind() != Kind::Multiplication
	) {
	
		AST* u1 = L->operand(0);
		AST* u2 = L->operand(1);


		if(isConstant(u1) && isConstant(u2)) {
			AST* P_ = mul({u1->copy(), u2->copy()});
			AST* P = reduceRNEAST(P_);
			
			delete P_;
		
			if(P->kind() == Kind::Integer && P->value() == 1) {
				delete P;
				return list({});
			}
			
			return list({ P });
		}

		if(u1->kind() == Kind::Infinity) {
			if(u2->kind() == Kind::Integer && u2->value() == 0) {
				return list({undefined()});
			}
			else if(u2->kind() == Kind::Integer && u2->value() == -1) {
				return list({new AST(Kind::MinusInfinity)});
			} else {
				return list({new AST(Kind::Infinity)});
			}
		} 

		if(u1->kind() == Kind::MinusInfinity) {
			if(u2->kind() == Kind::Integer && u2->value() == 0) {
				return list({undefined()});
			}
			else if(u2->kind() == Kind::Integer && u2->value() == -1) {
				return list({new AST(Kind::Infinity)});
			} else {
				return list({new AST(Kind::MinusInfinity)});
			}
		} 

		if(u2->kind() == Kind::Infinity) {
			if(u1->kind() == Kind::Integer && u1->value() == 0) {
				return list({undefined()});
			}
			else if(u1->kind() == Kind::Integer && u1->value() == -1) {
				return list({new AST(Kind::MinusInfinity)});
			} else {
				return list({new AST(Kind::Infinity)});
			}
		} 

		if(u2->kind() == Kind::MinusInfinity) {
			if(u1->kind() == Kind::Integer && u1->value() == 0) {
				return list({undefined()});
			}
			else if(u1->kind() == Kind::Integer && u1->value() == -1) {
				return list({new AST(Kind::Infinity)});
			} else {
				return list({new AST(Kind::MinusInfinity)});
			}
		} 

		if(u1->kind() == Kind::Integer && u1->value() == 1) {
			return list({u2->copy()});
		}

		if(u1->kind() == Kind::Integer && u1->value() == 0) {
			return list({integer(0)});
		}

		if(u2->kind() == Kind::Integer && u2->value() == 1) {
			return list({u1->copy()});
		}

		if(u2->kind() == Kind::Integer && u2->value() == 0) {
			return list({integer(0)});
		}

		AST* base_u1 = base(u1);
		AST* base_u2 = base(u2);
		
		if(base_u1->match(base_u2)) {
			AST* S_ = add({ expoent(u1), expoent(u2) });

			AST* P_ = power(base(u1), reduceAdditionAST(S_));
			AST* P = reducePowerAST(P_);
			
			delete S_;
			delete P_;
			delete base_u1;
			delete base_u2;

			if(P->kind() == Kind::Integer && P->value() == 1) {
				delete P;
				return list({});
			}
			
			return list({P});
		}
	
		delete base_u1;
		delete base_u2;

		if(orderRelation(u2, u1)) {
			return list({u2->copy(), u1->copy()});
			// AST* L_ = list({ u2->copy(), u1->copy() });
			// AST* R = simplifyProductRec(L_);

			// delete L_;

			// return R;
		}

		return list({u1->copy(), u2->copy()});
	}

	if(
		L->numberOfOperands() == 2 &&
		(
			L->operand(0)->kind() == Kind::Multiplication ||
			L->operand(1)->kind() == Kind::Multiplication
		)
	) {
		AST* u1 = L->operand(0);
		AST* u2 = L->operand(1);

		if(
			u1->kind() == Kind::Multiplication &&
			u2->kind() == Kind::Multiplication
		) {
			AST* U1 = new AST(Kind::List);
			AST* U2 = new AST(Kind::List);
			
			for(unsigned int i=0; i<u1->numberOfOperands(); i++)
				U1->includeOperand(u1->operand(i)->copy());
	
			for(unsigned int i=0; i<u2->numberOfOperands(); i++)
				U2->includeOperand(u2->operand(i)->copy());
			

			AST* L_ = mergeProducts(U1, U2);
			
			delete U1;
			delete U2;

			return L_;
		}
	
		if(u1->kind() == Kind::Multiplication) {
			AST* U1 = new AST(Kind::List);
			AST* U2 = new AST(Kind::List);
			
			for(unsigned int i=0; i<u1->numberOfOperands(); i++)
				U1->includeOperand(u1->operand(i)->copy());

			U2->includeOperand(u2->copy());
	
			AST* L_ = mergeProducts(U1, U2);
			
			delete U1;
			delete U2;
			return L_;
		}
	
		if(u2->kind() == Kind::Multiplication) {
			AST* U1 = new AST(Kind::List);
			AST* U2 = new AST(Kind::List);
			
			for(unsigned int i=0; i<u2->numberOfOperands(); i++)
				U2->includeOperand(u2->operand(i)->copy());

			U1->includeOperand(u1->copy());

			AST* L_ = mergeProducts(U1, U2);
			
			delete U1;
			delete U2;

			return L_;
		}
	}

	AST* u1 = L->operand(0);

	AST* restL = rest(L);
	
	AST* w = simplifyProductRec(restL);
	
	delete restL;

	if(u1->kind() == Kind::Multiplication) {
		AST* U1 = new AST(Kind::List);
		
		for(unsigned int i=0; i<u1->numberOfOperands(); i++)
			U1->includeOperand(u1->operand(i)->copy());


		AST* L_ = mergeProducts(U1, w);

		delete U1;
		delete w;

		return L_;
	}

	AST* U1 = new AST(Kind::List);
	U1->includeOperand(u1->copy());

	AST* L_ = mergeProducts(U1, w);
	
	delete w;
	delete U1;
	
	return L_;
}

AST* reduceMultiplicationAST(AST* u) {

	if(u->kind() == Kind::Undefined)
		return undefined();
	
	for(unsigned int i=0; i<u->numberOfOperands(); i++) {
		AST* o = u->operand(i);
		if(o->kind() == Kind::Integer && o->value() == 0)
			return integer(0);
	}

	if(u->numberOfOperands() == 1)
		return u->operand(0)->copy();

	AST* L = new AST(Kind::List);
	
	for(unsigned int i=0; i<u->numberOfOperands(); i++)
		L->includeOperand(u->operand(i)->copy());

	AST* R = simplifyProductRec(L);
	
	delete L;

	if(R->numberOfOperands() == 1) {
		AST* r = R->operand(0)->copy();
		delete R;
		return r;
	}

	if(R->numberOfOperands() == 0) {
		delete R;
		return integer(1);
	}
	
	AST* res = new AST(Kind::Multiplication);
	
	for(unsigned int i=0; i<R->numberOfOperands(); i++) {
		res->includeOperand(R->operand(i)->copy());
	}
	
	delete R;

	return res;
}

}

