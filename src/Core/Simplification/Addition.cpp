#include "Addition.hpp"
#include "Multiplication.hpp"
#include "Rationals.hpp"
#include "Power.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Algebra/List.hpp"

#include <vector>

using namespace ast;
using namespace expand;
using namespace algebra;

namespace simplification {

AST* simplifyAdditionRec(AST* L);

// std::vector<AST*> restAddition(std::vector<AST*> p, int from = 1) {
//     std::vector<AST*> a;

//     for(int i=from; i <= p.size() - 1; i++) 
//         a.push_back(p[i]->deepCopy());

//     return a;
// }

// std::vector<AST*> adjoinAdditions(AST* p, std::vector<AST*> q) {
// 	std::vector<AST*> tmp = std::vector<AST*>(0);
	
// 	if(q.size() == 0) {
// 		return {p->deepCopy()};
// 	}

// 	std::vector<AST*> u = simplifyAdditionRec({ p, q[0] });

// 	for(AST* k : u) {
// 		tmp.push_back(k);
// 	}

// 	for(int i=1; i<q.size(); i++) {
// 		tmp.push_back(q[i]->deepCopy());
// 	}

// 	return tmp;
// }


AST* mergeAdditions(AST* p, AST* q) {
	// return a copy of q

	if(p->numberOfOperands() == 0) {
		return q->deepCopy();
	}

	// return a copy of p
	if(q->numberOfOperands() == 0) {
		return p->deepCopy();
	}

	AST* L = list({ p->operand(0)->deepCopy(), q->operand(0)->deepCopy() });
	AST* H = simplifyAdditionRec(L);
	delete L;

	if(H->numberOfOperands() == 0) {
		AST* a = rest(p);
		AST* b = rest(q);
	
		AST* R = mergeAdditions(a, b);

		delete a;
		delete b;
		delete H;
	
		return R;
	}
	
	if(H->numberOfOperands() == 1) {
	
		AST* a = rest(p);
		AST* b = rest(q);

		AST* R = mergeAdditions(a, b);
	
		delete a;
		delete b;
	
		AST* R_ = R;

		R = adjoin(H->operand(0), R, simplifyAdditionRec);

		delete R_;
		delete H;

		return R;
	}

	if(H->operand(0)->match(p->operand(0))) {
		AST* restP = rest(p);
		AST* mer 	= mergeAdditions(restP, q);
		AST* res 	= adjoin(p->operand(0), mer, simplifyAdditionRec);
	
		delete mer;
		delete restP;
		delete H;

		return res;
	}

	AST* restQ = rest(q);
	AST* mer 	= mergeAdditions(p, restQ);
	AST* res 	= adjoin(q->operand(0), mer, simplifyAdditionRec);

	delete mer;
	delete restQ;
	delete H;

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
		return undefined();
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
			return integer(1);

		return a->deepCopy();
	}

	AST* res = new AST(Kind::Multiplication);
	
	for(int i=0; i<a->numberOfOperands(); i++)
		if(isConstant(a->operand(i)))
			res->includeOperand(a->operand(i)->deepCopy());

	if(res->numberOfOperands() == 0) {
		delete res;
		return integer(1);
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

AST* simplifyAdditionRec(AST* L) {
	if(
		L->numberOfOperands() == 2 &&
		L->operand(0)->kind() != Kind::Addition &&
		L->operand(1)->kind() != Kind::Addition
	) {
		AST* u1 = L->operand(0);
		AST* u2 = L->operand(1);

		if(isConstant(u1) && isConstant(u2)) {
			AST* P_ = add({u1->deepCopy(), u2->deepCopy()});
			AST* P = reduceRNEAST(P_);

			delete P_;
		
			if(P->kind() == Kind::Integer && P->value() == 0) {
				delete P;
				return list({});
			}
			
			return list({P});
		}

		if(u2->kind() == Kind::Infinity) {
			if(u1->kind() == Kind::MinusInfinity)
				return {undefined()};
			return list({new AST(Kind::Infinity)});
		} 

		if(u2->kind() == Kind::MinusInfinity) {
			if(u1->kind() == Kind::Infinity)
				return list({undefined()});
			return list({new AST(Kind::MinusInfinity)});
		} 

		if(u1->kind() == Kind::Integer && u1->value() == 0) {
			return list({u2->deepCopy()});
		}
	
		if(u2->kind() == Kind::Integer && u2->value() == 0) {
			return list({u1->deepCopy()});
		}

		AST* nc_u1 = nonConstantCoefficient(u1);
		AST* nc_u2 = nonConstantCoefficient(u2);
		
		if(nc_u1->match(nc_u2)) {
			AST* S_ = add({
				constantCoefficient(u1),
				constantCoefficient(u2)
			});

			AST* P_ = mul({ reduceAdditionAST(S_), nonConstantCoefficient(u1) });
			AST* P = reduceMultiplicationAST(P_);
			
			delete S_;
			delete P_;
			delete nc_u1;
			delete nc_u2;

			if(P->kind() == Kind::Integer && P->value() == 0) {
				delete P;
				return list({});
			}
			
			return list({ P });
		}
	
		delete nc_u1;
		delete nc_u2;


		if(orderRelation(u2, u1)) {
			return list({u2->deepCopy(), u1->deepCopy()});
			// AST* L_ = list({u2->deepCopy(), u1->deepCopy()});
			// AST* R = simplifyAdditionRec(L_);
			// delete L_;
			// return R;
		}

		// std::vector<AST*> L_;

		// for(AST* k : L)
		// 	L_.push_back(k->deepCopy());
		
		return list({u1->deepCopy(), u2->deepCopy()});
	}

	if(
		L->numberOfOperands() == 2 &&
		(
			L->operand(0)->kind() == Kind::Addition ||
			L->operand(1)->kind() == Kind::Addition
		)
	) {
		AST* u1 = L->operand(0);
		AST* u2 = L->operand(1);

		if(
			u1->kind() == Kind::Addition &&
			u2->kind() == Kind::Addition
		) {
			AST* U1 = new AST(Kind::List);
			AST* U2 = new AST(Kind::List);
			
			for(int i=0; i<u1->numberOfOperands(); i++)
				U1->includeOperand(u1->operand(i)->deepCopy());
	
			for(int i=0; i<u2->numberOfOperands(); i++)
				U2->includeOperand(u2->operand(i)->deepCopy());

			AST* L_ = mergeAdditions(U1, U2);
			
			delete U1;
			delete U2;

			return L_;
		}
	
		if(u1->kind() == Kind::Addition) {
			AST* U1 = new AST(Kind::List);
			AST* U2 = new AST(Kind::List);
			
			for(int i=0; i<u1->numberOfOperands(); i++)
				U1->includeOperand(u1->operand(i)->deepCopy());
	
			U2->includeOperand(u2->deepCopy());
			
			AST* L_ = mergeAdditions(U1, U2);
			
			delete U1;
			delete U2;

			return L_;
		}
	
		if(u2->kind() == Kind::Addition) {
			AST* U1 = new AST(Kind::List);
			AST* U2 = new AST(Kind::List);
			
			for(int i=0; i<u2->numberOfOperands(); i++)
				U2->includeOperand(u2->operand(i)->deepCopy());

			U1->includeOperand(u1->deepCopy());

			AST* L_ = mergeAdditions(U1, U2);

			delete U1;
			delete U2;

			return L_;
		}
	}

	AST* u1 = L->operand(0);

	AST* restL = rest(L);
	
	AST* w = simplifyAdditionRec(restL);

	delete restL;

	if(u1->kind() == Kind::Addition) {
		AST* U1;
		
		for(int i=0; i<u1->numberOfOperands(); i++)
			U1->includeOperand(u1->operand(i)->deepCopy());

		AST* L_ = mergeAdditions(U1, w);

		delete U1;
		delete w;
		
		return L_;
	}

	AST* U1 = new AST(Kind::List);

	U1->includeOperand(u1->deepCopy());

	AST* L_ = mergeAdditions(U1, w);
	
	delete U1;
	delete w;

	return L_;
}

AST* reduceAdditionAST(AST* u) {

	if(u->kind() == Kind::Undefined)
		return undefined();
	
	if(u->numberOfOperands() == 1)
		return u->operand(0)->deepCopy();

	AST* L = new AST(Kind::List);
	
	for(int i=0; i<u->numberOfOperands(); i++)
		L->includeOperand(u->operand(i)->deepCopy());

	AST* R = simplifyAdditionRec(L);
	
	delete L;

	if(R->numberOfOperands() == 0) {
		delete R;
		return integer(0);
	}

	if(R->numberOfOperands() == 1) {
		AST* r = R->operand(0)->deepCopy();
		delete R;
		return r;
	}

	AST* res = new AST(Kind::Addition);
	
	for(int i=0; i<R->numberOfOperands(); i++) {
		res->includeOperand(R->operand(i)->deepCopy());
	}
	
	delete R;

	return res;
}

}

