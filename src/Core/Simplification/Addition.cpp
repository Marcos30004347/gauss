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

AST* mergeAdditions(AST* p, AST* q) {
	// return a copy of q

	if(p->numberOfOperands() == 0) {
		return q->copy();
	}

	// return a copy of p
	if(q->numberOfOperands() == 0) {
		return p->copy();
	}

	AST* L = list({ p->operand(0)->copy(), q->operand(0)->copy() });
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
	if(a->kind() == Kind::FunctionCall) {
		bool non_constant = false;
		
		for(long i=0; i<a->numberOfOperands(); i++) {
			if(!isConstant(a->operand(i))) {
				non_constant = true;
				break;
			}
		}
	
		if(non_constant) {
			return a->copy();
		}
	
		return undefined();
	}

	if(a->kind() == Kind::Power) {
		if(!isConstant(a->operand(0)) || !isConstant(a->operand(1)))
			return a->copy();
	}

	AST* res = new AST(Kind::Multiplication);

	for(long i=0; i<a->numberOfOperands(); i++) {
		if(a->operand(i)->kind() == Kind::FunctionCall) {
			AST* k = nonConstantCoefficient(a->operand(i));
			if(k->kind() == Kind::Undefined) {
				delete k;
			} else {
				res->includeOperand(k);
			}
		} else if(!isConstant(a->operand(i))) {
			res->includeOperand(a->operand(i)->copy());
		}
	}

	if(res->numberOfOperands() == 0) {
		delete res;
		return undefined();
	}

	if(res->numberOfOperands() == 1) {
		AST* r = res->operand(0)->copy();
		delete res;
		return r;
	}

	return res;
}

AST* constantCoefficient(AST* a) {
	if(a->kind() == Kind::FunctionCall) {
		bool non_constant = false;
		
		for(long i=0; i<a->numberOfOperands(); i++) {
			if(!isConstant(a->operand(i))) {
				non_constant = true;
				break;
			}
		}
	
		if(non_constant) {
			return integer(1);
		}

		return a->copy();
	}

	if(a->kind() == Kind::Power) {
		if(!isConstant(a->operand(0)) || !isConstant(a->operand(1)))
			return integer(1);

		return a->copy();
	}

	AST* res = new AST(Kind::Multiplication);
	
	for(long i=0; i<a->numberOfOperands(); i++) {
		if(a->operand(i)->kind() == Kind::FunctionCall) {
			AST* k = constantCoefficient(a->operand(i));
			if(k->kind() == Kind::Integer && k->value() == 1) {
				delete k;
			} else {
				res->includeOperand(k);
			}
		} else if(isConstant(a->operand(i)))
			res->includeOperand(a->operand(i)->copy());
	}

	if(res->numberOfOperands() == 0) {
		delete res;
		return integer(1);
	}

	if(res->numberOfOperands() == 1) {
		AST* r = res->operand(0)->copy();
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
	) 
	{
		AST* u1 = L->operand(0);
		AST* u2 = L->operand(1);

		if(isConstant(u1) && isConstant(u2)) 
		{
			AST* P_ = add({u1->copy(), u2->copy()});

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
			{
				return { undefined() };
			}
	
			return list({new AST(Kind::Infinity)});
		} 

		if(u2->kind() == Kind::MinusInfinity) 
		{
			if(u1->kind() == Kind::Infinity)
			{
				return list({undefined()});
			}
	
			return list({new AST(Kind::MinusInfinity)});
		} 

		if(u1->kind() == Kind::Integer && u1->value() == 0) 
		{
			return list({u2->copy()});
		}
	
		if(u2->kind() == Kind::Integer && u2->value() == 0) 
		{
			return list({u1->copy()});
		}

		AST* nc_u1 = nonConstantCoefficient(u1);
		AST* nc_u2 = nonConstantCoefficient(u2);
	
		if(nc_u1->match(nc_u2)) 
		{
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

		if(orderRelation(u2, u1)) 
		{
			return list({u2->copy(), u1->copy()});
		}

		return list({u1->copy(), u2->copy()});
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
			
			for(long i=0; i<u1->numberOfOperands(); i++)
				U1->includeOperand(u1->operand(i)->copy());
	
			for(long i=0; i<u2->numberOfOperands(); i++)
				U2->includeOperand(u2->operand(i)->copy());

			AST* L_ = mergeAdditions(U1, U2);
			
			delete U1;
			delete U2;

			return L_;
		}
	
		if(u1->kind() == Kind::Addition) {
			AST* U1 = new AST(Kind::List);
			AST* U2 = new AST(Kind::List);
			
			for(long i=0; i<u1->numberOfOperands(); i++)
				U1->includeOperand(u1->operand(i)->copy());
	
			U2->includeOperand(u2->copy());
			
			AST* L_ = mergeAdditions(U1, U2);
			
			delete U1;
			delete U2;

			return L_;
		}
	
		if(u2->kind() == Kind::Addition) {
			AST* U1 = list({});
			AST* U2 = list({});
			
			for(long i=0; i<u2->numberOfOperands(); i++)
				U2->includeOperand(u2->operand(i)->copy());

			U1->includeOperand(u1->copy());

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
		AST* U1 = list({});
		
		for(long i=0; i<u1->numberOfOperands(); i++)
			U1->includeOperand(u1->operand(i)->copy());

		AST* L_ = mergeAdditions(U1, w);

		delete U1;
		delete w;
		
		return L_;
	}

	AST* U1 = list({});

	U1->includeOperand(u1->copy());

	AST* L_ = mergeAdditions(U1, w);
	
	delete U1;
	delete w;

	return L_;
}

AST* reduceAdditionAST(AST* u) {
	if(u->kind() == Kind::Undefined)
	{
		return undefined();
	}
	
	if(u->numberOfOperands() == 1)
	{
		return u->operand(0)->copy();
	}

	AST* L = list({});
	
	for(long i=0; i<u->numberOfOperands(); i++)
	{
		L->includeOperand(u->operand(i)->copy());
	}

	AST* R = simplifyAdditionRec(L);
	
	delete L;

	if(R->numberOfOperands() == 0) 
	{
		delete R;

		return integer(0);
	}

	if(R->numberOfOperands() == 1) 
	{
		AST* r = R->operand(0);
	
		R->removeOperand(0L);
	
		delete R;
		
		return r;
	}

	AST* res = new AST(Kind::Addition);
	
	while(R->numberOfOperands())
	{
		res->includeOperand(R->operand(0));
		R->removeOperand(0L);
	}
	
	delete R;

	return res;
}

}

