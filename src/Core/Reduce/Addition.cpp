#include "Addition.hpp"
#include <vector>
#include <assert.h>
#include "Rationals.hpp"
#include "Multiplication.hpp"

using namespace ast;
using namespace algebra;

namespace reduce {


std::vector<AST*> simplifySummationRec(std::vector<AST*> L);
std::vector<AST*> adjoinSummations(AST* p, std::vector<AST*> q);
std::vector<AST*> mergeSummations(std::vector<AST*> p, std::vector<AST*> q);

std::vector<AST*> transformSummation(AST* a, AST* b);
std::vector<AST*> transformSummationAssociative(std::vector<AST*> L);
std::vector<AST*> transformSummationExpand(std::vector<AST*> L);

std::vector<AST*> restAddition(std::vector<AST*> p, int from = 1) {
    std::vector<AST*> a;

    for(int i=from; i <= p.size() - 1; i++) 
        a.push_back(p[i]->deepCopy());

    return a;
}

std::vector<AST*> adjoinSummations(AST* p, std::vector<AST*> q) {
	std::vector<AST*> tmp = std::vector<AST*>(0);

	int i = 0;
	bool inc = false;

	while(i != q.size()) {
		if(!inc) {
			std::vector<AST*> t = transformSummation(p, q[i]);
		
			if(t.size() == 1) {
				tmp.push_back(t[0]->deepCopy());
				i++;
				inc = true;
			} else if(orderRelation(p, q[i])) {
				inc = true;
				tmp.push_back(p->deepCopy());
			} else {
				tmp.push_back(q[i++]->deepCopy());				
			}

			for(AST* k : t)
				delete k;
		
		} else {
			tmp.push_back(q[i++]->deepCopy());
		}
	}

	if(!q.size() || !inc) {
			tmp.push_back(p->deepCopy());
	}

	return tmp;
}

std::vector<AST*> mergeSummations(std::vector<AST*> p, std::vector<AST*> q) {
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

	std::vector<AST*> H = simplifySummationRec({ p[0], q[0]});
	
	std::vector<AST*> a = restAddition(p);
	std::vector<AST*> b = restAddition(q);

	std::vector<AST*> R = mergeSummations(a, b);

	for(AST* k: a) delete k;
	for(AST* k: b) delete k;

	for(int i=0; i<H.size(); i++) {
		std::vector<AST*> _R = R;

		R = adjoinSummations(H[i], R);

		for(AST* k : _R)
			delete k;

		delete H[i];
	}

	return (R);
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

std::vector<AST*> transformSummation(AST* a, AST* b) {
	// printf("a %s\n", a->toString().c_str());
	// printf("b %s\n", b->toString().c_str());

	if(isConstant(a) && isConstant(b)) {
		AST* v = add({ a->deepCopy(), b->deepCopy() });
		AST* res =	reduceRNEAST(v); 

		destroyASTs({v});

		return { res };
	}

	AST* nc_a = nonConstantCoefficient(a);
	AST* nc_b = nonConstantCoefficient(b);
	// printf("nc_a %s\n", nc_a->toString().c_str());
	// printf("nc_b %s\n", nc_b->toString().c_str());
	// printf("c_a %s\n", constantCoefficient(a)->toString().c_str());
	// printf("c_b %s\n", constantCoefficient(b)->toString().c_str());

	if(
		nc_a->kind() != Kind::Undefined &&
		nc_b->kind() != Kind::Undefined &&
		nc_a->match(nc_b)
	) {
	
		AST* k = add({
			constantCoefficient(a),
			constantCoefficient(b)
		});
		
		AST* m = mul({
			reduceRNEAST(k),
			nonConstantCoefficient(a)
		});
	
		AST* res = reduceMultiplicationAST(m);
	
		destroyASTs({ nc_a, nc_b, k, m });
	
		return { res };
	}

	destroyASTs({ nc_a, nc_b });

	if(orderRelation(a, b))
		return { a->deepCopy(), b->deepCopy() };

	return { b->deepCopy(), a->deepCopy() };
}


// simplify((a+b)+c) = simplify(a+b+c)
std::vector<AST*> transformSummationAssociative(std::vector<AST*> L) {
	std::vector<AST*> H;

	for(int j=0; j<L[0]->numberOfOperands(); j++) 
		H.push_back(L[0]->operand(j)->deepCopy());

	std::vector<AST*> restL 		= restAddition(L);
	std::vector<AST*> merged 		=	mergeSummations(H, restL);
	std::vector<AST*> response 	= simplifySummationRec(merged);

	for(AST* t : restL) 	delete t;
	for(AST* t : merged)	delete t;
	for(AST* t : H)  			delete t;

	return response;
}

std::vector<AST*> transformSummationExpand(std::vector<AST*> L) {
	if(L.size() == 2 && L[1]->kind() != Kind::Addition)
		return transformSummation(L[0], L[1]);

	std::vector<AST*> rest_L 	= restAddition(L);
	std::vector<AST*> simp 		= simplifySummationRec(rest_L);
	std::vector<AST*> res 		=	mergeSummations({ L[0] }, simp); 

	for(AST* t : rest_L) 	delete t;
	for(AST* t : simp) 		delete t;

	return res;
}

std::vector<AST*> simplifySummationRec(std::vector<AST*> L) {

	if(L.size() == 0 || (L.size() == 1 && L[0]->kind() != Kind::Addition)) {
		std::vector<AST*> l;

		for(AST* k: L) 
			l.push_back(k->deepCopy());

		return l;
	}



	// 0+a = a
	if(L[0]->kind() == Kind::Integer && L[0]->value() == 0) {
		std::vector<AST*> rest_L = restAddition(L);
		std::vector<AST*> res = simplifySummationRec(rest_L);
		
		for(AST* k : rest_L)
			delete k;
		
		return res;
	}


	// (a+b)+(c+d) = a+b+c+d
	if(L[0]->kind() == Kind::Addition)
		return transformSummationAssociative(L);

	return transformSummationExpand(L);       
}

AST* reduceAdditionAST(ast::AST* u) {

	if(u->kind() == Kind::Undefined)
		return new AST(Kind::Undefined);


	if(u->numberOfOperands() == 1)
		return u->operand(0)->deepCopy();

	std::vector<AST*> L;
	for(int i=0; i < u->numberOfOperands(); i++)
		L.push_back(u->operand(i)->deepCopy());


	std::vector<AST*> R = simplifySummationRec(L);
	
	for(AST* k : L)
		delete k;

	if(R.size() == 1)
		return R[0];

	if(R.size() == 0)
		return inte(0);

	AST* res = new AST(Kind::Addition);
	
	for(AST* t : R)
		res->includeOperand(t);
	
	return res;
}

}

