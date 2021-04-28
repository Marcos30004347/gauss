#include "Multiplication.hpp"
#include "Addition.hpp"
#include "Rationals.hpp"
#include "Power.hpp"
#include "Core/Expand/Expand.hpp"

#include <vector>

using namespace ast;
using namespace expand;
using namespace algebra;

namespace reduce {

std::vector<AST*> simplifyProductRec(std::vector<AST*> L);
std::vector<AST*> adjoinProducts(AST* p, std::vector<AST*> q);
std::vector<AST*> mergeProducts(std::vector<AST*> p, std::vector<AST*> q);

std::vector<AST*> transformProduct(AST* a, AST* b);
std::vector<AST*> transformProductSimplify(std::vector<AST*> L);
std::vector<AST*> transformProductDistributive(AST* a, AST* b);
std::vector<AST*> transformProductAssociative(std::vector<AST*> L);

std::vector<AST*> restMultiplication(std::vector<AST*> p, int from = 1) {
    std::vector<AST*> a;

    for(int i=from; i <= p.size() - 1; i++) 
        a.push_back(p[i]->deepCopy());

    return a;
}

// insert p in q so that q stills in order. This function assumes
// that q is in order when executing.
std::vector<AST*> adjoinProducts(AST* p, std::vector<AST*> q) {
	std::vector<AST*> tmp = std::vector<AST*>(0);

	int i = 0;
	bool inc = false;

	while(i != q.size()) {
		if(!inc) {
			std::vector<AST*> t = transformProduct(p, q[i]);
		
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

	std::vector<AST*> a = restMultiplication(p);
	std::vector<AST*> b = restMultiplication(q);

	std::vector<AST*> R = mergeProducts(a, b);
	
	for(AST* k: a) delete k;
	for(AST* k: b) delete k;

	for(int i=0; i<H.size(); i++) {
		std::vector<AST*> _R = R;

		R = adjoinProducts(H[i], R);

		for(AST* k : _R)
			delete k;

		delete H[i];
	}


	return R;
}


// std::vector<AST*> transformProductDistributive(AST* a, AST* b) {
// 	if(a->kind() == Kind::Addition || b->kind() == Kind::Addition) {
// 		AST* m = mul({a->deepCopy(), b->deepCopy()});
// 		AST* res = expandAST(m);
// 		delete m;
// 		return { res };
// 	}

// 	// if(a->kind() == Kind::Addition && b->kind() == Kind::Addition) {
// 	// 	AST* u = new AST(Kind::Addition);
	
// 	// 	for(int i=0; i<a->numberOfOperands(); i++) {
// 	// 		for(int j=0; j<b->numberOfOperands(); j++) {
				
// 	// 			std::vector<AST*> t =	mergeProducts(
// 	// 				{ a->operand(i) },
// 	// 				{ b->operand(j) }
// 	// 			);

// 	// 			if(t.size() == 1) {
// 	// 				u->includeOperand(t[0]);
// 	// 			} else {
// 	// 				AST* m = new AST(Kind::Multiplication);
					
// 	// 				for(AST* v : t)
// 	// 					m->includeOperand(v);
					
// 	// 				u->includeOperand(m);
// 	// 			}
// 	// 		}
// 	// 	}

	
// 	// 	AST* r = reduceAdditionAST(u);
		
// 	// 	delete u;
	
// 	// 	return { r };
// 	// }

// 	// if(a->kind() != Kind::Addition && b->kind() == Kind::Addition) {
// 	// 	// printf("ASDASDASD\n");
// 	// 	AST* u = new AST(Kind::Addition);
		
// 	// 	for(int j=0; j<b->numberOfOperands(); j++) {
// 	// 		std::vector<AST*> t =	mergeProducts(
// 	// 				{ b->operand(j) },
// 	// 				{ a }
// 	// 		);
// 	// 		if(t.size() == 1) {
// 	// 			u->includeOperand(t[0]);
// 	// 		} else {
// 	// 			AST* m = new AST(Kind::Multiplication);
				
// 	// 			for(AST* v : t)
// 	// 				m->includeOperand(v);
				
// 	// 			u->includeOperand(m);
// 	// 		}
// 	// 		// AST* m = new AST(Kind::Multiplication);
			
// 	// 		// for(AST* v : t)
// 	// 		// 	m->includeOperand(v);
			
// 	// 		// u->includeOperand(m);
// 	// 	}
// 	// 	// printf("U: ");
// 	// 	// u->print();
// 	// 	// printf("\n");
// 	// 	AST* r = reduceAdditionAST(u);
// 	// 	// r->print();
// 	// 	// printf("\n");
// 	// 	// printf("\n");
// 	// 	delete u;
// 	// 	return { r };
// 	// }

// 	// if(a->kind() == Kind::Addition && b->kind() != Kind::Addition) {
// 	// 	AST* u = new AST(Kind::Addition);
// 	// 	for(int j=0; j<a->numberOfOperands(); j++) {
// 	// 		std::vector<AST*> t =	mergeProducts(
// 	// 				{ b },
// 	// 				{ a->operand(j) }
// 	// 		);

// 	// 		if(t.size() == 1) {
// 	// 			u->includeOperand(t[0]);
// 	// 		} else {
// 	// 			AST* m = new AST(Kind::Multiplication);
				
// 	// 			for(AST* v : t)
// 	// 				m->includeOperand(v);
				
// 	// 			u->includeOperand(m);
// 	// 		}
// 	// 		// AST* m = new AST(Kind::Multiplication);
			
// 	// 		// for(AST* v : t)
// 	// 		// 	m->includeOperand(v);
			
// 	// 		// u->includeOperand(m);
// 	// 	}

// 	// 	AST* r = reduceAdditionAST(u);
// 	// 	delete u;
// 	// 	return { r };
// 	// }

// 	return { a->deepCopy(), b->deepCopy() };
// }


// simplify(a*a) = a^2, simplify(a*b) = a*b
std::vector<AST*> transformProduct(AST* a, AST* b) {
	if(b->kind() == Kind::Integer && b->value() == 1)
		return { a->deepCopy() };

	if(a->kind() == Kind::Integer && a->value() == 1)
		return { b->deepCopy() };

	if(isConstant(a) && isConstant(b)) {
		AST* t = mul({ a->deepCopy(), b->deepCopy() });
		AST* k = reduceRNEAST(t);
		delete t;
		return { k };
	}

	AST* base_a = base(a);
	AST* base_b = base(b);

	if(base_a->match(base_b)) {
	
		AST* e = add({ exp(a), exp(b) });
		AST* p = pow(base(a), reduceAdditionAST(e));
		AST* r = reducePowerAST(p);
	
		destroyASTs({ p, e, base_a, base_b });
	
		return { r };
	}

	destroyASTs({ base_a, base_b });

	if(orderRelation(a, b))
			return { a->deepCopy(), b->deepCopy() };

	return { b->deepCopy(), a->deepCopy() };
}

// simplify((a*b)*c) = simplify(a*b*c)
std::vector<AST*> transformProductAssociative(std::vector<AST*> L) {
	std::vector<AST*> H;

	for(int j=0; j < L[0]->numberOfOperands(); j++) 
		H.push_back(L[0]->operand(j)->deepCopy());

	std::vector<AST*> restL 		= restMultiplication(L);
	std::vector<AST*> merged 		=	mergeProducts(H, restL);
	std::vector<AST*> response 	= simplifyProductRec(merged);        
	// std::vector<AST*> simp 	= simplifyProductRec(restL);

	for(AST* t : restL) 	delete t;
	for(AST* t : merged)	delete t;
	for(AST* t : H)  			delete t;

	return response;
}

std::vector<AST*> transformProductSimplify(std::vector<AST*> L) {
	if(L.size() == 2) {
		if(
			L[0]->kind() == Kind::Addition ||
			L[1]->kind() == Kind::Addition
		) {
			AST* m = mul({L[0]->deepCopy(), L[1]->deepCopy()});
			AST* res = expandAST(m);
			delete m;
			return { res };
		} 
	
		if(
			L[0]->kind() != Kind::Multiplication &&
			L[1]->kind() != Kind::Multiplication
		) return transformProduct(L[0], L[1]);
	}

	std::vector<AST*> rest_L 	= restMultiplication(L);
	std::vector<AST*> simp 		= simplifyProductRec(rest_L);
	std::vector<AST*> res 		=	mergeProducts({ L[0] }, simp); 

	for(AST* t : rest_L) 	delete t;
	for(AST* t : simp) 		delete t;
	return res;
}



std::vector<AST*> simplifyProductRec(std::vector<AST*> L) {
	if(L.size() == 0)
		return {};

	if(L.size() == 0 || (L.size() == 1 && L[0]->kind() != Kind::Multiplication)) {
		std::vector<AST*> l;

		for(AST* k: L) 
			l.push_back(k->deepCopy());

		return l;
	}

	// 0 * a = 0
	if(L[0]->kind() == Kind::Integer && L[0]->value() == 0)
			return {};

	// if(L[0]->kind() == Kind::Integer && L[0]->value() == 1) {
	// 	std::vector<AST*> rest_L = restMultiplication(L);
	// 	std::vector<AST*> res = simplifyProductRec(rest_L);
	// 	for(AST* k : rest_L)
	// 		delete k;
	// 	return res;
	// }

	// 1 * a = a
	if(L[0]->kind() == Kind::Integer && L[0]->value() == 1) {
		std::vector<AST*> rest_L = restMultiplication(L);
		std::vector<AST*> simp_rest_L = simplifyProductRec(rest_L);

		for(AST* k: rest_L)
			delete k;

		return simp_rest_L;
	}

	// (a*b)*(c*d) = a*b*c*d
	if(L[0]->kind() == Kind::Multiplication)
		return transformProductAssociative(L);

	return transformProductSimplify(L);       
}

AST* reduceMultiplicationAST(AST* u) {
	if(u->kind() == Kind::Undefined)
		return new AST(Kind::Undefined);
	
	for(int i=0; i<u->numberOfOperands(); i++) {
		AST* o = u->operand(i);
		if(o->kind() == Kind::Integer && o->value() == 0)
			return inte(0);
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
		return inte(0);
	
	AST* res = new AST(Kind::Multiplication);
	
	for(AST* t : R)
		res->includeOperand(t);
	
	return res;
}

}

