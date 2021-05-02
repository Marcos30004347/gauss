#include "Factorization.hpp"
#include "Core/Expand/Expand.hpp"

using namespace ast;
using namespace expand;
using namespace algebra;

namespace polynomial {



std::vector<AST*> trueFactors(AST* u, AST* l, AST* x, AST* p, AST* k) {

	AST* U = u->deepCopy();
	AST* L = l->deepCopy();

	std::vector<AST*> factors = std::vector<AST*>(0);

	int m = 1;

	while(m < L->numberOfOperands()/2) {
		// TODO
		// All sets that contain m elements
		// of the set L, 
		// comb({a,b,c,d}, 2) -> {{a,b}, {a,c}, {a,d}, {b,c}, {b,d}, {c, d}}
		AST* C = comb(L, m); 

		while(C->kind() != Kind::Integer || C->value() == 0) {
			AST* t = C->operand(0);

			AST* T = new AST(Kind::Multiplication);
			for(int i=0; i<t->numberOfOperands(); i++)
				T->includeOperand(t->operand(i)->deepCopy());

			// TODO
			// Let m ≥ 2 be an integer, and let u = an xn + ··· + a0 be in Z[x].
			//
			// 1. For the non-negative representation of Zm, define
			// 		Tm(u) = irem(an, m) xn + ··· + irem(a0, m).
			//
			// 2. For the symmetric representation of Zm, define
			//		Tm(u) = Sm(irem(an, m))xn + ··· + Sm(irem(a0, m)).
			// 
			// Sm(b) =  b, 		if 0<=b<=iquot(m,2),
			//				  b-m, 	if iquot(m,2) < b < m
			T = Ts(expandAST(T), x->deepCopy(), pow(p->deepCopy(), k->deepCopy())); // todo

			std::pair<AST*, AST*> D = divideGPE(U,T,x);
			
			if(D.second->kind() == Kind::Integer && D.second->value() == 0) {
				factors.push_back(T->deepCopy());
				U = D.first;
				L = L ~ t; // set difference

				//TODO
				// cleanUp(C, t) receives C, a set of sets and t a set
				// and return a new set with the members s of C such that
				// s ∩ t = ∅
				// C = {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5}, {3, 4}, {3, 5}, {4, 5}}
				// and t = {1, 2}, then Clean up(C, t) → {{3, 4}, {3, 5}, {4, 5}}.
				C = cleanUp(C, t); // todo
			} else {
				C = C ~ {t}; // set difference
			}
		}
		m = m + 1;
	}

	if( U->kind() != Kind::Integer && U->value() != 1) {
		factors.push_back(U->deepCopy());
	}

	return factors;
}

}
