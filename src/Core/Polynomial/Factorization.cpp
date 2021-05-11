#include "Zp.hpp"
#include "Factorization.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"

#include <cmath>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace prime;

namespace polynomial {

AST* genExtendSigmaP(AST* V, AST* x, unsigned p) {
	if(V->numberOfOperands() == 2) {
		AST* k = extendedEuclideanAlgGPE_Sp(V->operand(0), V->operand(1), x, p);

		AST* A = k->operand(1)->deepCopy();
		AST* B = k->operand(2)->deepCopy();

		delete k;

		return list({ A, B });
	}

	AST* V_ = new AST(Kind::List);

	for(int i=0; i<V->numberOfOperands() - 1; i++)
		V_->includeOperand(V->operand(i)->deepCopy());

	AST* tal = genExtendSigmaP(V_, x, p);

	AST* gs_ = new AST(Kind::Multiplication);
	for(int j=0; j<V_->numberOfOperands(); j++) {
		gs_->includeOperand(V_->operand(j)->deepCopy());
	}

	AST* gs = Ts(gs_, x, p);
	
	AST* k = extendedEuclideanAlgGPE_Sp(V->operand(V->numberOfOperands() - 1), gs, x, p);

	delete gs;
	delete gs_;
	delete V_;

	AST* A = k->operand(1)->deepCopy();
	AST* B = k->operand(2)->deepCopy();
	
	delete k;
	
	AST* thetas = new AST(Kind::List);
	
	for(int i=0; i<tal->numberOfOperands(); i++) {
		AST* thetha = mul({ A->deepCopy(), tal->operand(i)->deepCopy() });
		thetas->includeOperand(Ts(thetha, x, p));
		delete thetha;
	}

	thetas->includeOperand(B);

	delete tal;

	return thetas;
}

AST* genExtendRP(AST* V, AST* S, AST* F, AST* x, unsigned p) {
	AST* Rs = new AST(Kind::List);

	for(int i=0; i<V->numberOfOperands(); i++) {
		AST* u_ = mul({ F, S->operand(i)->deepCopy() });
		AST* u = Ts(u_, x, p);

		AST* ri = remainderGPE_Sp(u, V->operand(i), x, p);
		Rs->includeOperand(ri);
	}

	return Rs;
}

// u is a polynomial in x, findPrime return
// and integer such tath p % leadCoeff(u,x) != 0
int findPrime(AST* u, AST* x) {
	AST* lc_ = leadingCoefficientGPE(u, x);
	AST* lc = expandAST(lc_);

	int p = nth_prime(0);

	// We are gonna look just for the first 5000000 primes
	for(int i=0; i < 5000000; i++) {
		if(lc->value() % nth_prime(i) != 0) {
			p = nth_prime(i);
			break;
		}
	}
	
	delete lc_, lc;
	return p;
}


unsigned long abs(signed long i) {
	if(i >= 0) return i;
	return -1*i;
}

// Get height of polynomial in Z[x]
AST* polynomialHeight_Z(AST* u, AST* x) {
	AST* u_ = expandAST(u);
	AST* d_ = degreeGPE(u_, x);
	
	unsigned long d = d_->value();
	unsigned long h = 0;

	for(int i=d; i>=0; i++) {
		AST* d = integer(i);
		AST* c = coefficientGPE(u_, x, d);
		
		unsigned long h_ = abs(c->value());
	
		if(h_ > h) 
			h = h_;
		
		delete d, c;
	}
	
	delete u_, d_;

	return integer(h);
}

unsigned long log(double base, int x) {
    return (unsigned long)(std::log(x) / std::log(base));
}

unsigned long findK(AST* u, AST* x, int p) {
	AST* h = polynomialHeight_Z(u, x);
	AST* n_ = degreeGPE(u, x);
	unsigned long n = n_->value();
	
	double B = std::pow(2, n) * std::sqrt(n+1) * h->value();
	
	delete h;
	
	return log((unsigned long)std::ceil(2*B), p);
}

AST* trueFactors(AST* u, AST* l, AST* x, int p, int k) {
	AST* U = u->deepCopy();
	AST* L = l->deepCopy();

	AST* factors = list({});

	int m = 1;

	while(m < L->numberOfOperands()/2) {
		AST* m_ = integer(m);
		AST* C = combination(L, m_); 
		delete m_;

		while(C->numberOfOperands() != 0) {
			AST* t = C->operand(0);

			AST* T_ = new AST(Kind::Multiplication);
			for(int i=0; i<t->numberOfOperands(); i++)
					T_->includeOperand(t->operand(i)->deepCopy());
	
			AST* T = Ts(T_, x, (int)std::pow(p, k));

			delete T_;
	
			AST* D = divideGPE(U,T,x);

			AST* Q = D->operand(0);	
			AST* R = D->operand(1);	

			if(R->kind() == Kind::Integer && R->value() == 0) {
				factors->includeOperand(T->deepCopy());
				U = Q;

				AST* L_ = L;
				L = remove(L, t);//  // L ~ t
				delete L_;
			
				AST* C_ = C;
				C = cleanUp(C, t);
				delete C_;
		
			} else {
				AST* T_ = new AST(Kind::Set, { t->deepCopy() });
				AST* C_ = C;
				
				C = difference(C, T_); // C ~ {t};
				
				delete T_;
				delete C_;
			}

			delete D;
		}
		m = m + 1;
	}

	if( U->kind() != Kind::Integer && U->value() != 1) {
		factors->includeOperand(U->deepCopy());
	}

	return factors;
}

}
