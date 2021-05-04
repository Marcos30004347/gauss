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

std::vector<AST*> genExtendSigmaP(std::vector<AST*> V, AST* x, unsigned p) {
	if(V.size() == 2) {
		std::vector<AST*> g;

		g.push_back(V[1]->deepCopy());
		g.push_back(V[0]->deepCopy());

		std::vector<AST*> k = extendedEuclideanAlgGPE_Sp(g[0], g[1], x, p);

		AST* gcd = k[0];
		// assert(gcd == 1)
		delete gcd;

		AST* A = k[1];
		AST* B = k[2];

		return { A, B };
	}

	std::vector<AST*> V_;
	for(int i=0; i<V.size() - 1; i++) {
		V_.push_back(V[i]->deepCopy());
	}

	std::vector<AST*> tal = genExtendSigmaP(V_, x, p);

	AST* g_ = new AST(Kind::Multiplication);
	for(int j=0; j<V_.size(); j++) {
		g_->includeOperand(V_[j]->deepCopy());
	}

	AST* g = Ts(g_, x, p);
	
	std::vector<AST*> k = extendedEuclideanAlgGPE_Sp(V[V.size() - 1], g, x, p);
	AST* gcd = k[0];

	delete gcd;

	AST* A = k[1];
	AST* B = k[2];

	std::vector<AST*> theta;
	
	for(int i=0; i<tal.size(); i++) {
		AST* t_ = mul({A->deepCopy(), tal[i]});
		theta.push_back(Ts(t_, x, p));
		delete t_;
	}

	theta.push_back(B);

	return theta;
}

std::vector<AST*> genExtendRP(std::vector<AST*> V, std::vector<AST*> S, AST* F, AST* x, unsigned p) {
	std::vector<AST*> rs;

	for(int i=0; i<V.size(); i++) {
		AST* u_ = mul({F, S[i]});
		AST* u = Ts(u_, x, p);

		AST* ri = remainderGPE_Sp(u, V[i], x, p);
		rs.push_back(ri);
	}

	return rs;
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
unsigned long polynomialHeight_Z(AST* u, AST* x) {
	// Todo 

	AST* u_ = expandAST(u);
	AST* d_ = degreeGPE(u_, x);
	
	unsigned long d = d_->value();
	unsigned long h = 0;

	for(int i=d; i>=0; i++) {
		AST* p = pow(x->deepCopy(), inte(i));
		AST* c = coefficientGPE(u_, p);
		
		unsigned long h_ = abs(c->value());
	
		if(h_ > h) 
			h = h_;
		
		delete p, c;
	}
	
	delete u_, d_;

	return h;
}

unsigned long log(double base, int x) {
    return (unsigned long)(std::log(x) / std::log(base));
}

unsigned long findK(AST* u, AST* x, int p) {
	unsigned long h = polynomialHeight_Z(u, x);
	AST* n_ = degreeGPE(u, x);
	unsigned long n = n_->value();
	
	double B = std::pow(2, n) * std::sqrt(n+1) * h;

	return log((unsigned long)std::ceil(2*B), p);
}

AST* trueFactors(AST* u, AST* l, AST* x, int p, int k) {

	AST* U = u->deepCopy();
	AST* L = l->deepCopy();

	AST* factors = list({});

	int m = 1;

	while(m < L->numberOfOperands()/2) {
		AST* m_ = inte(m);
		AST* C = combination(L, m_); 
		delete m_;

		while(C->numberOfOperands() != 0) {
			AST* t = C->operand(0);

			AST* T_ = new AST(Kind::Multiplication);
			for(int i=0; i<t->numberOfOperands(); i++)
					T_->includeOperand(t->operand(i)->deepCopy());
	
			AST* T = Ts(T_, x, (int)std::pow(p, k));

			delete T_;
	
			std::pair<AST*, AST*> D = divideGPE(U,T,x);

			AST* Q = D.first;	
			AST* R = D.second;	

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
		}
		m = m + 1;
	}

	if( U->kind() != Kind::Integer && U->value() != 1) {
		factors->includeOperand(U->deepCopy());
	}

	return factors;
}

}
