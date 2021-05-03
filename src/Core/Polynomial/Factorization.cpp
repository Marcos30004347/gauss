#include "Factorization.hpp"
#include "Zp.hpp"
#include "Core/Expand/Expand.hpp"

using namespace ast;
using namespace expand;
using namespace algebra;

namespace polynomial {




// genExtendSigmaP(V, x, p)
//
// V is a list of two or more relatively prime,
// monic polynomials in Zp[x]
// x is a symbol and p is a prime number and Zp is represented with the
// symmetric representation
//
// return list [σ0, ..., σi]
// where sum(i=1 to len(V)){ σi*gi } = 1
// v = V[1]*...*V[len(V)] 
// g[i] = v/v[i]
// 
// Ex:
//
// BASE CASE s = len(V) = 2
// v = V[1]*V[2], g[1]=v/V[1] = V[2], g[2] = v/V[1] = V[2]
// A*g[1] + B*g[2] = 1, A and B can be obtained with extended euclidean algorithm
//
// INDUCTION s = len(V) > 2
// For s' = s - 1 we have shown that for h[i] = (v[1]*...*v[s'])/(v[s]*v[i]) and
// note that h[i] is just the g[i] for the base case because we remove v[s] from v = v[0]*...*v[s']
// sum(i=1 to s' - 1){τ[i]*h[i]} = 1, note that τ[i] is the σ[i] from the base case.
//
// So we can define:
// v[s] = v[s]*sum(i=1 to s' - 1){τ[i]*h[i]} = sum(i=1 to s' - 1){τ[i]*g[i]}
// therefore 1 = A*sum(i=1 to s' - 1){τ[i]*g[i]} + B*g[s] = 1
//
// So we can obtain A and B using the extended euclidedan algorithm
// and σ[i] = 'A*τ[i]' for i < s and 'B' for i = s
std::vector<AST*> genExtendSigmaP(std::vector<AST*> V, AST* x, unsigned p) {
	// assert V.size() >= 2
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

	AST* g = Tss(g_, x, p);
	
	std::vector<AST*> k = extendedEuclideanAlgGPE_Sp(V[V.size() - 1], g, x, p);
	AST* gcd = k[0];
	// todo: assert(gcd == 1)
	delete gcd;

	AST* A = k[1];
	AST* B = k[2];

	std::vector<AST*> theta;
	
	for(int i=0; i<tal.size(); i++) {
		AST* t_ = mul({A->deepCopy(), tal[i]});
		theta.push_back(Tss(t_, x, p));
		delete t_;
	}

	theta.push_back(B);

	return theta;
}

std::vector<AST*> genExtendRP(std::vector<AST*> V, std::vector<AST*> S, AST* F, AST* x, unsigned p) {
	std::vector<AST*> rs;

	for(int i=0; i<V.size(); i++) {
		AST* u_ = mul({F, S[i]});
		AST* u = Tss(u_, x, s);

		AST* ri = remainderGPE_Sp(u, V[i], x, p);
		rs.push_back(ri);
	}

	return rs;
}

// #define MAX_SIZE 1000005
// void sieveOfEratosthenes(std::vector<int>& primes) {
// 	bool IsPrime[MAX_SIZE];
// 	memset(IsPrime, true, sizeof(IsPrime));

// 	for (int p = 2; p * p < MAX_SIZE; p++) {
// 		if (IsPrime[p] == true) {
// 			for (int i = p * p; i < MAX_SIZE; i += p)
// 				IsPrime[i] = false;
// 		}
// 	}

// 	for (int p = 2; p < MAX_SIZE; p++)
// 		if (IsPrime[p])
// 			primes.push_back(p);
// }

// u is a polynomial in x, findPrime return
// and integer such tath p % leadCoeff(u,x) != 0
int findPrime(AST* u, AST* x) {
	// TODO: calculate primes prior to this call

	// AST* lc_ = leadingCoefficientGPE(u, x);
	// AST* lc = expandAST(lc_);

	// int p = primes[0];

	// for(int i=0; i < 32768; i++) {
	// 	if(lc->value() % primes[i] != 0) {
	// 		p = primes[i];
	// 		break;
	// 	}
	// }
	
	// delete lc_, lc;
	// return p;
}

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
