#include "Zp.hpp"
#include "Factorization.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Primes/Primes.hpp"

#include <cmath>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace prime;

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

	AST* g = Ts(g_, x, p);
	
	std::vector<AST*> k = extendedEuclideanAlgGPE_Sp(V[V.size() - 1], g, x, p);
	AST* gcd = k[0];
	// todo: assert(gcd == 1)
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


// TODO
// All sets that contain m elements
// of the set L, 
// comb({a,b,c,d}, 2) -> {{a,b}, {a,c}, {a,d}, {b,c}, {b,d}, {c, d}}
void combUtil(std::vector<std::vector<AST*> >& ans, std::vector<AST*>& tmp, std::vector<AST*>& n, int left, int k) {
	if (k == 0) {
		ans.push_back(tmp);
		return;
	}

	for (int i = left; i < n.size(); ++i) {
		tmp.push_back(n[i]->deepCopy());
		combUtil(ans, tmp, n, i + 1, k - 1);
		tmp.pop_back();
	}
}

// Prints all combinations of size k of numbers
// from 1 to n.
std::vector<std::vector<AST*> > comb(std::vector<AST*>& n, int k) {
	std::vector<std::vector<AST*> > ans;
	std::vector<AST*> tmp;

	combUtil(ans, tmp, n, 0, k);

	return ans;
}



bool isListsEqual(std::vector<AST*> a, std::vector<AST*> b) {
	unsigned c = 0;

	for(int i=0; i<a.size(); i++) {
		for(int j=0; j<b.size(); j++) {
			if(a[i]->match(b[j])) {
				c++;
				break;
			}
		}
	}
	if(c >= a.size())
		return true;

	return false;
}

std::vector<AST*> list_diff(std::vector<AST*>& L, std::vector<AST*> N) {
	std::vector<AST*> L_;
	
	for(int i=0; i<L.size(); i++) {
		bool is_in_n = false;
		for(int j=0; j<N.size(); j++) {
			if(L[i]->match(N[j])) {
				is_in_n = true;
				break;
			}
			if(is_in_n) break;
		}
	
		if(!is_in_n) L_.push_back(L[i]->deepCopy());
	}

	return L_;
}

std::vector<std::vector<AST*>> set_difference(std::vector<std::vector<AST*>>& C, std::vector<std::vector<AST*>> N) {
	std::vector<std::vector<AST*>> C_;
	
	for(int i=0; i<C.size(); i++) {
		for(int t=0; t<N.size(); t++) {
			if(!isListsEqual(C[i], N[t])) {
				C_.push_back(std::vector<AST*>());

				for(int p=0; p<C[i].size(); p++)
					C_.back().push_back(C[i][p]->deepCopy());
			}
		}
	}

	return C_;
}


// cleanUp(C, t) receives C, a set of sets and t a set
// and return a new set with the members s of C such that
// s ∩ t = ∅
// C = {{1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5}, {3, 4}, {3, 5}, {4, 5}}
// and t = {1, 2}, then Clean up(C, t) → {{3, 4}, {3, 5}, {4, 5}}.
std::vector<std::vector<AST*>> cleanUp(std::vector<std::vector<AST*>>& C,std::vector<AST*>& n) {
	std::vector<std::vector<AST*>> C_;
	
	for(int i=0; i<C.size(); i++) {
		bool inc = false;
	
		for(int j=0; j<C[i].size(); j++) {
	
			for(int k=0; k<n.size(); k++) {
				if(n[i]->match(C[i][j])) {
					inc = true;
					break;
				}
			}
	
			if(inc) break;
		}
	
		if(inc) continue;
	
		C_.push_back(std::vector<AST*>());
		for(int p=0; p<C[i].size(); p++) {
			C_.back().push_back(C[i][p]->deepCopy());
		}
	}

	return C_;
}


std::vector<AST*> copyList(std::vector<AST*> l) {
	std::vector<AST*> l_;
	
	for(int i=0; i<l.size(); i++)
		l_.push_back(l[i]->deepCopy());
	
	return l_;
}

AST* build_product(std::vector<AST*> t) {
	AST* T = new AST(Kind::Multiplication);
		for(int i=0; i<t.size(); i++)
			T->includeOperand(t[i]->deepCopy());
	return T;
}

std::vector<AST*> trueFactors(AST* u, std::vector<AST*> l, AST* x, int p, int k) {

	AST* U = u->deepCopy();
	std::vector<AST*> L = copyList(l);

	std::vector<AST*> factors = std::vector<AST*>(0);

	int m = 1;

	while(m < L.size()/2) {

		std::vector<std::vector<AST*>> C = comb(L, m); 

		while(C.size() != 0) {
			std::vector<AST*> t = C[0];

			AST* T = build_product(t);
			T = Ts(T, x, (int)std::pow(p, k));
	
			std::pair<AST*, AST*> D = divideGPE(U,T,x);
			AST* Q = D.first;	
			AST* R = D.second;	
			if(R->kind() == Kind::Integer && R->value() == 0) {
				factors.push_back(T->deepCopy());
				U = Q;

				L = list_diff(L, t); // L ~ t
				C = cleanUp(C, t);
			} else {
				C = set_difference(C, {t}); // C ~ {t}; // set difference
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
