#include "Zp.hpp"
#include "Factorization.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Simplification/Simplification.hpp"

#include <cmath>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace rational;
using namespace simplification;

namespace polynomial {

int** R = nullptr;

int getRMatrixValue(int i, int j) {
	return R[i][j];
}

void destroyRMatrix(int n) {
	for(int i=0; i < n; i++) {
		delete []R[i];
	}

	delete []R;

	R = nullptr;
}

AST* genExtendSigmaP(AST* V, AST* x, unsigned p) {
	if(V->numberOfOperands() == 2) {
		AST* k = extendedEuclideanAlgGPE_Sp(V->operand(1), V->operand(0), x, p);

		AST* A = k->operand(1)->deepCopy();
		AST* B = k->operand(2)->deepCopy();

		delete k;

		return list({ A, B });
	}

	AST* V_ = new AST(Kind::List);

	for(unsigned int i=0; i<V->numberOfOperands() - 1; i++)
		V_->includeOperand(V->operand(i)->deepCopy());

	AST* tal = genExtendSigmaP(V_, x, p);

	AST* gs_ = construct(Kind::Multiplication, V_);
	AST* gs = Ts(gs_, x, p);

	AST* vs = V->operand(V->numberOfOperands() - 1);

	AST* k = extendedEuclideanAlgGPE_Sp(vs, gs, x, p);

	AST* A = k->operand(1);
	AST* B = k->operand(2)->deepCopy();

	AST* thetas = new AST(Kind::List);

	for(unsigned int i=0; i<tal->numberOfOperands(); i++) {
		AST* thetha = mul({ A->deepCopy(), tal->operand(i)->deepCopy() });
		
		thetas->includeOperand(Ts(thetha, x, p));
		
		delete thetha;
	}

	thetas->includeOperand(B);

	delete tal;
	delete k;
	delete gs;
	delete gs_;
	delete V_;

	return thetas;
}

AST* genExtendRP(AST* V, AST* S, AST* F, AST* x, unsigned p) {
	AST* Rs = new AST(Kind::List);
	for(unsigned int i=0; i<V->numberOfOperands(); i++) {
		AST* t = mul({ F->deepCopy(), S->operand(i)->deepCopy() });
		AST* u = Ts(t, x, p);

		Rs->includeOperand(
			remainderGPE_Sp(u, V->operand(i), x, p)
		);
	
		delete t;
		delete u;
	}

	return Rs;
}

// u is a polynomial in x, findPrime return
// and integer such tath p % leadCoeff(u,x) != 0
int findPrime(AST* u, AST* x) {
	AST* lc_ = leadingCoefficientGPE(u, x);
	AST* lc  = expandAST(lc_);

	int p = primes[0];

	// We are gonna look just for the first 5000000 primes
	for(unsigned int i=0; i < primes.count(); i++) {
		if(lc->value() % primes[i] != 0) {
			p = primes[i];
			break;
		}
	}
	
	delete lc_;
	delete lc;

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
		
		delete d;
		delete c;
	}
	
	delete u_;
	delete d_;

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

	unsigned int m = 1;

	while(m < L->numberOfOperands()/2) {
		AST* m_ = integer(m);
		AST* C = combination(L, m_); 
		delete m_;

		while(C->numberOfOperands() != 0) {
			AST* t = C->operand(0);

			AST* T_ = new AST(Kind::Multiplication);
			for(unsigned int i=0; i<t->numberOfOperands(); i++)
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

AST* henselLift(AST* u, AST* S, AST* x, int p, int k) {
	
	if(k==1)
	{
		return trueFactors(u, S, x, p, k);
	}

	AST* V = S->operandList();
	AST* G = genExtendSigmaP(V, x, p);

	for(int j=2; j<=k; j++) 
	{
	
		AST* Vp_ = construct(Kind::Multiplication, V);
		AST* Vp = algebraicExpand(Vp_);
	
		AST* E_ = sub({ u->deepCopy(), Vp->deepCopy() });
		AST* E = reduceAST(E_);

		if(E->kind() == Kind::Integer && E->value() == 0) {
			AST* r = construct(Kind::Set, V);	
			return r;
		}
	
		AST* E_TS = Ts(E, x, pow(p, j));
	
		AST* p_ = power(integer(p), sub({integer(j), integer(1)}));
		AST* f_ = div(E_TS, p_);
		AST* f = reduceAST(f_);

		AST* F = algebraicExpand(f);

		AST* R = genExtendRP(V, G, F, x, p);

		AST* Vnew = list({});

		for(unsigned int i=0; i<V->numberOfOperands(); i++) {
			AST* v_lift_ = add({
				V->operand(i)->deepCopy(),
				mul({
					power(integer(p), sub({
						integer(j),
						integer(1)
					})),
					R->operand(i)->deepCopy()
				})
			});

			AST* v_lift = algebraicExpand(v_lift_);

			AST* V__ = list({ v_lift });

			AST* Vnew_ = join(Vnew, V__);

			delete Vnew;
			Vnew = Vnew_;
		}

		V = Vnew;
	}

	AST* K = construct(Kind::Set, V);
	AST* r = trueFactors(u, K, x, p, k);

	return r;
}


void RMatrix(AST* u, AST* x, AST* n_, int p) {
	int n = n_->value();

	// if(R != nullptr) {
	// 	destroyRMatrix(n);
	// }
	
	R = new int*[n];
	for(int i=0; i<n; i++)
		R[i] = new int[n];

	AST* yk = integer(1);

	AST* n_min_one = integer(n-1);

	AST* v_ = integer(0);

	for(int i=0; i<n; i++) {
		AST* e = integer(i);
		AST* c = coefficientGPE(u, x, e);
		
		v_ = add({
			v_,
			mul({
				integer(
					mod(c->value(), p)
				),
				power(
					x->deepCopy(),
					e
				)
			})
		});
		
		delete c;
	}

	AST* v = reduceAST(v_);
	delete v_;

	for(int j=0; j<n; j++) {
		for(int i=0; i<n; i++) {
			AST* e = integer(i);
			AST* coeff_ = coefficientGPE(yk, x, e);
			AST* coeff = reduceAST(coeff_);
			
			if(i == j) {
				R[i][j] = mod(coeff->value() - 1,p);
			} else {
				R[i][j] = mod(coeff->value(),p);
			}
			
			delete e;
			delete coeff;
			delete coeff_;
		}
	
		if(j == n - 1) break;		
	
		for(int i = p*(j+1); i < p*(j+2); i++) {
			
			AST* ck = coefficientGPE(yk, x, n_min_one);
			AST* zk_ = integer(0);
			
			for(int i=n-2; i>=0; i--) {
				AST* e = integer(i);

				AST* c = coefficientGPE(yk, x, e);

				zk_ = add({
					zk_,
					mul({
						power(
							x->deepCopy(),
							e->deepCopy()
						),
						integer(mod(c->value(), p))
					})
				});
				delete e;
				delete c;
			}

			AST* zk = reduceAST(zk_);
			delete zk_;

			AST* yk_ = add({
				mul({ integer(-1), ck, v->deepCopy() }),
				mul({ x->deepCopy(), zk })
			});

			delete yk;
			yk = algebraicExpand(yk_);
			delete yk_;
			
			// project yk into Zp
			AST* yk_p = integer(0);
			AST* deg = degreeGPE(yk, x);

			for(int s=deg->value(); s>=0; s--) {
				AST* ex = integer(s);
				AST* coeff = coefficientGPE(yk, x, ex);
				
				yk_p = add({
					yk_p,
					mul({
						integer(mod(coeff->value(), p)),
						power(x->deepCopy(), ex)
					})
				});
				
				delete coeff;
			}
		
			delete yk;
			delete deg;
	
			yk = reduceAST(yk_p);
	
			delete yk_p;
		}
	}

	delete yk;
	delete n_min_one;
	delete v;
}

AST* auxiliaryBasis(AST* x, AST* n, int p) {
	int P[n->value()];
	
	for(int i=1; i<=n->value(); i++) {
		P[i-1] = 0;
	}

	AST* S = list({});

	for(int j=1; j<=n->value(); j++) {

		int i = 1;
		bool pivot_found = false;
		
		while(!pivot_found && i < n->value()) {
			if(R[i-1][j-1] != 0 && P[i-1] == 0) {
				pivot_found = true;
			} else {
				i = i+1;
			}
		}

		if(pivot_found) {
			P[i-1] = j;

			int a = modInverse_p(R[i-1][j-1], p);

			for(int l=1; l<=n->value(); l++) {
				R[i-1][l-1] = mod(a * R[i-1][l-1], p);
			}

			for(int k=1; k <= n->value(); k++) {
				if(k!=i) {
					int f = R[k-1][j-1];
					for(int l=1; l <= n->value(); l++) {
						R[k-1][l-1] = mod(R[k-1][l-1] - f * R[i-1][l-1], p);
					}
				}
			}

		} else if(!pivot_found) {
		
			AST* s = power(
				x->deepCopy(),
				sub({ integer(j), integer(1) })
			);
		
			for(int l=1; l <= j-1; l++) {
	
				int e = 0;
				int i = 1;
	
				while(e == 0 && i< n->value()) {
					if(l == P[i-1]) {
						e = i;
					} else {
						i = i+1;
					}
				}
				if(e > 0) {
					int c = mod(-1*R[e-1][j-1], p);
					s = add({ s, mul({ integer(c), power(x->deepCopy(), sub({ integer(l), integer(1) })) }) });
				}
			}
	
			AST* L = list({ algebraicExpand(s) });
			AST* S_ = join(S, L);
	
			delete s;
			delete S;
			delete L;
	
			S = S_;
		}
	}

	return S;
}

AST* findFactors(AST* u, AST* S, AST* x, int p) {
	signed long r = S->numberOfOperands();
	
	AST* factors = set({ u->deepCopy() });

	for(int k=2; k <= r; k++) {
		
		AST* b = S->operand(k - 1)->deepCopy();
	
		AST* old_factors = factors->deepCopy();
		
		for(unsigned int i = 0; i < old_factors->numberOfOperands(); i++) {
			AST* w = old_factors->operand(i)->deepCopy();
		
			int j = 0;
		
			while(j <= p - 1) {

				AST* b__ = add({
					b->deepCopy(),
					integer(mod(-1*j,p))
				});

				AST* b_ = reduceAST(b__);
	
				delete b__;	
			
				AST* g = gcdGPE_Zp(b_, w, x, p);
	
				delete b_;

				if(g->kind() == Kind::Integer && g->value() == 1) {
					j = j+1;
				} else if(g->match(w)) {
					j = p;
				} else {
					AST* factors_;
					AST* S0 = set({ w->deepCopy() });
					factors_ = difference(factors, S0);
					delete S0;
					delete factors;
					factors = factors_;

					AST* z = divideGPE_Zp(w, g, x, p);

					AST* q = z->operand(0)->deepCopy();

					delete z;

					AST* S1 = set({g->deepCopy(), q->deepCopy()});
					factors_ = unification(factors, S1);
					delete S1;
					delete factors;
					factors = factors_;
					

					if(factors->numberOfOperands() == r) {
						
						delete w;
						delete g;
						delete q;
						delete b;
						delete old_factors;
						
						return factors;
					} else {
						j = j + 1;
						
						delete w;
						
						w = q->deepCopy();
					}
				
					delete q;
				}

				delete g;
			}
			delete w;
		}
	
		delete b;
		delete old_factors;
	}

	return factors;
}

AST* berlekampFactor(AST* u, AST* x, int p) {
	AST* n = degreeGPE(u, x);
	if(
		(n->kind() == Kind::Integer && n->value() == 0) ||
		(n->kind() == Kind::Integer && n->value() == 1)
	) {

		delete n;

		return set({ u->deepCopy() });
	}

	RMatrix(u, x, n, p);
	
	AST* S = auxiliaryBasis(x, n, p);

	if(S->numberOfOperands() == 1) {
		
		delete S;
		delete n;
		
		destroyRMatrix(n->value());
		
		return set({u->deepCopy()});
	}

	AST* factors = findFactors(u, S, x, p);

	destroyRMatrix(n->value());

	delete S;
	delete n;

	return factors;
}

AST* irreducibleFactor(AST* u, AST* x, AST* y) {
	AST* n = degreeGPE(u, x);
	AST* l = leadingCoefficientGPE(u, x);

	AST* l_ = mul({
		power(l->deepCopy(), sub({n->deepCopy(), integer(1)})),
		u->deepCopy()
	});

	AST* x_ = div(y->deepCopy(), l->deepCopy());
	AST* V_ = deepReplace(l_, x, x_);
	AST* V = algebraicExpand(V_);

	int p = findPrime(V, y);

	AST* V_tnn = Tnn(V, y, p);

	AST* S = berlekampFactor(V_tnn, y, p);

	if(S->numberOfOperands() == 1) 
	{
		return u->deepCopy();
	}

	unsigned long k = findK(V, y, p);

	AST* k_ = mapAST(Ts, S, y, p);
	AST* W = henselLift(V, k_, y, p, k);
	AST* t = mul({ l->deepCopy(), x->deepCopy() });
	AST* W_ = deepReplace(W, y, t);

	delete W;
	W = W_;

	AST* M = integer(1);

	for(unsigned int i=0; i<W->numberOfOperands(); i++) {
		AST* w = W->operand(i);
		AST* L = list({});
		AST* Z = symbol("Z");

		AST* z_ = div(
			w->deepCopy(),
			polynomialContent(w, x, L, Z)
		);

		// AST* z = algebraicExpand(z_);
		AST* z = reduceAST(z_);

		M = mul({ M, z });
	}

	// return reduceAST(M);
	return M;
}


std::pair<AST*, AST*> getPolynomialInZ(AST* u, AST* x)
{
	AST* M = denominator(u->operand(0));

	for(int i=1; i<u->numberOfOperands(); i++)
	{
		AST* b = denominator(u->operand(i));

		M = leastCommomMultiple(M, b);

		delete b;
	}

	AST* v_ = mul({M->deepCopy(), u->deepCopy()});
	AST* v  = algebraicExpand(v_);

	delete v_;

	return {v, M};
}

}
