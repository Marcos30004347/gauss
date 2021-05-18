#include "Zp.hpp"
#include "Factorization.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Simplification/Simplification.hpp"

#include <cmath>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace prime;
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

AST* hanselLift(AST* u, AST* S, AST* x, int p, int k) {
	if(k==1) {
		return trueFactors(u, S, x, p, k);
	}

	AST* V = S->operandList();
	AST* G = genExtendSigmaP(V, x, p);

	for(int j=2; j<=k; j++) {
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

		for(int i=0; i<V->numberOfOperands(); i++) {
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

int _mod(int a, int b) {
	return (b + (a%b)) % b;
}

void RMatrix(AST* u, AST* x, AST* n_, int p) {
	int n = n_->value();

	if(R != nullptr) {
		destroyRMatrix(n);
	}
	
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
					_mod(c->value(), p)
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
				R[i][j] = _mod(coeff->value() - 1,p);
			} else {
				R[i][j] = _mod(coeff->value(),p);
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
						integer(_mod(c->value(), p))
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
						integer(_mod(coeff->value(), p)),
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
		// printf("%s\n", yk->toString().c_str());
	}

	delete yk;
	delete n_min_one;
	// delete v_;
	delete v;

	// for(int i=0; i<n; i++) {
	// 	for(int j=0; j<n; j++) {
	// 		printf("%i ", R[i][j]);
	// 	}
	// 	printf("\n");
	// }
}

AST* auxiliaryBasis(AST* x, AST* n, int p) {
	int P[n->value()];
	for(int i=0; i<n->value(); i++) {
		P[i] = 0;
	}
	
	AST* S = list({});
	for(int j=0; j<n->value(); j++) {
		int i = 0;
		bool pivot_found = false;
		while(!pivot_found && i < n->value()) {
			if(R[i][j] != 0 && P[i] == 0) {
				pivot_found = true;
			} else {
				i = i+1;
			}
		}
		if(pivot_found) {
			P[i] = j;
			int a = modInverse_p(R[i][j], p);
			for(int l=0; l<n->value(); l++) {
				R[i][l] = remainder(a*R[i][l], p);
			}
			for(int k=0; k<n->value(); k++) {
				if(k!=i) {
					int f = R[k][j];
					for(int l=0; l<n->value(); l++) {
						R[k][l] = remainder(R[k][l] - f * R[i][l], p);
					}
				}
			}
		} else if(!pivot_found) {
			AST* s = power(
				x->deepCopy(),
				sub({ integer(j), integer(1) })
			);
			for(int l=0; l<j-1; l++) {
				int e = 0;
				int i = 1;
				while(e == 0 && i< n->value()) {
					if(l == P[i]) {
						e = i;
					} else {
						i = i+1;
					}
				}
				if(e>0) {
					int c = remainder(-1*R[e][j], p);
					s = add({ s, mul({ integer(c), power(x->deepCopy(), sub({ integer(l), integer(1) })) }) });
				}
			}
			AST* L = list({ s->deepCopy() });
			AST* S_ = join(S, L);
			delete S;
			delete L;
			S = S_;
		}
	}

	return S;
}

AST* findFactors(AST* u, AST* S, AST* x, int p) {
	signed long r = S->numberOfOperands();
	
	AST* factors = set({u->deepCopy()});
	
	for(int k=2; k<=r; k++) {
		
		AST* b = S->operand(k);
		AST* old_factors = factors;
		
		for(int i=0; i<old_factors->numberOfOperands(); i++) {
			AST* w = old_factors->operand(i);
			int j = 0;
			while(j <= p-1) {
				AST* b__ = sub({
					b->deepCopy(),
					integer(j)
				});
				AST* b_ = reduceAST(b__);
				delete b__;
				
				AST* g = gcdGPE_Zp(b_, w, x, p);

				delete b_;

				if(g->kind() == Kind::Integer && g->value()==1) {
					j = j+1;
				} else if(g->match(w)) {
					j = p;
				} else {
					AST* W = set({w->deepCopy()});
					AST* factors_ = difference(factors, W);
					delete factors;
					delete W;
					factors = factors_;

					AST* q__ = divideGPE_Zp(w, g, x, p);
					AST* q = q__->operand(0)->deepCopy();
					delete q__;

					W = set({g->deepCopy(), q->deepCopy()});
					AST* factors__ = unification(factors, W);
					delete W;
					delete factors;
					factors = factors__;

					if(factors->numberOfOperands() == r) {
						return factors;
					} else {
						j = j+1;

						delete w;
						w = q;
					}
				}
			}
		}
	}

	return factors;
}

AST* berlekampFactor(AST* u, AST* x, int p) {
	AST* n = degreeGPE(u, x);
	if(
		n->kind() == Kind::Integer && n->value() == 0 ||
		n->kind() == Kind::Integer && n->value() == 1
	) {
		return set({u->deepCopy()});
	}

	RMatrix(u, x, n, p);
	AST* S = auxiliaryBasis(x, n, p);
	if(S->numberOfOperands() == 1) {
		return set({u->deepCopy()});
	}

	return findFactors(u, S, x, p);
}

AST* irreducibleFactor(AST* u, AST* x, AST* y) {
	AST* n = degreeGPE(u, x);
	AST* l = leadingCoefficientGPE(u, x);

	AST* l_ = mul({
		power(
			l->deepCopy(),
			sub({n->deepCopy(), integer(1)})
		),
		u->deepCopy()
	});

	AST* x_ = div(y->deepCopy(), l->deepCopy());
	AST* V_ = deepReplace(l_, x, x_);
	AST* V = algebraicExpand(V_);

	int p = findPrime(V, y);
	AST* V_tnn = Tnn(V, y, p);

	AST* S = berlekampFactor(V_tnn, y, p);

	if(S->numberOfOperands() == 1) {
		return u->deepCopy();
	}

	unsigned long k = findK(V, y, p);

	AST* k_ = mapAST(Ts, S, y, p);
	AST* W = hanselLift(V, k_, y, p, k);
	AST* t = mul({ l->deepCopy(), x->deepCopy() });
	AST* W_ = deepReplace(W, y, t);
	delete W;
	W = W_;

	AST* M = integer(1);
	for(int i=0; i<W->numberOfOperands(); i++) {
		AST* w = W->operand(i);
		AST* L = list({});
		AST* Z = symbol("Z");
		AST* z_ = div(
			w->deepCopy(),
			polynomialContent(w, x, L, Z)
		);

		AST* z = reduceAST(z_);

		AST* M = mul({ M, z });
	}

	return reduceAST(M);
}

}
