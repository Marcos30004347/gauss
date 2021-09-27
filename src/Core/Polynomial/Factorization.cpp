#include "Zp.hpp"
#include "Hensel.hpp"
#include "Factorization.hpp"
#include "Resultant.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Calculus/Calculus.hpp"

#include <cmath>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace rational;
using namespace calculus;
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
		AST* k = extendedEuclideanAlgGPE_sZp(V->operand(1), V->operand(0), x, p);

		AST* A = k->operand(1)->copy();
		AST* B = k->operand(2)->copy();

		delete k;

		return list({ A, B });
	}

	AST* V_ = new AST(Kind::List);

	for(unsigned int i=0; i<V->numberOfOperands() - 1; i++)
		V_->includeOperand(V->operand(i)->copy());

	AST* tal = genExtendSigmaP(V_, x, p);

	AST* gs_ = construct(Kind::Multiplication, V_);
	AST* gs = sZp(gs_, x, p);

	AST* vs = V->operand(V->numberOfOperands() - 1);

	AST* k = extendedEuclideanAlgGPE_sZp(vs, gs, x, p);

	AST* A = k->operand(1);
	AST* B = k->operand(2)->copy();

	AST* thetas = new AST(Kind::List);

	for(unsigned int i=0; i<tal->numberOfOperands(); i++) {
		AST* thetha = mul({ A->copy(), tal->operand(i)->copy() });

		thetas->includeOperand(sZp(thetha, x, p));

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
		AST* t = mul({ F->copy(), S->operand(i)->copy() });
		AST* u = sZp(t, x, p);

		Rs->includeOperand(
			remainderGPE_sZp(u, V->operand(i), x, p)
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

	for(int i=d; i>=0; i--) 
	{
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

AST* trueFactors(AST* u, AST* l, AST* x, int p, int k) 
{
	AST* U = u->copy();
	AST* L = l->copy();

	AST* factors = list({});

	unsigned int m = 1;

	while(m <= L->numberOfOperands()/2) 
	{
		AST* w = integer(m);

		AST* C = combination(L, w);

		delete w;

		while(C->numberOfOperands()) 
		{
			AST* t = C->operand(0);

			AST* Z = new AST(Kind::Multiplication);
			for(unsigned int i=0; i<t->numberOfOperands(); i++)
			{
				Z->includeOperand(t->operand(i)->copy());
			}

			AST* I = algebraicExpand(Z);

			AST* T = sZp(I, x, (long)std::pow(p, k));

			delete Z;
			delete I;

			AST* D = divideGPE(U, T, x);

			AST* Q = D->operand(0);
			AST* R = D->operand(1);

			if(R->is(0)) 
			{
				factors->includeOperand(T->copy());
		
				U = Q->copy();
	
				AST* R = difference(L, t);//  // L ~ t
				
				delete L;
				
				L = R->copy();
	
				AST* E = cleanUp(C, t);

				delete C;

				C = E;
			} 
			else 
			{
				AST* J = new AST(Kind::Set, { t->copy() });

				AST* T = difference(C, J); // C ~ {t};
		
				delete J;

				delete C;
				C = T;
			}

			delete D;
		}

		delete C;
	
		m = m + 1;
	}

	if(U->isNot(1))
	{
		factors->includeOperand(U->copy());
	}

	delete U;
	delete L;

	return factors;
}

AST* henselLift(AST* u, AST* S, AST* x, int p, int k) {

	if(k==1)
	{
		return trueFactors(u, S, x, p, k);
	}

	AST* V = S->operandList();
	AST* G = genExtendSigmaP(V, x, p);

	for(int j = 2; j<=k; j++)
	{

		AST* Vp_ = construct(Kind::Multiplication, V);
		AST* Vp = algebraicExpand(Vp_);

		AST* E_ = sub({ u->copy(), Vp->copy() });
		AST* E = reduceAST(E_);

		if(E->kind() == Kind::Integer && E->value() == 0) {
			AST* r = construct(Kind::Set, V);
			return r;
		}

		AST* E_TS = sZp(E, x, pow(p, j));

		AST* p_ = power(integer(p), sub({integer(j), integer(1)}));
		AST* f_ = div(E_TS, p_);
		AST* f = reduceAST(f_);

		AST* F = algebraicExpand(f);

		AST* R = genExtendRP(V, G, F, x, p);

		AST* Vnew = list({});

		for(unsigned int i=0; i<V->numberOfOperands(); i++) {
			AST* v_lift_ = add({
				V->operand(i)->copy(),
				mul({
					power(integer(p), sub({
						integer(j),
						integer(1)
					})),
					R->operand(i)->copy()
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

	printf("lift %s\n", K->toString().c_str());

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
					x->copy(),
					e
				)
			})
		});

		delete c;
	}

	AST* v = reduceAST(v_);
	delete v_;

	for(int j=0; j<n; j++)
	{
		for(int i=0; i<n; i++)
		{
			AST* e = integer(i);
			AST* coeff_ = coefficientGPE(yk, x, e);
			AST* coeff = reduceAST(coeff_);

			if(i == j)
			{
				R[i][j] = mod(coeff->value() - 1, p);
			} else
			{
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
							x->copy(),
							e->copy()
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
				mul({ integer(-1), ck, v->copy() }),
				mul({ x->copy(), zk })
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
						power(x->copy(), ex)
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
				x->copy(),
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
					s = add({ s, mul({ integer(c), power(x->copy(), sub({ integer(l), integer(1) })) }) });
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

	AST* factors = set({ u->copy() });

	for(int k=2; k <= r; k++) {

		AST* b = S->operand(k - 1)->copy();

		AST* old_factors = factors->copy();

		for(unsigned int i = 0; i < old_factors->numberOfOperands(); i++) {
			AST* w = old_factors->operand(i)->copy();

			int j = 0;

			while(j <= p - 1) {

				AST* b__ = add({
					b->copy(),
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
					AST* S0 = set({ w->copy() });
					factors_ = difference(factors, S0);
					delete S0;
					delete factors;
					factors = factors_;

					AST* z = divideGPE_Zp(w, g, x, p);

					AST* q = z->operand(0)->copy();

					delete z;

					AST* S1 = set({g->copy(), q->copy()});
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

						w = q->copy();
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

		return set({ u->copy() });
	}

	RMatrix(u, x, n, p);

	AST* S = auxiliaryBasis(x, n, p);

	if(S->numberOfOperands() == 1) {

		delete S;
		delete n;

		destroyRMatrix(n->value());

		return set({u->copy()});
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
		power(l->copy(), sub({ n->copy(), integer(1) })),
		u->copy()
	});

	AST* x_ = div(y->copy(), l->copy());
	AST* V_ = deepReplace(l_, x, x_);
	AST* V = algebraicExpand(V_);

	int p = findPrime(V, y);

	AST* V_tnn = Zp(V, y, p);

	AST* S = berlekampFactor(V_tnn, y, p);

	if(S->numberOfOperands() == 1)
	{
		return u->copy();
	}

	unsigned long k = findK(V, y, p);

	AST* k_ = mapAST(sZp, S, y, p);

	AST* W = henselLift(V, k_, y, p, k);
	AST* t = mul({ l->copy(), x->copy() });
	AST* W_ = deepReplace(W, y, t);

	delete W;
	W = W_;

	AST* M = integer(1);

	for(unsigned int i=0; i<W->numberOfOperands(); i++) 
	{
		AST* w = W->operand(i);
		AST* L = list({});
		AST* Z = symbol("Z");

		AST* z_ = div(w->copy(), polynomialContent(w, x, L, Z));

		AST* z = algebraicExpand(z_);
		// AST* z = reduceAST(z_);

		M = mul({ M, z });
	}

	// return reduceAST(M);
	return M;
}


std::pair<AST*, AST*> getPolynomialInZ(AST* u, AST* x)
{
	AST* j = integer(0), *b;
	AST* c = coefficientGPE(u->operand(0), x, j);
	AST* M = denominator(c);

	delete c;
	delete j;
	
	for(unsigned int i = 1; i<u->numberOfOperands(); i++)
	{

		j = degreeGPE(u->operand(i), x);
		c = coefficientGPE(u->operand(i), x, j);

		b = denominator(c);

		AST* m = leastCommomMultiple(M, b);

		delete M;
		M = m;
	
		delete c;
		delete b;
		delete j;

	}

	AST* k = mul({M->copy(), u->copy()});
	AST* v  = algebraicExpand(k);

	delete k;

	return { v, M };
}




ast::AST* squareFreeFactor(ast::AST* u, ast::AST* x)
{
	if(isEqZero(u))
	{
		return u->copy();
	}

	AST* c = leadingCoefficientGPE(u, x);
	AST* u_ = div(u, c);
	AST* U = algebraicExpand(u_);

	AST* P = integer(1);

	AST* Udx = derivate(U, x);

	AST* R = gcdGPE(U, Udx, x);
	AST* F = quotientGPE(U, R, x);
	AST* j = integer(1);

	while(R->kind() != Kind::Integer || (R->kind() == Kind::Integer && R->value() != 1))
	{
		AST* G = gcdGPE(R, F, x);
		AST* s = quotientGPE(F, G, x);

		P = mul({P, power(s, j)});

		R = quotientGPE(R, G, x);

		delete F;
		F = G;

		j = integer(j->value() + 1);
	}

	P = mul({P, power(F, j)});

	AST* ret = mul({c, P});

	return ret;
}


AST* squareFreeFactorization(AST* ax, AST* x)
{
	unsigned int i = 1;

	AST* out = integer(1);
	AST* bx = derivate(ax, x);
	AST* cx = gcdGPE(ax, bx, x);
	AST* wx = quotientGPE(ax, cx, x);

	while(
		cx->kind() != Kind::Integer ||
		(cx->kind() == Kind::Integer && cx->value() != 1)
	)
	{
		AST* yx = gcdGPE(wx, cx, x);
		AST* zx = quotientGPE(wx, yx, x);

		out = mul({ out, power(zx, integer(i)) });

		i = i + 1;

		delete wx;
		wx = yx;

		AST* qx = quotientGPE(cx, yx, x);

		delete cx;
		cx = qx;
	}

	delete bx;
	delete cx;

	out = mul({ out , power(wx, integer(i)) });

	AST* t = reduceAST(out);

	delete out;

	return t;
}

AST* squareFreeFactorization2(AST* ax, AST* x)
{
	unsigned int i = 1;
	AST* out = integer(1);

	AST* bx = derivate(ax, x);
	AST* cx = gcdGPE(ax, bx, x);
	AST* wx = nullptr;

	if(cx->is(1))
	{
		wx = ax->copy();
	}
	else
	{
		wx = quotientGPE(ax, cx, x);

		AST* yx = quotientGPE(bx, cx, x);
		AST* kx = sub({ yx, derivate(wx, x) });
		AST* zx = reduceAST(kx);

		delete kx;

		while(zx->isNot(0))
		{
			AST* gx = gcdGPE(wx, zx, x);
			out = mul({ out, power(gx->copy(), integer(i))});

			i = i + 1;

			AST* tx = quotientGPE(wx, gx, x);

			delete wx;
			wx = tx;

			yx = quotientGPE(zx, gx, x);

			AST* rx = sub({ yx, derivate(wx, x) });

			delete zx;
			zx = reduceAST(rx);

			delete rx;
			delete gx;
		}

		delete zx;
	}

	out = mul({out, power(wx, integer(i))});

	AST* tx = reduceAST(out);

	delete cx;
	delete bx;
	delete out;

	return tx;
}

AST* squareFreeFactorizationFiniteField(AST* ax, AST* x, AST* q)
{
	assert(
		q->kind() == Kind::Power,
		"p is not a order of a Galois Field, should be q = p^m"
	);

	AST* p = q->operand(0);

	unsigned int i = 1;

	AST* out = integer(1);
	AST* ux = derivate(ax, x);
	AST* bx = Zp(ux, x, p->value());

	delete ux;

	if(bx->isNot(0))
	{
		AST* cx = gcdGPE_Zp(ax, bx, x, p->value());
		AST* wx = quotientGPE_Zp(ax, cx, x, p->value());

		while(wx->isNot(1))
		{
			AST* yx = gcdGPE_Zp(wx, cx, x, p->value());
			AST* zx = quotientGPE_Zp(wx, yx, x, p->value());

			out = mul({ out, power(zx, integer(i))});

			i = i + 1;

			delete wx;
			wx = yx;

			AST* kx = quotientGPE_Zp(cx, yx, x, p->value());

			delete cx;
			cx = kx;
		}

		if(cx->isNot(1))
		{
			AST* kx = add({});
			AST* deg = degreeGPE(cx, x);

			for(unsigned int i = 0; i <= deg->value(); i++)
			{
				AST* j = integer(i);

				kx->includeOperand(mul({
					coefficientGPE(cx, x, j),
					power(x->copy(), integer(i/p->value()))
				}));

				delete j;
			}

			delete cx;
			cx = reduceAST(kx);

			delete deg;
			delete kx;

			AST* sx = squareFreeFactorizationFiniteField(cx, x, q);

			delete cx;
			cx = sx;

			out = mul({ out, power(cx->copy(), integer(p->value())) });
		}

		delete cx;
		delete wx;
	}
	else
	{
		AST* deg = degreeGPE(ax, x);
		AST* kx = add({});

		for(unsigned int i = 0; i <= deg->value(); i++)
		{
			AST* j = integer(i);

			kx->includeOperand(mul({
				coefficientGPE(ax, x, j),
				power(x->copy(), integer(i/p->value()))
			}));

			delete j;
		}

		delete deg;

		delete ax;
		ax = kx;

		AST* sx = squareFreeFactorizationFiniteField(ax, x, q);

		delete out;
		out = power(sx, integer(p->value()));
	}

	AST* tx = reduceAST(out);

	delete out;
	delete bx;

	return tx;
}

AST* formQRow(AST* ax, AST* x, unsigned int q, unsigned int n, AST* r)
{
	// TODO: Maybe there is an easy way to form the r vector
	// without computing the remainder every time

	AST* e = integer(0);

	AST* r0 = mul({
		integer(-1),
		r->operand(n - 1)->copy(),
		coefficientGPE(ax, x, e),
	});

	delete e;

	AST* ux = r0;

	for(unsigned int i = 1; i < n; i++)
	{
		AST* e = integer(i);

		AST* ri = sub({
			r->operand(i - 1)->copy(),
			mul({
				r->operand(n - 1)->copy(),
				coefficientGPE(ax, x, e)
			})
		});

		delete e;

		ux = add({ ux, mul({ri, power(x->copy(), integer(i))}) });
	}

	AST* kx = reduceAST(ux);
	delete ux;

	ux = remainderGPE_sZp(kx, ax, x, q);

	delete kx;

	AST* l = list({});

	for(unsigned int i=0; i < n; i++)
	{
		AST* e = integer(i);

		AST* ri = coefficientGPE(ux, x, e);

		AST* a = l;
		AST* b = list({ ri });

		l = append(a, b);

		delete e;
		delete a;
		delete b;
	}

	delete ux;

	return l;
}





AST* initialQRow(AST* n)
{
	AST* r = list({integer(1)});
	AST* z = list({integer(0)});

	for(unsigned int i = 1; i < n->value(); i++)
	{
		AST* k = append(r, z);
		delete r;
		r = k;
	}

	delete z;

	return r;
}

AST* formMatrixQ(AST* ax, AST* x, AST* q)
{
	assert(
		q->kind() == Kind::Integer,
		"q needs to be an integer"
	);

	AST* n = degreeGPE(ax, x);

	assert(
		n->kind() == Kind::Integer,
		"degree of the polynomial ax needs to be an integer"
	);

	unsigned int p = q->value();
	unsigned int e = n->value();

	AST* Q = matrix(n, n);

	AST* r = initialQRow(n);

	// Set Q row
	for(unsigned int i = 0; i < e; i++)
	{
		AST* ri = r->operand(i)->copy();

		Q->operand(0)->deleteOperand(i);
		Q->operand(0)->includeOperand(ri, i);
	}

	for(unsigned int m = 1; m <= (e - 1)*p; m++)
	{
		AST* j = formQRow(ax, x, p, e, r);

		delete r;
		r = j;

		if(m % p == 0)
		{
			// Set Q row
			for(unsigned int i = 0; i < e; i++)
			{
				AST* ri = r->operand(i)->copy();

				Q->operand(m / p)->deleteOperand(i);
				Q->operand(m / p)->includeOperand(ri, i);
			}
		}
	}

	delete n;
	delete r;

	return Q;
}



AST* polynomialMultiplication(AST* ax, AST* bx, AST* x)
{
	// TODO: override this with algebraic expand when it gets optimized

	AST* ux = add({});

	ax = reduceAST(ax);
	bx = reduceAST(bx);

	AST* da = degreeGPE(ax, x);
	AST* db = degreeGPE(bx, x);

	for(unsigned int i = 0; i <= da->value(); i++)
	{
		for(unsigned int j = 0; j <= db->value(); j++)
		{
			AST* ae = integer(i);
			AST* be = integer(j);

			AST* ca = coefficientGPE(ax, x, ae);
			AST* cb = coefficientGPE(bx, x, be);

			ux->includeOperand(
				mul({
					mul({ca, cb}),
					power(
						x->copy(),
						add({ ae, be })
					)
				})
			);
		}
	}

	AST* px = reduceAST(ux);

	delete da;
	delete db;
	
	delete ux;
	delete ax;
	delete bx;

	return px;
}

AST* formQRowBinaryExp(AST* ax, AST* x, signed long q, signed long n, signed long i, signed long** cache)
{
	// TODO: Maybe there is an easy way to form the r vector
	// without computing the remainder every time
	signed long t = i % 2;
	signed long j = (i - t)/2;


	AST* ux = nullptr;
	AST* px = nullptr;

	if(j >= 2*n)
	{
		AST* r = formQRowBinaryExp(ax, x, q, n, j, cache);

		ux = integer(0);

		for(signed long k=0; k<n; k++)
		{
			signed long int rk = r->operand(k)->value();

			ux = add({
				ux,
				mul({
					integer(rk),
					power(x->copy(), integer(k))
				})
			});
		}

		delete r;

		// ux = mul({ ux->copy(), ux });
		px = polynomialMultiplication(ux, ux, x);
		ux = remainderGPE_sZp(px, ax, x, q);
		
		delete px;
	
	}
	else
	{
		ux = integer(0);

		for(unsigned int k=0; k<n; k++)
		{
			signed long int rji = cache[j - 1][k];
			ux = add({
				ux,
				mul({
					integer(rji),
					power(x->copy(), integer(k))
				})
			});
		}

		AST* lx = polynomialMultiplication(ux, ux, x);
		delete ux;
	
		ux = remainderGPE_sZp(lx, ax, x, q);
		delete lx;
	}

	if(t == 1)
	{
		delete px;
		px = polynomialMultiplication(ux, x, x);

		delete ux;
		ux = remainderGPE_sZp(px, ax, x, q);
	}

	AST* l = list({});

	for(unsigned int i=0; i < n; i++)
	{
		AST* e = integer(i);

		AST* ri = coefficientGPE(ux, x, e);

		AST* a = l;
		AST* b = list({ ri });

		l = append(a, b);

		delete e;
		delete a;
		delete b;
	}

	delete px;
	delete ux;

	return l;
}




AST* formMatrixQBinary(AST* ax, AST* x, AST* q)
{
	assert(
		q->kind() == Kind::Integer,
		"q needs to be an integer"
	);

	AST* n = degreeGPE(ax, x);

	assert(
		n->kind() == Kind::Integer,
		"degree of the polynomial ax needs to be an integer"
	);

	unsigned int p = q->value();
	unsigned int e = n->value();

	AST* Q = matrix(n, n);

	AST* r = initialQRow(n);
	for(unsigned int i = 0; i < e; i++)
	{
		AST* ri = r->operand(i)->copy();

		Q->operand(0)->deleteOperand(i);
		Q->operand(0)->includeOperand(ri, i);
	}

	// This are the number of base steps that will
	// be available to the binary exponentiation
	// for small n, it maybe worthed to set a default
	// value
	unsigned long c = 2 * e < 20 ? 20 : 2 * e;

	signed long int** xn = new signed long int*[c];

	for(unsigned int m = 0; m < c; m++)
	{
		AST* j = formQRow(ax, x, p, e, r);

		delete r;
		r = j;

		xn[m] = new signed long[e];
		for(unsigned int t = 0; t < e; t++)
		{
			xn[m][t] = j->operand(t)->value();
		}
	}

	delete r;

	r = formQRowBinaryExp(ax, x, p, e, p, xn);

	AST* r0 = integer(0);

	for(unsigned int k=0; k<e; k++)
	{
		r0 = add({
			r0,
			mul({
				r->operand(k)->copy(),
				power(x->copy(), integer(k))
			})
		});
	}

	AST* tx = reduceAST(r0);

	AST* rx = tx;

	for(unsigned int i = 0; i < e; i++)
	{
		AST* e = integer(i);
		AST* ri = coefficientGPE(rx, x, e);

		Q->operand(1)->deleteOperand(i);
		Q->operand(1)->includeOperand(ri, i);

		delete e;
	}

	for(unsigned int m = 2; m < e; m++)
	{
		AST* zx = polynomialMultiplication(r0, rx, x);

		delete rx;

		rx = remainderGPE_sZp(zx, ax, x, p);

		for(unsigned int i = 0; i < e; i++)
		{
			AST* e = integer(i);
			AST* ri = coefficientGPE(rx, x, e);

			Q->operand(m)->deleteOperand(i);
			Q->operand(m)->includeOperand(ri, i);

			delete e;
		}

		delete zx;
	}

	delete rx;
	delete r0;
	delete n;
	delete r;

	for(unsigned int m = 0; m < c; m++)
	{
		delete[] xn[m];
	}

	delete[] xn;

	return Q;
}

AST* polyFromList(AST* l, AST* x)
{
	if(l->numberOfOperands() == 0)
	{
		return integer(0);
	}

	if(l->numberOfOperands() == 1)
	{
		return l->operand(0)->copy();
	}

	AST* px = add({});

	for(long i=0; i < l->numberOfOperands(); i++)
	{
		px->includeOperand(mul({ l->operand(i)->copy(), power(x->copy(), integer(i))}));
	}

	AST* ux = reduceAST(px);

	delete px;

	return ux;
}

AST* berlekamp(AST* ax, AST* x, AST* q)
{
	long p = q->value();

	AST* Q = formMatrixQBinary(ax, x, q);

	for(long i=0; i < Q->numberOfOperands(); i++)
	{
		long Qii = Q->operand(i)->operand(i)->value();
		Q->operand(i)->deleteOperand(i);
		Q->operand(i)->includeOperand(integer(sZp(Qii -1, p)), i);
	}

	AST* v = nullSpace_sZp(Q, p);

	AST* factors = list({ ax->copy() });

	long k = v->numberOfOperands();

	long r = 1;

	while(factors->numberOfOperands() < k)
	{
		for(long idx = 0; idx < factors->numberOfOperands(); idx++)
		{
			AST* ux = factors->operand(idx)->copy();
			AST* lx = polyFromList(v->operand(r), x);
			
			for(long s = 0; s < p; s++)
			{
				AST* kx = sub({ lx->copy(), integer(s) });

				AST* vx = reduceAST(kx);
				
				delete kx;

				AST* gx = gcdGPE_Zp(vx, ux, x, p);

				delete vx;
				
				if(gx->isNot(1) && !gx->match(ux))
				{
					factors->deleteOperand(idx);
					AST* tx = quotientGPE_Zp(ux, gx, x, p);

					delete ux;
					ux = tx;

					factors->includeOperand(ux->copy(), 0);
					factors->includeOperand(gx->copy(), 1);
				}
	
				delete gx;

				if(factors->numberOfOperands() == k)
				{
					delete Q;
					delete v;
					delete ux;
					delete lx;
					return factors;
				}
			}

			delete ux;
			delete lx;
		
			r = r + 1;
		}
	}

	delete Q;
	delete v;
	return factors;
}

long getPrime(AST* ux, AST* x)
{
	// this prime also needs to not divide the resultant between u(x) and u'(x)
	AST* lc = leadingCoefficientGPE(ux, x);

	int i = 1; // prime >= 3

	std::cout << primes[primes.count() - 1] << "\n";

	while(lc->value() % primes[i] == 0 && i < 50000000)
	{
		i++;
	}

	if(i == 50000000)
	{
		printf("polynomial leading coefficient needs to be smaller than %i", primes[i]);
		exit(1);
	}

	delete lc;
	
	return primes[i];
}

AST* magnitude(AST* ux, AST* x)
{
	AST* n = degreeGPE(ux, x);
	AST* c = leadingCoefficientGPE(ux, x);

	AST* px = sub({ux->copy(), mul({c->copy(), power(x->copy(), n->copy())})});
	
	ux = algebraicExpand(px);

	AST* u = c;
	delete px;
	delete n;

	while(ux->isNot(0))
	{

		n = degreeGPE(ux, x);
		c = leadingCoefficientGPE(ux, x);

		AST* t = max(c, u);
	
		delete u;

		u = t;
	
		px = sub({ ux->copy(), mul({c->copy(), power(x->copy(), n->copy())})});

		delete ux;

		ux = algebraicExpand(px);

		delete px;
		delete n;
		delete c;
	}

	delete ux;


	return u;
}

AST* findB(AST* ux, AST* x)
{
	AST* n = degreeGPE(ux, x);

	printf("----> %s\n", n->toString().c_str());

	AST* B = mul({
		power(integer(2), n->copy()),
		integer(ceil(sqrt(n->value() + 1))),
		magnitude(ux, x)
	});

	printf("----> %s\n", B->toString().c_str());
	AST* T = reduceAST(B);

	delete n;
	delete B;

	return T;
}

int findK(AST* p, AST* B)
{
	int k = 1;
	long q = p->value();

	while(q < 2 * B->value())
	{
		q = q * p->value();
		k++;
	}

	return k;
}

AST* irreductibleFactors(AST* ux, AST* x, AST* y)
{
	AST* t = nullptr;
	AST* g = nullptr;
	AST* j = nullptr;

	AST* n = degreeGPE(ux, x);
	AST* l = leadingCoefficientGPE(ux, x);

	t = mul({ power(l->copy(), sub({n->copy(), integer(1)})), ux->copy()});

	g = div(y->copy(), l->copy());


	printf("%s\n", t->toString().c_str());

	j = deepReplace(t, x, g);

	AST* V = algebraicExpand(j);

	printf("%s\n", V->toString().c_str());

	delete t;
	delete g;
	delete j;


	AST* p = integer(getPrime(V, y));
	printf("p = %s\n", p->toString().c_str());


	AST* vx = sZp(V, y, p->value());


	AST* S = berlekamp(vx, y, p);

	printf("%s\n", S->toString().c_str());

	// printf("%s\n", sZp(algebraicExpand(mul({S->operand(0)->copy(), S->operand(1)->copy()})), y, p->value())->toString().c_str());
	// printf("%s\n", sZp(vx, y, p->value())->toString().c_str());


	if(S->numberOfOperands() == 1)
	{
		return ux->copy();
	}

	AST* B = findB(V, y);

	int k = findK(p, B);


	AST* X = henselLift(V, S, y, p->value(), k);


	printf("hensel %s\n", X->toString().c_str());

	j = mul({l->copy(), x->copy()});

	AST* T = deepReplace(X, y, j);

	AST* W = algebraicExpand(T);

	AST* M = integer(1);

	for(unsigned int i=0; i<W->numberOfOperands(); i++)
	{
		AST* w = W->operand(i)->copy();
	
		M = mul({M, div(w, cont(w, x))});
	}

	AST* r = reduceAST(M);

	delete M;

	return r;
}



// AST* factors(AST* ux, AST* x)
// {
// 	std::pair<AST*, AST*> l = getPolynomialInZ(ux, x);
	
// 	AST* q = l.first;
// 	AST* M = l.second;
	
// 	AST* lc = leadingCoefficientGPE(q, x);

// 	long i = 0;

// 	while(lc->value() % primes[i] != 0 && i < 50000000)
// 	{
// 		i++;
// 	}

// 	if(i == 50000000)
// 	{
// 		printf("polynomial leading coefficient needs to be smaller than %li", primes[i]);
// 		exit(1);
// 	}

// 	AST* p = integer(primes[i]);

// 	AST* px = Zp(ux, x, p->value());

// 	AST* sx = squareFreeFactorization2(px, x);

// 	// AST* px = sZp(ux, x, p->value());

// 	// AST* f = berlekamp(px, x, p);

// 	// if(f->numberOfOperands() == 1)
// 	// {
// 	// 	return 
// 	// }

// }



AST* PRS_rec(AST* Gi2, AST* Gi1, AST* x, AST* hi2, int i)
{
	AST *Gi, *Hi, *hi1, *hi, *gi;
	AST *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *t9, *t10;

	t4 = pseudoRemainder(Gi2, Gi1, x);
	
	if(t4->is(0))
	{
		printf("-----> G[k] = %s\n", Gi1->toString().c_str());
		return nullptr;
	}

	t1 = sub({ degreeGPE(Gi2, x), degreeGPE(Gi1, x) });
	t2 = add({ t1, integer(1) });
	t3 = power(integer(-1), t2);

	t5 = mul({t3, t4});

	t6 = leadingCoefficientGPE(Gi2, x);
	t7 = power(hi2->copy(), sub({ degreeGPE(Gi2, x), degreeGPE(Gi1, x) }));
	t8 = mul({t6, t7});

	t9  = algebraicExpand(t5);
	t10 = algebraicExpand(t8);

	Gi = quotientGPE(t9, t10, x);

	printf("G[%i] = %s\n", i, Gi->toString().c_str());

	delete t5;
	delete t9;
	delete t8;
	delete t10;

	t1 = leadingCoefficientGPE(Gi1, x);
	t2 = sub({ degreeGPE(Gi2, x), degreeGPE(Gi1, x) });
	t3 = power(t1, t2);

	t4 = hi2->copy();
	t5 = sub({integer(1), t2->copy()});
	t6 = power(t4, t5);
	t7 = mul({t3, t6});

	hi1 = algebraicExpand(t7); // h4

	delete t7;

	t1 = leadingCoefficientGPE(Gi, x);
	t2 = sub({ degreeGPE(Gi1, x), degreeGPE(Gi, x) });
	t3 = power(t1, t2); // gi^d[i-1]

	t4 = sub({integer(1), sub({degreeGPE(Gi1, x), degreeGPE(Gi, x)})});
	t5 = power(hi1->copy(), t4);
	t6 = mul({t3, t5});

	hi = algebraicExpand(t6);

	delete t6;

	gi = leadingCoefficientGPE(Gi, x);
	
	t1 = mul({hi->copy(), Gi->copy()});
	t2 = algebraicExpand(t1);

	Hi = quotientGPE(t2, gi, x);

	delete t1;
	delete t2;

	// cleanup
	// printf("H[%i] = %s\n", i, Hi->toString().c_str());

	t8 = PRS_rec(Gi1, Gi, x, hi1 , i + 1);

	delete Gi;
	delete hi1;


	return t8;
}


AST* PRS(AST* F1, AST* F2, AST* x)
{
	AST *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8;
	AST *G1, *G2, *G3, *H3, *g3, *h2, *h3;

	G1 = F1;
	G2 = F2;

	printf("G[1] = %s\n", F1->toString().c_str());
	printf("G[2] = %s\n", F2->toString().c_str());

	// compute G[3]
	t1 = sub({ degreeGPE(F1, x), degreeGPE(F2, x) });
	t2 = add({ t1, integer(1) });

	t3 = power(integer(-1), t2);
	t4 = pseudoRemainder(G1, G2, x);
	t5 = mul({t3, t4});

	G3 = algebraicExpand(t5);

	delete t5;

	printf("G[3] = %s\n", G3->toString().c_str());

	// compute h[2]
	t6 = leadingCoefficientGPE(G2, x);
	t7 = sub({ degreeGPE(F1, x), degreeGPE(F2, x) });
	t8 = power(t6, t7);
	h2 = reduceAST(t8);

	delete t8;

	// compute h[3]
	g3 = leadingCoefficientGPE(G3, x);

	t1 = power(g3->copy(), sub({degreeGPE(G2, x), degreeGPE(G3, x)}));

	t2 = sub({degreeGPE(G2, x), degreeGPE(G3, x)});
	t3 = power(h2->copy(), sub({integer(1), t2}));

	t4 = mul({t1, t3});

	h3 = algebraicExpand(t4);

	delete t4;

	// compute H[3]
	t1 = mul({h3->copy(), G3->copy()});
	t2 = algebraicExpand(t1);

	H3 = quotientGPE(t2, g3, x);

	delete t1;
	delete t2;

	return PRS_rec(G2, G3, x, h2, 4);

}

AST* PPRS_rec(AST* Gi2, AST* Gi1, AST* x, int i)
{
	printf("***\n");

	AST *Pi, *bi;
	AST *Ri;

	Ri = pseudoRemainder(Gi2, Gi1, x);
	
	if(Ri->is(0))
	{
		printf("-----> G[k] = %s\n", Gi1->toString().c_str());
		return nullptr;
	}

	bi = cont(Ri, x);
	Pi = quotientGPE(Ri, bi, x);

	printf("%s\n", Pi->toString().c_str());
	printf("%s\n", bi->toString().c_str());
	return PPRS_rec(Gi1, Pi, x, i + 1);
}



AST* PPRS(AST* F1, AST* F2, AST* x)
{
	// AST* n = degreeGPE(F1, x);
	// AST* m = degreeGPE(F2, x);

	// assert(n->value() > m->value());

	// AST *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8; // tmp variables
	// AST *G1, *G2, *G3, *h2;

	// G1 = quotientGPE(F1, cont(F1, x));
	// G2 = quotientGPE(F2, cont(F2, x));

	// printf("G[1] = %s\n", F1->toString().c_str());
	// printf("G[2] = %s\n", F2->toString().c_str());

	// t1 = sub({ degreeGPE(F1, x), degreeGPE(F2, x) });
	// t2 = add({ t1, integer(1) });

	// t3 = power(integer(-1), t2);
	// t4 = pseudoRemainder(G1, G2, x);
	// t5 = mul({t3, t4});

	// G3 = algebraicExpand(t5);

	// printf("G[3] = %s\n", G3->toString().c_str());

	// t6 = leadingCoefficientGPE(G2, x);
	// t7 = sub({ degreeGPE(F1, x), degreeGPE(F2, x) });
	// t8 = power(t6, t7);
	// h2 = reduceAST(t8);

	return PPRS_rec(F1, F2, x, 4);

}


AST* res(AST* f, AST* g, AST* x)
{
	// printf("%s\n",
	// gcdGPE(
	// 	add({
	// 		power(symbol("x"), integer(8)),
	// 		power(symbol("x"), integer(6)),
	// 		mul({
	// 			integer(-3),
	// 			power(symbol("x"), integer(4)),
	// 		}),
	// 		mul({
	// 			integer(-3),
	// 			power(symbol("x"), integer(3)),
	// 		}),
	// 		mul({
	// 			integer(8),
	// 			power(symbol("x"), integer(2)),
	// 		}),
	// 		mul({
	// 			integer(2),
	// 			symbol("x")
	// 		}),
	// 		integer(-5)
	// 	}),
	// 	add({
	// 		mul({
	// 			integer(3),
	// 			power(symbol("x"), integer(6)),
	// 		}),
	// 		mul({
	// 			integer(5),
	// 			power(symbol("x"), integer(4)),
	// 		}),
	// 		mul({
	// 			integer(-4),
	// 			power(symbol("x"), integer(2)),
	// 		}),
	// 		mul({
	// 			integer(-9),
	// 			symbol("x")
	// 		}),
	// 		integer(21)
	// 	}),
	// 	symbol("x")
	// )->toString().c_str()	
	// );
	AST* u = add({
			power(symbol("x"), integer(8)),
			power(symbol("x"), integer(6)),
			mul({
				integer(-3),
				power(symbol("x"), integer(4)),
			}),
			mul({
				integer(-3),
				power(symbol("x"), integer(3)),
			}),
			mul({
				integer(8),
				power(symbol("x"), integer(2)),
			}),
			mul({
				integer(2),
				symbol("x")
			}),
			integer(-5)
		});

	AST* v = add({
			mul({
				integer(3),
				power(symbol("x"), integer(6)),
			}),
			mul({
				integer(5),
				power(symbol("x"), integer(4)),
			}),
			mul({
				integer(-4),
				power(symbol("x"), integer(2)),
			}),
			mul({
				integer(-9),
				symbol("x")
			}),
			integer(21)
		});
	
	// PPRS(u, v, x);
	PRS(u, v, x);

	// printf("******* = %s\n", multivariateResultant(u, v, list({x}), symbol("Z"))->toString().c_str());

	// PRS(f, g, x);
	return nullptr;


	// AST* p3 = pseudoRemainder(p1, p2, x);
	
	// printf("%s\n", p3->toString().c_str());
	// // printf("%s\n", quotientGPE(p1, p2, x)->toString().c_str());

	// AST* tmp = leadingCoefficientGPE(p2, x);

	// printf("%s\n", tmp->toString().c_str());
	// tmp = power(tmp, add({ sub({ degreeGPE(p1, x) , degreeGPE(p2, x) }), integer(1) }));
	// tmp = reduceAST(tmp);
	// printf("tmp = %s\n", tmp->toString().c_str());

	// AST* p4 = pseudoRemainder(p2, p3, x);
	
	// tmp = leadingCoefficientGPE(p3, x);

	// printf("%s\n", tmp->toString().c_str());
	// tmp = power(tmp, add({ sub({ degreeGPE(p2, x) , degreeGPE(p3, x) }), integer(1) }));
	// tmp = reduceAST(tmp);
	// printf("tmp = %s\n", tmp->toString().c_str());


	// printf("%s\n", p4->toString().c_str());
	// AST* p5 = pseudoRemainder(p3, p4, x);
	// printf("%s\n", p5->toString().c_str());


	// AST* t = extendedEuclideanAlgGPE(p1, p2, x);
	// AST* t1 = univariateResultant(f, g, x);

	// printf("%i - %s\n", i,  t1->toString().c_str());
	// AST* t2 = power(mul({integer(9), power(symbol("y"), integer(8))}), integer(-1));

	// AST* t3 = mul({t1, t2});
	// printf("%i - %s\n", i,  algebraicExpand(t3)->toString().c_str());

	// while(p3->isNot(0))
	// {
	// 	printf("%i - %s\n", i,  p3->toString().c_str());
		
	// 	delete p1;
		
	// 	p1 = p2;
	// 	p2 = p3;
		
	// 	p3 = remainderGPE(p1, p2, x);
		
	// 	i = i + 1;
	// }
	return nullptr;
}

// AST* algFactorization(AST* az, AST* z, AST* mx, AST* x, AST* alpha)
// {
// 	// az = a(z) is considere a bivariate polynomial in 'alpha' and 'z'
// 	AST* s = integer(0);
// 	AST* as = az->copy(); // as(alpha, z) = a(alpha, z)

// 	AST* bxz = deepReplace(as, alpha, x); // b(x, z) = as(x, z)
	
// 	AST* norm_as = univariateResultant(mx, bxz, x);
	
// 	delete bxz;

// 	AST* norm_as_ = derivate(norm_as, x);

// 	AST* g = gcdGPE(norm_as, norm_as_, x);
// 	AST* d = degreeGPE(norm_as, norm_as_);


// 	delete s;
// 	delete as;

// }

}
