#include "Berlekamp.hpp"

#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Polynomial/Zp.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;

namespace factorization
{

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

AST* berlekampFactors(AST* u, AST* x, int p) 
{
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


// AST* polyFromList(AST* l, AST* x)
// {
// 	if(l->numberOfOperands() == 0)
// 	{
// 		return integer(0);
// 	}

// 	if(l->numberOfOperands() == 1)
// 	{
// 		return l->operand(0)->copy();
// 	}

// 	AST* px = add({});

// 	for(long i=0; i < l->numberOfOperands(); i++)
// 	{
// 		px->includeOperand(mul({ l->operand(i)->copy(), power(x->copy(), integer(i))}));
// 	}

// 	AST* ux = reduceAST(px);

// 	delete px;

// 	return ux;
// }

// AST* polynomialMultiplication(AST* ax, AST* bx, AST* x)
// {
// 	// TODO: override this with algebraic expand when it gets optimized

// 	AST* ux = add({});

// 	ax = reduceAST(ax);
// 	bx = reduceAST(bx);

// 	AST* da = degreeGPE(ax, x);
// 	AST* db = degreeGPE(bx, x);

// 	for(unsigned int i = 0; i <= da->value(); i++)
// 	{
// 		for(unsigned int j = 0; j <= db->value(); j++)
// 		{
// 			AST* ae = integer(i);
// 			AST* be = integer(j);

// 			AST* ca = coefficientGPE(ax, x, ae);
// 			AST* cb = coefficientGPE(bx, x, be);

// 			ux->includeOperand(
// 				mul({
// 					mul({ca, cb}),
// 					power(
// 						x->copy(),
// 						add({ ae, be })
// 					)
// 				})
// 			);
// 		}
// 	}

// 	AST* px = reduceAST(ux);

// 	delete da;
// 	delete db;
	
// 	delete ux;
// 	delete ax;
// 	delete bx;

// 	return px;
// }


// AST* formQRowBinaryExp(AST* ax, AST* x, signed long q, signed long n, signed long i, signed long** cache)
// {
// 	// TODO: Maybe there is an easy way to form the r vector
// 	// without computing the remainder every time
// 	signed long t = i % 2;
// 	signed long j = (i - t)/2;


// 	AST* ux = nullptr;
// 	AST* px = nullptr;

// 	if(j >= 2*n)
// 	{
// 		AST* r = formQRowBinaryExp(ax, x, q, n, j, cache);

// 		ux = integer(0);

// 		for(signed long k=0; k<n; k++)
// 		{
// 			signed long int rk = r->operand(k)->value();

// 			ux = add({
// 				ux,
// 				mul({
// 					integer(rk),
// 					power(x->copy(), integer(k))
// 				})
// 			});
// 		}

// 		delete r;

// 		// ux = mul({ ux->copy(), ux });
// 		px = polynomialMultiplication(ux, ux, x);
// 		ux = remainderGPE_sZp(px, ax, x, q);
		
// 		delete px;
	
// 	}
// 	else
// 	{
// 		ux = integer(0);

// 		for(unsigned int k=0; k<n; k++)
// 		{
// 			signed long int rji = cache[j - 1][k];
// 			ux = add({
// 				ux,
// 				mul({
// 					integer(rji),
// 					power(x->copy(), integer(k))
// 				})
// 			});
// 		}

// 		AST* lx = polynomialMultiplication(ux, ux, x);
// 		delete ux;
	
// 		ux = remainderGPE_sZp(lx, ax, x, q);
// 		delete lx;
// 	}

// 	if(t == 1)
// 	{
// 		delete px;
// 		px = polynomialMultiplication(ux, x, x);

// 		delete ux;
// 		ux = remainderGPE_sZp(px, ax, x, q);
// 	}

// 	AST* l = list({});

// 	for(unsigned int i=0; i < n; i++)
// 	{
// 		AST* e = integer(i);

// 		AST* ri = coefficientGPE(ux, x, e);

// 		AST* a = l;
// 		AST* b = list({ ri });

// 		l = append(a, b);

// 		delete e;
// 		delete a;
// 		delete b;
// 	}

// 	delete px;
// 	delete ux;

// 	return l;
// }

// AST* formQRow(AST* ax, AST* x, unsigned int q, unsigned int n, AST* r)
// {
// 	// TODO: Maybe there is an easy way to form the r vector
// 	// without computing the remainder every time

// 	AST* e = integer(0);

// 	AST* r0 = mul({
// 		integer(-1),
// 		r->operand(n - 1)->copy(),
// 		coefficientGPE(ax, x, e),
// 	});

// 	delete e;

// 	AST* ux = r0;

// 	for(unsigned int i = 1; i < n; i++)
// 	{
// 		AST* e = integer(i);

// 		AST* ri = sub({
// 			r->operand(i - 1)->copy(),
// 			mul({
// 				r->operand(n - 1)->copy(),
// 				coefficientGPE(ax, x, e)
// 			})
// 		});

// 		delete e;

// 		ux = add({ ux, mul({ri, power(x->copy(), integer(i))}) });
// 	}

// 	AST* kx = reduceAST(ux);
// 	delete ux;

// 	ux = remainderGPE_sZp(kx, ax, x, q);

// 	delete kx;

// 	AST* l = list({});

// 	for(unsigned int i=0; i < n; i++)
// 	{
// 		AST* e = integer(i);

// 		AST* ri = coefficientGPE(ux, x, e);

// 		AST* a = l;
// 		AST* b = list({ ri });

// 		l = append(a, b);

// 		delete e;
// 		delete a;
// 		delete b;
// 	}

// 	delete ux;

// 	return l;
// }

// AST* initialQRow(AST* n)
// {
// 	AST* r = list({integer(1)});
// 	AST* z = list({integer(0)});

// 	for(unsigned int i = 1; i < n->value(); i++)
// 	{
// 		AST* k = append(r, z);
// 		delete r;
// 		r = k;
// 	}

// 	delete z;

// 	return r;
// }

// AST* formMatrixQ(AST* ax, AST* x, AST* q)
// {
// 	assert(
// 		q->kind() == Kind::Integer,
// 		"q needs to be an integer"
// 	);

// 	AST* n = degreeGPE(ax, x);

// 	assert(
// 		n->kind() == Kind::Integer,
// 		"degree of the polynomial ax needs to be an integer"
// 	);

// 	unsigned int p = q->value();
// 	unsigned int e = n->value();

// 	AST* Q = matrix(n, n);

// 	AST* r = initialQRow(n);

// 	// Set Q row
// 	for(unsigned int i = 0; i < e; i++)
// 	{
// 		AST* ri = r->operand(i)->copy();

// 		Q->operand(0)->deleteOperand(i);
// 		Q->operand(0)->includeOperand(ri, i);
// 	}

// 	for(unsigned int m = 1; m <= (e - 1)*p; m++)
// 	{
// 		AST* j = formQRow(ax, x, p, e, r);

// 		delete r;
// 		r = j;

// 		if(m % p == 0)
// 		{
// 			// Set Q row
// 			for(unsigned int i = 0; i < e; i++)
// 			{
// 				AST* ri = r->operand(i)->copy();

// 				Q->operand(m / p)->deleteOperand(i);
// 				Q->operand(m / p)->includeOperand(ri, i);
// 			}
// 		}
// 	}

// 	delete n;
// 	delete r;

// 	return Q;
// }


// AST* formMatrixQBinary(AST* ax, AST* x, AST* q)
// {
// 	assert(
// 		q->kind() == Kind::Integer,
// 		"q needs to be an integer"
// 	);

// 	AST* n = degreeGPE(ax, x);

// 	assert(
// 		n->kind() == Kind::Integer,
// 		"degree of the polynomial ax needs to be an integer"
// 	);

// 	unsigned int p = q->value();
// 	unsigned int e = n->value();

// 	AST* Q = matrix(n, n);
// 	AST* r = initialQRow(n);

// 	for(unsigned int i = 0; i < e; i++)
// 	{
// 		AST* ri = r->operand(i)->copy();

// 		Q->operand(0)->deleteOperand(i);
// 		Q->operand(0)->includeOperand(ri, i);
// 	}

// 	// This are the number of base steps that will
// 	// be available to the binary exponentiation
// 	// for small n, it maybe worthed to set a default
// 	// value
// 	unsigned long c = 2 * e < 20 ? 20 : 2 * e;

// 	signed long int** xn = new signed long int*[c];

// 	for(unsigned int m = 0; m < c; m++)
// 	{
// 		AST* j = formQRow(ax, x, p, e, r);

// 		delete r;
// 		r = j;

// 		xn[m] = new signed long[e];
// 		for(unsigned int t = 0; t < e; t++)
// 		{
// 			xn[m][t] = j->operand(t)->value();
// 		}
// 	}

// 	delete r;

// 	r = formQRowBinaryExp(ax, x, p, e, p, xn);

// 	AST* r0 = integer(0);

// 	for(unsigned int k=0; k<e; k++)
// 	{
// 		r0 = add({
// 			r0,
// 			mul({
// 				r->operand(k)->copy(),
// 				power(x->copy(), integer(k))
// 			})
// 		});
// 	}

// 	AST* tx = reduceAST(r0);

// 	AST* rx = tx;

// 	for(unsigned int i = 0; i < e; i++)
// 	{
// 		AST* e = integer(i);
// 		AST* ri = coefficientGPE(rx, x, e);

// 		Q->operand(1)->deleteOperand(i);
// 		Q->operand(1)->includeOperand(ri, i);

// 		delete e;
// 	}

// 	for(unsigned int m = 2; m < e; m++)
// 	{
// 		AST* zx = polynomialMultiplication(r0, rx, x);

// 		delete rx;

// 		rx = remainderGPE_sZp(zx, ax, x, p);

// 		for(unsigned int i = 0; i < e; i++)
// 		{
// 			AST* e = integer(i);
// 			AST* ri = coefficientGPE(rx, x, e);

// 			Q->operand(m)->deleteOperand(i);
// 			Q->operand(m)->includeOperand(ri, i);

// 			delete e;
// 		}

// 		delete zx;
// 	}

// 	delete rx;
// 	delete r0;
// 	delete n;
// 	delete r;

// 	for(unsigned int m = 0; m < c; m++)
// 	{
// 		delete[] xn[m];
// 	}

// 	delete[] xn;

// 	return Q;
// }

// AST* berlekamp(AST* ax, AST* x, AST* q)
// {
// 	long p = q->value();

// 	AST* Q = formMatrixQBinary(ax, x, q);

// 	printf("--------------> %s\n", ax->toString().c_str());
// 	printf("--------------> %s\n", Q->toString().c_str());

// 	for(long i=0; i < Q->numberOfOperands(); i++)
// 	{
// 		long Qii = Q->operand(i)->operand(i)->value();
// 		Q->operand(i)->deleteOperand(i);
// 		Q->operand(i)->includeOperand(integer(sZp(Qii -1, p)), i);
// 	}

// 	AST* v = nullSpace_sZp(Q, p);

// 	printf("%s\n", v->toString().c_str());

// 	AST* factors = list({ ax->copy() });

// 	long k = v->numberOfOperands();

// 	long r = 1;

// 	while(factors->numberOfOperands() < k)
// 	{
// 		for(long idx = 0; idx < factors->numberOfOperands(); idx++)
// 		{
// 			printf("A\n");
// 			printf("%s\n", ax->toString().c_str());
// 			printf("%s\n", factors->toString().c_str());

// 			AST* ux = factors->operand(idx)->copy();

// 			printf("%s\n", ux->toString().c_str());
// 			printf("%i\n", r);
// 			printf("%s\n", v->toString().c_str());

// 			AST* lx = polyFromList(v->operand(r), x);

// 			printf("B\n");

// 			for(long s = 0; s < p; s++)
// 			{
// 				AST* kx = sub({ lx->copy(), integer(s) });

// 				AST* vx = reduceAST(kx);
				
// 				delete kx;

// 				AST* gx = gcdGPE_Zp(vx, ux, x, p);

// 				delete vx;
				
// 				if(gx->isNot(1) && !gx->match(ux))
// 				{
// 					factors->deleteOperand(idx);
// 					AST* tx = quotientGPE_Zp(ux, gx, x, p);

// 					delete ux;
// 					ux = tx;

// 					factors->includeOperand(ux->copy(), 0);
// 					factors->includeOperand(gx->copy(), 1);
// 					idx = 0;
// 				}
	
// 				delete gx;

// 				if(factors->numberOfOperands() == k)
// 				{
// 					delete Q;
// 					delete v;
// 					delete ux;
// 					delete lx;
// 					return factors;
// 				}
// 			}
// 			printf("C\n");

// 			delete ux;
// 			delete lx;
		
// 			r = r + 1;
// 		}
// 	}

// 	delete Q;
// 	delete v;
// 	return factors;
// }

}

