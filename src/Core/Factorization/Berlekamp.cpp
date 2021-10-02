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
			AST* deg = degree(yk, x);

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

AST* auxiliaryBasis(AST* x, AST* n, int p) 
{
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

		if(pivot_found) 
		{
			P[i-1] = j;

			int a = modInverse_p(R[i-1][j-1], p);

			for(int l=1; l<=n->value(); l++) 
			{
				R[i-1][l-1] = mod(a * R[i-1][l-1], p);
			}

			for(int k=1; k <= n->value(); k++) 
			{
				if(k!=i) 
				{
					int f = R[k-1][j-1];
					for(int l=1; l <= n->value(); l++) 
					{
						R[k-1][l-1] = mod(R[k-1][l-1] - f * R[i-1][l-1], p);
					}
				}
			}

		} else if(!pivot_found) 
		{
			AST* s = power(
				x->copy(),
				sub({ integer(j), integer(1) })
			);

			for(int l=1; l <= j-1; l++) 
			{
				int e = 0;
				int i = 1;

				while(e == 0 && i< n->value()) 
				{
					if(l == P[i-1]) 
					{
						e = i;
					} else {
						i = i+1;
					}
				}
				if(e > 0) 
				{
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

bool isRowOfZeros(AST* M, int n, int j)
{
	for(int i=0; i<n; i++)
	{
		if(M->operand(j)->operand(i)->isNot(0))
			return false;
	}

	return true;
}

void swapColumn(AST* M, long j, long t)
{
	for(long i = 0; i < M->numberOfOperands(); i++)
	{
		AST* Mji = M->operand(j)->operand(i)->copy();
		AST* Mti = M->operand(t)->operand(i)->copy();

		M->operand(j)->deleteOperand(i);
		M->operand(t)->deleteOperand(i);
		
		M->operand(j)->includeOperand(Mti, i);
		M->operand(t)->includeOperand(Mji, i);
	}
}

AST* nullSpace_Zp(AST* A, signed long q)
{
	// M = M->copy();
	bool zeros, next_pivo;

	long k, i, j, n, t, p;

	n = A->numberOfOperands();

	// A = A->copy();
	// // M = M - I
	// for(i = 0; i < n; i++)
	// {
	// 	long tt = A->operand(i)->operand(i)->value() - 1;

	// 	A->operand(i)->deleteOperand(i);
	// 	A->operand(i)->includeOperand(integer(mod(tt, q)), i);
	// }
	
	AST* M = A->copy();
	// M = M - I
	for(i = 0; i < n; i++)
	{
		long tt = M->operand(i)->operand(i)->value() - 1;

		M->operand(i)->deleteOperand(i);
		M->operand(i)->includeOperand(integer(mod(tt, q)), i);
	}
	// for(i = 0; i < n; i++)
	// {
	// 	M->includeOperand(list({}));
	// }
	
	// for(i = 0; i < n; i++)
	// {
	// 	for(j = 0; j < n; j++)
	// 	{
	// 		M->operand(i)->includeOperand(A->operand(i)->operand(j)->copy(), j);
	// 	}
	// }

	printf("A %s\n", A->toString().c_str());
	printf("M %s\n", M->toString().c_str());

	t = n;

	printf("%s\n", M->toString().c_str());

	// perform gausian elimination on M
	for(j = 0; j < t; j++)
	{
		zeros = true;
		
		for(k = 0; k < n; k++)
		{
			if(M->operand(j)->operand(k)->isNot(0))
			{
				zeros = false;
			}
		
			if(!zeros) break;
		}
		
		if(zeros == true)
		{
			swapColumn(M, j, t - 1);

			t = t - 1;
		}
	}

	printf("%s\n", M->toString().c_str());

	p = 0;

	while(p < t && p < n)
	{
		next_pivo = true;

		j = 1;

		while(M->operand(p)->operand(p)->is(0))
		{
			if(p + j <= t)
			{
				p = p + 1;
				next_pivo = false;
			}
	
			if(!next_pivo) break;

			swapColumn(M, p, p + j);
			
			j = j + 1;
		}
	
		if(next_pivo)
		{
			for(j = 1; j < (t - p); j++)
			{
				AST* Mip = M->operand(p+j)->operand(p);
				AST* Mpp = M->operand(p)->operand(p);

				if(Mip->isNot(0))
				{
					long x = division_Zp(-Mip->value(), Mpp->value(), q); 
				
					printf("%li\n", x);
				
					for(k = p; k < n; k++)
					{
						long Mpk = M->operand(p)->operand(k)->value();
						long Mjk = M->operand(p + j)->operand(k)->value();

						// long v = Mpk * x + Mjk;
						long v = mod(mod(Mpk * x, q) + Mjk, q);
					
						M->operand(p + j)->deleteOperand(k);
						M->operand(p + j)->includeOperand(integer(v), k);
					}
				}
			}
			p = p + 1;
		}
	}

	printf("--> %s\n", M->toString().c_str());

	return list({});

	// for(k = 0; k < n; k++)
	// {
	// 	i = k;
	
	// 	printf("\n%s\n", M->toString().c_str());
	
	// 	while(i < n && M->operand(i)->operand(k)->is(0))
	// 	{
	// 		i = i + 1;
	// 	}

	// 	printf("(%li, %li)\n", i, k);

	// 	if(i < n)
	// 	{
	// 		// Normalized column i 
	// 		long d = M->operand(i)->operand(k)->value();
		
	// 		for(j = 0; j < n; j++)
	// 		{
	// 		 	long Mji = M->operand(i)->operand(j)->value();

	// 			long p = division_Zp(Mji, d, q);
			
	// 			M->operand(i)->deleteOperand(j);
	// 			M->operand(i)->includeOperand(integer(p), j);
	// 		}
	
	// 		printf("%s\n", M->toString().c_str());

	// 		// Switch column i with column k
	// 		for(j = 0; j < n; j++)
	// 		{
	// 			AST* Mji = M->operand(i)->operand(j)->copy();
	// 			AST* Mjk = M->operand(k)->operand(j)->copy();

	// 			M->operand(i)->deleteOperand(j);
	// 			M->operand(i)->includeOperand(Mjk, j);
			
	// 			M->operand(k)->deleteOperand(j);
	// 			M->operand(k)->includeOperand(Mji, j);
	// 		}

	// 		printf("%s\n", M->toString().c_str());

	// 		// Eliminate rest of row k via column operations
	// 		for(i = 0; i < n; i++)
	// 		{
	// 			if(i != k)
	// 			{
	// 				signed long Mki = M->operand(i)->operand(k)->value();

	// 				for(j = 0; j < n; j++)
	// 				{
	// 					signed long col_i = M->operand(i)->operand(j)->value();
	// 					signed long col_k = M->operand(k)->operand(j)->value();

	// 					signed long tt = mod(col_i - col_k * Mki, q);

	// 					M->operand(i)->deleteOperand(j);
	// 					M->operand(i)->includeOperand(integer(tt) ,j);
	// 				}
	// 			}
	// 		}

	// 		printf("%s\n", M->toString().c_str());
	// 	}
	// }

	// printf("%s\n", M->toString().c_str());
	// // Convert M to M-I
	// // for(i = 0; i < n; i++)
	// // {
	// // 	signed long tt = M->operand(i)->operand(i)->value() - 1;

	// // 	M->operand(i)->deleteOperand(i);
	// // 	M->operand(i)->includeOperand(integer(mod(tt, q)), i);
	// // }

	// i = 0;
	// j = 0;

	// printf("aaaaaa\n");
	// printf("%s\n", M->toString().c_str());

	// AST* v = list({});

	// while(j < n)
	// {
	// 	while(j < n && isRowOfZeros(M, n, j))
	// 	{
	// 		j = j + 1;
	// 	}

	// 	if(j < n)
	// 	{
	// 		i = i + 1;

	// 		AST* r = list({});

	// 		for(int k = 0; k < n; k++)
	// 		{
	// 			r->includeOperand(integer(-1 * M->operand(j)->operand(k)->value()));
	// 		}

	// 		v->includeOperand(r);
	// 	}

	// 	j = j + 1;
	// }

	// printf("aaaaaa\n");
	// delete M;
	// printf("aaaaaa\n");

	// return v;
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
	AST* n = degree(u, x);
	if(
		(n->kind() == Kind::Integer && n->value() == 0) ||
		(n->kind() == Kind::Integer && n->value() == 1)
	) {

		delete n;

		return set({ u->copy() });
	}

	RMatrix(u, x, n, p);

	printf("[");
	for(int i = 0; i < n->value(); i++)
	{
		printf("[");
		for(int j = 0; j < n->value(); j++)
		{
			printf("%i ", R[i][j]);
		}
		printf("]");
		if(i < n->value() - 1)
		{
			printf(", ");
		}
	}
	printf("]\n");

	AST* S = auxiliaryBasis(x, n, p);

	printf("%s\n", S->toString().c_str());

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

AST* initBerkelampBasisMatrix(AST* n)
{
	AST* Q = list({});

	Q->includeOperand(list({}));

	Q->operand(0)->includeOperand(integer(1), 0);

	for(long i = 1; i < n->value(); i++)
	{
		Q->operand(0)->includeOperand(integer(0));
	}

	return Q;
}

void addVectorToBerkelampBasisMatrixRow(AST* Q, AST* r, AST* x, long i, long n)
{
	AST *ex, *ri;

	Q->includeOperand(list({}), i);

	for(long k = 0; k < n; k++)
	{
		ex = integer(k);
		
		ri = coefficientGPE(r, x, ex);

		Q->operand(i)->includeOperand(ri);
	
		delete ex;
	}
}

// compute x^p mod a(x) by repeat squaring
AST* repeadSquaring(AST* x, AST* a, AST* p)
{
	AST *t1, *t2, *t3, *b[64];

	long v = p->value();

	long k = 0;

	while (v >>= 1) k++;
	
	b[k] = x->copy();

	for (long i = k - 1; i >= 0; i--)
	{
		t1 = mulPoly(b[i + 1], b[i + 1]);

		t2 = reduceAST(t1);
		
		delete t1;

		t1 = remainderGPE_Zp(t2, a, x, p->value());
	
		delete t2;

		if(p->value() & (1 << i))
		{
		
			t2 = mulPoly(t1, x);

			t3 = reduceAST(t2);
		
			b[i] = remainderGPE_Zp(t3, a, x, p->value());

			delete t2;
			delete t3;
		}
		else
		{
			b[i] = remainderGPE_Zp(t1, a, x, p->value());
		}

		delete t1;
	}

	for(int i = 1; i<=k; i++) delete b[i];
	
	return b[0];
}

AST* buildBerkelampBasisMatrix(AST* ax, AST* x, AST* p)
{
	AST *n, *Q, *r0, *rx, *zx, *qx;

	n = degree(ax, x);
	Q = initBerkelampBasisMatrix(n);

	// compute x^p mod a(x)
	r0 = repeadSquaring(x, ax, p);

	// add x^p mod a(x) to the first line of Q
	addVectorToBerkelampBasisMatrixRow(Q, r0, x, 1, n->value());

	rx = r0->copy();

	// compute and add x^(p*i) mod a(x) to Q[i]
	// x^(p*i) mod a(x) = x^(p*i-1) * x^p mod a(x)
	for(long m = 2; m < n->value(); m++)
	{
		zx = mulPoly(r0, rx);

		qx = reduceAST(zx);
	
		delete rx;

		// r(x) = x^(p*i) mod a(x)
		rx = remainderGPE_Zp(qx, ax, x, p->value());

		// add coefficients to Q[m][:]
		addVectorToBerkelampBasisMatrixRow(Q, rx, x, m, n->value());

		delete zx;
	}

	delete rx;
	delete r0;
	delete n;

	return Q;
}

AST* berlekamp(AST* ax, AST* x, AST* q)
{
	long p = q->value();

	AST* Q = buildBerkelampBasisMatrix(ax, x, q);

	// printf("--------------> %s\n", ax->toString().c_str());
	printf("%s\n", Q->toString().c_str());

	for(long i=0; i < Q->numberOfOperands(); i++)
	{
		long Qii = Q->operand(i)->operand(i)->value();
		Q->operand(i)->deleteOperand(i);
		Q->operand(i)->includeOperand(integer(sZp(Qii -1, p)), i);
	}

	AST* v = nullSpace_sZp(Q, p);
	printf("%s\n", v->toString().c_str());

	AST* factors = list({ ax->copy() });

	long k = v->numberOfOperands();

	long r = 1;

	while(factors->numberOfOperands() < k)
	{
		for(long idx = 0; idx < factors->numberOfOperands(); idx++)
		{
			// printf("A\n");
			// printf("%s\n", ax->toString().c_str());
			// printf("%s\n", factors->toString().c_str());

			AST* ux = factors->operand(idx)->copy();

			// printf("%s\n", ux->toString().c_str());
			// printf("%li\n", r);
			// printf("%s\n", v->toString().c_str());

			AST* lx = polyFromList(v->operand(r), x);

			printf("B\n");

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
					idx = 0;
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
			printf("C\n");

			delete ux;
			delete lx;
		
			r = r + 1;
		}
	}

	delete Q;
	delete v;
	return factors;
}

}

