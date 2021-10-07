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

void swapRows(AST* M, long j, long t)
{
	for(long i = 0; i < M->numberOfOperands(); i++)
	{
		AST* Mji = M->operand(j)->operand(i)->copy();
		AST* Mti = M->operand(t)->operand(i)->copy();
	
		M->operand(j)->deleteOperand(i);
		M->operand(j)->includeOperand(Mti, i);

		M->operand(t)->deleteOperand(i);
		M->operand(t)->includeOperand(Mji, i);
	}
}

AST* matGet(AST* M, long i, long j)
{
	return M->operand(i)->operand(j);
}

AST* matSet(AST* M, long i, long j, AST* Mij)
{
	M->operand(i)->deleteOperand(j);
	M->operand(i)->includeOperand(Mij, j);
	return M;
}

void addFreeVariableToBase(AST* v, long n, long var_idx)
{
	v->includeOperand(list({}));
	
	for(long k = 0; k < n; k++)
	{
		if(k == var_idx)
		{
			v->operand(v->numberOfOperands() - 1)->includeOperand(integer(1));
		}
		else
		{
			v->operand(v->numberOfOperands() - 1)->includeOperand(integer(0));
		}
	}
}

// from the matrix Q, compute the auxiliary basis
// that is the left nullspace of (M - I)
AST* buildBerlekampBasis(AST* A, AST* w)
{
	AST* M;
	long lead, row_count, col_count, r, i, n, j, k, x, q;

	q = w->value();
	
	n = A->numberOfOperands();

	M = list({});

	// M = (A - I)'
	for(i = 0; i < n; i++)
	{
		M->includeOperand(list({}));

		for(j = 0; j < n; j++)
		{
			if(i == j)
			{
				M->operand(i)->includeOperand(integer(mod(A->operand(j)->operand(i)->value() - 1, q)), j);
			}
			else
			{
				M->operand(i)->includeOperand(integer(A->operand(j)->operand(i)->value()), j);
			}
		}
	}

	printf("-> %s\n", M->toString().c_str());

	lead = 0;

	row_count = n;
	col_count = n;

	for(r = 0; r < row_count; r++)
	{
		printf("-> %s\n", M->toString().c_str());

		if(col_count <= lead)
		{
			break;
		}

		i = r;

		while(matGet(M, i, lead)->is(0))
		{
			i = i + 1;

			if(row_count == i)
			{
				i = r;
				lead = lead + 1;

				if(col_count == lead)
				{
					break;
				}
			}

		}
		if(col_count == lead)
		{
			break;
		}

		swapRows(M, i, r);

		if(matGet(M, r, lead)->isNot(0))
		{
			long x = matGet(M, r, lead)->value();
			
			for(j = 0; j < n; j++)
			{
				long v = matGet(M, r, j)->value();
				long Mrj = division_Zp(v, x, q);

				matSet(M, r, j, integer(Mrj));
			}
		}

		for(i = 0; i < row_count; i++)
		{
			if(i != r)
			{
				long x = matGet(M, i, lead)->value();
				
				for(j = 0; j < n; j++)
				{
					long v = matGet(M, r, j)->value();
					long t = matGet(M, i, j)->value();
					
					long Mij = mod(t - mod(x*v, q), q);

					matSet(M, i, j, integer(Mij));
				}
			}
		}
	}
	printf("-> %s\n", M->toString().c_str());

	AST* v = list({});

	k = 0;

	for(i = 0; i < n; i++)
	{
		addFreeVariableToBase(v, n, -1);
	}

	k = 0;

	for(i = 0; i < n; i++)
	{
		while(k < n && matGet(M, i, k)->is(0))
		{
			v->operand(k)->deleteOperand(k);
			v->operand(k)->includeOperand(integer(1), k);
			k++;
		}

		if(k < n)
		{
			for(j = 0; j < n; j++)
			{
				if(j != k)
				{
					x = mod(-1 * matGet(M, i, j)->value(), q);
					
					v->operand(j)->deleteOperand(k);
					v->operand(j)->includeOperand(integer(x), k);
				}
			}
		}

		k = k + 1;
	}


	delete M;

	for(i = 0; i < v->numberOfOperands(); i++)
	{
		bool rem = true;

		for(j = 0; j < n && rem; j++)
		{
			if(v->operand(i)->operand(j)->isNot(0))
			{
				rem = false;
			}
		}
	
		if(rem)
		{
			v->deleteOperand(i);
			i = i - 1;
		}
	}

	return v;
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
		
		ri = coeff(r, x, ex);

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

AST* buildBerkelampMatrix(AST* ax, AST* x, AST* p)
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
	
		delete zx;
	
		delete rx;

		// r(x) = x^(p*i) mod a(x)
		rx = remainderGPE_Zp(qx, ax, x, p->value());

		// add coefficients to Q[m][:]
		addVectorToBerkelampBasisMatrixRow(Q, rx, x, m, n->value());

		delete qx;
	}

	delete rx;
	delete r0;
	delete n;

	return Q;
}

AST* buildBerlekampBasisPolynomials(AST* B, AST* x, AST* n)
{
	long i, j;

	AST* basis = list({});

	// Build polynomial basis ignoring the '1' basis
	for(i = 1; i < B->numberOfOperands(); i++) 
	{
		AST* bx = add({});

		for(j = 0; j < B->operand(i)->numberOfOperands(); j++)
		{
			if(B->operand(i)->operand(j)->isNot(0))
			{
				bx->includeOperand(
					mul({
						B->operand(i)->operand(j)->copy(),
						power(
							x->copy(),
							integer(n->value() - j - 1)
						)
					})
				);		
			}
		}
	
		basis->includeOperand(reduceAST(bx));
		
		delete bx;
	}

	return basis;
}

AST* berlekampFactors(AST* sfx, AST* x, AST* p)
{
	long i, k, s, r;

	AST *Q, *B, *lc, *fx, *n, *H, *F, *h, *v, *g, *f;

	// f(x) = sf(x)/lc(sf(x))
	lc = leadCoeff(sfx, x);
	fx = quotientGPE_Zp(sfx, lc, x, p->value());

	n = degree(fx, x);

	printf("f(x) = %s\n", fx->toString().c_str());

	// Build the Berlekamp basis
	Q = buildBerkelampMatrix(fx, x, p);
	printf("Q = %s\n", Q->toString().c_str());

	B = buildBerlekampBasis(Q, p);
	printf("B = %s\n", B->toString().c_str());

	delete Q;

	if(B->numberOfOperands() == 1)
	{
		delete lc;
		delete fx;
	
		delete B;
	
		delete n;
	
		return sfx->copy();
	}

	r = B->numberOfOperands();

	H = buildBerlekampBasisPolynomials(B, x, n);
	printf("H = %s\n", H->toString().c_str());

	delete B;

	k = 0;
	i = 0;

	F = set({ fx });

	f = F->operand(0);
	
	for(k = 0; k < r; k++)
	{
		if(F->numberOfOperands() == r || f->is(1))
		{
			break;
		}

		f = F->operand(i);
		h = H->operand(k);

		for(s = 0; s < p->value(); s++)
		{
			g = sub({ h->copy(), integer(s) });
			
			v = Zp(g, x, p->value());
		
			delete g;
	
			g = gcdGPE_Zp(f, v, x, p->value());

			delete v;

			if(g->isNot(1) && !g->match(f))
			{
				v = quotientGPE_Zp(f, g, x, p->value());

				F->deleteOperand(i);
			
				F->includeOperand(g);
				F->includeOperand(v);
				
				i = F->numberOfOperands() - 1;
				
				f = F->operand(i);
			}
			else
			{
				delete g;
			}
		}
	}

	F->includeOperand(lc, 0);

	printf("%s\n", F->toString().c_str());

	delete n;
	delete H;

	return F;
}


// int** R = nullptr;

// int getRMatrixValue(int i, int j) {
// 	return R[i][j];
// }

// void destroyRMatrix(int n) {
// 	for(int i=0; i < n; i++) {
// 		delete []R[i];
// 	}

// 	delete []R;

// 	R = nullptr;
// }

// void RMatrix(AST* u, AST* x, AST* n_, int p) {
// 	int n = n_->value();

// 	// if(R != nullptr) {
// 	// 	destroyRMatrix(n);
// 	// }

// 	R = new int*[n];
// 	for(int i=0; i<n; i++)
// 		R[i] = new int[n];

// 	AST* yk = integer(1);

// 	AST* n_min_one = integer(n-1);

// 	AST* v_ = integer(0);

// 	for(int i=0; i<n; i++) {
// 		AST* e = integer(i);
// 		AST* c = coeff(u, x, e);

// 		v_ = add({
// 			v_,
// 			mul({
// 				integer(
// 					mod(c->value(), p)
// 				),
// 				power(
// 					x->copy(),
// 					e
// 				)
// 			})
// 		});

// 		delete c;
// 	}

// 	AST* v = reduceAST(v_);
// 	delete v_;

// 	for(int j=0; j<n; j++)
// 	{
// 		for(int i=0; i<n; i++)
// 		{
// 			AST* e = integer(i);
// 			AST* coeff_ = coeff(yk, x, e);
// 			AST* coeff = reduceAST(coeff_);

// 			if(i == j)
// 			{
// 				R[i][j] = mod(coeff->value() - 1, p);
// 			} else
// 			{
// 				R[i][j] = mod(coeff->value(),p);
// 			}

// 			delete e;
// 			delete coeff;
// 			delete coeff_;
// 		}

// 		if(j == n - 1) break;

// 		for(int i = p*(j+1); i < p*(j+2); i++) {

// 			AST* ck = coeff(yk, x, n_min_one);
// 			AST* zk_ = integer(0);

// 			for(int i=n-2; i>=0; i--) {
// 				AST* e = integer(i);

// 				AST* c = coeff(yk, x, e);

// 				zk_ = add({
// 					zk_,
// 					mul({
// 						power(
// 							x->copy(),
// 							e->copy()
// 						),
// 						integer(mod(c->value(), p))
// 					})
// 				});
// 				delete e;
// 				delete c;
// 			}

// 			AST* zk = reduceAST(zk_);
// 			delete zk_;

// 			AST* yk_ = add({
// 				mul({ integer(-1), ck, v->copy() }),
// 				mul({ x->copy(), zk })
// 			});

// 			delete yk;
// 			yk = algebraicExpand(yk_);
// 			delete yk_;

// 			// project yk into Zp
// 			AST* yk_p = integer(0);
// 			AST* deg = degree(yk, x);

// 			for(int s=deg->value(); s>=0; s--) {
// 				AST* ex = integer(s);
// 				AST* coeff = coeff(yk, x, ex);

// 				yk_p = add({
// 					yk_p,
// 					mul({
// 						integer(mod(coeff->value(), p)),
// 						power(x->copy(), ex)
// 					})
// 				});

// 				delete coeff;
// 			}

// 			delete yk;
// 			delete deg;

// 			yk = reduceAST(yk_p);

// 			delete yk_p;
// 		}
// 	}

// 	delete yk;
// 	delete n_min_one;
// 	delete v;
// }

// bool isRowOfZeros(AST* M, int n, int j)
// {
// 	for(int i=0; i<n; i++)
// 	{
// 		if(M->operand(j)->operand(i)->isNot(0))
// 			return false;
// 	}

// 	return true;
// }

// AST* findFactors(AST* u, AST* S, AST* x, int p) {
// 	signed long r = S->numberOfOperands();

// 	AST* factors = set({ u->copy() });

// 	for(int k=2; k <= r; k++) {

// 		AST* b = S->operand(k - 1)->copy();

// 		AST* old_factors = factors->copy();

// 		for(unsigned int i = 0; i < old_factors->numberOfOperands(); i++) {
// 			AST* w = old_factors->operand(i)->copy();

// 			int j = 0;

// 			while(j <= p - 1) {

// 				AST* b__ = add({
// 					b->copy(),
// 					integer(mod(-1*j,p))
// 				});

// 				AST* b_ = reduceAST(b__);

// 				delete b__;

// 				AST* g = gcdGPE_Zp(b_, w, x, p);

// 				delete b_;

// 				if(g->kind() == Kind::Integer && g->value() == 1) {
// 					j = j+1;
// 				} else if(g->match(w)) {
// 					j = p;
// 				} else {
// 					AST* factors_;
// 					AST* S0 = set({ w->copy() });
// 					factors_ = difference(factors, S0);
// 					delete S0;
// 					delete factors;
// 					factors = factors_;

// 					AST* z = divideGPE_Zp(w, g, x, p);

// 					AST* q = z->operand(0)->copy();

// 					delete z;

// 					AST* S1 = set({g->copy(), q->copy()});
// 					factors_ = unification(factors, S1);
// 					delete S1;
// 					delete factors;
// 					factors = factors_;


// 					if(factors->numberOfOperands() == r) {

// 						delete w;
// 						delete g;
// 						delete q;
// 						delete b;
// 						delete old_factors;

// 						return factors;
// 					} else {
// 						j = j + 1;

// 						delete w;

// 						w = q->copy();
// 					}

// 					delete q;
// 				}

// 				delete g;
// 			}
// 			delete w;
// 		}

// 		delete b;
// 		delete old_factors;
// 	}

// 	return factors;
// }

// void swapColumn(int** M, long j, long t, long n)
// {
// 	for(long i = 0; i < n; i++)
// 	{
// 		int Mji = M[j][i];
// 		int Mti = M[t][i];

// 		M[j][i] = Mti;
// 		M[t][i] = Mji;
// 	}
// }

// AST* buildBerlekampBasis(AST* x, AST* n, int p) {
// 	int P[n->value()];

// 	for(int i=0; i<n->value(); i++) 
// 	{
// 		P[i] = 0;
// 	}
// 	bool zeros = false;

// 	AST* S = list({});
// 	long t = n ->value();

// 	long j, k;

// 	for(j = 0; j < t; j++)
// 	{
// 		zeros = true;
		
// 		for(k = 0; k < n->value(); k++)
// 		{
// 			if(R[j][k] != 0)
// 			{
// 				zeros = false;
// 			}
		
// 			if(!zeros) break;
// 		}
		
// 		if(zeros == true)
// 		{
// 			swapColumn(R, j, t - 1, n->value());
// 			t = t - 1;
// 		}
// 	}

// 	for(j = 0; j < n->value(); j++)
// 	{
// 		for(k = 0; k < n->value(); k++)
// 		{
// 			printf("%i ", R[j][k]);
// 		}
// 		printf("\n");
// 	}

// 	for(int j=0; j<n->value(); j++) 
// 	{

// 		int i = 0;

// 		bool pivot_found = false;

// 		while(!pivot_found && i < n->value()) 
// 		{
// 			if(R[i][j] != 0 && P[i] == 0) 
// 			{
// 				pivot_found = true;
// 			} 
// 			else 
// 			{
// 				i = i + 1;
// 			}
// 		}

// 		if(pivot_found) 
// 		{
// 			P[i] = j;

// 			int a = modInverse_p(R[i][j], p);

// 			for(int l=0; l < n->value(); l++) 
// 			{
// 				R[i][l] = mod(a * R[i][l], p);
// 			}

// 			for(int k=0; k < n->value(); k++) 
// 			{
// 				if(k != i) 
// 				{
// 					int f = R[k][j];
// 					for(int l=0; l < n->value(); l++) 
// 					{
// 						R[k][l] = mod(R[k][l] - f * R[i][l], p);
// 					}
// 				}
// 			}
// 		} 
// 		else if(!pivot_found)
// 		{
// 			printf("***\n");
// 			AST* s = power(x->copy(), sub({ integer(j)}));

// 			for(int l = 0; l < j; l++) 
// 			{
// 				int e = -1;
// 				int i = 0;

// 				while(e < 0 && i < n->value()) 
// 				{
// 					if(l == P[i]) 
// 					{
// 						e = i;
// 					} else
// 					{
// 						i = i+1;
// 					}
// 				}
	
// 				if(e >= 0) 
// 				{
				
// 					int c = mod(-1*R[e][j], p);
				
// 					s = add({ s, mul({ integer(c), power(x->copy(), sub({ integer(n->value() - l) })) }) });
// 				}
// 			}

// 			AST* L = list({ algebraicExpand(s) });
// 			AST* S_ = join(S, L);

// 			delete s;
// 			delete S;
// 			delete L;

// 			S = S_;
// 		}
// 	}

// 	return S;
// }

// AST* berlekampFactors(AST* u, AST* x, int p) 
// {
// 	AST* n = degree(u, x);
// 	if(
// 		(n->kind() == Kind::Integer && n->value() == 0) ||
// 		(n->kind() == Kind::Integer && n->value() == 1)
// 	) {

// 		delete n;

// 		return set({ u->copy() });
// 	}

// 	RMatrix(u, x, n, p);

// 	printf("[");
// 	for(int i = 0; i < n->value(); i++)
// 	{
// 		printf("[");
// 		for(int j = 0; j < n->value(); j++)
// 		{
// 			printf("%i ", R[i][j]);
// 		}
// 		printf("]");
// 		if(i < n->value() - 1)
// 		{
// 			printf(", ");
// 		}
// 	}
// 	printf("]\n");

// 	AST* S = buildBerlekampBasis(x, n, p);

// 	printf("%s\n", S->toString().c_str());

// 	if(S->numberOfOperands() == 1) {

// 		delete S;
// 		delete n;

// 		destroyRMatrix(n->value());

// 		return set({u->copy()});
// 	}

// 	AST* factors = findFactors(u, S, x, p);

// 	destroyRMatrix(n->value());

// 	delete S;
// 	delete n;

// 	return factors;
// }

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


// AST* berlekamp(AST* ax, AST* x, AST* q)
// {
// 	long p = q->value();

// 	AST* Q = buildBerkelampMatrix(ax, x, q);

// 	// printf("--------------> %s\n", ax->toString().c_str());
// 	printf("%s\n", Q->toString().c_str());

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

	// while(factors->numberOfOperands() < k)
	// {
	// 	for(long idx = 0; idx < factors->numberOfOperands(); idx++)
	// 	{
	// 		AST* ux = factors->operand(idx)->copy();

	// 		AST* lx = polyFromList(v->operand(r), x);

	// 		printf("B\n");

	// 		for(long s = 0; s < p; s++)
	// 		{
	// 			AST* kx = sub({ lx->copy(), integer(s) });

	// 			AST* vx = reduceAST(kx);
				
	// 			delete kx;

	// 			AST* gx = gcdGPE_Zp(vx, ux, x, p);

	// 			delete vx;
				
	// 			if(gx->isNot(1) && !gx->match(ux))
	// 			{
	// 				factors->deleteOperand(idx);
	// 				AST* tx = quotientGPE_Zp(ux, gx, x, p);

	// 				delete ux;
	// 				ux = tx;

	// 				factors->includeOperand(ux->copy(), 0);
	// 				factors->includeOperand(gx->copy(), 1);
	// 				idx = 0;
	// 			}
	
	// 			delete gx;

	// 			if(factors->numberOfOperands() == k)
	// 			{
	// 				delete Q;
	// 				delete v;
	// 				delete ux;
	// 				delete lx;
	// 				return factors;
	// 			}
	// 		}
	// 		printf("C\n");

	// 		delete ux;
	// 		delete lx;
		
	// 		r = r + 1;
	// 	}
// 	}

// 	delete Q;
// 	delete v;
// 	return factors;
// }

}

