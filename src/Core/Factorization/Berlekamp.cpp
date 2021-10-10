#include "Berlekamp.hpp"

#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/GaloisField/GaloisField.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace galoisField;
using namespace simplification;

namespace factorization {

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

AST* buildBerlekampBasis(AST* A, AST* w, bool symmetric)
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

	lead = 0;

	row_count = n;
	col_count = n;

	for(r = 0; r < row_count; r++)
	{
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
				long Mrj = quoGf(v, x, q, symmetric);

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

void addVecToBasisMatrix(AST* Q, AST* r, AST* x, long i, long n)
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

AST* buildBerkelampMatrix(AST* ax, AST* x, AST* p, bool symmetric)
{
	AST *n, *Q, *r0, *rx, *zx;

	n = degree(ax, x);

	Q = initBerkelampBasisMatrix(n);

	r0 = powModPolyGf(x, ax, x, p->value(), p->value(), symmetric);

	addVecToBasisMatrix(Q, r0, x, 1, n->value());

	rx = r0->copy();

	for(long m = 2; m < n->value(); m++)
	{
		zx = mulPolyGf(r0, rx, x, p->value(), symmetric);

		delete rx;
	
		rx = remPolyGf(zx, ax, x, p->value(), symmetric);
	
		delete zx;

		addVecToBasisMatrix(Q, rx, x, m, n->value());
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

AST* berlekampFactors(AST* sfx, AST* x, AST* p, bool symmetric)
{
	long i, k, s, r;

	AST *Q, *B, *lc, *fx, *n, *H, *F, *h, *v, *g, *f;

	lc = leadCoeff(sfx, x);

	fx = quoPolyGf(sfx, lc, x, p->value(), symmetric);

	n = degree(fx, x);

	Q = buildBerkelampMatrix(fx, x, p, symmetric);
	B = buildBerlekampBasis(Q, p, symmetric);

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
			g = integer(s);
		
			v = subPolyGf(h, g, x, p->value(), symmetric);
		
			delete g;
	
			g = gcdPolyGf(f, v, x, p->value(), symmetric);

			delete v;

			if(g->isNot(1) && !g->match(f))
			{
				v = quoPolyGf(f, g, x, p->value(), symmetric);

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


	delete n;
	delete H;

	return F;
}

}

