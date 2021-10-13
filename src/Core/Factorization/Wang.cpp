#include "Wang.hpp"
#include "Utils.hpp"
#include "SquareFree.hpp"
#include "Berlekamp.hpp"
#include "Zassenhaus.hpp"

#include "Core/Algebra/Set.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/GaloisField/GaloisField.hpp"

#include <limits>
#include <random>

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace galoisField;
using namespace simplification;

namespace factorization {

long gcd(long  a, long  b) {
	if (a == 0)
	{
		return b;
	}

	return gcd(b % a, a);
}

AST* nondivisors(long G, AST* F, long d, AST* L, AST* K)
{
	assert(G != 0, "G needs to be different from zero!");
	assert(d != 0, "c needs to be different from zero!");

	long i, j, k, q, r;

	AST *Fi;
	
	k = F->numberOfOperands();

	long* x = new long[k + 1];

	x[0] = G * d;

	for(i = 1; i <= k; i++)
	{
		Fi = F->operand(i - 1);
	
		q = norm(Fi, L, K);
	
		for(j = i - 1; j >= 0; j--)
		{
			while(r != 1)
			{
				r = x[j];

				r = gcd(r, q);
				q = q / r;
			}

			if(q == 1)
			{
				return list({});
			}
		}
	
		x[i] = q;
	}	

	AST* p = list({});

	for(i = 1; i <= k; i++)
	{
		p->includeOperand(integer(x[i]));
	}

	delete x;
	
	return p;
}

AST* groundLeadCoeff(AST* f, AST* L)
{
	long i = 0;
	
	AST* p = f->copy();

	AST* t = nullptr;
	
	for(i = 0; i < L->numberOfOperands(); i++)
	{
		t = leadCoeff(p, L->operand(i));
		
		delete p;
		
		p = t;
	}

	return p;	
}


ast::AST* trialDivision(ast::AST* f, ast::AST* F, ast::AST* L, ast::AST* K)
{
	AST *v, *q, *r, *d, *t = list({});

	bool stop = false;

	f = f->copy();

	long i, k;

	for(i = 0; i < F->numberOfOperands(); i++)
	{
		k = 0;
		v = F->operand(i);

		while(!stop)
		{
			d =	recPolyDiv(f, v, L, K);

			q = d->operand(0);
			r = d->operand(1);

			if(r->is(0))
			{
				delete f;
			
				f = q;
			
				k = k + 1;
			}

			stop = !r->is(0);

			delete d;
		}

		t->includeOperand(list({ v->copy(), integer(k) }));
	}

	delete f;

	return t;
}


AST* univariateFactors(AST* f, AST* L, AST* K)
{
	assert(L->numberOfOperands() == 1, "L needs to have just one variable");

	AST *ct, *pr, *n, *x, *lc, *t1, *T;

	x = L->operand(0);

	ct = cont(f, L, K);

	pr = pp(f, ct, L, K);

	n = degree(pr, x);
	
	lc = leadCoeff(f, x);

	if(lc->value() < 0)
	{
		t1 = mul({integer(-1), ct});
		ct = reduceAST(t1);
		
		delete t1;
		
		t1 = mul({integer(-1), pr});
		pr = reduceAST(t1);
		
		delete t1;
	}

	if(n->value() <= 0)
	{
		return list({ ct, list({ integer(1), integer(1) }) });
	}

	if(n->value() == 1)
	{
		return list({ ct, list({ pr, integer(1) }) });
	}

	t1 = pr;

	pr = squareFreePart(t1, L, K);

	delete t1;

	T = list({});
}

AST* sqf_factors(AST* f, AST* x, AST* K)
{
	AST *n, *cn, *pr, *lc, *L, *t1, *F;

	L = list({x->copy()});

	cn = cont(f, L, K);
	pr = pp(f, cn, L, K);

	lc = leadCoeff(pr, x);
	
	if(lc->value() < 0)
	{
		t1 = invert(cn);
		delete cn;
		cn = t1;

		t1 = invert(pr);
		delete pr;
		pr = t1;
	}

	delete lc;

	n = degree(pr, x);

	if(n->value() <= 0)
	{
		return list({cn, list({})});
	}

	if(n->value() == 1)
	{
		return list({cn, list({pr})});
	}

	F = zassenhaus(pr, x);

	delete n;
	delete pr;
	delete L;

	return list({cn, F});
}


AST* factors(AST* f, AST* L, AST* K)
{
	AST *x, *F, *c, *p, *T, *lc, *t1, *t2, *G, *g, *s, *n, *H, *R, *b, *e;

	x = L->operand(0);
	R = rest(L);

	if(L->numberOfOperands() == 1)
	{

		c = cont(f, L, K);
		p = pp(f, c, L, K);
	
		F = zassenhaus(p, x);
		T = trialDivision(p, F, x, K);

		return list({c, T});
	}

	if(f->is(0))
	{
		return list({integer(0), list({integer(0), integer(1)})});
	}

	c = cont(f, L, K);
	p = pp(f, c, L, K);

	lc = groundLeadCoeff(p, L);
	
	if(lc->value() < 0)
	{
		t1 = integer(-1);
		t2 = mulPoly(c, t1);
		delete c;
		c = reduceAST(t2);
		delete t2;
	
		t2 = mulPoly(p, t1);
		delete p;
		p = reduceAST(t2);
		delete t2;
	
		delete t1;
	}


	G = cont(p, L, K);
	g = pp(p, G, K, K);

	F = list({});

	n = degree(g, x);

	if(n->value() > 0)
	{
		s = squareFreePart(g, L, K);

		H = factorsWang(g, L, K);
	
		delete F;
	
		F = trialDivision(f, H, L, K);
	}

	t1 = factors(G, R, K);
	
	while(t1->numberOfOperands())
	{
		b = t1->operand(0);
		e = t1->operand(1);
	
		F->includeOperand(list({b, e}), 0L);
	}

	return list({c, F});
}

AST* invert(AST* p)
{
	AST *t1, *t2, *t3;
	
	t1 = integer(-1);
	
	t2 = mulPoly(p, t1);
	
	delete t1;

	t3 = reduceAST(t2);

	delete t2;

	return t2;
}

AST* eval(AST* f, AST* L, AST* a, long j)
{
	AST* k;

	AST* p = deepReplace(f, L->operand(j), a->operand(0));

	for(long i = j + 1; i < L->numberOfOperands(); i++)
	{
		k = deepReplace(p, L->operand(i), a->operand(i - j));
		delete p;
		p = k;
	}

	k = reduceAST(p);

	delete p;

	return k;
}

AST* testEvaluationPoints(AST* U, AST* G, AST* F, AST* a, AST* L, AST* K)
{
	assert(G->kind() == Kind::Integer, "Gamma parameter needs to be an integer");

	long i;

	AST *x, *V, *U0, *g, *delta, *pr, *lc, *t1, *t2, *R, *E, *d;
	
	x = L->operand(0);
	
	R = rest(L);

	assert(R->numberOfOperands() == a->numberOfOperands(), "Wrong numbers of test points");

	// 	Test Wang condition 1: V(a1, ..., ar) = lc(f(x))(a1,...,ar) != 0
	V = leadCoeff(U, x);
	g = eval(V, R, a, 1);
	
	if(g->is(0))
	{
		delete g;
		delete V;
		delete R;

		return fail();
	}
	
	delete g;
	delete V;

	// Test Wang condition 3: U0(x) = U(x, a1, ..., at) is square free
	U0 = eval(U, L, a, 1);
	
	if(!isSquareFree(U0, x, K))
	{
		delete U0;
		delete R;

		return fail();
	}

	// Test Wang condition 2: For each F[i], E[i] = F[i](a1, ..., ar) 
	// has at least one prime division p[i] which does not divide 
	// any E[j] j < i, Gamma, or the content of U0
	delta = cont(U0, L, K);
	pr = pp(U0, delta, L, K);

	lc = groundLeadCoeff(pr, L);

	if(lc->value() < 0)
	{
		t1 = invert(delta);
		delete delta;
		delta = t1;

		t1 = invert(pr);
		delete pr;
		pr = t1;
	}

	E = list({});

	for(i = 0; i < F->numberOfOperands(); i++)
	{
		E->includeOperand(eval(F->operand(i), R, a, 0));
	}
	
	d = nondivisors(G->value(), E, delta->value(), R, K);

	if(d->numberOfOperands() == 0)
	{
		delete d;

		return fail();
	}

	// return the content of U0, the primitive part of U0 and the E[i] = F[i](a1, ... ar)
	// paper defines delta = cn, and u(x) = pp(U0)
	return list({ delta, pr, E });
}

AST* getEvaluationPoints(ast::AST* f, ast::AST* G, ast::AST* F, ast::AST* L, ast::AST* K, long p)
{
	long i, t, j, t, r, r_;

	r_ = -1;
	r  = -1;

	AST *ux, *R, *x, *cn, *pr, *g, *lc, *t1, *t2, *Fi, *E, *d, *a, *s, *t3, *delta, *pr_u0;

	t = L->numberOfOperands() - 1;

	a = list({});

	for(i = 0; i < t; i++)
	{
		a->includeOperand(integer(0));
	}

	AST* c = set({});

	x = L->operand(0);

	while(c->numberOfOperands() < 3)
	{
		for(t = 0; t < 5; t++)
		{
			for(i = 0; i < t; i++)
			{
				a->deleteOperand(0L);
				a->includeOperand(integer(mod(random(), p, true)), 0L);
			}

			s = testEvaluationPoints(f, G, F, a, L, K);
			
			if(s->kind() == Kind::Fail)
			{
				delete s;
				continue;
			}

			delta = s->operand(0);
			pr_u0 = s->operand(1);
			E     = s->operand(2);
			
			ux = sqf_factors(pr_u0, x, K);
			
			cn = ux->operand(0);
			pr = ux->operand(1);

			// Verify that the sets a[i] are
			// given same low r value
			r_ = pr->numberOfOperands();
	
			if(r == -1)
			{
				r = r_;
			

				delete c;
	
				c = set({
					list({ 
						delta->copy(), // paper delta
						pr_u0->copy(), // paper pr(U0)
						E->copy(), 	   // paper ~F[i]
						pr->copy(),		 // paper u[i](x)*...*u[r](x)
						a->copy()			 // paper a[i]
					})
				});

				continue;
			}


			// If this config leads to a smaller number of factors
			// than the current configurations, erase current configs
			if(r_ < r)
			{
				r = r_;

				delete c;
				c = set({});
			}

			// If this config leads to the same number of factors
			// save it
			if(r_ == r)
			{
				t2 = list({ 
					delta->copy(), // paper delta
					pr_u0->copy(), // paper pr(U0)
					E->copy(), 	   // paper ~F[i]
					pr->copy(),		 // paper u[i](x)*...*u[r](x)
					a->copy()			 // paper a[i]
				});
		
				t1 = set({ t2 });
				t3 = unification(c, t1);

				delete t1;
				
				c = t3;
			}

			if(c->numberOfOperands() < 3)
			{
				break;
			}
		}

		p = p + 1;
	}

	return c;
}

AST* wangLeadingCoeff(AST* U, AST* delta, AST* u, AST* F, AST* sF, AST* a, AST* L, AST* K)
{

	/**
	 * From the Wang's paper:
	 * 
	 * If none of u[1](x),...,u[r](x) is extraneous, then U factors into r distinct 
	 * irreductible polynomials U = prod i = 1 to r G[i](x[0], ..., x[t]).
	 * 
	 * Let C[i](x2, ..., x[t]) = lc(G[i]), ~C[i] = C[i](a[1], ..., a[t - 1]),
	 * and G[i](x, a[1], ..., a[t - 1]) = delta[i] * u[i] where delta[i]
	 * is some divisor of delta.
	 * 
	 * Lemma: If there are no extraneous factors, then for all i and m, F[k]^m
	 * divides C[i], then ~C[i] = ~F[1]^s1 * ... * F[k]^s[k]*w where w | G, s[i] >= 0
	 * and s[k] < m. Thus p[k]^m dows not divided ~C[i], which implies that ~F[k]^m 
	 * does not divide lc(u[i])*delta 
	 * 
	 * This lemma enables one to distribute all F[k] first, then all F[k-1], etc.
	 * Thus D[i](x[2], ... ,x[t]) can be determined as products of powers of F[i],
	 * Now let ~D[i] = D[i](a2, ..., a[t]). If delta == 1, then C[i] = lc(u[i]/~D[i])*D[i].
	 * Otherwise, if delta != 1, the following steps are carried out for all i = 1, ..., r to
	 * correctly distribute the factors of delta.
	 * 
	 * 	1. Let d = gcd(lc(u[i]), ~D[i]) and C[i] = D[i] * lc(u[i]) / d
	 *  2. Let u[i] = (~D[i] / d) * u[i]
	 *  3. Let delta = delta / (D[i] / d)
	 *	
	 * The process ends when delta = 1, Otherwise, let u[i] = delta * u[i], C[i] = delta * C[i], u = delta^(r - 1) * U.
	 * In this case, when the true factors over Z of U are found, they may have integer contents which should be removed.
	 * 
	 * In the above process, if any factors of V[n] is not distributed, then there are extraneous factors, and the program
	 * goes back for different substitutions that lower r. 
	 */

	bool extraneous = false;

	long i, d, j, m, k;

	AST *x, *R, *Ci, *C, *ui, *sFk, *Fk, *lc_ui;
	
	x = L->operand(0);
	R = rest(L);

	// 1. Distribute C[i]
	C = list({});
	
	bool* was_set = new bool[sF->numberOfOperands()];

	for(i = 0; i < sF->numberOfOperands(); i++)
	{
		was_set[i] = false;
	}

	for(i = 0; i < u->numberOfOperands(); i++)
	{
		ui = u->operand(i);
		lc_ui = leadCoeff(ui, x);
		
		Ci = integer(1);

		/**
		 * Aplying lemma: It there are no extraneous factors, then,
		 * for all i and m, F[k]^m divides C[i] if and only if ~F[k]^m
		 * divides lc(u[i])*delta
		 * 
		 */
		d  = lc_ui->value() * delta->value();
		
		delete lc_ui;

		// Distribute F[k], then k - 1, ...
		for(k = sF->numberOfOperands() - 1; k >= 0; k--)
		{
			m   = 0;

			sFk = sF->operand(k);
			
			// find expoent m
			while(d % sFk->value() == 0)
			{
				d = d / sFk->value();
				m = m + 1;
			}

			Fk  = F->operand(k)->operand(0);

			if(m != 0)
			{
				Ci = mul({ Ci, power(Fk->copy(), integer(m)) });
				was_set[k] = true;
			}
		}
	
		lc_ui = reduceAST(Ci);
	
		delete Ci;
	
		Ci = lc_ui;

		C->includeOperand(Ci);
	}

	extraneous = false;

	for(int i = 0; i < sF->numberOfOperands(); i++)
	{
		if(!was_set[i])
		{
			// Extraneous factors found, stop and try again for another set of a[i]'s
			extraneous = true;
			break;
		}
	}

	delete was_set;

	if(extraneous)
	{
		return fail();
	}	


}

AST* factorsWangRec(AST* f, AST* L, AST* K, long mod)
{
	AST* x, *lc, *R, *H, *G, *Vn;
	
	// First step: factor lc(f)
	x  = L->operand(0);
	lc = leadCoeff(f, x);
	
	R = rest(L);

	H = factors(lc, R, K);

	G  = H->operand(0L);
	Vn = H->operand(1L);

	H->removeOperand(0L);
	H->removeOperand(0L);

	delete H;

	// Second step: find integers a1, ..., ar
	AST* S = getEvaluationPoints(f, G, Vn, L, K, mod);
	
	// todo: set c to config where pp_u0 have largest max norm
	// or tests all configs on wang leading coeff and verify if
	// one works
	AST* c = S->operand(0);

	AST* delta = c->operand(0);
	AST* pp_u0 = c->operand(1);
	AST* sF    = c->operand(2);
	AST* u   	 = c->operand(3);
	AST* a     = c->operand(4);
	
	AST* k = wangLeadingCoeff(f, delta, u, Vn, sF, a, L, K);
	
	
	return nullptr;
}

AST* factorsWang(AST* f, AST* L, AST* K)
{
	return factorsWangRec(f, L, K, 3);
}

}
