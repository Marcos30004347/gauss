
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

AST* invert(AST* p)
{
	AST *t1, *t2, *t3;

	t1 = integer(-1);

	t2 = mulPoly(p, t1);

	delete t1;

	t3 = reduceAST(t2);

	delete t2;

	return t3;
}

AST* nondivisors(Int G, AST* F, Int d, AST* L, AST* K)
{
	assert(G != 0, "G needs to be different from zero!");
	assert(d != 0, "d needs to be different from zero!");

	long i, j, k;

	Int q, r;

	AST *Fi;

	k = F->numberOfOperands();

	Int* x = new Int[k + 1];

	x[0] = G * d;
	for(i = 1; i <= k; i++)
	{
		Fi = F->operand(i - 1);

		q = norm(Fi, L, K);

		for(j = i - 1; j >= 0; j--)
		{
			r = x[j];

			while(abs(r) != 1)
			{
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

	delete[] x;

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

Int groundContRec(AST* f, AST* L, AST* K)
{
	if(f->kind() == Kind::Integer)
	{
		return f->value();
	}

	AST* p = f->copy();

	Int g = 0;

	AST *r, *u, *e, *t, *x, *R;

	x = L->operand(0);

	R = rest(L);
	
	AST* d = degree(f, x);

	while(!p->is(0))
	{
		t = leadCoeff(p, x);

		g = gcd(g, groundContRec(t, R, K));

		e = power(x->copy(), degree(p, x));
		u = mulPoly(t, e);
		t = subPoly(p, u);

		delete e;
		delete u;
		delete p;

		p = t;	
	}

	delete R;
	delete d;

	return g;
}

AST* groundCont(AST* f, AST* L, AST* K)
{
	return integer(groundContRec(f, L, K));
}

AST* groundPP(AST* f, AST* L, AST* K)
{
	AST* c = groundCont(f, L, K);
	AST* p = recQuotient(f, c, L, K);

	delete c;

	return p;
}

AST* groundPP(AST* f, AST* c, AST* L, AST* K)
{
	AST* p = recQuotient(f, c, L, K);

	return p;
}

AST* trialDivision(AST* f, AST* F, AST* L, AST* K)
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

AST* sqfFactors(AST* f, AST* x, AST* K)
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
		return list({cn, list({ pr })});
	}

	F = zassenhaus(pr, x);

	delete n;
	delete pr;
	delete L;

	return list({cn, F});
}

AST* factors(AST* f, AST* L, AST* K)
{
	if(f->is(0))
	{
		return list({integer(0), list({integer(0), integer(1)})});
	}

	AST *x, *F, *c, *p, *T, *lc, *t1, *t2, *G, *g, *s, *n, *H, *R, *b, *e;

	x = L->operand(0);
	R = rest(L);

	AST* cnt = groundCont(f, L, K);
	AST* prp = groundPP(f, cnt, L, K);

	printf("%s\n", cnt->toString().c_str());
	printf("%s\n", prp->toString().c_str());

	lc = groundLeadCoeff(prp, L);

	if(lc->value() < 0)
	{
		t1 = integer(-1);

		t2 = mulPoly(cnt, t1);
		delete cnt;
		cnt = t2;

		t2 = mulPoly(prp, t1);
		delete prp;
		prp = t2;

		delete t1;
	}

	bool is_const = true;

	for(Int i = 0; i < L->numberOfOperands() && is_const; i++)
	{
		AST* d = degree(prp, L->operand(i));
	
		if(!d->is(0))
		{
			is_const = false;
		}
	
		delete d;
	}

	if(is_const)
	{
		return list({ cnt, list({integer(0), integer(1)}) });
	}

	// if(L->numberOfOperands() == 1)
	// {
	// 	c = cont(f, L, K);
	// 	p = pp(f, c, L, K);

	// 	F = zassenhaus(p, x);
	// 	T = trialDivision(p, F, x, K);

	// 	return list({c, T});
	// }

	// c = cont(f, L, K);
	// p = pp(f, c, L, K);

	G = cont(prp, x);
	g = pp(prp, G, x);

	F = list({});

	n = degree(g, x);

	if(n->value() > 0)
	{
		s = squareFreePart(g, L, K);

		H = factorsWang(s, L, K);

		delete F;

		F = trialDivision(f, H, L, K);

		delete s;
	}

	t1 = factors(G, R, K);

	while(t1->numberOfOperands())
	{
		b = t1->operand(0);
		e = t1->operand(1);

		F->includeOperand(list({b, e}), 0L);
	}

	return list({cnt, F});
}



AST* eval(AST* f, AST* L, AST* a, long j = 0)
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

	AST *x, *V, *U0, *g, *delta, *pr, *lc, *t1, *R, *E, *d;

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

	if(delta->value() == 0)
	{
		return fail();
	}

	printf("NON\n");

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

Int degreeSum(AST* f, AST* L)
{
	AST* n;

	Int s = 0;

	for(long i=0; i < L->numberOfOperands(); i++)
	{
		n = degree(f, L->operand(i));

		if(n->kind() == Kind::Integer)
		{
			s += n->value();
		}

		delete n;
	}

	return s;
}

Int mignotteBound(AST* f, AST* L, AST* K)
{
	AST* l = groundLeadCoeff(f, L);

	Int a = norm(f, L, K);
	Int b = l->value();
	Int n = degreeSum(f, L);

	return Int(std::sqrt(n.longValue() + 1) * std::pow(2, n.longValue()) * a.longValue() * b.longValue());
}

// return l, such that p^l is a bound to the coefficients of the factors of f in K[L]
Int mignoteExpoent(AST* f, AST* L, AST* K, Int p)
{
	return std::ceil(std::log(2*mignotteBound(f, L, K).longValue() + 1) / std::log(p.longValue()));
}

AST* getEvaluationPoints(AST* f, AST* G, AST* F, AST* L, AST* K, Int p)
{
	Int i, t, r, r_;

	r_ = -1;
	r  = -1;

	AST *ux, *x, *pr, *t1, *t2, *E, *a, *s, *t3, *delta, *pr_u0;

	t = L->numberOfOperands() - 1;

	AST* c = set({});

	x = L->operand(0);

	while(c->numberOfOperands() < 3)
	{
		for(t = 0; t < 5; t++)
		{
			a = list({});

			for(i = 0; i < L->numberOfOperands() - 1; i++)
			{
				a->includeOperand(integer(mod(random(), p, true)));
			}

			// for(i = 0; i < t; i++)
			// {
			// 	printf("AQUI\n");
			// 	a->deleteOperand(0L);
			// 	a->includeOperand(integer(mod(random(), p, true)));
			// }
	
			printf("a = %s\n", a->toString().c_str());
	
			s = testEvaluationPoints(f, G, F, a, L, K);

			if(s->kind() == Kind::Fail)
			{
				delete s;
				delete a;
	
				continue;
			}

			delta = s->operand(0);
			pr_u0 = s->operand(1);
			E     = s->operand(2);

			ux = sqfFactors(pr_u0, x, K);

			// cn = ux->operand(0);
			pr = ux->operand(1);

			ux->removeOperand(1);
			delete ux;

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
			
				delete a;
			
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
				printf("UNI\n");
				delete t1;

				c = t3;
			}

			delete a;

			if(c->numberOfOperands() < 3)
			{
				break;
			}
		}

		p = p + 3;
	}

	return c;
}

AST* wangLeadingCoeff(AST* f, AST* delta, AST* u, AST* F, AST* sF, AST* a, AST* L, AST* K)
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

	Int i, m, k;

	Int d, di, dt;

	AST *x, *R, *Di, *D, *ui, *sFk, *C, *ci, *Fk, *lc, *sDi, *U, *t1, *t2;

	x = L->operand(0);
	R = rest(L);

	D = list({});
	C = list({});
	U = list({});

	// 1. Distribute D[i]
	bool* was_set = new bool[sF->numberOfOperands()];
	printf("A\n");
	for(i = 0; i < sF->numberOfOperands(); i++)
	{
		was_set[i.longValue()] = false;
	}

	for(i = 0; i < u->numberOfOperands(); i++)
	{
		printf("B\n");

		ui = u->operand(i);
		printf("ui = %s\n", ui->toString().c_str());
		lc = leadCoeff(ui, x);

		Di = integer(1);

		/**
		 * Aplying lemma: It there are no extraneous factors, then,
		 * for all i and m, F[k]^m divides C[i] if and only if ~F[k]^m
		 * divides lc(u[i])*delta
		 *
		 */
		d  = lc->value() * delta->value();

		printf("lc = %s\n", lc->toString().c_str());
		printf("delt = %s\n", delta->toString().c_str());
		delete lc;
		printf("d = %s\n", d.to_string().c_str());

		// Distribute F[k], then k - 1, ...
		for(k = sF->numberOfOperands() - 1; k >= 0; k--)
		{
			printf("Z\n");

			m   = 0;

			sFk = sF->operand(k);
	
			printf("%s\n", sFk->toString().c_str());
			printf("%s\n", d.to_string().c_str());

			// find expoent m
			while(d % sFk->value() == 0)
			{
				d = d / sFk->value();
				m = m + 1;
			}

			Fk  = F->operand(k)->operand(0);

			if(m != 0)
			{
				Di = mul({ Di, power(Fk->copy(), integer(m)) });
				was_set[k.longValue()] = true;
			}
		}

		lc = reduceAST(Di);

		delete Di;

		Di = lc;

		D->includeOperand(Di);
	}
	printf("C\n");

	extraneous = false;

	for(i = 0; i < sF->numberOfOperands(); i++)
	{
		if(!was_set[i.longValue()])
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

	printf("D\n");

	dt = delta->value();

	// otherwise, if delta != 1, the following steps
	// are carried out, for all i = 1, ..., r, to
	// correctly distribute the factors of delta
	for(i = 0; i < D->numberOfOperands(); i++)
	{
		ui = u->operand(i)->copy();
		Di = D->operand(i)->copy();

		lc = leadCoeff(ui, x);
		sDi = eval(Di, R, a, 0);
		// assert(sDi->kind() == Kind:::Integer);
		// assert(lc->kind() == Kind:::Integer);

		di = sDi->value();
		
		printf("-->1\n");

		// if delta == 1, Then Ci = (lc(ui) / sD[i])*D[i]
		if(delta->is(1))
		{
			printf("KKKKKKK\n");
			printf("%s\n", lc->toString().c_str());
			printf("%s\n", Di->toString().c_str());

			ci = integer(lc->value() / sDi->value());
		}
		else
		{
	  	// * 	1. Let d = gcd(lc(u[i]), ~D[i]) and C[i] = D[i] * lc(u[i]) / d
			printf("-->AAA\n");
		
			d = gcd(lc->value(), di);

			ci = integer(lc->value() / d);
			printf("-->BBB\n");

			di = di / d;

			// *  2. Let u[i] = (~D[i] / d) * u[i]
			t1 = integer(di);
			t2 = mulPoly(ui, t1);
			printf("-->CCC\n");

			delete ui;

			ui = t2;

			delete t1;

	 	 	// *  3. Let delta = delta / (D[i] / d)
			dt = dt / di;
		}
		printf("-->2\n");

		// Ci = (lc(ui)/d) * Di = ci * Di
		t1 = mulPoly(ci, Di);
		delete ci;
		ci = t1;

		C->includeOperand(ci);
		U->includeOperand(ui);
	}

	printf("E\n");

	// Now if dt == 1, the process ends
	if(dt == 1)
	{
		return list({f->copy(), U, C});
	}

	// otherwise, let ui = delta*ui, Ci = delta*Ci, U = delta^(r - 1) * U
	t1 = integer(dt);
	printf("F\n");

	for(i = 0; i < C->numberOfOperands(); i++)
	{
		ui = U->operand(i);
		ci = C->operand(i);

		U->removeOperand(i);
		U->includeOperand(mulPoly(ui, t1), i);

		C->removeOperand(i);
		C->includeOperand(mulPoly(ci, t1), i);

		delete ui;
		delete ci;
	}

	dt = pow(dt, u->numberOfOperands() - 1);
	t1 = integer(dt);
	t2 = mulPoly(f, t1);

	return list({ t2, U, C  });
}


AST* EEAlift(AST* a, AST* b, AST* x, Int p, Int k)
{
	Int j, modulus;
	AST *t1, *t2, *t3, *t4, *t5, *_sig, *_tal, *tal, *sig;

	AST *amodp, *bmodp, *smodp, *tmodp, *s, *t, *G, *e, *c, *mod, *q;

	amodp = gf(a, x, p);
	bmodp = gf(b, x, p);

	G = extendedEuclidGf(amodp, bmodp, x, p);

	s = G->operand(1);
	t = G->operand(2);

	G->removeOperand(2);
	G->removeOperand(1);

	delete G;

	smodp = s->copy();
	tmodp = t->copy();

	modulus = p;

	for(j = 1; j <= k - 1; j++)
	{
		t1 = integer(1);
		t2 = mulPoly(s, a);
		t3 = mulPoly(t, b);
		t4 = subPoly(t1, t2);
		t5 = subPoly(t4, t3);

		e = reduceAST(t5);

		delete t1;
		delete t2;
		delete t3;
		delete t4;
		delete t5;

		mod = integer(modulus);

		c = quoPolyGf(e, mod, x, p, true);

		_sig = mulPoly(smodp, c);
		_tal = mulPoly(tmodp, c);

		t3 = divideGPE(_sig, bmodp, x);

		q = t3->operand(0);

		sig = t3->operand(1);

		t3->removeOperand(0L);
		t3->removeOperand(0L);

		delete t3;

		t3 = mulPoly(q, amodp);
		t5 = addPoly(_tal, t3);

		tal = gf(t5, x, p, true);

		delete t3;
		delete t5;

		t3 = mulPoly(sig, mod);
		t5 = addPoly(s, t3);
		delete t3;

		delete s;
		s = t5;

		t3 = mulPoly(tal, mod);
		t5 = addPoly(t, t3);
		delete t3;
		delete t;

		t = t5;

		delete mod;

		modulus = modulus * p;

		delete _sig;
		delete _tal;
		delete sig;
		delete tal;
		delete mod;
		delete c;
	}

	delete smodp;
	delete tmodp;
	delete amodp;
	delete bmodp;

	printf("s, t = [%s, %s]\n", s->toString().c_str(), t->toString().c_str());

	return list({s, t});
}

AST* multiTermEEAlift(AST* a, AST* L, Int p, Int k)
{
	long j, r;

	AST *q, *t1, *t2, *s, *bet, *sig;

	r = a->numberOfOperands();

	q = list({ a->operand(r - 1)->copy() });

	printf("G = %s\n", q->toString().c_str());
	printf("a = %s\n", a->toString().c_str());

	for(j = r - 2; j >= 1; j--)
	{
		t1 = mulPoly(a->operand(j), q->operand(0L));
		q->includeOperand(reduceAST(t1), 0L);
		delete t1;
	}

	bet = list({ integer(1) });

	t2 = list({});
	s = list({});

	printf("G = %s\n", q->toString().c_str());

	for(j = 0; j < r - 1; j++)
	{
		printf("================================\n");
		t1 = list({ q->operand(j)->copy(), a->operand(j)->copy() });

		printf("list = %s\n", t1->toString().c_str());
		printf("%s, %s\n", a->toString().c_str(), q->toString().c_str());

		sig = multivariateDiophant(t1, bet->operand(bet->numberOfOperands() - 1), L, t2, 0, p, k);

		bet->includeOperand(sig->operand(0));
		s->includeOperand(sig->operand(1));

		sig->removeOperand(1L);
		sig->removeOperand(0L);

		delete t1;
		delete sig;
	}

	s->includeOperand(bet->operand(r - 1)->copy());

	delete q;
	delete t2;
	delete bet;

	printf("aaaa\n");
	return s;
}

AST* replaceAndReduce(AST* f, AST* x, AST* a)
{
	AST* g = deepReplace(f, x, a);
	AST* r = reduceAST(g);

	delete g;

	return r;
}

AST* diff(AST* f, Int j, AST* x)
{
	AST *t, *g = f->copy();

	for(int i = 0; i < j; i++)
	{
		t = derivate(g, x);
		delete g;
		g = t;
	}

	return g;
}

/**
 * @brief Find the coefficient of the taylor expansion of f in the variable L[j] at a.
 *
 * @param f A polynomial expresiion in Z[L]
 * @param m The order of the derivative
 * @param j The index of the variable in L
 * @param L The list of symbols in f
 * @param a Value that taylor should be taken.
 * @return The coefficient of f in the Taylor expansion of e about L[j] = a
 */
AST* taylorExpansionCoeffAt(AST* f, Int m, Int j, AST* L, AST* a)
{
	AST* g = diff(f, m, L->operand(j));
	AST* t = replaceAndReduce(g, L->operand(j), a);


	AST* n = integer(fact(m));
	AST* K = symbol("Z");
	printf("m = %li\n", m.longValue());
	printf("j = %li\n", j.longValue());
	printf("f = %s\n", f->toString().c_str());
	printf("f' = %s\n", g->toString().c_str());
	printf("C = %s\n", t->toString().c_str());
	AST* q = recQuotient(t, n, L, K);

	delete g;
	delete n;
	delete t;
	delete K;

	return q;
}

AST* multivariateDiophant(AST* a, AST* c, AST* L, AST* I, Int d, Int p, Int k)
{
	long long i, j;

	Int m, v, r;

	AST *K, *x1, *ds, *monomial, *cm, *e, *sig, *R, *xv, *av, *A, *t1, *t2, *t3, *b, *anew, *Inew, *cnew;

	// 1. Initialization
	r = 		a->numberOfOperands();
	v = 1 + I->numberOfOperands();

	if(v > 1)
	{
		K = symbol("Z");

		printf("FIRST OPTION\n");

		xv = L->operand(L->numberOfOperands() - 1);
		av = I->operand(I->numberOfOperands() - 1);

		// 2.1. Multivariate case
		A = integer(1);

		for(i = 0; i < r; i++)
		{
			t1 = mulPoly(A, a->operand(i));
			delete A;
			A = t1;
		}

		t1 = reduceAST(A);

		delete A;

		A = t1;

		b = list({});
		anew = list({});

		for(j = 0; j < r; j++)
		{
			b->includeOperand(recQuotient(A, a->operand(j), L, K));
		}

		for(j = 0; j < a->numberOfOperands(); j++)
		{
			anew->includeOperand(replaceAndReduce(a->operand(j), xv, av));
		}

		cnew = replaceAndReduce(c, xv, av);

		printf("%s\n", c->toString().c_str());

		Inew = I->copy();
		Inew->deleteOperand(Inew->numberOfOperands() - 1);

		R = L->copy();
		R->deleteOperand(R->numberOfOperands() - 1);

		printf("CALLING MULTIVARIATEDIOPHANT:\n");
		printf("	%s\n", anew->toString().c_str());
		printf("	%s\n", cnew->toString().c_str());
		printf("	%s\n", Inew->toString().c_str());

		sig = multivariateDiophant(anew, cnew, R, Inew, d, p, k);

		printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
		printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");
		printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n");

		printf("S = %s\n", sig->toString().c_str());

		delete R;

		t1 = integer(0);

		for(j = 0; j < sig->numberOfOperands(); j++)
		{
			t2 = mulPoly(sig->operand(j), b->operand(j));
			t3 = addPoly(t1, t2);

			delete t2;
			delete t1;

			t1 = t3;
		}

		t2 = subPoly(c, t1);

		delete t1;

		e = gf(t2, pow(p, k), true);

		delete t2;

		monomial = integer(1);

		for(m = 1; m < d; m++)
		{
			if(e->is(0)) break;

			t1 = subPoly(xv, av);
			t2 = mulPoly(monomial, t1);

			delete t1;
			delete monomial;

			monomial = t2;

			cm = taylorExpansionCoeffAt(e, m, v - 1, L, av);

			printf("###############################\n");
			printf("cm = %s\n", cm->toString().c_str());

			if(cm->isNot(0))
			{
				ds = multivariateDiophant(anew, cm, L, Inew, d, p, k);

				for(j = 0; j < ds->numberOfOperands(); j++)
				{
					t1 = ds->operand(j);
					ds->removeOperand(j);

					t2 = mulPoly(t1, monomial);

					ds->includeOperand(t2, j);

					delete t1;
				}

				for(j = 0; j < ds->numberOfOperands(); j++)
				{
					t1 = ds->operand(j);

					t2 = sig->operand(j);

					sig->removeOperand(j);

					t3 = addPoly(t1, t2);

					sig->includeOperand(t3, j);

					delete t2;
				}

				t1 = integer(0);
				for(j = 0; j < r; j++)
				{
					t2 = mulPoly(ds->operand(j), b->operand(j));
					t3 = addPoly(t1, t2);

					delete t1;
					delete t2;

					t1 = t3;
				}

				t2 = subPoly(e, t1);

				delete t1;

				delete e;

				e = gf(t2, pow(p, k), true);

				delete t2;
				delete ds;
			}

			delete cm;
		}

		delete e;
		delete monomial;

		delete K;
		delete A;
		delete b;
		delete anew;
		delete cnew;
		delete Inew;
	}
	else
	{
		printf("SECOND OPTION\n");

		x1 = L->operand(0);

		sig = list({});
		for(j = 0; j < r; j++)
		{
			sig->includeOperand(integer(0));
		}

		printf("**************************\n");
		printf("%s\n", c->toString().c_str());

		AST* C = c->copy();

		while(C->isNot(0))
		{
			t1 = degree(C, x1);
			m = t1->value();
			cm = leadCoeff(C, x1);

			delete t1;
			printf("c = %s\n", c->toString().c_str());
			printf("z = %s\n", cm->toString().c_str());

			printf("****** deg, cm = %li %s\n", m, cm->toString().c_str());
			ds = univariateDiophant(a, L, m, p, k);
			printf("%s\n", sig->toString().c_str());

			for(i = 0; i < ds->numberOfOperands(); i++)
			{
				printf("--> (%s) * (%s)\n", ds->operand(i)->toString().c_str(), cm->toString().c_str());

				t2 = mulPoly(ds->operand(i), cm);

				t3 = addPolyGf(sig->operand(i), t2, x1, pow(p, k), true);

				printf("--> (%s) + (%s)\n", sig->operand(i)->toString().c_str(), t2->toString().c_str());

				delete t2;

				sig->deleteOperand(i);
				sig->includeOperand(t3, i);
			}

			t1 = mul({
				cm->copy(),
			 	power(x1->copy(), integer(m))
			});

			t2 = subPoly(C, t1);

			delete C;

			C = reduceAST(t2);

			delete t1;
			delete t2;
			delete cm;
			delete ds;
		}

		delete C;
	}

	for(j = 0; j < sig->numberOfOperands(); j++)
	{
		t2 = sig->operand(j);

		sig->removeOperand(j);

		sig->includeOperand(gf(t2, pow(p, k), true), j);

		delete t2;
	}

	return sig;
}

AST* univariateDiophant(AST* a, AST* L, Int m, Int p, Int k)
{
	AST *x, *s, *t1, *t2, *t3, *t4, *result, *u, *v;
	printf("UNIVARIATE DIOPHANTINE\n");

	x = L->operand(0);

	long long r, j;

	r = a->numberOfOperands();

	result = list({});

	if(r > 2)
	{
		printf("r > 2\n");
		s = multiTermEEAlift(a, L, p, k);

		printf("S = %s\n", s->toString().c_str());
		printf("F = %s\n", a->toString().c_str());
		printf("%li\n", r);

		for(j = 0; j < r; j++)
		{
			t1 = power(x->copy(), integer(m));
			t2 = mulPoly(s->operand(j), t1);

			result->includeOperand(remPolyGf(t2, a->operand(j), x, pow(p, k), true));

			delete t1;
			delete t2;
		}

		printf("RESULT = %s\n", result->toString().c_str());
		delete s;
	}
	else
	{
		printf("len == 2\n");
		printf("EEAlift [%s, %s]\n", a->operand(1)->toString().c_str(),a->operand(0)->toString().c_str() );

		s = EEAlift(a->operand(1), a->operand(0), x, p, k);

		t1 = power(x->copy(), integer(m));

		t2 = mulPoly(s->operand(0), t1);

		t3 = divPolyGf(t2, a->operand(0), x, pow(p, k), true);

		u = t3->operand(0);
		v = t3->operand(1);

		t3->removeOperand(0L);
		t3->removeOperand(0L);

		delete t1;
		delete t2;
		delete t3;

		t1 = power(x->copy(), integer(m));
		t2 = mulPoly(s->operand(1), t1);

		t3 = mulPoly(u, a->operand(1));
		t4 = addPolyGf(t2, t3, x, pow(p, k));

		result->includeOperand(v);
		result->includeOperand(t4);

		delete u;

		delete t1;
		delete t2;
		delete t3;

		delete s;

		printf("result = %s\n", result->toString().c_str());
	}

	return result;
}

AST* factorsWangRec(AST* f, AST* L, AST* K, Int mod)
{
	long long i, j;

	Int nrm1 = std::numeric_limits<long long>::min(), nrm2 = std::numeric_limits<long long>::min();
	
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

	j = 0;

	for(i = 0; i < S->numberOfOperands(); i++)
	{
		AST* pp_u0 = S->operand(i)->operand(1);

		nrm2 = norm(pp_u0, x);
	
		if(nrm2 > nrm1)
		{
			nrm1 = nrm2;
			j = i;
		}
	}
	// todo: set c to config where pp_u0 have largest max norm
	// or tests all configs on wang leading coeff and verify if
	// one works
	AST* c = S->operand(j);

	AST* delta = c->operand(0);
	AST* pp_u0 = c->operand(1);
	AST* sF    = c->operand(2);
	AST* u   	 = c->operand(3);
	AST* a     = c->operand(4);

	AST* k = wangLeadingCoeff(f, delta, u, Vn, sF, a, L, K);

	if(k->kind() == Kind::Fail)
	{
		// try again
		return factorsWangRec(f, L, K, mod + 1);
	}

	return nullptr;
}

AST* factorsWang(AST* f, AST* L, AST* K)
{
	return factorsWangRec(f, L, K, 3);
}

}
