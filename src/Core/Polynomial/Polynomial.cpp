#include "Polynomial.hpp"
#include "Resultant.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Calculus/Calculus.hpp"

#include <numeric>

using namespace ast;
using namespace expand;
using namespace simplification;
using namespace algebra;
using namespace calculus;


namespace polynomial {

void includeVariable(std::vector<Expr>& vars, Expr u) {
	bool included = false;

	for(Expr k : vars) {
		if(k==(u)){
			included = true;
			break;
		}
	}

	if(!included) {
		vars.push_back(u);
	}
}

bool isGeneralMonomial(Expr u, Expr v) {
	Expr S;
	if(v.kind() != Kind::Set) {
		S = set({v});
	} else {
		S = v;
	}

	if(exists(S, u)) {
		
		return true;
	} else if(u.kind() == Kind::Power){
		Expr b = u[0];
		Expr e = u[1];

		if(exists(S, b) && e.kind() == Kind::Integer && e.value() > 1) {
			
			return true;
		}
	} else if(u.kind() == Kind::Multiplication) {
		for(unsigned int i=0; i<u.size(); i++) {
			if(isGeneralMonomial(u[i], S) == false) {
				
				return false;
			}
		}
		
		return true;
	}

	return u.freeOfElementsInSet(S);
}

bool isGerenalPolynomial(Expr u, Expr v) {
	Expr S;

	if(v.kind() != Kind::Set) {
		S = set({ v });
	} else {
		S = v;
	}

	if(u.kind() != Kind::Addition && u.kind() != Kind::Subtraction) {
		bool r = isGeneralMonomial(u, S);
		
		return r;
	}

	if(exists(S, u)) {
		
		return true;
	}

	for(unsigned int i=0; i<u.size(); i++) {
		if(isGeneralMonomial(u[i], S) == false) {
			
			return false;
		}
	}

	
	return true;
}

Expr coeffVarMonomial(Expr u, Expr S) {
	if(!isGeneralMonomial(u, S))
		return undefined();

	if(isConstant(u))
		return list({ u, integer(1) });

	if(exists(S, u))
		return list({ integer(1), u });
	
	if(u.kind() == Kind::Power && exists(S, u[0]))
		return list({ integer(1), u });

	if(u.kind() == Kind::Multiplication) {
		Expr C = list({});
		Expr V = list({});
	
		for(unsigned int i=0; i<u.size(); i++) {
			Expr L = coeffVarMonomial(u[i], S);
		
			Expr CL = list({L[0]});
			Expr VL = list({L[1]});
	
			Expr C_ = join(C, CL);
			Expr V_ = join(V, VL);
			
			
			
			
			C = C_;
			V = V_;
		
			
			
			
		}

		Expr coefs = mul({});
		Expr vars  = mul({});
		
		for(unsigned int i=0; i<C.size(); i++) {
			if(
				C[i].kind() == Kind::Integer &&
				C[i].value() == 1
			) continue;
				
				coefs.insert(C[i]);
		}
		
		for(unsigned int i=0; i<V.size(); i++) {
			if(
				V[i].kind() == Kind::Integer &&
				V[i].value() == 1
			) continue;
			vars.insert(V[i]);
		}
		
		
		
	
		if(coefs.size() == 0) {
			
			coefs = integer(1);
		} else if(coefs.size() == 1) {
			Expr coefs_ = coefs[0];
			
			coefs = coefs_;
		}
	
		if(vars.size() == 0) {
			
			vars = integer(1);
		} else if(vars.size() == 1) {
			Expr vars_ = vars[0];
			
			vars = vars_;
		}
	
		return list({ coefs, vars });
	}

	return list({ u, integer(1) });
}

Expr collectTerms(Expr u, Expr S) {
	if(u.kind() != Kind::Addition) {
		Expr L = coeffVarMonomial(u, S);
		if(L.kind() == Kind::Undefined) {
			
			return undefined();
		}

		return u;
	}

	if(exists(S, u)) {
		return u;
	}

	int N = 0;

	Expr T = list({});

	for(unsigned int i=0; i<u.size(); i++) {
		Expr f = coeffVarMonomial(u[i], S);
		
		if(f.kind() == Kind::Undefined) {
			
			return undefined();
		}
	
		int j = 1;
		bool combined = false;
		
		while(!combined && j <= N) {
			int j_ = j - 1;

			if(f[1]==(T[j_][1])) {
				
				Expr Tj = list({
					add({ 
						T[j_][0],
						f[0]
					}),
					f[1]
				}); 

				Expr Tj_ = T[j_];
				T.remove(j_);
				
				
				
				T.insert(Tj, j_);
				
				combined = true;
			}

			j = j+1;
		}

		if(!combined) {
			T.insert(f, N);
			N = N + 1;
		}
	
		
	}

	Expr v = add({});

	for(int j=0; j<N; j++) {
		if(
			T[j][1].kind() == Kind::Integer &&
			T[j][1].value() == 1
		) {
			v.insert(T[j][0]);
		} else {
			v.insert(mul({
				T[j][0],
				T[j][1],
			}));
		}
	}

	

	if(v.size() == 0) {
		
		return integer(0);
	}

	if(v.size() == 1) {
		Expr v_ = v[0];
		
		v = v_;
	}

	return v;
}

Expr degreeGME(Expr u, Expr v) {
	if(u.kind() == Kind::Integer && u.value() == 0)
		return Expr(Kind::MinusInfinity);
	
	if(isConstant(u))
		return integer(0);

	Expr S;
	if(v.kind() != Kind::Set) {
		S = set({ v });
	} else {
		S = v;
	}

	if(exists(S, u)) {

		
		return integer(1);
	} else if(u.kind() == Kind::Power){
		Expr b = u[0];
		Expr e = u[1];

		if(exists(S, b) && isConstant(e)) {
			
			return e;
		}

	} else if(u.kind() == Kind::Multiplication) {
		Expr deg = integer(0);
		for(unsigned int i=0; i<u.size(); i++) {
			Expr deg_ = degreeGME(u[i], S);
			if(deg_.value() > deg.value()) {
				
				deg = deg_;
			} else {
				
			}
		}
		
		return deg;
	}
	
	
	return integer(0);
}

Expr degree(Expr u, Expr v) {
	Expr S;

	if(u.kind() == Kind::Integer && u.value() == 0) {
		return Expr(Kind::MinusInfinity);
	}
	
	if(v.kind() != Kind::Set) {
		S = set({v});
	} else {
		S = v;
	}

	if(u.kind() != Kind::Addition && u.kind() != Kind::Subtraction) {
		Expr r = degreeGME(u, S);
		
		return r;
	}

	if(exists(S, u)) {
		
		return integer(1);
	}

	Expr deg = integer(0);

	for(unsigned int i=0; i<u.size(); i++) {
		Expr deg_ = degreeGME(u[i], S);

		if(deg_.value() > deg.value()) {
			
			deg = deg_;
		} else {
			
		}
	}

	

	return deg;
}

Expr variables(Expr u) {
	if(
		u.kind() == Kind::Integer ||
		u.kind() == Kind::Fraction
	)	return set({});

	if(u.kind() == Kind::Power) {
		Expr b = u[0];
		Expr e = u[1];

		if(e.kind() == Kind::Integer && e.value() > 1)
			return set({ b });

		return set({ u });
	}

	if(
		u.kind() == Kind::Addition ||
		u.kind() == Kind::Subtraction ||
		u.kind() == Kind::Multiplication
	) {
		Expr S = set({});

		for(unsigned int i=0; i<u.size(); i++) {
			Expr S_ = variables(u[i]);
			Expr S__ = unification(S, S_);
			
			
			S = S__;
		}

		return S;
	}

	return set({ u });
}

Expr coefficientGME(Expr u, Expr x) 
{
	if(u==(x)) 
	{
		return list({ integer(1), integer(1) });
	}

	if(u.kind() == Kind::Power) 
	{
		Expr b = u[0];
		Expr e = u[1];

		if(b==(x) && e.kind() == Kind::Integer && e.value() > 0) 
		{
			return list({ integer(1), e });
		}

	} 
	else if(u.kind() == Kind::Multiplication) 
	{
		Expr m = integer(0);
		Expr c = u;

		for(unsigned int i=0; i<u.size(); i++) 
		{
			Expr f =	coefficientGME(u[i], x);
			
			if(f.kind() == Kind::Undefined) 
			{
				
				
				
		
				return undefined();
			} 
			if(
				f[1].kind() != Kind::Integer ||
				f[1].value() != 0
			) 
			{
				
				
	
				m = f[1];
				Expr c_ = div(u, power(x, m));
				c = algebraicExpand(c_);
				
			}
	
			
		}

		return list({ c, m });
	}

	if(u.freeOf(x)) 
	{
		return list({ u, integer(0) });
	}

	return undefined();
}

Expr coeff(Expr u, Expr x, Expr j) 
{
	
	if(u.kind() != Kind::Addition && u.kind() != Kind::Subtraction) {
		Expr f = coefficientGME(u, x);
		
		if(f.kind() == Kind::Undefined) return f;
		
		if(j==(f[1])) {
			Expr k = f[0];
			
			

			return k;
		}

		

		return integer(0);
	}

	if(x==(u)) {
		if(j.kind() == Kind::Integer && j.value() == 1) {
			return integer(1);
		}

		return integer(0);
	}

	Expr c = integer(0);

	for(unsigned int i=0; i<u.size(); i++) 
	{
		Expr f = coefficientGME(u[i], x);

		if(f.kind() == Kind::Undefined) return f;
		
		if(j==(f[1])) {
			Expr k = f[0];

			if(c.kind() == Kind::Integer && c.value() == 0) {
				
				c = Expr(u.kind());
				c.insert(k);
			} else {
				c.insert(k);
			}
		}

		
	}

	if(c.kind()!= Kind::Integer && c.size() == 1) {
		Expr l = c[0];
		
		return l;
	}
	
	return c;
}

Expr leadCoeff(Expr u, Expr x) {
	Expr d = degree(u, x);
	Expr lc = coeff(u, x, d);

	return lc;
}

Expr divideGPE(Expr u, Expr v, Expr x) {
	Expr t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;

	Expr q = integer(0);
	Expr r = u;

	Expr m = degree(r, x);
	Expr n = degree(v, x);

	Expr lcv = leadCoeff(v, x);

	while(
		m.kind() != Kind::MinusInfinity &&
		(m.kind() == Kind::Integer && n.kind() == Kind::Integer &&
		m.value() >= n.value())
	) {
		Expr lcr = leadCoeff(r, x);

		Expr s = div(lcr, lcv);

		t1 = power(x, m - n);
		t2 = mulPoly(s, t1);
		t3 = addPoly(q, t2);

		q = reduceAST(t3);
	
		t1 = power(x, m);
		t2 = mulPoly(lcr, t1);
		t3 = subPoly(r, t2);

		t4 = power(x, n);
		t5 = mulPoly(lcv, t4);
		t6 = subPoly(v, t5);

		t7 = mulPoly(t6, s);
		t8 = power(x, sub({m, n}));
		t9 = mulPoly(t7, t8);
		t10 = subPoly(t3, t9);
	
		r = reduceAST(t10);
		
		m = degree(r, x);

	}

	Expr res = list({ algebraicExpand(q), algebraicExpand(r) });
	
	return res;
}

Expr quotientGPE(Expr u, Expr v, Expr x) {
	return divideGPE(u,v,x)[0];
}

Expr remainderGPE(Expr u, Expr v, Expr x) {
	return divideGPE(u,v,x)[1];
}

Expr expandGPE(Expr u, Expr v, Expr x, Expr t) {
	if(u == 0)
		return 0;

	Expr d = divideGPE(u, v, x);

	Expr q = d[0];
	Expr r = d[1];

	Expr expoent = add({
		mul({
			t,
			expandGPE(q, v, x, t)
		}),
		r
	});

	return algebraicExpand(expoent);
}


Expr gcdGPE(Expr u, Expr v, Expr x) {
	if(
		u.kind() == Kind::Integer && u.value() == 0 &&
		v.kind() == Kind::Integer && v.value() == 0
	) {
		return integer(0);
	}

	Expr U = u;
	Expr V = v;

	while (V.kind() != Kind::Integer || (V.kind() == Kind::Integer && V.value() != 0)) {
		Expr R = remainderGPE(U, V, x);

		U = V;
		V = R;
	}

	Expr e = mul({div(integer(1), leadCoeff(U,x)), U});
	Expr res = algebraicExpand(e);

	return res;
}

Expr extendedEuclideanAlgGPE(Expr u, Expr v, Expr x) {
	if(
		u.kind() == Kind::Integer && u.value() == 0 &&
		v.kind() == Kind::Integer && v.value() == 0
	) {
		return list({ integer(0), integer(0), integer(0) });
	}

	Expr U 		= u;
	Expr V 		= v;

	Expr App 	= 1, Ap = 0, Bpp = 0, Bp = 1;

	while (V != 0) {
		Expr d = divideGPE(U,V,x);

		Expr q = d[0];
		Expr r = d[1];

		Expr A_ = sub({
			App,
			mul({
				q,
				Ap
			})
		});

		Expr B_ = sub({
			Bpp,
			mul({
				q,
				Bp
			})
		});
	
		Expr A = algebraicExpand(A_);
		Expr B = algebraicExpand(B_);
		
		App = Ap;
		
		Ap 	= A;

		Bpp = Bp;

		Bp 	= B;
		
		U = V;
		
		V = r;
		
	}

	Expr c = leadCoeff(U, x);

	Expr App_ = quotientGPE(App, c, x);
	
	App = App_;

	Expr Bpp_ = quotientGPE(Bpp, c, x);
	
	Bpp = Bpp_;

	Expr U_ = quotientGPE(U, c, x);
	
	U = U_;

	return list({ U, App, Bpp });
}

Expr mulPolyRec(Expr p1, Expr p2)
{
	if(p1.kind() == Kind::Addition)
	{
		Expr res = add({});

		for(unsigned int i = 0; i < p1.size(); i++)
		{
			res.insert(mulPoly(p1[i], p2));
		}

		return res;
	}

	if(p2.kind() == Kind::Addition)
	{
		Expr res = add({});
		
		for(unsigned int i = 0; i < p2.size(); i++)
		{
			res.insert(mulPoly(p2[i], p1));
		}

		return res;
	}

	return mul({ p1, p2 });
}

Expr mulPoly(Expr p1, Expr p2)
{
	Expr t1 = mulPolyRec(p1, p2);
	Expr t2 = reduceAST(t1);
	
	

	return t2;
}

Expr subPoly(Expr p1, Expr p2)
{
	Expr s = sub({ p1, p2 });
	Expr p = reduceAST(s);
	
	return p;
}

Expr addPoly(Expr p1, Expr p2)
{
	Expr s = add({ p1, p2 });
	Expr p = reduceAST(s);
	
	return p;
}

Expr recPolyDiv(Expr u, Expr v, Expr L, Expr K) {
	assert(
		K.identifier() == "Z" || K.identifier() == "Q",
		"Field needs to be Z or Q"
	);

	if(L.size() == 0) 
	{
		Expr k = div(u, v);
	
		Expr d = algebraicExpand(k);

		if(K.identifier() == "Z") 
		{
			if(d.kind() == Kind::Integer) 
			{
				return list({ d, integer(0) });
			}
			
			return list({ integer(0), u });
		}
	
		return list({ d, integer(0) });
	}

	Expr x = first(L);
	Expr r = u;

	Expr m = degree(r, x);
	Expr n = degree(v, x);
	
	Expr q = integer(0);
	Expr lcv = leadCoeff(v, x);

	while(m.kind() != Kind::MinusInfinity && m.value() >= n.value()) 
	{
		Expr lcr = leadCoeff(r, x);
	
		Expr R = rest(L);

		Expr d = recPolyDiv(lcr, lcv, R, K);

		if(d[1] != 0) 
		{
			return list({ algebraicExpand(q), r });
		}

		Expr j = power(x, sub({ m, n }));

		q = add({q, mul({ d[0], j})});

		Expr t1 = mulPoly(v, d[0]);
		Expr t2 = mulPoly(t1, j);
		Expr t3 = subPoly(r, t2);

		r = reduceAST(t3);
	
		m = degree(r, x);
	}

	return list({ algebraicExpand(q), r });
}

Expr recQuotient(Expr u, Expr v, Expr L, Expr K) {
	Expr r = recPolyDiv(u, v, L, K);
	Expr q = r[0];
	
	return q;
}

Expr recRemainder(Expr u, Expr v, Expr L, Expr K) {
	Expr r = recPolyDiv(u, v, L, K);
	Expr q = r[1];
	
	return q;
}

Expr pdiv(Expr f, Expr g, Expr x)
{
	assert(g != 0, "Division by zero!");

	Expr lg, k, q, r, t, m, n, j;
	Expr t1, t2, t3, t4, t5, t6;

	m = degree(f, x);
	n = degree(g, x);

	if(m.value() < n.value())
	{
		return list({ 0, f });
	}

	if(g == 1)
	{
		return list({ f, 0});
	}

	q = 0;
	r = f;
	t = m;

	k = m - n + 1;

	lg = leadCoeff(g, x);

	// printf("df\n");
	// printf("%s\n", m->toString().c_str());
 
	// printf("dg\n");
	// printf("%s\n", n->toString().c_str());

	// printf("lc_g\n");
	// printf("%s\n", lg->toString().c_str());

	while(true)
	{
		t1 = leadCoeff(r, x);
		j = sub({ t, n });
		k = sub({ k, integer(1) });
		t3 = power(x, j);

		t2 = mulPoly(q, lg); //mul({ q, lg });
		t4 = mulPoly(t1, t3); 
		q  = addPoly(t2, t4);

		
		
		// q = add({ t2, mul({ t1, t3}) });
		
		t4 = mulPoly(r, lg);// mul({ r, lg });
		t5 = mulPoly(g, t1); //mul({ g, t1, t3 });
		t6 = mulPoly(t5, t3);
		r  = subPoly(t4, t6);

		// t5 = sub({ t4, t5 });
		// r = algebraicExpand(t6);

		
		
		
		

		t = degree(r, x);

		if(t.kind() == Kind::MinusInfinity || t.value() < n.value())
		{
			break;
		}

		
	}

	q = mul({q, power(lg, k)});
	r = mul({r, power(lg, k)});

	
	
	

	t1 = algebraicExpand(q);
	t2 = algebraicExpand(r);
	
	
	
	
	return list({t1, t2});
}

Expr pseudoDivision(Expr u, Expr v, Expr x) 
{
	Expr p = integer(0);
	Expr s = u;

	Expr m = degree(s, x);
	Expr n = degree(v, x);

	Expr delta = integer(max(m.value() - n.value() + 1, Int(0)));

	Expr lcv = leadCoeff(v, x);

	long tal = 0;

	while(m.kind() != Kind::MinusInfinity && m.value() >= n.value()) 
	{
		Expr lcs = leadCoeff(s, x);

		Expr j = power(x, sub({ m, n }));
		
		Expr t1 = mulPoly(lcv, p);
		Expr t2 = mulPoly(lcs, j);
		Expr t3 = addPoly(t1, t2);
		
		
		

		
		p = reduceAST(t3);
		
		

		Expr t4 = mulPoly(lcv, s);
		Expr t5 = mulPoly(lcs, v);
		Expr t6 = mulPoly(t5, j);

		
	
		s = subPoly(t4, t6);
	
		
		
		

		tal = tal + 1;

		
	
		m = degree(s, x);
	
		
	
		
	}

	Expr k = power(lcv, integer(delta.value() - tal));

	Expr A = mulPoly(k, p);
	Expr B = mulPoly(k, s);

	Expr Q = reduceAST(A);
	Expr R = reduceAST(B);

	
	
	
	
	
	
	
	

	
	
	
	return list({ Q, R });
}

Expr pseudoQuotient(Expr u, Expr v, Expr x) {
	Expr r = pseudoDivision(u, v, x);
	Expr q = r[0];
	
	return q;
}

Expr pseudoRemainder(Expr u, Expr v, Expr x) {
	Expr r = pseudoDivision(u, v, x);
	Expr q = r[1];
	
	return q;
}


Expr getNormalizationFactor(Expr u, Expr L, Expr K) {
	assert(
		K.identifier() == "Z" || K.identifier() == "Q", 
		"field must be Z or Q"
	);

	if(u==(0))
	{
		return integer(0);
	}

	if(isConstant(u)) 
	{
		if(u > 0) 
		{
			if(K.identifier() == "Z") 
			{
				return integer(1);
			}
	
			return power(u, integer(-1));
		}
		else 
		{
			if(K.identifier() == "Z") 
			{
				return -1;
			}

			return -1 * power(u, -1);
		}
	}

	if(L.size() == 0)
	{
		return undefined();
	}
	
	Expr lc = leadCoeff(u, L[0]);

	Expr rL = rest(L);

	Expr cf = getNormalizationFactor(lc, rL, K);

	
	

	return cf;
}

Expr normalizePoly(Expr u, Expr L, Expr K) 
{
	if(u.kind() == Kind::Integer && u.value() == 0)
	{
		return integer(0);
	}

	Expr u__ = mul({ getNormalizationFactor(u, L, K), u });

	Expr u_ = algebraicExpand(u__);

	
	
	return u_;
}


Expr unitNormal(Expr v, Expr K)
{
	assert(
		K.identifier() == "Z" || K.identifier() == "Q", 
		"field must be Z or Q"
	);

	if(K.identifier() == "Z")
	{
		if(v < 0)
		{
			return -1;
		}
	
		return 1;
	}

	if(K.identifier() == "Q")
	{
		if(v < 0)
		{
			return mul({ integer(-1), power(v, integer(-1)) });
		}

		return power(v, integer(-1));
	}

	return integer(1);
}

// Finds the content of u with respect to x using
// the auxiliary variables R with coeff domain K,
// with is Z or Q
Expr polynomialContent(Expr u, Expr x, Expr R, Expr K) 
{	
	if(u==(0))
	{
		return integer(0);
	}

	Expr n = degree(u, x);

	Expr g = coeff(u, x, n);

	Expr k = sub({ u, mul({g, power(x, n)}) });

	Expr v = algebraicExpand(k);

	

	if(v==(0))
	{
		Expr un = unitNormal(g, K);
		Expr t = mul({un, g});
	
		
	
		g = reduceAST(t);
		
		
	}
	else
	{
		while(v != 0)
		{
			Expr d = degree(v, x);
			Expr c = leadCoeff(v, x);
		
			Expr t = mvPolyGCD(g, c, R, K);

			

			g = t;

			k = sub({v, mul({c, power(x, d)})});
			v = algebraicExpand(k);

			
			

			
		}
	}

	
	

	return g;
}


// Finds the content of u with respect to x using
// the auxiliary variables R with coeff domain K,
// with is Z or Q
Expr polynomialContentSubResultant(Expr u, Expr x, Expr R, Expr K) 
{
	if(u==(0))
	{
		return integer(0);
	}

	Expr n = degree(u, x);

	Expr g = coeff(u, x, n);

	Expr k = sub({u, mul({g, power(x, n)})});

	Expr v = algebraicExpand(k);

	

	if(v==(0))
	{
		Expr un = unitNormal(g, K);
		Expr t = mul({un, g});
	
		
	
		g = reduceAST(t);
		
		
	}
	else
	{
		while(v != 0)
		{
			Expr d = degree(v, x);
			Expr c = leadCoeff(v, x);
		
			Expr t = mvSubResultantGCD(g, c, R, K);

			

			g = t;

			k = sub({v, mul({c, power(x, d)})});
			v = algebraicExpand(k);

			
			

			
		}
	}


	
	
	
	return g;
}

Expr subResultantGCDRec(Expr u, Expr v, Expr L, Expr K)
{
	if(L.size() == 0) 
	{
		if(K.identifier() == "Z") 
		{ 
			return integerGCD(u, v); 
		}
	
		if(K.identifier() == "Q") 
		{ 
			return integer(1); 
		}
	}

	Expr x = first(L);

	Expr du = degree(u, x);
	Expr dv = degree(v, x);
	
	Expr U = undefined();
	Expr V = undefined();

	if(du.value() >= dv.value())
	{
		U = u;
		V = v;
	}
	else
	{
		U = v;
		V = u;
	}

	
	

	Expr R = rest(L);
	
	Expr contU = polynomialContentSubResultant(U, x, R, K);
	Expr contV = polynomialContentSubResultant(V, x, R, K);

	Expr d = subResultantGCDRec(contU, contV, R, K);

	Expr tmp1 = recQuotient(U, contU, L, K);
	Expr tmp2 = recQuotient(V, contV, L, K);

	
	U = tmp1;
	
	
	V = tmp2;

	Expr tmp3 = leadCoeff(U, x);
	Expr tmp4 = leadCoeff(V, x);

	Expr g = subResultantGCDRec(tmp3, tmp4, R, K);

	
	

	int i = 1;

	Expr delta = undefined();
	Expr y = undefined();
	Expr b = undefined();
	Expr dp = undefined();

	while (V != 0)
	{
		Expr r = pseudoRemainder(U, V, x);
	
		if(r != 0)
		{
			if(i == 1)
			{
	
				Expr tmp3 = add({
					degree(U, x),
					mul({integer(-1), degree(V, x) }),
					integer(1)
				});

				delta = algebraicExpand(tmp3);

				

				y = integer(-1);
				
				Expr tmp4 = power(integer(-1), delta);
				
				b = algebraicExpand(tmp4);
				
				
			}
			else
			{
				dp = delta;
			
				Expr tmp3 = add({
					degree(U, x),
					mul({integer(-1), degree(V, x)}),
					integer(1)
				});

				
				delta = algebraicExpand(tmp3);

				

				Expr f = leadCoeff(U, x);

				Expr tmp4 = power(mul({integer(-1), f}), sub({dp, integer(1)}));
				Expr tmp5 = power(y, sub({dp, integer(2)}));
				
				Expr tmp6 = algebraicExpand(tmp4);
				Expr tmp7 = algebraicExpand(tmp5);
				
				
				
				
				y = recQuotient(tmp6, tmp7, R, K);

				
				

				Expr tmp8 = mul({
					integer(-1),
					f,
					power(y, sub({ delta, integer(1) }))
				});
				
				b = algebraicExpand(tmp8);
				
				
			}
			
			
			U = V;

			

			V = recQuotient(r, b, L, K);

			i = i + 1;
		}
		else
		{
			
			U = V;

			
			V = r;
		}

		
	}

	
	
	

	Expr tmp5 = leadCoeff(U, x);

	Expr s = recQuotient(tmp5, g, R, K);

	

	Expr W = recQuotient(U, s, L, K);

	

	Expr contW = polynomialContentSubResultant(W, x, R, K);
	Expr ppW = recQuotient(W, contW, L, K);
		
	Expr tmp6 = mul({d, ppW});
	Expr res = algebraicExpand(tmp6);
	
	
	
	

	
	
	

	
	

	
	

	
	

	return res;
}

Expr mvSubResultantGCD(Expr u, Expr v, Expr L, Expr K) {
	if(u==(0)) {
		return normalizePoly(v, L, K);
	}

	if(v==(0)) {
		return normalizePoly(u, L, K);
	}

	Expr gcd = subResultantGCDRec(u, v, L, K);

	Expr r = normalizePoly(gcd, L, K);

	
	
	return r;
}

Expr mvPolyGCDRec(Expr u, Expr v, Expr L, Expr K) 
{
	if(L.size() == 0) 
	{
		if(K.identifier() == "Z") 
		{ 
			return integerGCD(u, v); 
		}
	
		if(K.identifier() == "Q") 
		{ 
			return integer(1); 
		}
	}

	Expr x = first(L);
	Expr R = rest(L);
	Expr cont_u = polynomialContent(u, x, R, K);

	Expr cont_v = polynomialContent(v, x, R, K);

	Expr d = mvPolyGCDRec(cont_u, cont_v, R, K);

	Expr pp_u = recQuotient(u, cont_u, L, K);
	Expr pp_v = recQuotient(v, cont_v, L, K);

	while(pp_v != 0) 
	{
		Expr r = pseudoRemainder(pp_u, pp_v, x);
	
		Expr pp_r = undefined();

		if(r==(0)) {
			pp_r = integer(0);
		} 
		else 
		{
			Expr cont_r = polynomialContent(r, x, R, K);
			pp_r = recQuotient(r, cont_r, L, K);
			
			
		}

		
	
		
	
		pp_u = pp_v;
		pp_v = pp_r;
	}
	
	Expr k = mul({ d, pp_u });
	Expr result = algebraicExpand(k);

	

	
	
	
	
	
	
	

	return result;
}


Expr mvPolyGCD(Expr u, Expr v, Expr L, Expr K) {
	if(u==(0)) {
		return normalizePoly(v, L, K);
	}

	if(v==(0)) {
		return normalizePoly(u, L, K);
	}

	Expr gcd = mvPolyGCDRec(u, v, L, K);
	Expr r = normalizePoly(gcd, L, K);

	
	
	return r;
}

Expr leadMonomial(Expr u, Expr L) {
	if(L.size() == 0) {
		return u;
	}

	Expr x = first(L);
	Expr m = degree(u, x);

	Expr c = coeff(u, x, m);

	Expr restL = rest(L);

	Expr r_ = mul({
		power(x, m),
		leadMonomial(c, restL)
	});

	
	
	
	

	Expr r = algebraicExpand(r_);

	

	return r;
}

bool wasSimplified(Expr u) {
	if(u.kind() == Kind::Symbol)
		return true;

	if(u.kind() == Kind::Integer)
		return true;

	if(u.kind() == Kind::Division)
		return false;
	
	if(u.kind() == Kind::Multiplication) {
		
		for(unsigned int i=0; i<u.size(); i++) {
			if(u[i].kind() == Kind::Fraction)
			{
				return false;
			}

			if(
				u[i].kind() == Kind::Power &&
				u[i][1].kind() == Kind::Integer &&
				u[i][1].value() < 0
			)	
			{
				return false;
			}
		}

		return true;
	}

	return false;
}


/**
 * Return summation(u[i]/v) if v divides u[i]
 */
Expr G(Expr u, Expr v) {
	if(u.kind() == Kind::Addition || u.kind() == Kind::Subtraction) {
		Expr k = Expr(u.kind());

		for(unsigned int i=0; i<u.size(); i++) 
		{
			Expr z_ = div(u[i], v);
			Expr z  = algebraicExpand(z_);

			

			if(wasSimplified(z))
			{
				k.insert(z);
			}
			else
			{
				
			}
		}

		if(k.size() == 0) {
			
			return integer(0);
		}
	
		return k;
	}

	Expr z_ = div(u, v);

	Expr z = algebraicExpand(z_);

	

	if(wasSimplified(z)) 
	{
		return z;
	}

	

	return integer(0);

	// printf("%i\n", wasSimplified(z));
	// printf("%i\n", k.size());

	// if(k.size() == 0) {
	// 	
	// 	return integer(0);
	// }

	// Expr r = algebraicExpand(k);

	// 
	
	// return r;
}

Expr monomialPolyDiv(Expr u, Expr v, Expr L) {
	Expr q = integer(0);
	Expr r = u;
	Expr vt = leadMonomial(v, L);

	Expr f = G(r, vt);

	// printf("-> %s\n", r->toString().c_str());
	// printf("-> %s\n", vt->toString().c_str());
	// printf("-> %s\n", f->toString().c_str());
	// printf("-> %i\n", f.kind());

	while(f.kind() != Kind::Integer || f.value() != 0) {
		// printf("-> %s\n", f->toString().c_str());
		// printf("-> %i\n", f.kind());

		q = add({ q, f });

		Expr r_ = sub({ r, mul({ f, v }) });
		// printf("-> %s\n", r_->toString().c_str());
	
		r = algebraicExpand(r_);
		// printf("-> %s\n", r_->toString().c_str());

		

		f = G(r, vt);
		// printf("\n");
	}
	
	Expr l = list({ reduceAST(q), reduceAST(r) });
	
	
	
	
	

	return l;
}


// TODO

// monomialBasedPolyExpansion(a^2*b + 2*a*b^2 + b^3 + 2*a + 2*b + 3, a+b, [a, b], t) -> b*t^2 + 2*t + 3
Expr monomialBasedPolyExpansion(Expr u, Expr v, Expr L, Expr t) {
	if(u.kind() == Kind::Integer && u.value() == 0) {
		return integer(0);
	}

	Expr d = monomialPolyDiv(u, v, L);
	Expr q = d[0];
	Expr r = d[1];

	Expr k = add({
		mul({
			t,
			monomialBasedPolyExpansion(q, v, L, t),
		}),
		r
	});

	Expr x = algebraicExpand(k);

	
	
	

	return x;
}

// monomialPolyRem can be used for simplification, for example
// monomialPolyRem(a*i^3 + b*i^2 + c*i + d, i^2 + 1, [i]) -> -a*i - b + c*i + d
// simplification when i^2 + 1 = 0
// also
//monomialPolyRem(sin^4(x)+sin^3(x)+2*sin^2(x)cos^2(x)+cos^4(x), sin^2(x)+cos^2(x)-1, [cos(x), sin(x)]) -> 1 + sin^3(x)
Expr monomialPolyRem(Expr u, Expr v, Expr L) {
	Expr d = monomialPolyDiv(u,v,L);
	Expr r = d[1];
	
	return r;
}

Expr monomialPolyQuo(Expr u, Expr v, Expr L)
{
	Expr d = monomialPolyDiv(u,v,L);
	Expr r = d[0];
	
	return r;
}

Expr algebraicExpandRec(Expr u);

Expr expandProduct(Expr r, Expr s) 
{
	if(r==(0)) return integer(0);
	if(s==(0)) return integer(0);

	if(r.kind() == Kind::Addition && r.size() == 0) return integer(0);
	if(s.kind() == Kind::Addition && s.size() == 0) return integer(0);
	
	if(r.kind() == Kind::Addition) 
	{
		Expr f = r[0];

		Expr k = r;
	
		k.remove(0);

		Expr z = add({ expandProduct(f, s), expandProduct(k, s) });

		

		Expr y = reduceAST(z);

		

		return y;
	}
	else if(s.kind() == Kind::Addition) 
	{
		return expandProduct(s, r);
	}

	Expr a = algebraicExpandRec(r);
	Expr b = algebraicExpandRec(s);

	if(a.kind() == Kind::Addition || b.kind() == Kind::Addition) 
	{
		Expr t = expandProduct(a, b);

		
		
		
		return t;
	}

	return mul({ a, b });
}


Expr expandPower(Expr u, Expr n)
{
	if(u.kind() == Kind::Addition) 
	{
		Expr f = u[0];
	
		Expr o = sub({ u, f });

		Expr s = integer(0);

		for(long k = 0; k <= n.value(); k++) 
		{
			Int a = n.value();
		
			Int d = fact(a) / (fact(k) * fact(a - k));
	
			Expr z = mul({
				integer(d),
				power(f, integer(a - k))
			});

			Expr b = integer(k);
			Expr t = expandPower(o, b);
		
			s = add({ s, expandProduct(z, t) });

			
			
			
		}
	
		

		return s;
	}

	Expr v = power(u, n);

	return v;
}


Expr expandProductRoot(Expr r, Expr s) 
{
	if(r.kind() == Kind::Addition) 
	{
		Expr f = r[0];

		Expr k = r;
	
		k.remove(0);
	
		return f * s + k * s;
	}

	if(s.kind() == Kind::Addition) 
	{
		return expandProductRoot(s, r);
	}

	return r * s;
}

Expr expandPowerRoot(Expr u, Expr n) 
{
	if(u.kind() == Kind::Addition) {
		Expr f = u[0];

		Expr r_ = sub({ u, f });

		Expr r = reduceAST(r_);
		
		Expr s = integer(0);
	
		for(int k_ = 0; k_ <= n.value(); k_++) 
		{
			Expr k = integer(k_);

			Expr c_ = div(
				integer(fact(n.value())),
				integer(fact(k.value()) * fact(n.value() - k.value()))
			);
	
			Expr c = reduceAST(c_);
	
	
			Expr z_ = mul({
				c,
				power(
					f,
					integer(n.value() - k.value())
				)
			});

			Expr z = reduceAST(z_);


			Expr t = expandPowerRoot(r, k);
		
			s = s + expandProductRoot(z, t);
		}
		
		return s;
	}
	
	Expr v = power(u, n);

	return reduceAST(v);
}

Expr algebraicExpandRoot(Expr u) 
{
	if(u.isTerminal())
		return reduceAST(u);

	Expr u_ = reduceAST(u);

	if(u_.kind() == Kind::Addition)
	{
		Expr v = u_[0];
	
		Expr a = sub({
			u_,
			v
		});

		Expr k = reduceAST(a);
	
		Expr t = add({
			algebraicExpandRoot(v),
			algebraicExpandRoot(k)
		});

		
		u_ = reduceAST(t);
		
		
		
		
	}

	if(u_.kind() == Kind::Multiplication)
	{
		Expr v = u_[0];
		Expr e = div(
			u_,
			v
		);
		
		Expr t = reduceAST(e);

		Expr z = expandProductRoot(t, v);

		
		u_ = reduceAST(z);

		
		
		

	}

	if(u_.kind() == Kind::Power) {

		Expr b = u_[0];
		Expr e = u_[1];

		if(e.kind() == Kind::Integer && e.value() >= 2) {
			Expr t = expandPowerRoot(b, e);

			
			u_ = reduceAST(t);
			
		}
	
		if(e.kind() == Kind::Integer && e.value() <= -2) {
			Expr p_ = power(u_, integer(-1));
			Expr p = reduceAST(p_);

			
			
			Expr b_ = p[0];
			Expr e_ = p[1];
	
			Expr t = expandPowerRoot(b_, e_);
			
			
			
			
			u_ = power(t, integer(-1));
		}

		
		
	}
	
	Expr k = reduceAST(u_);

	

	return k;
}


Expr algebraicExpandRec(Expr u) 
{
	if(u.isTerminal()) return u;
	
	// Those kinds needs to be reduced to other
	// kinds before expansion, 
	// subtraction -> addition
	// division    -> multiplication by inverses
	// factorial   -> possible to an integer
	Expr z = undefined();

	if(
		u.kind() == Kind::Subtraction ||
		u.kind() == Kind::Division    ||
		u.kind() == Kind::Factorial 
	)
	{
		z = reduceAST(u);
	}
	else
	{
		z = u;
	}


	if(z.kind() == Kind::Power) 
	{
		// printf("pow -> : %s\n", z->toString().c_str());
		Expr b = z[0];
		Expr e = z[1];

		if(e.kind() == Kind::Integer)
		{
			Expr t = expandPower(b, e);

			

			z = t;
		}
		// printf("pow <- : %s\n", z->toString().c_str());
	}

	if(z.kind() == Kind::Multiplication) 
	{
		// printf("mul -> : %s\n", z->toString().c_str());

		Expr v = z[0];

		z.remove(0L);

		if(z.size() == 0)
		{
			
			z = algebraicExpandRec(v);

			
		}
		else
		{
			Expr q = expandProduct(z, v);

			
		
			z = q;

			
		}
		// printf("mul <- : %s\n", z->toString().c_str());

	}

	if(z.kind() == Kind::Addition) 
	{
		// printf("add -> : %s\n", z->toString().c_str());

		Expr v = z[0];
	
		z.remove(0L);

		if(z.size() == 0)
		{
			

			z = algebraicExpandRec(v);

			
		}
		else
		{
			Expr t = add({ algebraicExpandRec(v), algebraicExpandRec(z) });

			

			

			z = t;
		}
	
		// printf("add <- : %s\n", z->toString().c_str());
	} 

	Expr k = reduceAST(z);

	

	return k;
}

Expr algebraicExpand(Expr u)
{
	Expr t = algebraicExpandRec(u);
	
	Expr k = reduceAST(t);
	
	
	return k;
} 


Expr cont(Expr u, Expr x)
{
	Expr n, c, c1, c2, tmp;

	u = algebraicExpand(u);

	if(u==(0))
	{
		

		return integer(0);
	}

	if(u.size() >= 2)
	{
		c1 = coeff(u, x, degree(u, x));
		u = algebraicExpand(u - c1 * power(x, n));
		
		c2 = coeff(u, x, degree(u, x));
		u = algebraicExpand(u - c2 * power(x, n));

		c = integerGCD(c1, c2);

		while(u != 0)
		{
			c1 = coeff(u, x, degree(u, x));

			tmp = u - c1 * power(x, n);
		
			u = algebraicExpand(tmp);

			c2 = integerGCD(c, c1);

			c = c2;
		}

		return c;
	}

	if(u.size() == 1)
	{
		return coeff(u, x, degree(u, x));
	}

	return 0;
}

Expr cont(Expr u, Expr L, Expr K)
{
	Expr R = rest(L);
	Expr C = polynomialContentSubResultant(u, L[0], R, K);

	

	return C;
}

Expr pp(Expr u, Expr L, Expr K)
{
	Expr c = cont(u, L, K);

	Expr p = recQuotient(u, c, L, K);
	
	

	return p;
}

Expr pp(Expr u, Expr c, Expr L, Expr K)
{
	Expr R = rest(L);

	Expr p = recQuotient(u, c, L, K);
	
	

	return p;
}


}
