#include "Core/AST/AST.hpp"
#include "Utils.hpp"
#include "Hensel.hpp"
#include "Berlekamp.hpp"
#include "Zassenhaus.hpp"
#include "SquareFree.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Debug/Assert.hpp"
#include "Core/Primes/Primes.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Calculus/Calculus.hpp"
#include "Core/GaloisField/GaloisField.hpp"

#include <cmath>

using namespace ast;
using namespace algebra;
using namespace calculus;
using namespace polynomial;
using namespace galoisField;
using namespace simplification;

namespace factorization {

void subsetsRec(Expr& arr, Expr& data, Expr& s, Int start, Int end, Int index, Int r)
{
	Expr c;

	long i, j;

	if (index == r)
	{
		s.insert(set({}));

		c = s[s.size() - 1];

		for(j = 0; j < r; j++)
		{
			c.insert(data[j]);
		}
	}
	else
	{
		for (i = start.longValue(); i <= end && end.longValue() - i + 1 >= r - index; i++)
		{
			data.insert(arr[i]);
			subsetsRec(arr, data, s, i+1, end, index+1, r);
			data.remove(data.size() - 1);
		}
	}
}

Expr subset(Expr s, Int r)
{
	long n = s.size();

	Expr d, res;

	d = set({});

	res = set({});

	subsetsRec(s, d, res, 0, n - 1, 0, r);

	return res;
}

// from Algorithms for Computer Algebra Geddes
Expr cantorZassenhausDDF(Expr v, Expr x, Int p)
{
	long i;

	Expr h, f, t, G, g, n;

	i = 1;
	h = x;
	f = v;

	G = list({});

	g = 1;

	n = degree(f, x);

	while(n.value() >= 2*i)
	{
		t = powModPolyGf(h, f, x, p, p, true);

		h = t;

		t = subPolyGf(h, x, x, p, true);

		g = gcdPolyGf(t, f, x, p, true);

		if(g != 1)
		{
			G.insert(list({ g, integer(i) }));

			t = quoPolyGf(f, g, x, p, true);

			f = t;

			t = remPolyGf(h, f, x, p, true);

			h = t;
		}

		n = degree(f, x);

		i = i + 1;
	};

	if(f != 1)
	{
		G.insert(list({ f, degree(f, x) }));
	}

	return G;
}

// from Algorithms for Computer Algebra Geddes
Expr cantorZassenhausEDF(Expr a, Expr x, Int n, Int p)
{
	Int m, i;

	Expr g, da, F, v, h, k, f1, f2, t;

	da = degree(a, x);

	if(da.value() <= n)
	{
		return list({ a });
	}

	m = da.value() / n;

	F = list({ a });

	while(F.size() < m)
	{
		v = randPolyGf(2*n - 1, x, p);

		if(p == 2)
		{
			t = v;

			for(i = 0; i < pow(2, n * m - 1); i++)
			{
				h = powModPolyGf(t, a, x, 2, p, true);
				k = addPolyGf(v, h, x, p, true);

				t = h;
				v = k;
			}
		}
		else
		{
			h = powModPolyGf(v, a, x, (pow(p, n) - 1) / 2, p, true);

			v = h;
			k = 1;

			h = subPolyGf(v, k, x, p, true);

			v = h;
		}

		g = gcdPolyGf(a, v, x, p, true);

		if(g != 1 && g != a)
		{
			k = quoPolyGf(a, g, x, p, true);

			f1 = cantorZassenhausEDF(g, x, n, p);
			f2 = cantorZassenhausEDF(k, x, n, p);

			F = append(f1, f2);

		}
	}

	return F;
}

Expr cantorZassenhaus(Expr u, Expr x, Int m)
{
	// [lc, f]= monicGf(f, p)
	// if  deg(f) < 1: return [lc, []]

	Expr F = cantorZassenhausDDF(u, x, m);

	Expr f = list({});

	for(long i = 0; i < F.size(); i++)
	{
		Expr k = F[i][0];
		Int n = F[i][1].value();

		Expr T = cantorZassenhausEDF(k, x, n, m);

		while(T.size())
		{
			f.insert(T[0]);
			T.remove(0);
		}

	}

	return f;
}


Expr squareFreeFactoringGf(Expr u, Expr x, Int m)
{
	Expr T = monicPolyGf(u, x, m, false);

	Expr lc = T[0];
	Expr f = T[1];

	// printf("%s\n", T->toString().c_str());

	Expr n = degree(f, x);

	if(n.value() == 1)
	{
		return list({lc, list({})});
	}

	Expr F = cantorZassenhaus(f, x, m);

	return list({lc, F});
}

// From modern computer algebra by Gathen
Expr zassenhaus(Expr f, Expr x)
{
	bool stop = false;

	Int s, i, j, l, p, A, B, C, gamma, gcd;

	Expr g, n, b, F, D, E, H, Z, G, T, S, M, u, v, gi, L, K, I;

	L = list({ x });

	K = Expr("Z");

	n = degree(f, x);

	if(n == 1)
	{
		return list({ f });
	}

	A = norm(f, x);

	b = leadCoeff(f, x);

	B = Int(std::abs(std::sqrt(n.value().longValue() + 1))) * Int(pow(2, n.value())) * A * b.value();

	C = pow(n.value() + 1, 2 * n.value()) * pow(A, 2 * n.value() - 1);

	gamma = std::ceil(2 * (2 * n.value().longValue() * log2(n.value().longValue()+1) + (2 * n.value().longValue() - 1) * log2(A.longValue())));

	// choose a prime number p such that f be square free in Zp[x]
	// and such that p dont divide lc(f)

	for(i = 1; primes[i.longValue()] <= 2 * gamma.longValue() * std::log(gamma.longValue()); i++)
	{
		p = primes[i.longValue()];

		if(b.value() % p == 0)
		{
			continue;
		}

		F = gf(f, p, true);

		D = derivate(F, x);

		E = gf(D, p, true);

		D = gcdPolyGf(F, E, x, p, false);

		gcd = D.value();

		if(b.value() % p > 0 && gcd == 1)
		{
			break;
		}
	}

	l = std::ceil(std::log(2*B.longValue() + 1) / std::log(p.longValue()));

	I = squareFreeFactoringGf(f, x, p);

	Z = I[1];

	I.remove(1);

	g = multifactorHenselLifting(f, Z, x, p, l);

	T = set({});

	for(i = 0; i < g.size(); i++)
	{
		T.insert(integer(i));
	}

	F = list({});

	s = 1;

	while(2*s <= T.size())
	{
		stop = false;

		M = subset(T, s);

		for(j = 0; j < M.size(); j++)
		{
			S = M[j];

			H = mul({b});
			G = mul({b});

			Z = difference(T, S);

			for(i = 0; i < S.size(); i++)
			{
				gi = g[S[i].value()];

				G.insert(gi);
			}

			for(i = 0; i < Z.size(); i++)
			{
				gi = g[Z[i].value()];

				H.insert(gi);
			}

			u = gf(G, pow(p, l), true);
			v = gf(H, pow(p, l), true);

			G = u;
			H = v;

			if(norm(G, x) > pow(p, l) / 2)
			{
				continue;
			}

			if(norm(H, x) > pow(p, l) / 2)
			{
				continue;
			}

			if(l1norm(G, x) * l1norm(H, x) <= B)
			{
				T = Z;

				F.insert(pp(G, L, K));

				f = pp(H, L, K);

				b = leadCoeff(f, x);

				stop = true;
			}

			if(stop) break;
		}

		if(!stop) s = s + 1;
	}

	F.insert(f);

	return F;
}


}
