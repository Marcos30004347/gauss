#include "Subtraction.hpp"
#include "Addition.hpp"
#include "MathSystem/AST/AST.hpp"
#include "Multiplication.hpp"

#include "MathSystem/Expand/Expand.hpp"
#include <cstddef>
#include <cstdio>
#include <ratio>
#include <tuple>

using namespace ast;
using namespace algebra;

namespace simplification {

// negate[sub[add[a, b, c], sub[d, e]]] 						-> add[a + b + c + -d + -e]
// negate[sub[add[a, b, c], sub[d, e], sub[f, g]]]  -> add[a + b + c + -d + -e + -f + -g]
// negate[sub[sub[a, b, c], sub[d, e], sub[f, g]]]	-> add[a + -b + -c -d + e + -f + g]
Expr subRec(Expr u)
{
	Expr v = undefined();

	if(u[0].kind() == Kind::Subtraction)
	{
		v = subRec(u[0]);
	}
	else if(u[0].kind() == Kind::Addition)
	{
		v = reduceAdditionAST(u[0]);
	}
	else
	{
		v = u[0];
	}

	Expr r = add({});

	for(unsigned int j = 1; j < u.size(); j++)
	{
		if(u[j].kind() == Kind::Subtraction)
		{
			r.insert(subRec(u[j]));
		}
		else if(u[j].kind() == Kind::Addition)
		{
			r.insert(reduceAdditionAST(u[j]));
		}
		else
		{
			r.insert(u[j]);
		}
	}

	if(r.size() > 0 && v.kind() != Kind::Addition)
	{
		v = add({ v });
	}

	for(unsigned int j = 0; j < r.size(); j++)
	{
		Expr rj = r[j];

		if(rj.kind() == Kind::Addition)
		{
			for(unsigned int i = 0; i < rj.size(); i++)
			{
				v.insert(reduceMultiplicationAST(-1 * rj[i]));
			}
		}
		else
		{
			v.insert(reduceMultiplicationAST(-1 * rj));
		}
	}

	return v;
}

Expr reduceSubtractionAST(Expr u)
{
	if(u.kind() != Kind::Subtraction)
	{
		return u;
	}

	return reduceAdditionAST(subRec(u));
}

void flatSubtraction(Expr &u, std::vector<Expr> &L, bool invert = false) {

	// (a - b) - (c - d) - (e + f) =	a, -b, -c, +d, -e, -f

	if(u.kind() == Kind::Subtraction) {
		if(u.size() == 0) return;

		flatSubtraction(u[0], L, invert);

		invert = !invert;

		for (size_t i = 1; i < u.size(); i++) {
			flatSubtraction(u[i], L, invert);
		}
		return;
	}

	if(u.kind() == Kind::Addition) {
		for (size_t i = 0; i < u.size(); i++) {
			flatSubtraction(u[i], L, invert);
		}
		return;
	}

	L.push_back(invert ? reduceMultiplicationExpr(-u) : Expr(u));
}

Expr reduceSubtractionExpr(Expr &&u) {
	if(u.kind() != Kind::Subtraction) return Expr(u);
	if(u.size() == 0) return 0;
	if(u.size() == 1) return u[0];

	std::vector<Expr> L;

	flatSubtraction(u, L);

	Expr f = L[0];
	L[0] = 0;

	sort(L);

	Expr S = Expr(Kind::Addition, { f });

	for(size_t i = 0; i < L.size(); i++) {
		//S.insert(reduceMultiplicationExpr(-L[i]));
		S.insert(L[i]);
	}

	return reduceAdditionExpr(S);
}

Expr reduceSubtractionExpr(Expr &u) {
	if(u.kind() != Kind::Subtraction) return Expr(u);
	return reduceSubtractionExpr(std::forward<Expr>(u));
}



}
