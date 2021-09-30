#include "Subtraction.hpp"
#include "Addition.hpp"
#include "Multiplication.hpp"

#include "Core/Expand/Expand.hpp"

using namespace ast;
using namespace algebra;

namespace simplification {

// negate[sub[add[a, b, c], sub[d, e]]] 						-> add[a + b + c + -d + -e]
// negate[sub[add[a, b, c], sub[d, e], sub[f, g]]]  -> add[a + b + c + -d + -e + -f + -g]
// negate[sub[sub[a, b, c], sub[d, e], sub[f, g]]]	-> add[a + -b + -c -d + e + -f + g]
AST* subRec(AST* u)
{
	AST* v = nullptr;

	if(u->operand(0)->kind() == Kind::Subtraction)
	{
		v = subRec(u->operand(0));
	}
	else if(u->operand(0)->kind() == Kind::Addition)
	{
		v = reduceAdditionAST(u->operand(0));
	}
	else
	{
		v = u->operand(0)->copy();
	}

	AST* r = add({});

	for(unsigned int j = 1; j < u->numberOfOperands(); j++)
	{
		if(u->operand(j)->kind() == Kind::Subtraction)
		{
			r->includeOperand(subRec(u->operand(j)));
		}
		else if(u->operand(j)->kind() == Kind::Addition)
		{
			r->includeOperand(reduceAdditionAST(u->operand(j)));
		}
		else
		{
			r->includeOperand(u->operand(j)->copy());
		}
	}

	if(r->numberOfOperands() > 0 && v->kind() != Kind::Addition)
	{
		v = add({ v });
	}

	for(unsigned int j = 0; j < r->numberOfOperands(); j++)
	{
		AST* rj = r->operand(j);

		if(rj->kind() == Kind::Addition)
		{
			for(unsigned int i = 0; i < rj->numberOfOperands(); i++)
			{
				AST* k = mul({integer(-1), rj->operand(i)->copy()});
				v->includeOperand(reduceMultiplicationAST(k));
				delete k;
			}
		}
		else
		{

			AST* k = mul({integer(-1), rj->copy()});
			v->includeOperand(reduceMultiplicationAST(k));
			delete k;
		}
	}

	delete r;

	return v;
}

AST* reduceSubtractionAST(AST* u) 
{
	if(u->kind() != Kind::Subtraction)
	{
		return u->copy();
	}

	AST* s = subRec(u);

	// AST* s = new AST(Kind::Addition);
	
	// s->includeOperand(u->operand(0)->copy());

	// for(unsigned int i=1; i<u->numberOfOperands(); i++) 
	// {
	// 	AST* t = mul({integer(-1), u->operand(i)->copy()});

	// 	s->includeOperand(reduceMultiplicationAST(t));

	// 	delete t;
	// }

	AST* r = reduceAdditionAST(s);

	delete s;

	return r;
}

}
