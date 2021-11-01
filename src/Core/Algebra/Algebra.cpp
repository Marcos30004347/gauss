#include <assert.h>
#include <algorithm>
#include <string.h>
#include <cstdio>
#include "Algebra.hpp"

#include "Core/Simplification/Rationals.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Rational/Rational.hpp"
// #include "Core/Debug/Assert.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/List.hpp"

using namespace ast;
using namespace polynomial;
using namespace rational;
using namespace simplification;

namespace algebra {

AST* integer(Int val) {
	return new AST(Kind::Integer, val);
}

AST* symbol(const char* identifier) {
	return new AST(Kind::Symbol, identifier);
}

AST* fraction(Int n, Int d) {
	return new AST(Kind::Fraction, {
		new AST(Kind::Integer, n),
		new AST(Kind::Integer, d)
	});
}

AST* fraction(AST* n, AST* d) {
	assert(isConstant(n));
	assert(isConstant(d));
	return new AST(Kind::Fraction, { n, d });
}

AST* add(std::vector<AST*> terms) {
	return new AST(Kind::Addition, terms);
}

AST* sub(std::vector<AST*> terms) {
	return new AST(Kind::Subtraction, terms);
}

AST* mul(std::vector<AST*> terms) {
	return new AST(Kind::Multiplication, terms);
}

AST* div(AST* numerator, AST* denominator) {
	return new AST(Kind::Division, { numerator, denominator });
}

AST* power(AST* bas, AST* expoent) {
	return new AST(Kind::Power, { bas, expoent });
}

AST* factorial(AST* u) {
	return new AST(Kind::Factorial, {
		u,
	});
}

bool isConstant(AST* u) {
	if(
		u->kind() == Kind::Symbol ||
		u->kind() == Kind::Infinity ||
		u->kind() == Kind::MinusInfinity ||
		u->kind() == Kind::Undefined
	) return false;
	
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction
	) return true;
	
	for(unsigned int i=0; i<u->numberOfOperands(); i++) {
		if(isConstant(u->operand(i)))
			return false;
	}

	return false; 
}

long gcd(long a, long b) {
	if (b == 0)
		return a;

	return gcd(b, a % b);
}



AST* base(AST* u) {
	if(u->kind() == Kind::Power) 
		return u->operand(0)->copy();
	
	return u->copy();
}

AST* expoent(AST* u) {
	if(u->kind() == Kind::Power) 
		return u->operand(1)->copy();

	return integer(1);
}


bool isRNE(AST* u) {
	if(u->kind() == Kind::Integer)
		return true;

	if(u->kind() == Kind::Fraction)
		return isConstant(u->operand(0)) && isConstant(u->operand(1));

	if(u->kind() == Kind::Addition && u->numberOfOperands() <= 2) {
		for(unsigned int i=0; i < u->numberOfOperands(); i++) {
	
			if(!isRNE(u->operand(i)))
				return false;
		}
		return true;
	}

	if(u->kind() == Kind::Subtraction && u->numberOfOperands() <= 2) {
		for(unsigned int i=0; i < u->numberOfOperands(); i++)
			if(!isRNE(u->operand(i)))
				return false;
		return true;
	}

	if(u->kind() == Kind::Multiplication && u->numberOfOperands() == 2) {
		for(unsigned int i=0; i < u->numberOfOperands(); i++)
			if(!isRNE(u->operand(i)))
				return false;
		return true;
	}

	if(u->kind() == Kind::Division) {
		for(unsigned int i=0; i < u->numberOfOperands(); i++)
			if(!isRNE(u->operand(i)))
				return false;
		return true;
	}

	if(u->kind() == Kind::Power) {
		AST* b = base(u);
		AST* e = expoent(u);
	
		bool is_rne = isRNE(b) && e->kind() == Kind::Integer;

		destroyASTs({b, e});
	
		return is_rne;
	}

	return false;
}

bool compareSymbols(std::string a, std::string b) {
    return std::lexicographical_compare(a.c_str(), a.c_str() + a.length(), b.c_str(), b.c_str() + b.length());
}

bool compareConstants(AST* u, AST* v) {
	if(u->kind() == Kind::Integer && v->kind() == Kind::Integer)
			return u->value() < v->value();

	AST* d = integerGCD(u, v);
	AST* num_u = numerator(u);
	AST* num_v = numerator(v);

	if(
			d->kind() == Kind::Integer &&
			num_u->kind() == Kind::Integer &&
			num_v->kind() == Kind::Integer
	) {
		AST* o_e = mul({d->copy(), num_u->copy()});
		AST* o_f = mul({d->copy(), num_v->copy()});
	
		AST* e = reduceRNEAST(o_e);
		AST* f = reduceRNEAST(o_f);
		
		bool res = e->value() < f->value();
		
		delete d;
		delete num_u;
		delete num_v;
		delete o_e;
		delete o_f;
		delete e;
		delete f;
		
		return  res;
	}

	delete d;
	delete num_u;
	delete num_v;

	return false;
}

bool compareProductsAndSummations(AST* u, AST* v) {
	unsigned int m = u->numberOfOperands() - 1;
	unsigned int n = v->numberOfOperands() - 1;

	for(unsigned int k=0; k <= std::min(m, n); k++) {
		if(!u->operand(m - k)->match(v->operand(n - k))) {
			return orderRelation(u->operand(m - k), v->operand(n - k));
		}
	}
	
	return m < n;
}

bool comparePowers(AST* u, AST* v) {
	AST* b_u = base(u);
	AST* b_v = base(v);

	if(!b_u->match(b_v)) {
		bool res = orderRelation(b_u, b_v);
		destroyASTs({b_u, b_v});
		return res; 
	}
	
	AST* e_u = expoent(u);
	AST* e_v = expoent(v);

	bool res = orderRelation(e_u, e_v);
	destroyASTs({e_u, e_v, b_u, b_v});
	return res;
}

bool compareFactorials(AST* u, AST* v) {
    return orderRelation(u->operand(0), v->operand(0));
}

bool compareFunctions(AST* u, AST* v) {
	if(u->funName() != v->funName())
		return std::lexicographical_compare(
			u->funName().c_str(), u->funName().c_str() + u->funName().length(),
			v->funName().c_str(), v->funName().c_str() + v->funName().length()
		);

	AST* argsu = u->operand(0);
	AST* argsv = v->operand(0);

	if(argsu->numberOfOperands() >= 1 && argsv->numberOfOperands() >= 1) {
		unsigned int m = argsu->numberOfOperands() - 1;
		unsigned int n = argsv->numberOfOperands() - 1;

		for(unsigned int k=0; k <= std::min(m, n); k++) {
			if(!argsu->operand(m - k)->match(argsv->operand(n - k))) {
				bool res = orderRelation(argsu->operand(m - k), argsv->operand(n - k));
				return res;
			}
		}
		
		return m < n;
	}

	// destroyASTs({ argsu, argsv });
	return true;
}

bool orderRelation(AST* u, AST* v) {

	if(u->kind() == Kind::Infinity)
		return true;
	if(v->kind() == Kind::Infinity)
		return false;
	if(u->kind() == Kind::MinusInfinity)
		return true;
	if(v->kind() == Kind::MinusInfinity)
		return false;


 	if(isConstant(u) && isConstant(v))
		return compareConstants(u, v);

	if(u->kind() == Kind::Symbol && v->kind() == Kind::Symbol)
			return compareSymbols(u->identifier(), v->identifier());

	if(u->kind() == Kind::Addition && v->kind() == Kind::Addition)
			return compareProductsAndSummations(u, v);

	if(u->kind() == Kind::Multiplication && v->kind() == Kind::Multiplication)
			return compareProductsAndSummations(u, v);
	
	if(u->kind() == Kind::Power && v->kind() == Kind::Power)
			return comparePowers(u, v);

	if(u->kind() == Kind::Factorial && v->kind() == Kind::Factorial)
			return compareFactorials(u, v);

	if(u->kind() == Kind::FunctionCall && v->kind() == Kind::FunctionCall)
			return compareFunctions(u, v);

	if(isConstant(u))
			return true;

	// if(IsConstant(v))
	//     return false;

	if(
		u->kind() == Kind::Multiplication && (
		v->kind() == Kind::Power ||
		v->kind() == Kind::Addition ||
		v->kind() == Kind::Factorial ||
		v->kind() == Kind::FunctionCall ||
		v->kind() == Kind::Symbol
	)) {
		AST* m = mul({v->copy()});
		bool res = orderRelation(u, m);
		destroyASTs({ m });
		return res;
	} 

	if(
		u->kind() == Kind::Power && (
		v->kind() == Kind::Addition ||
		v->kind() == Kind::Factorial ||
		v->kind() == Kind::FunctionCall ||
		v->kind() == Kind::Symbol
	)) {

		AST* m = power(v->copy(), integer(1));
		bool res = orderRelation(u, m);
		delete m;
		return res;
	} 

	if(
		u->kind() == Kind::Addition && (
		v->kind() == Kind::Factorial ||
		v->kind() == Kind::FunctionCall ||
		v->kind() == Kind::Symbol
	)) {
		AST* m = add({v->copy()});
		bool res = orderRelation(u, m);
		destroyASTs({ m });
		return res;
	}

	if(
		u->kind() == Kind::Factorial && (
		v->kind() == Kind::FunctionCall ||
		v->kind() == Kind::Symbol
	)) {
		if(u->operand(0)->match(v)) {
			return false;
		} else {
			AST* m = factorial(v->copy());
			bool res = orderRelation(u, m);
			destroyASTs({m});
			return res;
		}
	}

	if(u->kind() == Kind::FunctionCall && v->kind() == Kind::Symbol) {
		if(u->operand(0)->identifier() == v->identifier()) {
			return false;
		} else {
			return orderRelation(u->operand(0), v);
		}
	}

	return !orderRelation(v, u);
}

AST* binomial(Int n, std::vector<Int> ks) {
	AST* p = new AST(Kind::Multiplication);
	for(Int k : ks)
		p->includeOperand(factorial(integer(k)));
	return div(factorial(integer(n)), p);
}

AST* funCall(const char* id, std::vector<AST*> args) {
	AST* f = new AST(Kind::FunctionCall);
	f->includeOperand(symbol(id));
	for(AST* a : args)
		f->includeOperand(a);
	return f;
}

AST* integerGCD(AST*  a, AST*  b) {
	if (a->value() == 0)
		return b->copy();
	
	AST* b_ = integer(b->value() % a->value());

	AST* gcd = integerGCD(b_, a);

	delete b_;

	return gcd;
}



AST* min(AST* a, AST* b) {
	if(a->kind() != Kind::Integer || b->kind() != Kind::Integer)
		return undefined();

	if(a->value() > b->value())
		return b->copy();
	return a->copy();
}

AST* max(AST* a, AST* b) {
	if(a->kind() != Kind::Integer || b->kind() != Kind::Integer)
		return undefined();
	
	if(a->value() > b->value())
		return a->copy();
	return b->copy();
}

AST* undefined() {
	return new AST(Kind::Undefined);
}
AST* fail() {
	return new AST(Kind::Fail);
}
bool isGreaterZero(AST* u) {
	AST* t = algebraicExpand(u);
	
	bool r = false;
	
	if(t->kind() == Kind::Integer) {
		r = t->value() > 0;
	} else
	if(t->kind() == Kind::Fraction) {
		AST* n = t->operand(0);
		AST* d = t->operand(1);

		r = (n->value() > 0 && d->value() > 0) ||
				(n->value() < 0 && d->value() < 0);
	}

	delete t;

	return r;
}

bool isLessEqZero(ast::AST* u) {
	return isEqZero(u) || isLessZero(u);
}

bool isLessZero(ast::AST* u) {
	AST* t = algebraicExpand(u);
	
	bool r = false;

	if(t->kind() == Kind::Integer) 
	{
		r = t->value() < 0;
	} 
	else
	if(t->kind() == Kind::Fraction) 
	{
		AST* n = t->operand(0);
		AST* d = t->operand(1);

		r = (n->value() < 0 && d->value() > 0) ||
				(n->value() > 0 && d->value() < 0);
	}

	delete t;

	return r;
}

bool isEqZero(ast::AST* u) {
	AST* t = algebraicExpand(u);
	
	bool r = false;
	
	if(t->kind() == Kind::Integer)
		r = t->value() == 0;

	delete t;

	return r;
}


bool isGreaterEqZero(AST* u) {
	return isEqZero(u) || isGreaterZero(u);
}


AST* completeSubExpressions(AST* u) {
	if(u->isTerminal())
		return set({ u->copy() });

	AST* S = set({u->copy()});
	for(unsigned int i=0; i<u->numberOfOperands(); i++) {
		AST* S_ = unification(S, completeSubExpressions(u->operand(i)));
		delete S;
		S = S_;
	}

	return S;
}
bool isDivisionByZero(AST* k) {
	AST* d = denominator(k);
	
	if(d->kind() == Kind::Integer && d->value() == 0) {
		delete d;
		return true;
	}

	delete d;
	return false;
}

int mod(int a, int b) {
	return (b + (a % b)) % b;
}

AST* leastCommomMultiple(AST* a, AST* b)
{
	return integer(
		abs(
			a->value() * b->value()
		).abs() / gcd(
			a->value(), b->value()
		)
	);
}

AST* leastCommomMultiple(AST* l)
{
	assert(l->kind() == Kind::List);

	if(l->numberOfOperands() == 2)
	{
		assert(l->operand(0)->kind() == Kind::Integer);
		assert(l->operand(0)->value() != 0);
		assert(l->operand(1)->kind() == Kind::Integer);
		assert(l->operand(1)->value() != 0);
	
		return leastCommomMultiple(l->operand(0), l->operand(1));
	}

	// lcm(b0, ... bn) = lcm(lcm(b0, ..., bn-1), bn)
	AST* j = rest(l);
	AST* a = first(l);
	AST* b = leastCommomMultiple(j);
	
	AST* lcm = leastCommomMultiple(a, b);

	delete j;
	delete a;
	delete b;

	return lcm;
}


std::pair<ast::AST*, ast::AST*> linearForm(ast::AST* u, ast::AST* x)
{
	if(u->match(x))
	{
		return {integer(1), integer(0)};
	}
	
	if(
		u->kind() == Kind::Symbol  ||
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction
	)
	{
		return {integer(0), u->copy()};
	}

	if(u->kind() == Kind::Multiplication)
	{
		if(u->freeOf(x))
		{
			return {integer(0), u->copy()};
		}

		AST* t = div(u->copy(), x->copy());
		AST* k = algebraicExpand(t);
	
		delete t;
		
		if(k->freeOf(x))
		{
			return { k, integer(0) };
		}
	
		delete k;
	
		return { nullptr, nullptr };
	}

	if(u->kind() == Kind::Addition)
	{
		std::pair<AST*, AST*> f = linearForm(u->operand(0), x);

		if(f.first == nullptr && f.second == nullptr)
		{
			return { nullptr, nullptr };
		}
	
		AST* t = sub({u->copy(), u->operand(0)->copy()});
		AST* k = algebraicExpand(t);
	
		std::pair<AST*, AST*> r = linearForm(k, x);

		delete t;
		delete k;
	
		if(r.first == nullptr && r.second == nullptr)
		{
			return { nullptr, nullptr };
		}
		
		AST* l = add({ f.first->copy(), r.first->copy() });
		AST* p = add({ f.second->copy(), r.second->copy() });
		
		AST* s = reduceAST(l);
		AST* z = reduceAST(p);
	
		delete l;
		delete p;
	
		return { s, z };
	}

	if(u->freeOf(x))
	{
		return { integer(0), u->copy() };
	}

	return { nullptr, nullptr };
}


AST* sinh(AST* x)
{
	return funCall("sinh", { x->copy() });
}


AST* cosh(AST* x)
{
	return funCall("cosh", { x->copy() });
}

AST* tanh(AST* x)
{
	return funCall("tanh", { x->copy() });
}

AST* exp(AST* x)
{
	return funCall("exp", { x->copy() });
}

AST* cos(AST* x)
{
	return funCall("cos", { x->copy() });
}

AST* sin(AST* x)
{
	return funCall("sin", { x->copy() });
}

AST* tan(AST* x)
{
	return funCall("tan", { x->copy() });
}

AST* csc(AST* x)
{
	return funCall("csc", { x->copy() });
}

AST* cot(AST* x)
{
	return funCall("cot", { x->copy() });
}

AST* log(AST* x)
{
	return funCall("log", { x->copy() });
}

AST* ln(AST* x)
{
	return funCall("ln", { x->copy() });
}

AST* sec(AST* x)
{
	return funCall("sec", { x->copy() });
}

AST* coth(AST* x)
{
	return funCall("coth", { x->copy() });
}

AST* sech(AST* x)
{
	return funCall("sech", { x->copy() });
}

AST* csch(AST* x)
{
	return funCall("csch", { x->copy() });
}

ast::AST* abs(ast::AST* x)
{
	return funCall("abs", { x->copy() });
}

AST* arccos(AST* x)
{
	return funCall("arccos", { x->copy() });
}

AST* arcsin(AST* x)
{
	return funCall("arcsin", { x->copy() });
}

AST* arctan(AST* x)
{
	return funCall("arctan", { x->copy() });
}

AST* arccot(AST* x)
{
	return funCall("arccot", { x->copy() });
}

AST* arcsec(AST* x)
{
	return funCall("arcsec", { x->copy() });
}

AST* arccsc(AST* x)
{
	return funCall("arccsc", { x->copy() });
}

AST* arccosh(AST* x)
{
	return funCall("arccosh", { x->copy() });
}

AST* arctanh(AST* x)
{
	return funCall("arctanh", { x->copy() });
}

AST* matrix(AST* rows, AST* cols)
{
	// assert(
	// 	rows->kind() == Kind::Integer, 
	// 	"matrix rows needs to be an integer"
	// );

	// assert(
	// 	cols->kind() == Kind::Integer, 
	// 	"matrix cols needs to be an integer"
	// );

	AST* m = new AST(Kind::Matrix);

	for (unsigned int i = 0; i < rows->value(); i++)
	{
		std::vector<AST*> r;

		for (unsigned int j = 0; j < cols->value(); j++)
		{
			r.push_back(integer(0));
		}
		
		m->includeOperand(list(r));
	}

	return m;
}

AST* matrix(std::vector<AST*> t)
{
	return new AST(Kind::Matrix, t);
}

bool isGreatherThan(ast::AST* a, ast::AST* b)
{
	if(
		a->kind() == Kind::Undefined ||
		a->kind() == Kind::Fail
	)
	{
		printf("Comparison with 'Undefined' or 'Fail' is illegal!\n");
		abort();
	}

	AST* u = algebraicExpand(a);
	AST* v = algebraicExpand(b);

	if(u->kind() == Kind::Infinity)
	{
		delete u;
		delete v;

		return true;
	}

	if(u->kind() == Kind::MinusInfinity)
	{
		delete u;
		delete v;

		return false;
	}

	if(
		a->kind() == Kind::Symbol && 
		b->kind() == Kind::Symbol
	)
	{
		return false;
	}

	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction 
	) {

		if(
			b->kind() != Kind::Integer || 
			b->kind() != Kind::Fraction
		) {
			delete u;
			delete v;
			
			return false;
		}

		AST* nu = numerator(u);
		AST* du = denominator(u);

		AST* nv = numerator(v);
		AST* dv = denominator(v);
	
		Int t0 = nu->value();
		Int t1 = du->value();
		Int t2 = nv->value();
		Int t3 = dv->value();

		Int Y = t0 * t3 - t1 * t2;
	
		delete nu;
		delete du;
		delete nv;
		delete dv;
	
		return Y > 0;	
	}

	if(u->kind() == Kind::Division)
	{
		AST* a = numerator(u);
		AST* b = denominator(u);

		AST* c = numerator(v);
		AST* d = denominator(v);

		AST* r = mulPoly(a, d);
		AST* k = mulPoly(c, b);

		return isGreatherThan(r, k);
	}
	
	// Addition,
	// Subtraction,
	// Multiplication,
	// Power,
	// Factorial,

	// FunctionCall,

	// Integral,
	// Derivative,

	// List,
	// Set,



	if(u->kind() == Kind::Addition)
	{
	}

	// AST* big_deg_u = list({});
	// AST* big_deg_v = list({});

	// if(u->kind() == Kind::Division)
	// {
	// 	AST* t0 = u->operand(0)->copy();
	// 	AST* t1 = u->operand(1)->copy();

	// 	AST* t2 = nullptr;
	// 	AST* t3 = nullptr;
	
	// 	if(v->kind() == Kind::Division)
	// 	{
	// 		t2 = v->operand(0)->copy();
	// 		t3 = v->operand(1)->copy();
	// 	}
	// 	else
	// 	{
	// 		t2 = v->copy();
	// 		t3 = integer(1);
	// 	}
	
	// 	AST* k = sub({ mul({ t0, t3 }), mul({ t1, t2 }) });
	
	// 	AST* z = integer(0);
	
	// 	bool res = isGreatherThan(k, z);

	// 	delete u;
	// 	delete v;
	// 	delete k;
	// 	delete z;
	
	// 	return res;		
	// }

	// if(
	// 	u->kind() == Kind::Symbol && (
	// 		v->kind() == Kind::Integer  ||
	// 		v->kind() == Kind::Fraction
	// ))
	// {
	// 	return true;
	// }

	// if(a->kind() == Kind::Power)
	// {
		
	// }

	// TODO: Tensor,
	// TODO: Matrix,



	return !isGreatherThan(b, a);

}

bool isLessThan(ast::AST* a, ast::AST* b);
bool isGreatherOrEqualThan(ast::AST* a, ast::AST* b);
bool isLessOrEqualThan(ast::AST* a, ast::AST* b);



ast::AST* getSymbols(ast::AST* u)
{
	if(u->kind() == Kind::Symbol)
	{
		return new AST(Kind::Set, { u->copy() });
	}

	AST* syms = set({});

	if(
		u->kind() == Kind::Addition       ||
		u->kind() == Kind::Subtraction    ||
		u->kind() == Kind::Power          ||
		u->kind() == Kind::Division		   ||
		u->kind() == Kind::Multiplication ||
		u->kind() == Kind::Matrix			   ||
		u->kind() == Kind::Set			   		 ||
		u->kind() == Kind::List
	)
	{
		for(unsigned int i = 0; i < u->numberOfOperands(); i++)
		{
			AST* s = getSymbols(u->operand(i));

			if(s->numberOfOperands() > 0)
			{
				for(unsigned int k = 0; k < s->numberOfOperands(); k++)
				{
					AST* t = unification(syms, s);
					delete syms;
					syms = t;
				}
			}

			delete s;
		}
	}

	if(
		u->kind() == Kind::Derivative ||
		u->kind() == Kind::Integral   ||
		u->kind() == Kind::Factorial
	)
	{
		AST* s = getSymbols(u->operand(0));
		if(s->numberOfOperands() > 0)
		{
			for(unsigned int k = 0; k < s->numberOfOperands(); k++)
			{
				AST* t = unification(syms, s);
				delete syms;
				syms = t;
			}
		}

		delete s;
	}

	return syms;
}

long fat(long a)
{
	long f = 1;

	while(a > 1)
	{
		f = f * a;
		a = a - 1;
	}

	return f;
}


} // algebra
