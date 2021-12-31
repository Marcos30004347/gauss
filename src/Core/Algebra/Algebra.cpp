#include "Algebra.hpp"
#include <algorithm>
#include <assert.h>
#include <cstdio>
#include <string.h>
#include <tuple>
#include <utility>
#include <vector>

#include "Core/AST/AST.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Rational/Rational.hpp"
#include "Core/Simplification/Rationals.hpp"
#include "Core/Simplification/Simplification.hpp"

using namespace ast;
using namespace polynomial;
using namespace rational;
using namespace simplification;

namespace algebra {

Expr gcdConstants(Expr a, Expr b) {
	assert(a.kind() == Kind::Integer || a.kind() == Kind::Fraction);
	assert(b.kind() == Kind::Integer || b.kind() == Kind::Fraction);

	if(a.kind() == Kind::Integer && b.kind() == Kind::Integer) {
		return abs(gcd(a.value(), b.value()));
	}

	Expr A = numerator(a);
	Expr C = numerator(b);
	Expr B = denominator(a);
	Expr D = denominator(b);

	assert(A.kind() == Kind::Integer);
	assert(B.kind() == Kind::Integer);
	assert(C.kind() == Kind::Integer);
	assert(D.kind() == Kind::Integer);

	Int t = abs(gcd(A.value(), C.value()));
	Int k = abs(lcm(B.value(), D.value()));

	if(k == 1) return t;

	Int j = gcd(t, k);

	if(k / j == 1) return t;

	return fraction(t / j, k / j);
}


Expr lcmConstants(Expr a, Expr b) {
	assert(a.kind() == Kind::Integer || a.kind() == Kind::Fraction);
	assert(b.kind() == Kind::Integer || b.kind() == Kind::Fraction);

	if(a.kind() == Kind::Integer && b.kind() == Kind::Integer) {
		return abs(lcm(a.value(), b.value()));
	}

	Expr A = numerator(a);
	Expr C = numerator(b);
	Expr B = denominator(a);
	Expr D = denominator(b);

	assert(A.kind() == Kind::Integer);
	assert(B.kind() == Kind::Integer);
	assert(C.kind() == Kind::Integer);
	assert(D.kind() == Kind::Integer);

	Int t = abs(lcm(A.value(), C.value()));
	Int k = abs(gcd(B.value(), D.value()));

	if(k == 1) return t;

	Int j = gcd(t, k);

	if(k / j == 1) return t;

	return fraction(t / j, k / j);
}



Expr integer(Int val) { return Expr(Kind::Integer, val); }

Expr symbol(const char *identifier) { return Expr(Kind::Symbol, identifier); }

Expr fraction(Expr n, Expr d) {
  assert(isConstant(n));
  assert(isConstant(d));
  return Expr(Kind::Fraction, {n, d});
}

Expr add(std::vector<Expr> terms) { return Expr(Kind::Addition, terms); }

Expr sub(std::vector<Expr> terms) { return Expr(Kind::Subtraction, terms); }

Expr mul(std::vector<Expr> terms) { return Expr(Kind::Multiplication, terms); }

Expr div(Expr numerator, Expr denominator) {
  return Expr(Kind::Division, {numerator, denominator});
}

Expr power(Expr bas, Expr expoent) { return Expr(Kind::Power, {bas, expoent}); }

Expr factorial(Expr u) {
  return Expr(Kind::Factorial, {
                                   u,
                               });
}

bool isConstant(Expr u) {
  if (u.kind() == Kind::Symbol || u.kind() == Kind::Infinity ||
      u.kind() == Kind::MinusInfinity || u.kind() == Kind::Undefined)
    return false;

  if (u.kind() == Kind::Integer || u.kind() == Kind::Fraction)
    return true;

  if (u.kind() == Kind::Division) {
    return isConstant(u[0]) && isConstant(u[1]);
  }

  for (size_t i = 0; i < u.size(); i++) {
    if (isConstant(u[i]))
      return false;
  }

  return false;
}

long gcd(long a, long b) {
  if (b == 0)
    return a;

  return gcd(b, a % b);
}

Expr base(Expr u) {
  if (u.kind() == Kind::Power)
    return u[0];

  return u;
}

Expr expoent(Expr u) {
  if (u.kind() == Kind::Power)
    return u[1];

  return integer(1);
}

bool isRNE(Expr u) {
  if (u.kind() == Kind::Integer)
    return true;

  if (u.kind() == Kind::Fraction)
    return isConstant(u[0]) && isConstant(u[1]);

  if (u.kind() == Kind::Addition && u.size() <= 2) {
    for (unsigned int i = 0; i < u.size(); i++) {

      if (!isRNE(u[i]))
        return false;
    }
    return true;
  }

  if (u.kind() == Kind::Subtraction && u.size() <= 2) {
    for (unsigned int i = 0; i < u.size(); i++)
      if (!isRNE(u[i]))
        return false;
    return true;
  }

  if (u.kind() == Kind::Multiplication && u.size() == 2) {
    for (unsigned int i = 0; i < u.size(); i++)
      if (!isRNE(u[i]))
        return false;
    return true;
  }

  if (u.kind() == Kind::Division) {
    for (unsigned int i = 0; i < u.size(); i++)
      if (!isRNE(u[i]))
        return false;
    return true;
  }

  if (u.kind() == Kind::Power) {
    Expr b = base(u);
    Expr e = expoent(u);

    bool is_rne = isRNE(b) && e.kind() == Kind::Integer;

    return is_rne;
  }

  return false;
}

bool compareSymbols(std::string a, std::string b) {
  return std::lexicographical_compare(a.c_str(), a.c_str() + a.length(),
                                      b.c_str(), b.c_str() + b.length());
}

bool compareConstants(Expr &u, Expr &v) {
  if (u.kind() == Kind::Integer && v.kind() == Kind::Integer)
    return u.value() < v.value();
  Expr d = gcd(u.value(), v.value());
  Expr num_u = numerator(u);
  Expr num_v = numerator(v);

  if (d.kind() == Kind::Integer && num_u.kind() == Kind::Integer &&
      num_v.kind() == Kind::Integer) {
    Expr e = reduceRNEAST(d * num_u);
    Expr f = reduceRNEAST(d * num_v);

    bool res = e.value() < f.value();

    return res;
  }

  return false;
}

bool compareProductsAndSummations(Expr &u, Expr &v) {
  unsigned int m = u.size() - 1;
  unsigned int n = v.size() - 1;

  for (unsigned int k = 0; k <= std::min(m, n); k++) {
    if (u[m - k] != v[n - k]) {
      return orderRelation(u[m - k], v[n - k]);
    }
  }

  return m < n;
}

bool comparePowers(Expr &u, Expr &v) {
  if (u[0] != v[0]) {
    return orderRelation(u[0], v[0]);
  }

  return orderRelation(u[1], v[1]);
}

bool compareFactorials(Expr &u, Expr &v) { return orderRelation(u[0], v[0]); }

bool compareFunctions(Expr &u, Expr &v) {
  if (u.funName() != v.funName())
    return std::lexicographical_compare(
        u.funName().c_str(), u.funName().c_str() + u.funName().length(),
        v.funName().c_str(), v.funName().c_str() + v.funName().length());

  Expr argsu = u[0];
  Expr argsv = v[0];

  if (argsu.size() >= 1 && argsv.size() >= 1) {
    unsigned int m = argsu.size() - 1;
    unsigned int n = argsv.size() - 1;

    for (unsigned int k = 0; k <= std::min(m, n); k++) {
      if (argsu[m - k] != argsv[n - k]) {
        bool res = orderRelation(argsu[m - k], argsv[n - k]);
        return res;
      }
    }

    return m < n;
  }

  return true;
}

bool orderRelation(Expr &u, Expr &v) {
  if (u.kind() == Kind::Infinity)
    return true;
  if (v.kind() == Kind::Infinity)
    return false;
  if (u.kind() == Kind::MinusInfinity)
    return true;
  if (v.kind() == Kind::MinusInfinity)
    return false;

  if (isConstant(u) && isConstant(v))
    return compareConstants(u, v);

  if (u.kind() == Kind::Symbol && v.kind() == Kind::Symbol)
    return compareSymbols(u.identifier().c_str(), v.identifier().c_str());

  if (u.kind() == Kind::Addition && v.kind() == Kind::Addition)
    return compareProductsAndSummations(u, v);

  if (u.kind() == Kind::Multiplication && v.kind() == Kind::Multiplication)
    return compareProductsAndSummations(u, v);

  if (u.kind() == Kind::Power && v.kind() == Kind::Power)
    return comparePowers(u, v);

  if (u.kind() == Kind::Factorial && v.kind() == Kind::Factorial)
    return compareFactorials(u, v);

  if (u.kind() == Kind::FunctionCall && v.kind() == Kind::FunctionCall)
    return compareFunctions(u, v);

  if (isConstant(u))
    return false;

  if (u.kind() == Kind::Multiplication &&
      (v.kind() == Kind::Power || v.kind() == Kind::Addition ||
       v.kind() == Kind::Factorial || v.kind() == Kind::FunctionCall ||
       v.kind() == Kind::Symbol)) {
    Expr m = mul({v});
    return orderRelation(u, m);
  }

  if (u.kind() == Kind::Power &&
      (v.kind() == Kind::Addition || v.kind() == Kind::Factorial ||
       v.kind() == Kind::FunctionCall || v.kind() == Kind::Symbol)) {

    Expr m = power(v, 1);
    bool res = orderRelation(u, m);

    return res;
  }

  if (u.kind() == Kind::Addition &&
      (v.kind() == Kind::Factorial || v.kind() == Kind::FunctionCall ||
       v.kind() == Kind::Symbol)) {
    Expr m = add({v});
    bool res = orderRelation(u, m);
    return res;
  }

  if (u.kind() == Kind::Factorial &&
      (v.kind() == Kind::FunctionCall || v.kind() == Kind::Symbol)) {
    if (u[0].match(v)) {
      return false;
    } else {
      Expr m = factorial(v);
      bool res = orderRelation(u, m);
      return res;
    }
  }

  if (u.kind() == Kind::FunctionCall && v.kind() == Kind::Symbol) {
    if (u[0].identifier() == v.identifier()) {
      return false;
    } else {
      return orderRelation(u[0], v);
    }
  }

  return !orderRelation(v, u);
}

Expr binomial(Int n, std::vector<Int> ks) {
  Expr p = Expr(Kind::Multiplication);
  for (Int k : ks)
    p.insert(factorial(integer(k)));
  return div(factorial(integer(n)), p);
}

Expr funCall(const char *id, std::vector<Expr> args) {
  Expr f = Expr(Kind::FunctionCall);
  f.insert(symbol(id));
  for (Expr a : args) {
    f.insert(a);
  }
  return f;
}

Expr min(Expr a, Expr b) {
  if (a.kind() != Kind::Integer || b.kind() != Kind::Integer)
    return undefined();

  if (a.value() > b.value())
    return b;
  return a;
}

Expr max(Expr a, Expr b) {
  if (a.kind() != Kind::Integer || b.kind() != Kind::Integer)
    return undefined();

  if (a.value() > b.value())
    return a;
  return b;
}

bool isGreaterZero(Expr u) {
  Expr t = algebraicExpand(u);

  bool r = false;

  if (t.kind() == Kind::Integer) {
    r = t.value() > 0;
  } else if (t.kind() == Kind::Fraction) {
    Expr n = t[0];
    Expr d = t[1];

    r = (n.value() > 0 && d.value() > 0) || (n.value() < 0 && d.value() < 0);
  }

  return r;
}

Expr completeSubExpressions(Expr u) {
  if (u.isTerminal())
    return set({u});

  Expr S = set({u});
  for (unsigned int i = 0; i < u.size(); i++) {
    Expr S_ = unification(S, completeSubExpressions(u[i]));

    S = S_;
  }

  return S;
}
bool isDivisionByZero(Expr k) {
  Expr d = denominator(k);

  if (d.kind() == Kind::Integer && d.value() == 0) {

    return true;
  }

  return false;
}

int mod(int a, int b) { return (b + (a % b)) % b; }

std::pair<ast::Expr, ast::Expr> linearForm(ast::Expr u, ast::Expr x) {
  if (u.match(x)) {
    return {integer(1), integer(0)};
  }

  if (u.kind() == Kind::Symbol || u.kind() == Kind::Integer ||
      u.kind() == Kind::Fraction) {
    return {integer(0), u};
  }

  if (u.kind() == Kind::Multiplication) {
    if (u.freeOf(x)) {
      return {integer(0), u};
    }

    Expr t = div(u, x);
    Expr k = algebraicExpand(t);

    if (k.freeOf(x)) {
      return {k, integer(0)};
    }

    return {undefined(), undefined()};
  }

  if (u.kind() == Kind::Addition) {
    std::pair<Expr, Expr> f = linearForm(u[0], x);

    if (f.first == undefined() && f.second == undefined()) {
      return {undefined(), undefined()};
    }

    Expr t = sub({u, u[0]});
    Expr k = algebraicExpand(t);

    std::pair<Expr, Expr> r = linearForm(k, x);

    if (r.first == undefined() && r.second == undefined()) {
      return {undefined(), undefined()};
    }

    Expr l = add({f.first, r.first});
    Expr p = add({f.second, r.second});

    Expr s = reduceAST(l);
    Expr z = reduceAST(p);

    return {s, z};
  }

  if (u.freeOf(x)) {
    return {integer(0), u};
  }

  return {undefined(), undefined()};
}

Expr sinh(Expr x) { return funCall("sinh", {x}); }

Expr cosh(Expr x) { return funCall("cosh", {x}); }

Expr tanh(Expr x) { return funCall("tanh", {x}); }

Expr exp(Expr x) { return funCall("exp", {x}); }

Expr cos(Expr x) { return funCall("cos", {x}); }

Expr sin(Expr x) { return funCall("sin", {x}); }

Expr tan(Expr x) { return funCall("tan", {x}); }

Expr csc(Expr x) { return funCall("csc", {x}); }

Expr cot(Expr x) { return funCall("cot", {x}); }

Expr log(Expr x) { return funCall("log", {x}); }

Expr ln(Expr x) { return funCall("ln", {x}); }

Expr sec(Expr x) { return funCall("sec", {x}); }

Expr coth(Expr x) { return funCall("coth", {x}); }

Expr sech(Expr x) { return funCall("sech", {x}); }

Expr csch(Expr x) { return funCall("csch", {x}); }

ast::Expr abs(ast::Expr x) { return funCall("abs", {x}); }

Expr arccos(Expr x) { return funCall("arccos", {x}); }

Expr arcsin(Expr x) { return funCall("arcsin", {x}); }

Expr arctan(Expr x) { return funCall("arctan", {x}); }

Expr arccot(Expr x) { return funCall("arccot", {x}); }

Expr arcsec(Expr x) { return funCall("arcsec", {x}); }

Expr arccsc(Expr x) { return funCall("arccsc", {x}); }

Expr arccosh(Expr x) { return funCall("arccosh", {x}); }

Expr arctanh(Expr x) { return funCall("arctanh", {x}); }

Expr matrix(Expr rows, Expr cols) {
  Expr m = Expr(Kind::Matrix);

  for (unsigned int i = 0; i < rows.value(); i++) {
    std::vector<Expr> r;

    for (unsigned int j = 0; j < cols.value(); j++) {
      r.push_back(integer(0));
    }

    m.insert(list(r));
  }

  return m;
}

Expr matrix(std::vector<Expr> t) { return Expr(Kind::Matrix, t); }

ast::Expr getSymbols(ast::Expr u) {
  if (u.kind() == Kind::Symbol) {
    return Expr(Kind::Set, {u});
  }

  Expr syms = set({});

  if (u.kind() == Kind::Addition || u.kind() == Kind::Subtraction ||
      u.kind() == Kind::Power || u.kind() == Kind::Division ||
      u.kind() == Kind::Multiplication || u.kind() == Kind::Matrix ||
      u.kind() == Kind::Set || u.kind() == Kind::List) {
    for (unsigned int i = 0; i < u.size(); i++) {
      Expr s = getSymbols(u[i]);

      if (s.size() > 0) {
        for (unsigned int k = 0; k < s.size(); k++) {
          Expr t = unification(syms, s);

          syms = t;
        }
      }
    }
  }

  if (u.kind() == Kind::Derivative || u.kind() == Kind::Integral ||
      u.kind() == Kind::Factorial) {
    Expr s = getSymbols(u[0]);
    if (s.size() > 0) {
      for (unsigned int k = 0; k < s.size(); k++) {
        Expr t = unification(syms, s);

        syms = t;
      }
    }
  }

  return syms;
}
void mergeRec(std::vector<Expr> &L, std::vector<Expr> &temp, long l, long m,
              long &r, bool reversed) {

  // printf("\n******\n***** merging %li %li %li\n******\n", l ,m , r);
  size_t left_pos = l;
  size_t left_end = m;

  size_t temp_pos = l;

  size_t righ_end = r;
  size_t righ_pos = m + 1;

  while (left_pos <= left_end && righ_pos <= righ_end) {
    bool unordered = orderRelation(L[left_pos], L[righ_pos]);
    if (reversed)
      unordered = !unordered;

    if (unordered) {
      temp[temp_pos++] = std::move(L[left_pos++]);
    } else {
      temp[temp_pos++] = std::move(L[righ_pos++]);
    }
  }

  while (left_pos <= left_end) {
    temp[temp_pos++] = std::move(L[left_pos++]);
  }

  while (righ_pos <= righ_end) {
    temp[temp_pos++] = std::move(L[righ_pos++]);
  }

  size_t num = r - l + 1;

  for (size_t i = 0; i < num; i++) {
    L[righ_end] = std::move(temp[righ_end]);
    righ_end--;
  }
}

void sortRec(std::vector<Expr> &L, std::vector<Expr> &tmp, long l, long r,
             bool reversed) {
  if (l < r) {
    long m = l + (r - l) / 2;

    sortRec(L, tmp, l, m, reversed);
    sortRec(L, tmp, m + 1, r, reversed);
    mergeRec(L, tmp, l, m, r, reversed);
  }
}

void sort(std::vector<ast::Expr> &L, bool reversed) {
  std::vector<Expr> tmp(L.size(), 0);
  sortRec(L, tmp, 0, L.size() - 1, reversed);
}

} // namespace algebra
