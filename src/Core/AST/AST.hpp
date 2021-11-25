#ifndef MATH_ABSTRACT_SYNTAX_TREE_AST_H
#define MATH_ABSTRACT_SYNTAX_TREE_AST_H

#include <cstddef>
#include <string>
#include <vector>

#include "Integer.hpp"

namespace ast {

enum Kind {
  Undefined = 0,
  Integer,
  Symbol,
  Infinity,
  MinusInfinity,
  Fraction,
  Tensor,
  Matrix,
  Fail,

  Addition,
  Subtraction,
  Multiplication,
  Division,
  Power,
  Factorial,

  FunctionCall,

  Integral,
  Derivative,

  List,
  Set,

};

class Expr {
public:
  Expr();
  Expr(Kind kind);
  Expr(Kind kind, Int value);
  Expr(Kind kind, const char *identifier);
  Expr(Kind kind, std::vector<Expr> operands);
	Expr(Int v);
	Expr(int v);
	Expr(long long v);
	Expr(std::string v);

  Kind kind() const;

  std::string toString();

  bool match(Expr other);
  bool freeOf(Expr other);
  bool freeOfElementsInSet(Expr const set);
  bool isTerminal();

  // bool is(Int i);
  // bool isNot(Int i);

  Expr operandList();

  Expr symbols();

  bool insert(Expr expr);
  bool insert(Expr expr, signed long i);
  bool insert(Expr expr, Int i);

  //bool remove(Expr expr);
  bool remove(signed long i);
  bool remove(Int i);

  unsigned size() const;
  Int value() const;

  const std::string identifier();
  const std::string funName();

  Expr(const Expr &other);
  Expr(const Expr &&other);

  Expr operator+(Expr &&other);
  Expr operator+(Expr &other);
  Expr operator-(Expr &&other);
  Expr operator-(Expr &other);
  Expr operator*(Expr &&other);
  Expr operator*(Expr &other);
  Expr operator/(Expr &&other);
  Expr operator/(Expr &other);
  Expr operator[](size_t idx);
  Expr operator[](Int idx);
  Expr operator-();
  Expr operator+();

  Expr operator+=(Expr &other);
  Expr operator+=(Expr &&other);
  Expr operator-=(Expr &other);
  Expr operator-=(Expr &&other);

  bool operator==(Expr &other);
  bool operator==(Expr &&other);
  //bool operator==(std::string other);

	bool operator!=(Expr &other);
  bool operator!=(Expr &&other);
  //bool operator!=(const Int other);
  //bool operator!=(std::string other);

	bool operator>(const Int other);
	bool operator>=(const Int other);
	bool operator<(const Int other);
	bool operator<=(const Int other);

	//Expr operator=(Int i);
  //Expr operator=(std::string);
  Expr operator=(const Expr &i);
  Expr operator=(const Expr &&i);

	std::vector<Expr> operands() const;
private:
  std::vector<Expr> _operands;
  Kind _kind;
  Int _value;
  std::string _identifier;

  Expr operand(unsigned long i);
  Expr operand(Int i);

  Expr(Kind kind, const Int value, const std::string identifier);
};

Expr operator*(Int i, Expr &&other);
Expr operator*(Int i, Expr &other);
Expr operator+(Int i, Expr &&other);
Expr operator+(Int i, Expr &other);
Expr operator-(Int i, Expr &&other);
Expr operator-(Int i, Expr &other);
Expr operator/(Int i, Expr &&other);
Expr operator/(Int i, Expr &other);
Expr operator*(int i, Expr &&other);
Expr operator*(int i, Expr &other);
Expr operator+(int i, Expr &&other);
Expr operator+(int i, Expr &other);
Expr operator-(int i, Expr &&other);
Expr operator-(int i, Expr &other);
Expr operator/(int i, Expr &&other);
Expr operator/(int i, Expr &other);

// void destroyExprs(std::vector<Expr>);
Expr mapBinaryAST(Expr a, Expr n, Expr (*)(Expr, Expr));
Expr mapUnaryAST(Expr u, Expr (*f)(Expr));
Expr deepReplace(Expr tree, Expr subtree, Expr v);
Expr construct(Kind kind, Expr L);

template <typename... types>
Expr mapExpr(Expr (*f)(Expr, types... args), Expr u, types... params) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Fraction ||
      u.kind() == Kind::Symbol || u.kind() == Kind::Infinity ||
      u.kind() == Kind::MinusInfinity) {
    return f(u, params...);
  }

  if (u.size() == 0) {
    return f(u, params...);
  }

  Expr t = Expr(u.kind());

  if (u.kind() == Kind::FunctionCall) {
    t.insert(Expr(Kind::Symbol, u.funName().c_str()));
  }

  for (unsigned int i = 0; i < u.size(); i++) {
    t.insert(f(u[i]), params...);
  }

  return t;
}

Expr undefined();
Expr fail();
Expr inf();
Expr minInf();

} // namespace ast

#endif
