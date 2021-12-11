#include "AST.hpp"
#include <assert.h>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

namespace ast {
Expr::Expr() { this->_kind = Kind::Undefined; }

Expr::Expr(Kind kind) {
  this->_kind = kind;
  this->_operands = std::vector<Expr>(0);
  this->_identifier = "";
  this->_value = 0;
}

Expr::Expr(Kind kind, Int value) {
  this->_kind = kind;
  this->_operands = std::vector<Expr>(0);
  this->_value = value;
  this->_identifier = "";
}

Expr::Expr(Kind kind, const char *identifier) {
  this->_kind = kind;
  this->_operands = std::vector<Expr>(0);
  this->_identifier = identifier;
  this->_value = 0;
}

Expr::Expr(Kind kind, std::vector<Expr> operands) {
  this->_kind = kind;
  this->_operands = operands;
  this->_identifier = "";
  this->_value = 0;
}

Expr::Expr(Kind kind, const Int value, const std::string identifier) {
  this->_kind = kind;
  this->_operands = std::vector<Expr>(0);
  this->_identifier = identifier;
  this->_value = value;
}

Expr::Expr(Int v) {
  this->_kind = Kind::Integer;
  this->_value = v;
}
Expr::Expr(int v) {
  this->_kind = Kind::Integer;
  this->_value = v;
}

Expr::Expr(long int v) {
  this->_kind = Kind::Integer;
  this->_value = v;
}

Expr::Expr(long long v) {
  this->_kind = Kind::Integer;
  this->_value = v;
}

Expr::Expr(std::string v) {
  this->_kind = Kind::Symbol;
  this->_identifier = v;
}
Kind Expr::kind() const { return this->_kind; }

	Expr &Expr::operand(Int i) { return this->operand(i.longValue()); }

Expr &Expr::operand(unsigned long i) {
  if (this->kind() == Kind::Integer || this->kind() == Kind::Symbol ||
      this->kind() == Kind::Infinity || this->kind() == Kind::MinusInfinity) {
    return *this;
  }

  if (this->kind() == Kind::FunctionCall) {
    i = i + 1;
  }

  return this->_operands[i];
}

bool Expr::insert(Expr &&expr) {
  if (this->kind() == Kind::Set) {
    for (unsigned int i = 0; i < this->size(); i++) {
      if (this->operand(i) == expr) {
        return false;
      }
    }
  }

  this->_operands.push_back(std::move(expr));

  return true;
}

bool Expr::insert(Expr &expr) {
  if (this->kind() == Kind::Set) {
    for (unsigned int i = 0; i < this->size(); i++) {
      if (this->operand(i) == expr) {
        return false;
      }
    }
  }

  this->_operands.push_back(expr);

  return true;
}

bool Expr::insert(Expr &expr, Int i) {
  return this->insert(expr, i.longValue());
}

bool Expr::insert(Expr &&expr, Int i) {
  return this->insert(expr, i.longValue());
}

std::vector<Expr> Expr::operands() const { return this->_operands; }

bool Expr::insert(Expr &&expr, signed long i) {
  if (expr.kind() == Kind::Set) {
    for (unsigned int i = 0; i < this->size(); i++) {
      if (this->operand(i) == expr)
        return false;
    }
  }

  std::vector<Expr>::iterator it = this->_operands.begin();
  std::advance(it, i);
  this->_operands.insert(it, std::move(expr));
  return true;
}

bool Expr::insert(Expr &expr, signed long i) {
  if (expr.kind() == Kind::Set) {
    for (unsigned int i = 0; i < this->size(); i++) {
      if (this->operand(i) == expr)
        return false;
    }
  }

  std::vector<Expr>::iterator it = this->_operands.begin();
  std::advance(it, i);
  this->_operands.insert(it, expr);
  return true;
}

bool Expr::remove(Int i) { return this->remove(i.longValue()); }

bool Expr::remove(signed long i) {
  this->_operands.erase(this->_operands.begin() + i);
  return true;
}

unsigned Expr::size() const {
  switch (this->kind()) {
  case Kind::Integer:
  case Kind::Fraction:
  case Kind::Symbol:
  case Kind::Infinity:
  case Kind::MinusInfinity:
    return 1;

  case Kind::FunctionCall:
    return this->_operands.size() - 1;

  default:
    return this->_operands.size();
  }
}

Int Expr::value() const { return this->_value; }

const std::string Expr::identifier() { return this->_identifier; }

bool Expr::match(Expr &other) { return this->match(std::forward<Expr>(other)); }

bool Expr::match(Expr &&other) {
  if (this->kind() != other.kind()) {
		//printf("A\n");
    return false;
  }

  if (this->size() != other.size()) {
		//printf("B\n");
		//printf("%s\n", this->operand(0).toString().c_str());
		//printf("%u\n", this->size());
		//printf("%u\n", other.size());

    return false;
  }

  if (this->kind() == Kind::FunctionCall) {
    if (this->funName() != other.funName()) {
      return false;
    }

		// for(Int i = 0; i < size(); i++) {
		// 	if(!operand(i).match(other.operand(i))) {
		// 		return false;
		// 	}
		// }
		// return true;
    // TODO: match arguments
  }

  if (this->kind() == Kind::Fraction) {
		//printf("B\n");
    return this->operand(0).match(other[0]) && this->operand(1).match(other[1]);
  }

  if (this->kind() == Kind::Integer) {
		//printf("C\n");
    return this->value() == other.value();
  }

  if (this->kind() == Kind::Symbol) {
		//printf("D\n");
    return this->identifier() == other.identifier();
  }

  if (this->kind() == Kind::Undefined) {
		//printf("E\n");
    return this->value() == other.value();
  }

  if (this->kind() == Kind::Factorial) {
		//printf("F\n");
    return this->operand(0).match(other[0]);
  }

  if (this->kind() == Kind::Division) {
		//printf("G\n");
    return this->operand(0).match(other[0]) && this->operand(1).match(other[1]);
  }

  if (this->kind() == Kind::Infinity || this->kind() == Kind::MinusInfinity) {
		//printf("H\n");
    return this->kind() == other.kind();
  }

  if (this->kind() == Kind::Subtraction) {
    unsigned int matches = 0;

    if (!this->operand(0).match(other[0])) {
      return false;
    }

    matches++;

    for (unsigned int i = 1; i < this->size(); i++) {
      for (unsigned int j = 1; j < other.size(); j++) {
        if (this->operand(i).match(other[j])) {
          matches++;
          break;
        }
      }
    }

		//printf("I\n");
    return matches == this->size();
  }

  if (this->kind() == Kind::Addition || this->kind() == Kind::Multiplication ||
      this->kind() == Kind::Set) {

    unsigned int matches = 0;

		std::vector<unsigned int> m;
		std::vector<unsigned int> o;

		for (unsigned int i = 0; i < this->size(); i++) {
      for (unsigned int j = 0; j < other.size(); j++) {
				if (this->operand(i).match(other[j])) {
					m.push_back(i);
					matches = matches + 1;
          break;
        } else {
					o.push_back(i);
					o.push_back(j);
				}
			}
    }

		//for(unsigned int k : m) printf("%u ", k);
		//printf("\n");
		//for (unsigned int k = 0; k < o.size(); k += 2) {
		//if (o[k] == o[k + 1]) {
				//printf("----> %u\n", o[k]);
				//printf("----> %s\n", operand(o[k]).toString().c_str());
				//printf("----> %u\n", operand(o[k]).kind());
				//	printf("----> %s\n", other.operand(o[k]).toString().c_str());
				//printf("----> %u\n", other.operand(o[k]).kind());
				//printf("%i\n", this->operand(o[k]).match(other.operand(o[k])));
		//}

		//}

		//printf("A(%u) = %s\n", this->size(), this->toString().c_str());
		//printf("B(%u) = %s\n", other.size(), other.toString().c_str());

		//printf("H %u\n", matches);
		//printf("H %u\n", this->size());

		return matches == this->size();
  }

  // order of the operators does matter
  for (unsigned int i = 0; i < this->size(); i++) {
    if (!this->operand(i).match(other[i])) {
			//printf("J\n");
      return false;
    }
  }

	//printf("K\n");
  return true;
}

bool Expr::isTerminal() {
  if (this->kind() == Kind::Integer || this->kind() == Kind::Fraction ||
      this->kind() == Kind::Infinity || this->kind() == Kind::MinusInfinity ||
      this->kind() == Kind::Symbol || this->kind() == Kind::Tensor)
    return true;

  return false;
}

bool Expr::freeOf(Expr &other) {
  return this->freeOf(std::forward<Expr>(other));
}

bool Expr::freeOf(Expr &&other) {
  if (this->match(other))
    return false;

  if (this->kind() == Kind::Integer || this->kind() == Kind::Fraction ||
      this->kind() == Kind::Symbol)
    return true;

  for (unsigned int i = 0; i < this->size(); i++) {
    if (!this->operand(i).freeOf(other))
      return false;
  }

  return true;
}

bool Expr::freeOfElementsInSet(Expr &S) {
  return this->freeOfElementsInSet(std::forward<Expr>(S));
}

bool Expr::freeOfElementsInSet(Expr &&S) {
  for (unsigned int i = 0; i < S.size(); i++) {
    if (!this->freeOf(S[i])) {
      return false;
    }
  }

  return true;
}
std::string Expr::toString() {
  std::string res = "";

  switch (this->kind()) {
  case Kind::Fail:
    res += "Fail";
    break;

  case Kind::Addition:
    // res += "(";
    for (unsigned int i = 0; i < this->size(); i++) {
      res += this->operand(i).toString();
      if (i != this->size() - 1)
        res += " + ";
    }
    // res += ")";
    break;

  case Kind::Subtraction:
    // res += "(";
    for (unsigned int i = 0; i < this->size(); i++) {
      res += this->operand(i).toString();
      if (i != this->size() - 1)
        res += " - ";
    }
    // res += ")";
    break;

  case Kind::Power:
    // res += "(";
    if ((this->operand(0).kind() == Kind::Integer &&
         this->operand(0).value() < 0) ||
        this->operand(0).size() > 1) {
      res += "(";
    }
    res += this->operand(0).toString();
    if ((this->operand(0).kind() == Kind::Integer &&
         this->operand(0).value() < 0) ||
        this->operand(0).size() > 1) {
      res += ")";
    }
    if (this->operand(1).size() > 1 ||
        this->operand(1).kind() == Kind::Fraction) {
      res += "^(";
    } else {
      res += "^";
    }
    res += this->operand(1).toString();
    if (this->operand(1).size() > 1 ||
        this->operand(1).kind() == Kind::Fraction) {
      res += ")";
    }
    // res += ")";
    break;

  case Kind::Multiplication:
    // res += "(";
    for (unsigned int i = 0; i < this->size(); i++) {
      // if(i == 0 && this->operand(i).kind() == Kind::Integer &&
      // this->operand(i).value() == -1) { 	res += "-"; 	continue;
      // }
      if (this->operand(i).kind() == Kind::Addition ||
          this->operand(i).kind() == Kind::Subtraction ||
          this->operand(i).kind() == Kind::Power ||
          this->operand(i).kind() == Kind::Division ||
          this->operand(i).kind() == Kind::Fraction)
        res += "(";

      res += this->operand(i).toString();

      if (this->operand(i).kind() == Kind::Addition ||
          this->operand(i).kind() == Kind::Subtraction ||
          this->operand(i).kind() == Kind::Power ||
          this->operand(i).kind() == Kind::Division ||
          this->operand(i).kind() == Kind::Fraction)
        res += ")";

      if (i != this->size() - 1)
        res += "*";
    }
    // res += ")";
    break;

  case Kind::Division:
    res += "(";
    res += this->operand(0).toString();
    res += ")";
    res += "/";
    res += "(";
    res += this->operand(1).toString();
    res += ")";
    // res += ")";
    break;
  case Kind::Fraction:
    // res += "(";
    res += this->operand(0).toString();
    res += "/";
    res += this->operand(1).toString();
    // res += ")";
    break;

  case Kind::Factorial:
    res += "!";
    // res += "(";
    res += this->operand(0).toString();
    // res += ")";
    break;

  case Kind::Integer:
    res += this->value().to_string();
    break;
  case Kind::Infinity:
    res += "∞";
    break;
  case Kind::MinusInfinity:
    res += "-∞";
    break;
  case Kind::Symbol:
    res += this->identifier();
    break;

  case Kind::FunctionCall:
    res += this->funName();
    res += "(";
    for (unsigned int i = 0; i < this->size(); i++) {
      res += this->operand(i).toString();
      if (i != this->size() - 1)
        res += ", ";
    }
    res += ")";
    break;
  case Kind::List:
    res += "[";
    for (unsigned int i = 0; i < this->size(); i++) {
      res += this->operand(i).toString();
      if (i != this->size() - 1)
        res += ", ";
    }
    res += "]";
    break;
  case Kind::Set:
    res += "{";
    for (unsigned int i = 0; i < this->size(); i++) {
      res += this->operand(i).toString();
      if (i != this->size() - 1)
        res += ", ";
    }
    res += "}";
    break;
  case Kind::Undefined:
    res += "Undefined";
    break;

  case Kind::Derivative:
    res += "diff(";
    res += this->operand(0).toString();
    res += ", ";
    res += this->operand(1).toString();
    res += ")";
    break;

  case Kind::Integral:
    res += "inte(";
    res += this->operand(0).toString();
    res += ", ";
    res += this->operand(1).toString();
    res += ")";
    break;

  case Kind::Matrix:
    res += "[";
    for (unsigned int i = 0; i < this->size(); i++) {
      res += "[";
      for (unsigned int j = 0; j < this->size(); j++) {
        res += this[i][j].toString();
        if (j < this->operand(i).size() - 1) {
          res += ",";
        }
      }
      res += "]";

      if (i < this->size() - 1) {
        res += ",";
      }
    }
    res += "]";
    break;

  default:
    res += "Not implemented(" + std::to_string(this->kind()) + ")";
    break;
  }

  return res;
}

Expr Expr::operandList() {
  Expr L = Expr(Kind::List);

  for (unsigned int i = 0; i < this->size(); i++) {
    L.insert(this->operand(i));
  }

  return L;
}

const std::string Expr::funName() { return this->_operands[0].identifier(); }

Expr mapUnaryAST(Expr &u, Expr (*f)(Expr)) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Fraction ||
      u.kind() == Kind::Symbol || u.kind() == Kind::Infinity ||
      u.kind() == Kind::MinusInfinity) {
    return f(u);
  }

  if (u.size() == 0) {
    return f(u);
  }

  Expr t = Expr(u.kind());

  if (u.kind() == Kind::FunctionCall) {
    t.insert(Expr(Kind::Symbol, u.funName().c_str()));
  }

  for (unsigned int i = 0; i < u.size(); i++) {
    t.insert(f(u[i]));
  }

  return t;
}

Expr mapUnaryAST(Expr &u, Expr (*f)(Expr &)) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Fraction ||
      u.kind() == Kind::Symbol || u.kind() == Kind::Infinity ||
      u.kind() == Kind::MinusInfinity) {
    return f(u);
  }

  if (u.size() == 0) {
    return f(u);
  }

  Expr t = Expr(u.kind());

  if (u.kind() == Kind::FunctionCall) {
    t.insert(Expr(Kind::Symbol, u.funName().c_str()));
  }

  for (unsigned int i = 0; i < u.size(); i++) {
    t.insert(f(u[i]));
  }

  return t;
}

Expr Expr::symbols() {
  if (this->kind() == Kind::Symbol) {
    return Expr(Kind::Set, {*this});
  }

  Expr syms = Expr(Kind::Set);

  if (this->kind() == Kind::Addition || this->kind() == Kind::Subtraction ||
      this->kind() == Kind::Power || this->kind() == Kind::Division ||
      this->kind() == Kind::Multiplication || this->kind() == Kind::Matrix ||
      this->kind() == Kind::Set || this->kind() == Kind::List) {
    for (unsigned int i = 0; i < this->size(); i++) {
      Expr s = this->operand(i).symbols();

      if (s.size() > 0) {
        for (unsigned int k = 0; k < s.size(); k++) {
          syms.insert(s[k]);
        }
      }
    }
  }

  if (this->kind() == Kind::Derivative || this->kind() == Kind::Integral ||
      this->kind() == Kind::Factorial) {
    Expr s = this->operand(0).symbols();
    if (s.size() > 0) {
      for (unsigned int k = 0; k < s.size(); k++) {
        syms.insert(s[k]);
      }
    }
  }

  if (this->kind() == Kind::Derivative || this->kind() == Kind::Integral ||
      this->kind() == Kind::Factorial) {
    Expr s = this->operand(0).symbols();

    if (s.size() > 0) {
      for (unsigned int k = 0; k < s.size(); k++) {
        syms.insert(s[k]);
      }
    }
  }

  return syms;
}
Expr mapBinaryAST(Expr &u, Expr &v, Expr (*f)(Expr &, Expr &)) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Symbol ||
      u.kind() == Kind::Infinity || u.kind() == Kind::MinusInfinity)
    return f(u, v);

  if (u.size() == 0)
    return f(u, v);

  Expr t = Expr(u.kind());

  if (u.kind() == Kind::FunctionCall) {
    t.insert(Expr(Kind::Symbol, u.funName().c_str()));
  }

  for (unsigned int i = 0; i < u.size(); i++)
    t.insert(f(u[i], v));

  return t;
}

Expr mapBinaryAST(Expr &u, Expr &v, Expr (*f)(Expr, Expr)) {
  if (u.kind() == Kind::Integer || u.kind() == Kind::Symbol ||
      u.kind() == Kind::Infinity || u.kind() == Kind::MinusInfinity)
    return f(std::forward<Expr>(u), std::forward<Expr>(v));

  if (u.size() == 0)
    return f(std::forward<Expr>(u), std::forward<Expr>(v));

  Expr t = Expr(u.kind());

  if (u.kind() == Kind::FunctionCall) {
    t.insert(Expr(Kind::Symbol, u.funName().c_str()));
  }

  for (unsigned int i = 0; i < u.size(); i++)
    t.insert(f(std::forward<Expr>(u[i]), std::forward<Expr>(v)));

  return t;
}

Expr deepReplace(Expr &tree, Expr &subtree, Expr &v) {
  return deepReplace(tree, std::forward<Expr>(subtree), std::forward<Expr>(v));
}

Expr deepReplace(Expr &tree, Expr &&subtree, Expr &&v) {
  if (tree.kind() == subtree.kind()) {
    if (tree.match(subtree)) {
      return Expr(v);
    }
  }

  if (tree.size() > 1) {
    Expr t = Expr(tree.kind());

    for (unsigned int i = 0; i < tree.size(); i++) {
      t.insert(deepReplace(tree[i], subtree, v));
    }

    return t;
  }

  return tree;
}

Expr construct(Kind kind, Expr L) {
  Expr u = Expr(kind);

  for (unsigned int i = 0; i < L.size(); i++) {
    u.insert(L[i]);
  }

  return u;
}

// Expr::Expr(const Expr &&other) {
//   this->_identifier = std::move(other._identifier);
//   this->_operands = std::move(other._operands);
//   this->_kind = std::move(other._kind);
//   this->_value = std::move(other._value);
// }
// Expr::Expr(const Expr &other) {
//   this->_identifier = other._identifier;
//   this->_operands = other._operands;
//   this->_kind = other._kind;
//   this->_value = other._value;
// }
Expr Expr::operator+(Expr &&other) {
  if (this->kind() == Kind::Addition) {
    Expr a = *this;
    a.insert(other);
    return a;
  }

  return Expr(Kind::Addition, {*this, other});
}
Expr Expr::operator+(Expr &other) {

  if (this->kind() == Kind::Addition) {
    Expr a = *this;
    a.insert(other);
    return a;
  }

  return Expr(Kind::Addition, {*this, other});
}
Expr Expr::operator-(Expr &&other) {

  if (this->kind() == Kind::Subtraction) {
    Expr a = *this;
    a.insert(other);
    return a;
  }

  return Expr(Kind::Subtraction, {*this, other});
}
Expr Expr::operator-(Expr &other) {

  if (this->kind() == Kind::Subtraction) {
    Expr a = *this;
    a.insert(other);
    return a;
  }

  return Expr(Kind::Subtraction, {*this, other});
}

Expr Expr::operator*(Expr &&other) {
  if (this->kind() == Kind::Multiplication) {
    Expr a = *this;
    a.insert(std::move(other));
    return a;
  }

  return Expr(Kind::Multiplication, {*this, other});
}

Expr Expr::operator*(Expr &other) {
  if (this->kind() == Kind::Multiplication) {
    Expr a = *this;
    a.insert(other);
    return a;
  }

  return Expr(Kind::Multiplication, {*this, other});
}

Expr Expr::operator/(Expr &&other) {
  if (this->kind() == Kind::Integer && other.kind() == Kind::Integer) {
    return Expr(Kind::Fraction, {*this, other});
  }

  return Expr(Kind::Division, {*this, other});
}

Expr Expr::operator/(Expr &other) {
  if (this->kind() == Kind::Integer && other.kind() == Kind::Integer) {
    return Expr(Kind::Fraction, {*this, other});
  }

  return Expr(Kind::Division, {*this, other});
}

bool Expr::operator==(Expr &&other) { return this->match(other); }
bool Expr::operator==(Expr &other) { return this->match(other); }

// bool Expr::operator==(long long other) {
//  return this->kind() == Kind::Integer && this->value() == other;
//}

// bool Expr::operator==(const char *other) {
//	printf("MAYBE\n");
//  return this->kind() == Kind::Symbol &&
//         this->identifier() == std::string(other);
//}

/*
bool Expr::operator==(const std::string other) {
return this->kind() == Kind::Symbol && this->identifier() == other;
}
*/

bool Expr::operator!=(Expr &&other) { return !this->match(other); }
bool Expr::operator!=(Expr &other) { return !this->match(other); }

/*
bool Expr::operator!=(const std::string other) {
 return this->kind() != Kind::Symbol || this->identifier() != other;
}
*/
bool Expr::operator>(const Int other) {
  if (this->kind() == Kind::MinusInfinity) {
    return false;
  }
  if (this->kind() == Kind::Infinity) {
    return true;
  }
  return this->kind() == Kind::Integer && this->value() > other;
}

bool Expr::operator>=(const Int other) {
  if (this->kind() == Kind::MinusInfinity) {
    return false;
  }
  if (this->kind() == Kind::Infinity) {
    return true;
  }
  return this->kind() == Kind::Integer && this->value() >= other;
}
bool Expr::operator<(const Int other) {
  if (this->kind() == Kind::MinusInfinity) {
    return true;
  }

  if (this->kind() == Kind::Infinity) {
    return false;
  }
  return this->kind() == Kind::Integer && this->value() < other;
}
bool Expr::operator<=(const Int other) {
  if (this->kind() == Kind::MinusInfinity) {
    return true;
  }
  if (this->kind() == Kind::Infinity) {
    return false;
  }
  return this->kind() == Kind::Integer && this->value() <= other;
}
Expr Expr::operator+() { return *this; }

Expr Expr::operator-() {
  if (this->kind() == Kind::Infinity)
    return minInf();
  if (this->kind() == Kind::MinusInfinity)
    return inf();

  return -1 * (*this);
}

Expr &Expr::operator[](size_t idx) { return this->operand(idx); }
Expr &Expr::operator[](Int idx) { return this->operand(idx); }

Expr operator*(Int i, Expr &&other) { return Expr(Kind::Integer, i) * other; }

Expr operator*(Int i, Expr &other) { return Expr(Kind::Integer, i) * other; }

Expr operator+(Int i, Expr &&other) { return Expr(Kind::Integer, i) + other; }

Expr operator+(Int i, Expr &other) { return Expr(Kind::Integer, i) + other; }

Expr operator-(Int i, Expr &&other) { return Expr(Kind::Integer, i) - other; }

Expr operator-(Int i, Expr &other) { return Expr(Kind::Integer, i) - other; }

Expr operator/(Int i, Expr &&other) { return Expr(Kind::Integer, i) / other; }

Expr operator/(Int i, Expr &other) { return Expr(Kind::Integer, i) / other; }

Expr operator*(int i, Expr &&other) { return Expr(Kind::Integer, i) * other; }

Expr operator*(int i, Expr &other) { return Expr(Kind::Integer, i) * other; }

Expr operator+(int i, Expr &&other) { return Expr(Kind::Integer, i) + other; }

Expr operator+(int i, Expr &other) { return Expr(Kind::Integer, i) + other; }

Expr operator-(int i, Expr &&other) { return Expr(Kind::Integer, i) - other; }

Expr operator-(int i, Expr &other) { return Expr(Kind::Integer, i) - other; }

Expr operator/(int i, Expr &&other) { return Expr(Kind::Integer, i) / other; }

Expr operator/(int i, Expr &other) { return Expr(Kind::Integer, i) / other; }

Expr operator*(long i, Expr &&other) { return Expr(Kind::Integer, i) * other; }

Expr operator*(long i, Expr &other) { return Expr(Kind::Integer, i) * other; }

Expr operator+(long i, Expr &&other) { return Expr(Kind::Integer, i) + other; }

Expr operator+(long i, Expr &other) { return Expr(Kind::Integer, i) + other; }

Expr operator-(long i, Expr &&other) { return Expr(Kind::Integer, i) - other; }

Expr operator-(long i, Expr &other) { return Expr(Kind::Integer, i) - other; }

Expr operator/(long i, Expr &&other) { return Expr(Kind::Integer, i) / other; }

Expr operator/(long i, Expr &other) { return Expr(Kind::Integer, i) / other; }

Expr operator*(long long i, Expr &&other) {
  return Expr(Kind::Integer, i) * other;
}

Expr operator*(long long i, Expr &other) {
  return Expr(Kind::Integer, i) * other;
}

Expr operator+(long long i, Expr &&other) {
  return Expr(Kind::Integer, i) + other;
}

Expr operator+(long long i, Expr &other) {
  return Expr(Kind::Integer, i) + other;
}

Expr operator-(long long i, Expr &&other) {
  return Expr(Kind::Integer, i) - other;
}

Expr operator-(long long i, Expr &other) {
  return Expr(Kind::Integer, i) - other;
}

Expr operator/(long long i, Expr &&other) {
  return Expr(Kind::Integer, i) / other;
}

Expr operator/(long long i, Expr &other) {
  return Expr(Kind::Integer, i) / other;
}

Expr undefined() { return Expr(Kind::Undefined); }
Expr fail() { return Expr(Kind::Fail); }

Expr inf() { return Expr(Kind::Infinity); }
Expr minInf() { return Expr(Kind::MinusInfinity); }

// Expr& Expr::operator=(const Expr &other) {
//   this->_operands = other._operands;
//   this->_value = other._value;
//   this->_identifier = other._identifier;
//   this->_kind = other._kind;
//   return *this;
// }

// Expr& Expr::operator=(const Expr &&other) {
//   this->_operands = std::move(other._operands);
//   this->_value = std::move(other._value);
//   this->_identifier = std::move(other._identifier);
//   this->_kind = std::move(other._kind);
//   return *this;
// }
} // namespace ast
