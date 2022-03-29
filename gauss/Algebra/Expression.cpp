#include "Expression.hpp"

#include "Matrix.hpp"
#include "Utils.hpp"
#include "Sorting.hpp"
#include "Reduction.hpp"

#include <iostream>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <initializer_list>
#include <limits>
#include <math.h>
#include <string>
#include <utility>
#include <vector>

namespace alg {

using namespace utils;

Int rows(expr a) {
  assert(is(&a, kind::MAT));

  return a.expr_mat->lines();
}

Int columns(expr a) {
  assert(is(&a, kind::MAT));
  return a.expr_mat->columns();
}


expr::expr(expr &&other) {
  kind_of = other.kind_of;
  expr_info = other.expr_info;

  switch (kind_of) {
	// case kind::ERROR: {
	// 	expr_sym = other.expr_sym;
	// 	other.expr_sym = 0;
	// 	return;
	// }
  case kind::SYM: {
    expr_sym = other.expr_sym;
    other.expr_sym = 0;
    return;
  }
  case kind::INT: {
    expr_int = other.expr_int;
    other.expr_int = 0;
    return;
  }

  case kind::LIST: {
    expr_list = other.expr_list;
    other.expr_list = 0;
    return;
  }

  case kind::SET: {
    expr_set = other.expr_set;
    other.expr_set = 0;
    return;
  }

  case kind::MAT: {
    expr_mat = other.expr_mat;
    other.expr_mat = 0;
    return;
  }

  case kind::FUNC: {
    expr_sym = other.expr_sym;
    other.expr_sym = 0;
    expr_childs = std::move(other.expr_childs);

    return;
  }
  case kind::FACT: {
    expr_childs = std::move(other.expr_childs);
    return;
  }
  case kind::POW: {
    expr_childs = std::move(other.expr_childs);
    return;
  }
  case kind::MUL: {
    expr_childs = std::move(other.expr_childs);
    return;
  }
  case kind::ADD: {
    expr_childs = std::move(other.expr_childs);
    return;
  }
  case kind::DIV: {
    expr_childs = std::move(other.expr_childs);
    return;
  }
  case kind::ROOT: {
    expr_childs = std::move(other.expr_childs);
    return;
  }
  case kind::SUB: {
    expr_childs = std::move(other.expr_childs);
    return;
  }

  case kind::FRAC: {
    expr_childs = std::move(other.expr_childs);
    return;
  }

  case kind::INF:
    return;
  case kind::UNDEF:
    return;
  case kind::FAIL:
    return;
  case kind::CONST:
    return;
  case kind::TERMINAL:
    return;
  case kind::SUMMABLE:
    return;
  case kind::NON_CONSTANT:
    return;
  case kind::MULTIPLICABLE:
    return;
  case kind::ORDERED:
    return;
  }
}

expr::expr(const expr &other) {
  kind_of = other.kind_of;
  expr_info = other.expr_info;

  switch (kind_of) {
	// case kind::ERROR: {
	// 	expr_sym = strdup(other.expr_sym);
	// 	return;
	// }
  case kind::SYM: {
    expr_sym = strdup(other.expr_sym);
    return;
  }
  case kind::INT: {
    Int &a = *other.expr_int;

    expr_int = new Int(a);

    return;
  }

  case kind::LIST: {
    expr_list = new list(*other.expr_list);
    return;
  }

  case kind::SET: {
    expr_set = new set(*other.expr_set);
    return;
  }

  case kind::FUNC: {
    expr_sym = strdup(other.expr_sym);

    expr_childs = other.expr_childs;

    return;
  }

  case kind::MAT: {
    expr_mat = matrix::copy(other.expr_mat);
    return;
  }

  case kind::FACT: {
    expr_childs = (other.expr_childs);
    return;
  }
  case kind::POW: {
    expr_childs = (other.expr_childs);
    return;
  }
  case kind::MUL: {
    expr_childs = (other.expr_childs);
    return;
  }
  case kind::ADD: {
    expr_childs = (other.expr_childs);
    return;
  }
  case kind::DIV: {
    expr_childs = (other.expr_childs);
    return;
  }
  case kind::ROOT: {
    expr_childs = (other.expr_childs);
    return;
  }
  case kind::SUB: {
    expr_childs = (other.expr_childs);
    return;
  }

  case kind::FRAC: {
    expr_childs = (other.expr_childs);
    return;
  }

  case kind::INF:
    return;
  case kind::UNDEF:
    return;
  case kind::FAIL:
    return;
  case kind::CONST:
    return;
  case kind::TERMINAL:
    return;
  case kind::SUMMABLE:
    return;
  case kind::NON_CONSTANT:
    return;
  case kind::MULTIPLICABLE:
    return;
  case kind::ORDERED:
    return;
  }
}

expr &expr::operator=(const expr &other) {
  expr_set_kind(this, other.kind());

  expr_info = other.expr_info;

  switch (kind_of) {
	// case kind::ERROR: {
  //   expr_sym = strdup(other.expr_sym);
	// 	expr_childs.clear();
	// 	return *this;
	// }
  case kind::SYM: {
    expr_sym = strdup(other.expr_sym);
    expr_childs.clear();
    return *this;
  }

  case kind::INT: {
    expr_int = new Int(*other.expr_int);
    expr_childs.clear();
    return *this;
  }

  case kind::LIST: {
    expr_list = new list(*other.expr_list);
    expr_childs.clear();
    return *this;
  }

  case kind::SET: {
    expr_set = new set(*other.expr_set);
    expr_childs.clear();
    return *this;
  }

  case kind::FUNC: {
    expr_sym = strdup(other.expr_sym);
    expr_childs = other.expr_childs;
    return *this;
  }

  case kind::MAT: {
    expr_mat = matrix::copy(other.expr_mat);
    return *this;
  }

  case kind::FACT: {
    expr_childs = (other.expr_childs);
    return *this;
  }

  case kind::POW: {
    expr_childs = (other.expr_childs);
    return *this;
  }

  case kind::MUL: {
    expr_childs = (other.expr_childs);
    return *this;
  }

  case kind::ADD: {
    expr_childs = (other.expr_childs);
    return *this;
  }

  case kind::DIV: {
    expr_childs = (other.expr_childs);
    return *this;
  }

  case kind::ROOT: {
    expr_childs = (other.expr_childs);
    return *this;
  }

  case kind::SUB: {
    expr_childs = (other.expr_childs);
    return *this;
  }

  case kind::FRAC: {
    expr_childs = (other.expr_childs);
    return *this;
  }

  case kind::INF:
    return *this;
  case kind::UNDEF:
    return *this;
  case kind::FAIL:
    return *this;
  case kind::CONST:
    return *this;
  case kind::TERMINAL:
    return *this;
  case kind::SUMMABLE:
    return *this;
  case kind::NON_CONSTANT:
    return *this;
  case kind::MULTIPLICABLE:
    return *this;
  case kind::ORDERED:
    return *this;
  }

  return *this;
}

expr &expr::operator=(expr &&other) {
  expr_set_kind(this, other.kind());

  expr_info = other.expr_info;

  switch (kind_of) {
	// case kind::ERROR: {
  //   expr_sym = other.expr_sym;
  //   other.expr_sym = 0;
  //   return *this;
	// }
  case kind::SYM: {
    expr_sym = other.expr_sym;
    other.expr_sym = 0;
    return *this;
  }

  case kind::INT: {
    expr_int = other.expr_int;
    other.expr_int = 0;
    return *this;
  }

  case kind::LIST: {
    expr_list = other.expr_list;
    other.expr_list = 0;
    return *this;
  }

  case kind::SET: {
    expr_set = other.expr_set;
    other.expr_set = 0;
    return *this;
  }

  case kind::FUNC: {
    expr_sym = other.expr_sym;
    other.expr_sym = 0;
    expr_childs = std::move(other.expr_childs);

    return *this;
  }

  case kind::MAT: {
    expr_mat = other.expr_mat;
    other.expr_mat = 0;
    return *this;
  }

  case kind::FACT: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }

  case kind::POW: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }

  case kind::MUL: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }

  case kind::ADD: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }

  case kind::DIV: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }

  case kind::ROOT: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }

  case kind::SUB: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }

  case kind::FRAC: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }

  case kind::INF:
    return *this;
  case kind::UNDEF:
    return *this;
  case kind::FAIL:
    return *this;
  case kind::CONST:
    return *this;
  case kind::TERMINAL:
    return *this;
  case kind::SUMMABLE:
    return *this;
  case kind::NON_CONSTANT:
    return *this;
  case kind::MULTIPLICABLE:
    return *this;
  case kind::ORDERED:
    return *this;
  }

  return *this;
}

expr::expr(enum kind k) {
  expr_info = info::UNKNOWN;
  kind_of = k;
}

enum kind expr::kind() const { return kind_of; }

expr::expr(list &s) {
  kind_of = kind::LIST;

  expr_info = info::UNKNOWN;

  expr_list = new list(s);
}

expr::expr(list &&s) {
  kind_of = kind::LIST;

  expr_info = info::UNKNOWN;

  expr_list = new list(s);
}

expr::expr(set &s) {
  kind_of = kind::SET;

  expr_info = info::UNKNOWN;

  expr_set = new set(s);
}

expr::expr(set &&s) {
  kind_of = kind::SET;

  expr_info = info::UNKNOWN;

  expr_set = new set(s);
}

expr::expr(enum kind k, std::initializer_list<expr> &&a) {
  expr_info = info::UNKNOWN;

  kind_of = k;

  if (k == kind::LIST) {
    this->expr_list = new list(a);
    return;
  }

  if (k == kind::SET) {
    this->expr_set = new set(a);
    return;
  }

  expr_childs = std::move(a);
}

size_t expr::size() { return size_of(this); }

expr::expr(Int v) {
  kind_of = kind::INT;

  expr_info = info::EXPANDED | info::REDUCED | info::SORTED;

  this->expr_int = new Int(v);
}

expr::expr(int v) {
	kind_of = kind::INT;
  expr_info = info::EXPANDED | info::REDUCED | info::SORTED;
  this->expr_int = new Int(v);
}

expr::expr(long int v) {
	kind_of = kind::INT;
  expr_info = info::EXPANDED | info::REDUCED | info::SORTED;
  this->expr_int = new Int(v);
}

expr::expr(long long v) {
  kind_of = kind::INT;
  expr_info = info::EXPANDED | info::REDUCED | info::SORTED;
  this->expr_int = new Int(v);
}

expr::expr(std::string v) {
  kind_of = kind::SYM;
  expr_info = info::EXPANDED | info::REDUCED | info::SORTED;
  this->expr_sym = strdup(v.c_str());
}

expr::expr() {
  kind_of = kind::UNDEF;
  expr_info = info::UNKNOWN;
}

// expr error(const char* msg) {
// 	expr r = expr(kind::ERROR);
// 	r.expr_sym = strdup(msg);
// 	return r;
// }

expr fromDouble(double v, Int max_den) {
	int sign = 1;

	if(v < 0) {
		sign = -1;
		v = -v;
	}

	double integral, fractional;

  fractional = std::modf(v, &integral);

  if(fractional > std::numeric_limits<double>::epsilon()) {
		Int n, d;

		alg::decimalToFraction(fractional, max_den, n, d);

		alg::expr r = sign * (Int(integral) + alg::fraction(n, d));

		reduce(&r);

		return r;
	}

	return Int(integral);
}

expr mat(unsigned int l, unsigned int c) {
  expr m = expr(kind::MAT);

  m.expr_mat = new matrix(l, c);

  return m;
}

expr mat(unsigned int l, unsigned int c, std::initializer_list<double> data) {
  expr m = expr(kind::MAT);

  m.expr_mat = new matrix(l, c, data);

  return m;
}

expr identity_matrix(unsigned int l, unsigned int c) {
  expr m = expr(kind::MAT);

  m.expr_mat = new matrix(l, c);

  for (unsigned i = 0; i < std::min(l, c); i++) {
    m.expr_mat->set(i, i, 1);
  }

  return m;
}

expr mat_get(expr &a, unsigned i, unsigned j) {
  assert(is(&a, kind::MAT));

  return fromDouble(a.expr_mat->get(i, j));
}

void mat_set(expr &a, unsigned i, unsigned j, expr v) {
  assert(is(&v, kind::INT | kind::FRAC));
  if (is(&v, kind::INT))
    a.expr_mat->set(i, j, v.expr_int->doubleValue());

  a.expr_mat->set(i, j,
                  numerator(v).expr_int->doubleValue() /
                      denominator(v).expr_int->doubleValue());
}

expr::~expr() {
  switch (kind_of) {
	// case kind::ERROR: {
	// 	if(expr_sym) {
	// 		free(expr_sym);
	// 	}
	// 	break;
	// }
  case kind::INT: {
    if (expr_int)
      delete expr_int;
    break;
  }

  case kind::SYM: {
    if (expr_sym) {
      free(expr_sym);
    }
    break;
  }

  case kind::LIST: {
    if (expr_list)
      delete expr_list;
    break;
  }

  case kind::SET: {
    if (expr_set)
      delete expr_set;
    break;
  }

  case kind::FUNC: {
    if (expr_sym) {
      delete[] expr_sym;
    }
    break;
  }

  case kind::MAT: {
    if (expr_mat) {
      delete expr_mat;
    }
    break;
  }

  default:
    break;
  }
}

expr &expr::operator[](size_t idx) {
  if (is(this, kind::LIST)) {
    return expr_list->members[idx];
  }

  if (is(this, kind::SET)) {
    return expr_set->members[idx];
  }

  return expr_childs[idx];
}

expr &expr::operator[](Int idx) {
  long long r = idx.longValue();

  if (is(this, kind::LIST)) {
    return expr_list->members[r];
  }

  if (is(this, kind::SET)) {
    return expr_set->members[r];
  }

  return expr_childs[r];
}

expr create(kind kind) { return expr(kind); }

expr func_call(const char *id, std::initializer_list<expr> &&l) {
  expr f = create(kind::FUNC);

  f.expr_sym = strdup(id);

  f.expr_childs = std::move(l);

  return f;
}

expr create(kind kind, std::initializer_list<expr> &&l) {
  // TODO return expr(kind, l)
  expr u(kind);
  u.expr_childs = l;
  return u;
}


inline bool is_sorted(expr *a, enum kind k) {
  return a->expr_info & info::SORTED && a->sort_kind == k;
}

inline bool is_expanded(expr *a) { return a->expr_info & info::EXPANDED; }

expr symbol(const char *id) {
  expr a = create(kind::SYM);

  a.expr_sym = strdup(id);

  return a;
}

expr integer(Int value) {
  expr a = create(kind::INT);

  a.expr_int = new Int(value);

  return a;
}

expr fraction(Int num, Int den) {
  return create(kind::FRAC, {integer(num), integer(den)});
}

void expr::insert(const expr &b, size_t idx) {
  assert(!is(this, kind::SET));

  this->expr_info = info::UNKNOWN;

  if (is(this, kind::LIST)) {
    return this->expr_list->insert(b, idx);
  }

  this->expr_childs.insert(this->expr_childs.begin() + idx, b);
}

void expr::insert(expr &&b, size_t idx) {
  assert(!is(this, kind::SET));

  this->expr_info = info::UNKNOWN;

  if (is(this, kind::LIST)) {
    return this->expr_list->insert(b, idx);
  }

  this->expr_childs.insert(this->expr_childs.begin() + idx, std::move(b));
}

void expr::insert(const expr &b) {
  this->expr_info = info::UNKNOWN;

  if (is(this, kind::LIST)) {
    this->expr_list->insert(b);
  }

  if (is(this, kind::SET)) {
    this->expr_set->insert(b);
  }

  this->expr_childs.push_back(b);
}

void expr::insert(expr &&b) {
  this->expr_info = info::UNKNOWN;

  if (is(this, kind::LIST)) {
    this->expr_list->insert(b);
  }

  if (is(this, kind::SET)) {
    this->expr_set->insert(b);
  }

  this->expr_childs.push_back((b));
}

void expr::remove(list &l) {
  assert(is(this, kind::LIST));
  expr_list->remove(l);
}

void expr::remove(list &&l) {
  assert(is(this, kind::LIST));
  expr_list->remove((l));
}

void expr::remove(size_t idx) {
  if (is(this, kind::LIST)) {
    this->expr_list->members.erase(this->expr_list->members.begin() + idx);
    return;
  }

  if (is(this, kind::SET)) {
    this->expr_set->members.erase(this->expr_set->members.begin() + idx);
    return;
  }

  this->expr_childs.erase(this->expr_childs.begin() + idx);
}

void expr::remove() {
  if (is(this, kind::LIST)) {
    this->expr_list->members.pop_back();
    return;
  }

  if (is(this, kind::SET)) {
    this->expr_set->members.pop_back();
    return;
  }

  this->expr_childs.pop_back();
}


// const char* error_message(expr e) {
// 	assert(kind_of(&e) == kind::ERROR);

// 	return e.expr_sym;
// }

// const char* error_message(expr* e) {
// 	assert(kind_of(e) == kind::ERROR);

// 	return e->expr_sym;
// }

std::string kind_of_id(expr *a) {
  switch (kind_of(a)) {

  case kind::INT: {
    return "integer";
  }
  case kind::SYM: {
    return "symbol";
  }
  case kind::FUNC: {
    return "funcall";
  }
  case kind::FACT: {
    return "fact";
  }
  case kind::POW: {
    return "pow";
  }
  case kind::MUL: {
    return "mul";
  }
  case kind::ADD: {
    return "add";
  }
  case kind::SUB: {
    return "div";
  }
  case kind::ROOT: {
    return "sqrt";
  }
  case kind::INF: {
    return "infinity";
  }

  case kind::UNDEF: {
    return "undefined";
  }

  case kind::FAIL: {
    return "fail";
  }

  case kind::FRAC: {
    return "fraction";
  }

  case kind::DIV: {
    return "div";
  }

  case kind::LIST: {
    return "list";
  }
  case kind::SET: {
    return "set";
  }
	// case kind::ERROR: {
  //   return "error";
  // }

  default:
    return "kind id not implemented";
  }
}



std::string to_string(expr *tree) {
  if (!tree) {
    return "null";
	}

	// if (is(tree, kind::ERROR)) {
	// 	printf("error(%s)\n", error_message(tree));
	// }

  if (is(tree, kind::MAT)) {
    return matrixToString(tree->expr_mat);
  }

  if (is(tree, kind::INT)) {
    return tree->expr_int->to_string();
  }

  if (is(tree, kind::SYM)) {
    return std::string(tree->expr_sym);
  }

  if (is(tree, kind::UNDEF)) {
    return "undefined";
  }

  if (is(tree, kind::FAIL)) {
    return "fail";
  }

  if (is(tree, kind::INF)) {
    return "inf";
  }

  if (is(tree, kind::FRAC)) {
    return to_string(operand(tree, 0)) + "/" + to_string(operand(tree, 1));
  }

  if (is(tree, kind::ROOT)) {
    return "sqrt(" + to_string(operand(tree, 0)) + "," +  to_string(operand(tree, 1)) + ")";
  }

  if (is(tree, kind::FUNC)) {
    std::string r = std::string(get_func_id(tree)) + "(";

    for (size_t i = 0; i < size_of(tree); i++) {
      r += to_string(operand(tree, i));

      if (i < size_of(tree) - 1) {
        r += ", ";
      }
    }

    r += ")";

    return r;
  }

  if (is(tree, kind::POW)) {
    std::string r = "";

    if (operand(tree, 0) &&
        is(operand(tree, 0), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += "(";
    }

    r += to_string(operand(tree, 0));

    if (operand(tree, 0) &&
        is(operand(tree, 0), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += ")";
    }

    r += "^";

    if (operand(tree, 1) &&
        is(operand(tree, 1), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += "(";
    }

    r += to_string(operand(tree, 1));

    if (operand(tree, 1) &&
        is(operand(tree, 1), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += ")";
    }

    return r;
  }

  if (is(tree, kind::DIV)) {
    std::string r = "";

    if (is(operand(tree, 0), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += "(";
    }

    r += to_string(operand(tree, 0));

    if (is(operand(tree, 0), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += ")";
    }

    r += " ÷ ";

    if (is(operand(tree, 1), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += "(";
    }

    r += to_string(operand(tree, 1));

    if (is(operand(tree, 1), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += ")";
    }

    return r;
  }

  if (is(tree, kind::ADD)) {
    std::string r = "";

    for (size_t i = 0; i < size_of(tree); i++) {
      if (operand(tree, i) && is(operand(tree, i), kind::SUB | kind::ADD)) {

        r += "(";
      }

      r += to_string(operand(tree, i));

      if (operand(tree, i) && is(operand(tree, i), kind::SUB | kind::ADD)) {
        r += ")";
      }

      if (i < size_of(tree) - 1) {
        r += " + ";
      }
    }

    return r;
  }

  if (is(tree, kind::SUB)) {
    std::string r = "";

    for (size_t i = 0; i < size_of(tree); i++) {
      if (operand(tree, i) && is(operand(tree, i), kind::SUB | kind::ADD)) {

        r += "(";
      }

      r += to_string(operand(tree, i));

      if (operand(tree, i) && is(operand(tree, i), kind::SUB | kind::ADD)) {
        r += ")";
      }

      if (i < size_of(tree) - 1) {
        r += " - ";
      }
    }

    return r;
  }

  if (is(tree, kind::MUL)) {
    std::string r = "";

    for (size_t i = 0; i < size_of(tree); i++) {

      if (operand(tree, i) == nullptr) {
        continue;
      }

      if (is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL | kind::FRAC)) {
        r += "(";
      }

      r += to_string(operand(tree, i));

      if (is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL | kind::FRAC)) {
        r += ")";
      }

      if (i < size_of(tree) - 1) {
        r += "⋅";
      }
    }

    return r;
  }

  if (is(tree, kind::FACT)) {
    return to_string(operand(tree, 0)) + "!";
  }

  if (is(tree, kind::LIST)) {
    return to_string(tree->expr_list);
	}

	if (is(tree, kind::SET)) {
    return to_string(tree->expr_set);
	}
  return "to string not implemented for kind " + kind_of_id(tree);
}



std::string to_latex(expr *tree, bool fractions, unsigned long max_den) {
  if (!tree)
    return "null";

  // if (is(tree, kind::ERROR)) {
	// 	printf("error(%s)\n", error_message(tree));
	// }

  if (is(tree, kind::MAT)) {
    return matrixToLatex(tree->expr_mat, fractions, max_den);
  }

  if (is(tree, kind::INT)) {
    return tree->expr_int->to_string();
  }

  if (is(tree, kind::SYM)) {
    return std::string(tree->expr_sym);
  }

  if (is(tree, kind::UNDEF)) {
    return "undefined";
  }

  if (is(tree, kind::FAIL)) {
    return "fail";
  }

  if (is(tree, kind::INF)) {
    return "\\infty";
  }

  if (is(tree, kind::FRAC)) {
    return "\\frac{" + to_latex(operand(tree, 0), fractions, max_den) + "}{" + to_latex(operand(tree, 1), fractions, max_den) + "}";
  }

  if (is(tree, kind::ROOT)) {
    return "\\sqrt[" +  to_latex(operand(tree, 1), fractions, max_den) + "]{" + to_latex(operand(tree, 0), fractions, max_den) + "}";
  }

  if (is(tree, kind::FUNC)) {
		if(strcmp(get_func_id(tree), "ln")) {
			return "\\ln " + to_latex(operand(tree, 0), fractions, max_den);
		}

		if(strcmp(get_func_id(tree), "root")) {
			return "\\log_{" + to_latex(operand(tree, 1), fractions, max_den) + "}" + to_latex(operand(tree, 0), fractions, max_den);
		}

		if(strcmp(get_func_id(tree), "log")) {
			return "\\log_{" + to_latex(operand(tree, 1), fractions, max_den) + "}" + to_latex(operand(tree, 0), fractions, max_den);
		}

		if(strcmp(get_func_id(tree), "derivative") == 0) {
			return "\\dv{" + to_latex(operand(tree, 1)) + "}" + to_latex(operand(tree, 0), fractions, max_den);
		}

    std::string r = std::string(get_func_id(tree)) + "(";

    for (size_t i = 0; i < size_of(tree); i++) {
      r += to_latex(operand(tree, i), fractions, max_den);

      if (i < size_of(tree) - 1) {
        r += ", ";
      }
    }

    r += ")";

    return r;
  }

  if (is(tree, kind::POW)) {
    std::string r = "";

    if (operand(tree, 0) && is(operand(tree, 0), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += "(";
    }

    r += to_latex(operand(tree, 0), fractions, max_den);

    if (operand(tree, 0) && is(operand(tree, 0), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += ")";
    }

    r += "^";

    if (operand(tree, 1) && is(operand(tree, 1), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += "{";
    }

    r += to_latex(operand(tree, 1), fractions, max_den);

    if (operand(tree, 1) && is(operand(tree, 1), kind::SUB | kind::ADD | kind::MUL | kind::DIV)) {
      r += "}";
    }

    return r;
  }

  if (is(tree, kind::DIV)) {
    std::string r = "\\frac{";

    r += to_latex(operand(tree, 0), fractions, max_den);

    r += "}{";

    r += to_latex(operand(tree, 1), fractions, max_den);

		r += "}";

    return r;
  }

  if (is(tree, kind::ADD)) {
    std::string r = "";

    for (size_t i = 0; i < size_of(tree); i++) {
      if (operand(tree, i) && is(operand(tree, i), kind::SUB | kind::ADD)) {

        r += "(";
      }

      r += to_latex(operand(tree, i), fractions, max_den);

      if (operand(tree, i) && is(operand(tree, i), kind::SUB | kind::ADD)) {
        r += ")";
      }

      if (i < size_of(tree) - 1) {
        r += " + ";
      }
    }

    return r;
  }

  if (is(tree, kind::SUB)) {
    std::string r = "";

    for (size_t i = 0; i < size_of(tree); i++) {
      if (operand(tree, i) && is(operand(tree, i), kind::SUB | kind::ADD)) {

        r += "(";
      }

      r += to_latex(operand(tree, i), fractions, max_den);

      if (operand(tree, i) && is(operand(tree, i), kind::SUB | kind::ADD)) {
        r += ")";
      }

      if (i < size_of(tree) - 1) {
        r += " - ";
      }
    }

    return r;
  }

  if (is(tree, kind::MUL)) {
		int start = 0;

		std::string r = "";
		if(is(operand(tree, 0), kind::INT) && *operand(tree, 0)->expr_int == -1) {
			start = 1;
			r += "-";
		}

    for (size_t i = start; i < size_of(tree); i++) {

      if (operand(tree, i) == nullptr) {
        continue;
      }

      if (is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL)) {
        r += "(";
      }

      r += to_latex(operand(tree, i), fractions, max_den);

      if (is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL)) {
        r += ")";
      }

      if (i < size_of(tree) - 1) {
        if (kind_of(operand(tree, i)) == kind_of(operand(tree, i + 1)) ||
            (is(operand(tree, 0), kind::FUNC | kind::SYM) &&
             is(operand(tree, 1), kind::FUNC | kind::SYM))) {
          r += "\\cdot";
        }
      }
    }

    return r;
  }

  if (is(tree, kind::FACT)) {
    return to_string(operand(tree, 0)) + " \\endspace !";
  }

  if (is(tree, kind::LIST)) {
    return to_latex(tree->expr_list, fractions, max_den);

	}
	if (is(tree, kind::SET)) {
    return to_latex(tree->expr_set, fractions, max_den);
	}

  return "to string not implemented for kind " + kind_of_id(tree);
}


std::string to_latex(expr tree, bool fractions, unsigned long max_den) {
	return to_latex(&tree, fractions, max_den);
}


// void expr_print(expr *a, int tabs) {
//   printf("%*c<expr ", tabs, ' ');
//   printf("address=\"%p\" ", a);
//   printf("kind=\"%s\"", kind_of_id(a).c_str());

//   if (kind_of(a) == kind::INT) {
//     printf(" value=\"%s\"", get_val(a).to_string().c_str());
//   }

//   if (kind_of(a) == kind::SYM || kind_of(a) == kind::FUNC) {
//     printf(" id=\"%s\"", get_id(a));
//   }

//   if (is_reduced(a)) {
//     printf(" reduced=\"true\"");
//   } else {
//     printf(" reduced=\"false\"");
//   }
//   if (is_expanded(a)) {
//     printf(" expanded=\"true\"");
//   } else {
//     printf(" expanded=\"false\"");
//   }

//   if (size_of(a)) {
//     printf(" size=\"%li\" >\n", size_of(a));

//     for (size_t i = 0; i < size_of(a); i++) {
//       expr_print(operand(a, i), tabs + 3);
//     }

//     printf("%*c</expr>\n", tabs, ' ');
//   } else {
//     printf(">\n");
//   }
// }


expr expand_mul(expr *a, size_t i, expr *b, size_t j) {
  expr *r = operand(a, i);
  expr *s = operand(b, j);

  if (is(r, kind::ADD) && is(s, kind::ADD)) {
    expr u = create(kind::ADD);

    for (size_t k = 0; k < size_of(r); k++) {
      for (size_t t = 0; t < size_of(s); t++) {
        u.insert(create(kind::MUL,
                        {expr(r->expr_childs[k]), expr(s->expr_childs[t])}));
      }
    }

    return u;
  }

  if (is(r, kind::ADD)) {
    expr u = create(kind::ADD);

    for (size_t k = 0; k < size_of(r); k++) {
      u.insert(create(kind::MUL,
                      {expr(r->expr_childs[k]), expr(b->expr_childs[j])}));
    }

    return u;
  }

  if (is(s, kind::ADD)) {
    expr u = create(kind::ADD);

    for (size_t k = 0; k < size_of(s); k++) {
      u.insert(create(kind::MUL,
                      {expr(a->expr_childs[i]), expr(s->expr_childs[k])}));
    }

    return u;
  }

  return create(kind::MUL, {expr(a->expr_childs[i]), expr(b->expr_childs[j])});
}


expr expand_mul(expr *r, expr *s) {
	// printf("expanding %s * %s = ", to_string(r).c_str(), to_string(s).c_str());

  if (is(r, kind::ADD) && is(s, kind::ADD)) {
    expr u = create(kind::ADD);

    for (size_t k = 0; k < size_of(r); k++) {
      for (size_t t = 0; t < size_of(s); t++) {
        u.insert(create(kind::MUL,
                        {expr(r->expr_childs[k]), expr(s->expr_childs[t])}));
      }
    }

		// printf("%s\n", to_string(u).c_str());

		return u;
  }

  if (is(r, kind::ADD)) {
    expr u = create(kind::ADD);

    for (size_t k = 0; k < size_of(r); k++) {
      u.insert(create(kind::MUL,
                      {expr(r->expr_childs[k]), *s}));
    }

		// printf("%s\n", to_string(u).c_str());
    return u;
  }

  if (is(s, kind::ADD)) {
    expr u = create(kind::ADD);

    for (size_t k = 0; k < size_of(s); k++) {
      u.insert(create(kind::MUL,
                      {*r, expr(s->expr_childs[k])}));
    }

    return u;
  }

	// printf("%s\n", to_string(create(kind::MUL, {*r, *s})).c_str());

  return create(kind::MUL, {*r, *s});
}


bool expand_pow(expr u, Int n, expr *a) {
	if (n == 1) {
    *a = u;
    return true;
  }

  if (is(&u, kind::TERMINAL)) {
    *a = pow(u, n);
    return true;
  }

  if (n == 0) {
    *a = integer(1);
    return true;
  }

  // if (is(&u, kind::ADD)) {
  //   Int c = fact(n);

  //   expr o = u;

  //   expr f = o[0];

  //   o.remove(0);

  //   if (size_of(&o) == 0) {
  //     expr_set_to_int(&o, 0);
  //   }

  //   if (size_of(&o) == 1) {
  //     expr_raise_to_first_op(&o);
  //   }

  //   expr s = create(kind::ADD);

  //   expr t;

  //   for (Int k = 0; k <= n; k++) {
  //     expr z = (c / (fact(k) * fact(n - k))) * pow(f, n - k);

  //     t = expand_pow(o, k, &t) ? z * t : z;

  //     s.insert(t);
  //   }
  //   *a = s;
  //   return true;
  // }

	if(is(&u, kind::ADD)) {
		// TODO: naive aproach with logarithmic exponentiation,
		// optmize this.

		expr g = 1;

		expr x = u;

		while (n) {
			if (n % 2 == 1) {
				g = expand_mul(&g, &x);
			}

			n = n / 2;

			x = expand_mul(&x, &x);
			x = reduce(x);
		}

		g = reduce(g);

		*a = g;

		return true;
	}

  return false;
}

void expand(expr *a) {
  if (is_expanded(a)) {
    return;
  }

  if (is(a, kind::TERMINAL)) {
    set_to_expanded(a);
    return;
  }

  if (is(a, kind::SUB | kind::DIV | kind::FACT | kind::FRAC | kind::ROOT)) {
    reduce(a);
  }

  if (is(a, kind::POW)) {
    expand(operand(a, 0));
    expand(operand(a, 1));

    if (is(operand(a, 1), kind::INT)) {
      expand_pow(a->expr_childs[0], get_val(operand(a, 1)), a);
    }

    set_to_unreduced(a);

    reduce(a);
  }

  if (is(a, kind::MUL)) {

    while (size_of(a) > 1) {
      expand(operand(a, 0));
      expand(operand(a, 1));

      expr t = expand_mul(a, 0, a, 1);

      a->insert(t, 0);

      a->remove(1);
      a->remove(1);
    }

    expr_raise_to_first_op(a);

    set_to_unreduced(a);

    reduce(a);
  }

  if (is(a, kind::ADD)) {

    for (size_t i = 0; i < size_of(a); i++) {
      expand(operand(a, i));
    }

    set_to_unreduced(a);

    reduce(a);
  }

  set_to_expanded(a);
}

expr &expr::operator+=(const expr &a) {
  if (is(this, kind::ADD)) {
    this->insert(a);
  } else {
    *this = create(kind::ADD, {*this, a});
  }

  return *this;
}

expr &expr::operator+=(expr &&a) {
  if (is(this, kind::ADD)) {
    this->insert(a);
  } else {
    *this = create(kind::ADD, {*this, a});
  }

  return *this;
}

expr &expr::operator-=(const expr &a) {
  if (is(this, kind::ADD)) {
    this->insert(create(kind::MUL, {integer(-1), a}));
  } else {
    *this = create(kind::ADD, {*this, create(kind::MUL, {integer(-1), a})});
  }

  return *this;
}

expr &expr::operator-=(expr &&a) {
  if (is(this, kind::ADD)) {
    this->insert(create(kind::MUL, {integer(-1), a}));
  } else {
    *this = create(kind::ADD, {*this, create(kind::MUL, {integer(-1), a})});
  }

  return *this;
}

expr expr::operator+(const expr &a) {
  if (is(this, kind::ADD)) {
    expr t = *this;

    t.insert(a);

    return t;
  }

  return create(kind::ADD, {*this, a});
}

expr expr::operator+(expr &&a) {
  if (is(this, kind::ADD)) {
    expr t = *this;
    t.insert(a);

    return t;
  }

  return create(kind::ADD, {*this, a});
}

expr expr::operator-(const expr &a) {
  if (is(this, kind::SUB)) {
    expr t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::SUB, {*this, a});
}

expr expr::operator-(expr &&a) {
  if (is(this, kind::SUB)) {
    expr t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::SUB, {*this, a});
}

expr expr::operator*(const expr &a) {
  if (is(this, kind::MUL)) {
    expr t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::MUL, {*this, a});
}

expr expr::operator*(expr &&a) {
  if (is(this, kind::MUL)) {
    expr t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::MUL, {*this, a});
}

expr expr::operator/(const expr &a) {
  if (is(this, kind::INT) && is(&a, kind::INT)) {
    return create(kind::FRAC, {*this, a});
  }

  if (is(this, kind::DIV)) {
    expr t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::DIV, {*this, a});
}

expr expr::operator/(expr &&a) {
  if (is(this, kind::INT) && is(&a, kind::INT)) {
    return create(kind::FRAC, {*this, a});
  }

  if (is(this, kind::DIV)) {
    expr t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::DIV, {*this, a});
}

bool expr::match(expr *other) {
  expr a = *other;
  expr b = *this;

  sort(&a, kind::UNDEF);
  sort(&b, kind::UNDEF);

  return compare(&a, &b, kind::UNDEF) == 0;
}

bool expr::operator==(const expr &other) {
  if (this->kind() != other.kind()) {
    return false;
  }

  if (is(this, kind::ADD | kind::SUB | kind::MUL) &&
      size_of(this) != size_of(&other)) {
    return false;
  }

  if (is(this, kind::MAT) && is(&other, kind::MAT)) {
    if (rows(this) != rows(other) || columns(this) != columns(other)) {
      return false;
    }

    for (unsigned i = 0; i < this->expr_mat->lines(); i++) {
      for (unsigned j = 0; j < this->expr_mat->columns(); j++) {
        if (this->expr_mat->get(i, j) != other.expr_mat->get(i, j)) {
          return false;
        }
      }
    }

    return true;
  }

  expr a = other;
  expr b = *this;

  sort(&a, kind::UNDEF);
  sort(&b, kind::UNDEF);

  return compare(&a, &b, kind::UNDEF) == 0;
}

bool expr::operator==(expr &&a) {
  if (this->kind() != a.kind()) {
    return false;
  }

  if (is(this, kind::ADD | kind::SUB | kind::MUL) &&
      size_of(this) != size_of(&a)) {
    return false;
  }

  if (is(this, kind::MAT) && is(&a, kind::MAT)) {
    if (rows(this) != rows(a) || columns(this) != columns(a)) {
      return false;
    }

    for (unsigned i = 0; i < this->expr_mat->lines(); i++) {
      for (unsigned j = 0; j < this->expr_mat->columns(); j++) {
        if (this->expr_mat->get(i, j) != a.expr_mat->get(i, j)) {
          return false;
        }
      }
    }

    return true;
  }

  expr b = *this;

  sort(&a, kind::UNDEF);

  sort(&b, kind::UNDEF);

  return compare(&a, &b, kind::UNDEF) == 0;
}

bool expr::operator!=(const expr &other) {
  if (this->kind() != other.kind()) {
    return true;
  }

  if (is(this, kind::ADD | kind::SUB | kind::MUL) &&
      size_of(this) != size_of(&other)) {
    return true;
  }
  if (is(this, kind::MAT) && is(&other, kind::MAT)) {
    if (rows(this) != rows(other) || columns(this) != columns(other)) {
      return true;
    }

    for (unsigned i = 0; i < this->expr_mat->lines(); i++) {
      for (unsigned j = 0; j < this->expr_mat->columns(); j++) {
        if (this->expr_mat->get(i, j) != other.expr_mat->get(i, j)) {
          return true;
        }
      }
    }

    return false;
  }

  expr a = other;
  expr b = *this;

  sort(&a, kind::UNDEF);
  sort(&b, kind::UNDEF);

  return compare(&a, &b, kind::UNDEF) != 0;
}

bool expr::operator!=(expr &&a) {
  if (this->kind() != a.kind()) {
    return true;
  }

  if (is(this, kind::ADD | kind::SUB | kind::MUL) &&
      size_of(this) != size_of(&a)) {
    return true;
  }

  if (is(this, kind::MAT) && is(&a, kind::MAT)) {
    if (rows(this) != rows(a) || columns(this) != columns(a)) {
      return true;
    }

    for (unsigned i = 0; i < this->expr_mat->lines(); i++) {
      for (unsigned j = 0; j < this->expr_mat->columns(); j++) {
        if (this->expr_mat->get(i, j) != a.expr_mat->get(i, j)) {
          return true;
        }
      }
    }

    return false;
  }

  expr b = *this;

  sort(&a, kind::UNDEF);
  sort(&b, kind::UNDEF);

  return compare(&a, &b, kind::UNDEF) != 0;
}

expr expr::operator+() { return *this; }

expr expr::operator-() { return create(kind::MUL, {integer(-1), *this}); }

expr operator*(Int i, expr &&other) { return integer(i) * other; }

expr operator*(Int i, expr &other) { return integer(i) * other; }

expr operator+(Int i, expr &&other) { return integer(i) + other; }

expr operator+(Int i, expr &other) { return integer(i) + other; }

expr operator-(Int i, expr &&other) { return integer(i) - other; }

expr operator-(Int i, expr &other) { return integer(i) - other; }

expr operator/(Int i, expr &&other) { return integer(i) / other; }

expr operator/(Int i, expr &other) { return integer(i) / other; }

expr operator*(int i, expr &&other) { return integer(i) * other; }

expr operator*(int i, expr &other) { return integer(i) * other; }

expr operator+(int i, expr &&other) { return integer(i) + other; }

expr operator+(int i, expr &other) { return integer(i) + other; }

expr operator-(int i, expr &&other) { return integer(i) - other; }

expr operator-(int i, expr &other) { return integer(i) - other; }

expr operator/(int i, expr &&other) { return integer(i) / other; }

expr operator/(int i, expr &other) { return integer(i) / other; }

expr operator*(long i, expr &&other) { return integer(i) * other; }

expr operator*(long i, expr &other) { return integer(i) * other; }

expr operator+(long i, expr &&other) { return integer(i) + other; }

expr operator+(long i, expr &other) { return integer(i) + other; }

expr operator-(long i, expr &&other) { return integer(i) - other; }

expr operator-(long i, expr &other) { return integer(i) - other; }

expr operator/(long i, expr &&other) { return integer(i) / other; }

expr operator/(long i, expr &other) { return integer(i) / other; }

expr operator*(long long i, expr &&other) { return integer(i) * other; }

expr operator*(long long i, expr &other) { return integer(i) * other; }

expr operator+(long long i, expr &&other) { return integer(i) + other; }

expr operator+(long long i, expr &other) { return integer(i) + other; }

expr operator-(long long i, expr &&other) { return integer(i) - other; }

expr operator-(long long i, expr &other) { return integer(i) - other; }

expr operator/(long long i, expr &&other) { return integer(i) / other; }

expr operator/(long long i, expr &other) { return integer(i) / other; }

expr pow(const expr &a, const expr &b) { return create(kind::POW, {a, b}); }

expr pow(expr &&a, expr &&b) { return create(kind::POW, {a, b}); }

expr pow(const expr &a, expr &&b) { return create(kind::POW, {a, b}); }

expr pow(expr &&a, const expr &b) { return create(kind::POW, {a, b}); }

expr sqrt(const expr &a, expr b) { return create(kind::ROOT, {a, b}); }

expr sqrt(expr &&a, expr b) { return create(kind::ROOT, {a, b}); }

expr fact(const expr &a) { return create(kind::FACT, {a}); }

expr fact(expr &&a) { return create(kind::FACT, {a}); }

expr undefined() { return create(kind::UNDEF); }

expr fail() { return create(kind::FAIL); }

expr inf() { return create(kind::INF); }



expr expand(expr &a) {
  expr b = a;

  expand(&b);

  return b;
}

expr expand(expr &&a) {
  expr b = a;

  expand(&b);

  return b;
}

std::string to_string(expr &a) { return to_string(&a); }

std::string to_string(expr &&a) { return to_string(&a); }

expr first(expr &a) {
  if (is(&a, kind::LIST)) {
    return a.expr_list->members[0];
  }

  if (is(&a, kind::SET)) {
    return a.expr_set->members[0];
  }

  assert(!is(&a, kind::TERMINAL));

  return a[0];
}

expr rest(expr &a) {
  if (is(&a, kind::LIST)) {
    list L = *a.expr_list;
    return expr(rest(L));
  }

  if (is(&a, kind::SET)) {
    set L = *a.expr_set;
    return expr(rest(L));
  }

  assert(!is(&a, kind::TERMINAL));

  return a[0];
}

expr append(const expr &a, const expr &b) {
  assert(is(&a, kind::LIST) && is(&b, kind::LIST));

  list L = *a.expr_list;
  list M = *b.expr_list;

  return expr(append(L, M));
}

expr append(const expr &a, expr &&b) {
  assert(is(&a, kind::LIST) && is(&b, kind::LIST));

  list L = *a.expr_list;
  list M = *b.expr_list;

  return expr(append(L, M));
}

expr join(const expr &a, const expr &b) {
  assert(is(&a, kind::LIST) && is(&b, kind::LIST));

  list L = *a.expr_list;
  list M = *b.expr_list;

  return expr(join(L, M));
}

expr join(const expr &a, expr &&b) {
  assert(is(&a, kind::LIST) && is(&b, kind::LIST));

  list L = *a.expr_list;
  list M = *b.expr_list;

  return expr(join(L, M));
}

expr difference(const expr &a, expr &&b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.expr_set;
  set M = *b.expr_set;

  return expr(difference(L, M));
}

expr difference(const expr &a, const expr &b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.expr_set;
  set M = *b.expr_set;

  return expr(difference(L, M));
}

expr unification(const expr &a, expr &&b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.expr_set;
  set M = *b.expr_set;

  return expr(unification(L, M));
}

expr unification(const expr &a, const expr &b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));
  set L = *a.expr_set;
  set M = *b.expr_set;

  return expr(unification(L, M));
}
expr intersection(const expr &a, expr &&b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.expr_set;
  set M = *b.expr_set;

  return expr(intersection(L, M));
}

expr intersection(const expr &a, const expr &b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.expr_set;
  set M = *b.expr_set;

  return expr(intersection(L, M));
}

bool exists(const expr &a, expr &b) {
  assert(is(&a, kind::SET));

  set L = *a.expr_set;

  return set_exists(L, b);
}

bool replace_rec(expr *a, expr *b, expr *c) {
  if (a->match(b)) {
    expr_replace_with(a, c);

    set_to_unreduced(a);
    set_to_unsorted(a);

    return true;
  }

  if (is(a, kind::TERMINAL))
    return false;

  bool replaced = false;

  for (size_t i = 0; i < size_of(a); i++) {

    if (replace_rec(operand(a, i), b, c)) {
      replaced = true;

      set_to_unreduced(a);
      set_to_unsorted(a);
    }
  }

  return replaced;
}

expr replace(expr &a, expr &b, expr &c) {
  expr d = a;

  replace_rec(&d, &b, &c);

  return d;
}

expr replace(expr &a, expr &&b, expr &&c) {
  expr d = a;

  replace_rec(&d, &b, &c);

  return d;
}

expr replace(expr &a, expr &&b, expr &c) {
  expr d = a;

  replace_rec(&d, &b, &c);

  return d;
}

expr replace(expr &a, expr &b, expr &&c) {
  expr d = a;

  replace_rec(&d, &b, &c);

  return d;
}

expr map(expr &u, expr &v, expr (*f)(expr &, expr &)) {
  if (is(&u, kind::TERMINAL)) {
    return f(u, v);
  }

  if (size_of(&u) == 0) {
    return f(u, v);
  }

  expr t = create(kind_of(&u));

  if (is(&u, kind::FUNC)) {
    u.expr_sym = strdup(u.expr_sym);
  }

  for (size_t i = 0; i < size_of(&u); i++)
    t.insert(f(u[i], v));

  return t;
}

expr map(expr &u, expr (*f)(expr &)) {
  if (is(&u, kind::TERMINAL)) {
    return f(u);
  }

  if (size_of(&u) == 0) {
    return f(u);
  }

  expr t = create(kind_of(&u));

  if (is(&u, kind::FUNC)) {
    u.expr_sym = strdup(u.expr_sym);
  }

  for (size_t i = 0; i < size_of(&u); i++) {
    t.insert(f(u[i]));
  }

  return t;
}

// TODO: free_of_rec is currently sorting a and b
// on each recursive call, this can be avoided by
//  sorting on the non recursive methods and just
// calling compare here
bool free_of_rec(expr *a, expr *b) {
  if (a->match(b)) {
    return false;
  }

  if (is(a, kind::TERMINAL)) {
    return true;
  }

  for (size_t i = 0; i < size_of(a); i++) {
    if (!free_of_rec(operand(a, i), b)) {
      return false;
    }
  }

  return true;
}

bool expr::freeOf(expr &a) { return free_of_rec(this, &a); }
bool expr::freeOf(expr &&a) { return free_of_rec(this, &a); }

list::list(std::initializer_list<expr> &&a) { members = a; }

list::list(std::vector<expr> &&a) { members = std::move(a); }

list::list(std::vector<expr> &a) { members = a; }

void list::append(list &&a) {
  for (size_t i = 0; i < a.size(); i++) {
    members.push_back(a[i]);
  }
}

void list::append(list &a) {
  for (size_t i = 0; i < a.size(); i++) {
    members.push_back(a[i]);
  }
}

void list::append(list *a) {
  for (size_t i = 0; i < a->members.size(); i++) {
    members.push_back(a->members[i]);
  }
}

void list::insert(const expr &a, size_t idx) {
  members.insert(members.begin() + idx, a);
}

void list::insert(expr &&a, size_t idx) {
  members.insert(members.begin() + idx, a);
}

void list::insert(const expr &a) { members.push_back(a); }

void list::insert(expr &&a) { members.push_back(a); }

list append(list &a, list &b) {
  list L = a;
  L.append(b);
  return L;
}

list append(list &a, list &&b) {
  list L = a;
  L.append(b);
  return L;
}

list insert(list &a, const expr &b) {
  list L = a;
  L.insert(b);
  return L;
}

list insert(list &a, expr &&b) {
  list L = a;
  L.insert(b);
  return L;
}

list append(list &a, list *b) {
  list L = a;
  L.append(b);
  return L;
}

list remove(list &a, list &b) {
  list L = a;
  L.remove(b);
  return L;
}

list remove(list &a, list &&b) {
  list L = a;
  L.remove(b);
  return L;
}

void list::remove(list &M) {
  std::vector<expr> L;
  size_t p = 0;

  for (size_t i = 0; i < size(); i++) {
    bool inc = false;

    for (size_t j = p; j < M.size(); j++) {
      if (members[i] == M[j]) {
        p = j + 1;
        inc = true;
        break;
      }
    }

    if (inc)
      continue;

    L.push_back(members[i]);
  }

  members = std::move(L);
}

void list::remove(list &&M) {
  std::vector<expr> L;

  size_t p = 0;

  for (size_t i = 0; i < size(); i++) {
    bool inc = false;

    for (size_t j = p; j < M.size(); j++) {
      if (members[i] == M[j]) {
        p = j + 1;
        inc = true;
        break;
      }
    }

    if (inc)
      continue;

    L.push_back(members[i]);
  }

  members = std::move(L);
}

list list::rest(size_t from) {
  std::vector<expr> l;

  for (size_t i = from; i < size(); i++) {
    l.push_back(members[i]);
  }

  return l;
}

list rest(list &a, size_t from) {
  std::vector<expr> l;

  for (size_t i = from; i < a.size(); i++) {
    l.push_back(a[i]);
  }

  return l;
}

expr fist(list &l) { return l[0]; }

void list::join(list &a) {
  for (size_t i = 0; i < a.size(); i++) {
    members.push_back(a[i]);
  }
}

void list::join(list &&a) {
  for (size_t i = 0; i < a.size(); i++) {
    members.push_back(a[i]);
  }
}

list join(list &a, list &b) {
  list L = a;
  L.join(b);
  return L;
}

void list::remove(size_t i) { members.erase(members.begin() + i); }

void list::remove() { members.pop_back(); }

bool list::match(list *a) {
  if (size() != a->size()) {
    return (long)size() - (long)a->size();
  }

  for (size_t i = 0; i < size(); i++) {
    bool found = false;

    for (size_t j = 0; j < size(); j++) {
      sort(&members[i], kind::UNDEF);
      sort(&a->members[i], kind::UNDEF);

      long cmp = compare(&members[i], &a->members[j], kind::UNDEF);

      if (cmp == 0) {
        found = true;
        break;
      }
    }

    if (found == false) {
      return -1;
    }
  }

  return 0;
}

list remove(list &a, size_t i) {
  list L = a;
  L.remove(i);
  return L;
}

list join(list &a, list &&b) {
  list L = a;
  L.join(b);
  return L;
}

std::string to_string(list &a) {
  std::string s = "{";
  for (size_t i = 0; i < a.size(); i++) {
    s += to_string(&a.members[i]);
    if (i < a.size() - 1) {
      s += ", ";
    }
  }

  s += "}";

  return s;
}

std::string to_string(list *a) {
  std::string s = "{";
  for (size_t i = 0; i < a->size(); i++) {
    s += to_string(&a->members[i]);

    if (i < a->size() - 1) {
      s += ", ";
    }
  }

  s += "}";

  return s;
}

std::string to_latex(list *a, bool f, unsigned long max_den) {
  std::string s = "\\left[";

  for (size_t i = 0; i < a->size(); i++) {
    s += to_latex(&a->members[i], f, max_den);

    if (i < a->size() - 1) {
      s += ", ";
    }
  }

  s += "//right]";

  return s;
}



bool list::operator==(list &a) { return this->match(&a) == 0; }

bool list::operator==(list &&a) { return this->match(&a) == 0; }

bool list::operator!=(list &a) { return this->match(&a) != 0; }

bool list::operator!=(list &&a) { return this->match(&a) != 0; }

long int sort_split(std::vector<expr> &a, kind k, long l, long r) {
  long int i = l - 1;

  expr &p = a[r];

  for (long int j = l; j < r; j++) {
    if (compare(&a[j], &p, k) < 0) {
      std::swap(a[++i], a[j]);
    }
  }

  i = i + 1;

  std::swap(a[i], a[r]);

  return i;
}

void sort_vec(std::vector<expr> &a, kind k, long int l, long int r) {
  if (l < r) {
    long int m = sort_split(a, k, l, r);

    sort_vec(a, k, l, m - 1);
    sort_vec(a, k, m + 1, r);
  }
}

std::vector<expr> expr::operands() {
  if (is(this, kind::LIST)) {
    return expr_list->members;
  }

  if (is(this, kind::SET)) {
    return expr_set->members;
  }

  return this->expr_childs;
}

set::set(std::initializer_list<expr> &&a) {
  members = std::move(a);
  trim(this);
}

set::set(std::vector<expr> &&a) {
  members = std::move(a);
  trim(this);
}

set::set(std::vector<expr> &a) {
  members = a;
  trim(this);
}

set difference(set &L, set &M) {
  set t = {};

  size_t j = 0;
  size_t i = 0;

  while (true) {
    if (j >= M.size() || i >= L.size()) {
      break;
    }

    int cmp = compare(&L[i], &M[j], kind::UNDEF);

    if (cmp < 0) {
      t.members.push_back(L[i++]);
    }

    if (cmp > 0) {
      t.members.push_back(M[j++]);
    }

    if (cmp == 0) {
      j = j + 1;
      i = i + 1;
    }
  }

  for (; i < L.size(); i++) {
    t.members.push_back(L[i]);
  }

  for (; j < M.size(); j++) {
    t.members.push_back(M[j]);
  }

  return t;
}

set unification(set &L, set &M) {
  set t = {};

  size_t j = 0;
  size_t i = 0;

  while (true) {
    if (j >= M.size() || i >= L.size()) {
      break;
    }

    int cmp = compare(&L[i], &M[j], kind::UNDEF);

    if (cmp > 0) {
      t.members.push_back(M[j++]);
    }

    if (cmp < 0) {
      t.members.push_back(L[i++]);
    }

    if (cmp == 0) {
      t.members.push_back(M[j++]);
      i = i + 1;
    }
  }

  for (; i < L.size(); i++) {
    t.members.push_back(L[i]);
  }

  for (; j < M.size(); j++) {
    t.members.push_back(M[j]);
  }

  return t;
}

set intersection(set &L, set &M) {
  set t = {};

  size_t j = 0;
  size_t i = 0;

  while (true) {
    if (j >= M.size() || i >= L.size()) {
      break;
    }

    int cmp = compare(&L[i], &M[j], kind::UNDEF);

    if (cmp < 0) {
      i = i + 1;
    }

    if (cmp > 0) {
      j = j + 1;
    }

    if (cmp == 0) {
      t.members.push_back(L[i++]);
      j = j + 1;
    }
  }

  return t;
}

int search(std::vector<expr> &a, expr &x, int l, int r) {
  if (r >= l) {
    int m = l + (r - l) / 2;

    int cmp = compare(&a[m], &x, kind::UNDEF);

    if (cmp == 0) {
      return m;
    }

    if (cmp > 1) {
      return search(a, x, l, m - 1);
    }

    return search(a, x, m + 1, r);
  }

  return -1;
}

bool set_exists(set &L, expr &e) {
  return search(L.members, e, 0, L.size() - 1) >= 0;
}

void set::sort() { sort_vec(members, kind::UNDEF, 0, size() - 1); }

void set::sort(enum kind k) { sort_vec(members, k, 0, size() - 1); }

set rest(set &a, size_t from) {
  std::vector<expr> l;

  for (size_t i = from; i < a.size(); i++) {
    l.push_back(a[i]);
  }

  return l;
}

expr fist(set &l) { return l[0]; }

void trim(set *s) {
  s->sort();

  for (long i = 0; i < ((long)s->size()) - 1; i++) {
    int cmp = compare(&s->members[i], &s->members[i + 1], kind::UNDEF);

    if (cmp == 0) {
      s->members.erase(s->members.begin() + i);
      i = i - 1;
    }
  }
}

long set::match(set *a) {
  if (size() != a->size()) {
    return (long)size() - (long)a->size();
  }
  for (size_t i = 0; i < size(); i++) {
    long cmp = compare(&members[i], &a->members[i], kind::UNDEF);
    if (cmp != 0) {
      return cmp;
    }
  }

  return 0;
}

bool set::operator==(set &a) { return this->match(&a) == 0; }

bool set::operator==(set &&a) { return this->match(&a) == 0; }

bool set::operator!=(set &a) { return this->match(&a) != 0; }

bool set::operator!=(set &&a) { return this->match(&a) != 0; }

std::string to_string(set *a) {
  std::string s = "{";

  for (size_t i = 0; i < a->size(); i++) {
    s += to_string(&a->members[i]);

    if (i < a->size() - 1) {
      s += ", ";
    }
  }

  s += "}";

  return s;
}



std::string to_string(set &a) {
  std::string s = "{";

  for (size_t i = 0; i < a.size(); i++) {
    s += to_string(&a.members[i]);

    if (i < a.size() - 1) {
      s += ", ";
    }
  }

  s += "}";

  return s;
}

std::string to_latex(set *a, bool f, unsigned long max_den) {
  std::string s = "\\left\\{";

  for (size_t i = 0; i < a->size(); i++) {
    s += to_latex(&a->members[i], f, max_den);

    if (i < a->size() - 1) {
      s += ", ";
    }
  }

  s += "//right\\}";

  return s;
}


expr gcd(expr &a, expr &b) {
  assert(a.kind() == kind::INT || a.kind() == kind::FRAC);
  assert(b.kind() == kind::INT || b.kind() == kind::FRAC);

  if (a.kind() == kind::INT && b.kind() == kind::INT) {
    return abs(gcd(a.value(), b.value()));
  }

  expr A = numerator(a);
  expr C = numerator(b);
  expr B = denominator(a);
  expr D = denominator(b);

  assert(A.kind() == kind::INT);
  assert(B.kind() == kind::INT);
  assert(C.kind() == kind::INT);
  assert(D.kind() == kind::INT);

  Int t = abs(gcd(A.value(), C.value()));
  Int k = abs(lcm(B.value(), D.value()));

  if (k == 1)
    return t;

  Int j = gcd(t, k);

  if (k / j == 1)
    return t;

  return fraction(t / j, k / j);
}

expr lcm(expr &a, expr &b) {
  assert(is(&a, kind::INT | kind::FRAC));
  assert(is(&b, kind::INT | kind::FRAC));

  if (a.kind() == kind::INT && b.kind() == kind::INT) {
    return abs(lcm(a.value(), b.value()));
  }

  expr A = numerator(a);
  expr C = numerator(b);
  expr B = denominator(a);
  expr D = denominator(b);

  assert(is(&A, kind::INT));
  assert(is(&B, kind::INT));
  assert(is(&C, kind::INT));
  assert(is(&D, kind::INT));

  Int t = abs(lcm(A.value(), C.value()));
  Int k = abs(gcd(B.value(), D.value()));

  if (k == 1)
    return t;

  Int j = gcd(t, k);

  if (k / j == 1)
    return t;

  return fraction(t / j, k / j);
}

expr binomial(Int n, std::vector<Int> &ks) {
  expr p = expr(kind::MUL);

  for (Int k : ks) {
    p.insert(fact(k));
  }

  return fact(n) / p;
}

expr binomial(Int n, std::vector<Int> &&ks) {
  expr p = expr(kind::MUL);

  for (Int k : ks) {
    p.insert(fact(k));
  }

  return fact(n) / p;
}

expr &base(expr &u) {
  assert(is(&u, kind::POW));
  return u[0];
}

expr &degree(expr &u) {
  assert(is(&u, kind::POW));
  return u[1];
}

expr &numerator(expr &u) {
  assert(is(&u, kind::FRAC));
  return u[0];
}

expr &denominator(expr &u) {
  assert(is(&u, kind::FRAC));
  return u[1];
}


expr abs(expr x) { return func_call("abs", {x}); }

expr diff(expr f, expr dx) { return func_call("diff", {f, dx}); }

void list::sortMembers() {
  long s = size();

  sort_vec(this->members, kind::UNDEF, 0, s - 1);
}

bool set::insert(const expr &a) {
  size_t i = 0;

  for (; i < size(); i++) {
    int cmp = compare((expr *)&a, &members[i], kind::UNDEF);

    if (cmp < 0) {
      members.insert(members.begin() + i, a);
      return true;
    }

    if (cmp == 0) {
      return false;
    }
  }

  members.push_back(a);
  return true;
}

bool set::insert(expr &&a) {
  size_t i = 0;

  for (; i < size(); i++) {
    int cmp = compare(&a, &members[i], kind::UNDEF);

    if (cmp < 0) {
      members.insert(members.begin() + i, std::move(a));
      return true;
    }

    if (cmp == 0) {
      return false;
    }
  }

  members.push_back(std::move(a));

  return true;
}

void set::remove(expr &a) {
  long idx = search(members, a, 0, members.size() - 1);

  if (idx == -1)
    return;

  remove(idx);
}

void set::remove(expr &&a) {
  long idx = search(members, a, 0, members.size() - 1);

  if (idx == -1)
    return;

  remove(idx);
}

void set::remove(size_t idx) { members.erase(members.begin() + idx); }

void decimalToFraction(double input, Int maxden, Int &n, Int &d) {
	assert(input < 1 && input >= 0);

  Int m[2][2];
  double x, startx;

  unsigned long long ai;

  startx = x = input;

  m[0][0] = m[1][1] = 1;
  m[0][1] = m[1][0] = 0;

  while (m[1][0] * (ai = (long)x) + m[1][1] <= maxden) {
    Int t = m[0][0] * ai + m[0][1];

    m[0][1] = m[0][0];
    m[0][0] = t;

    t = m[1][0] * ai + m[1][1];

    m[1][1] = m[1][0];
    m[1][0] = t;

    if (x == (double)ai) {
      break;
    }

    x = 1 / (x - (double)ai);

    if (x > (double)0x7FFFFFFF) {
      break;
    }
  }

  Int n0 = m[0][0];
  Int d0 = m[1][0];

  double err0 = startx - ((double)m[0][0].doubleValue() / (double)m[1][0].doubleValue());

  ai = (maxden.doubleValue() - m[1][1].doubleValue()) / m[1][0].doubleValue();

  m[0][0] = m[0][0] * ai + m[0][1];
  m[1][0] = m[1][0] * ai + m[1][1];

  Int n1 = m[0][0];
  Int d1 = m[1][0];

  double err1 = startx - (m[0][0].doubleValue() / m[1][0].doubleValue());

  if (fabs(err0) > fabs(err1)) {
    n = n1;
    d = d1;
  } else {
    n = n0;
    d = d0;
  }
}

set freeVariablesRec(expr &a) {
  if (is(&a, kind::SYM)) {
    return set({a});
  }

  if (is(&a, kind::TERMINAL)) {
    return {};
  }

  set r = {};

  for (size_t i = 0; i < size_of(&a); i++) {
    set t = freeVariablesRec(a[i]);

    r = unification(r, t);
  }

  return r;
}

list freeVariables(expr &a) {
  set t = freeVariablesRec(a);

  list r = {};

  for (size_t i = 0; i < t.size(); i++) {
    r.insert(t[i]);
  }

  return r;
}


expr inverse_matrix(expr A) {
	assert(is(&A, kind::MAT));
	double eps = std::numeric_limits<double>::epsilon();

	if (rows(A) == columns(A) && matrix::det_ptr(A.expr_mat) > eps) {
		// TODO: if condition above are not true a error should be
		//  emited, maybe? ...
		matrix *t = matrix::inv_ptr(A.expr_mat);

		expr r = expr(kind::MAT);

		r.expr_mat = t;

		return r;
	}

	return fail();
}

expr svd_matrix(expr A) {
	assert(is(&A, kind::MAT));
	std::tuple<matrix*, matrix*, matrix*> t = matrix::svd_ptr(A.expr_mat);

	expr U  = expr(kind::MAT);
	expr D  = expr(kind::MAT);
	expr VT = expr(kind::MAT);

	U.expr_mat  = std::get<0>(t);
	D.expr_mat  = std::get<1>(t);
	VT.expr_mat = std::get<2>(t);

	return list({ U, D, VT });
}

expr transpose_matrix(expr A) {
	assert(is(&A, kind::MAT));
	expr U  = expr(kind::MAT);

	U.expr_mat = matrix::transp_ptr(A.expr_mat);

	return U;
}

expr determinant_matrix(expr A) {
	assert(is(&A, kind::MAT));
	return fromDouble(matrix::det_ptr(A.expr_mat));
}

expr solve_linear_system(expr A, expr b) {
	assert(is(&A, kind::MAT));
	assert(is(&b, kind::MAT));

	expr x = expr(kind::MAT);

	x.expr_mat = matrix::solve_ptr(A.expr_mat, b.expr_mat);

	return x;
}

expr exp(expr x) { return func_call("exp", {x}); }

expr log(expr x, expr base) { return func_call("log", {x, base}); }

expr ln(expr x) { return func_call("ln", {x}); }


double doubleFromExpr(expr a) {
	assert(is(&a, kind::INT | kind::FRAC));

	if(is(&a, kind::INT)) return get_val(&a).doubleValue();
	if(is(&a, kind::FRAC)) {
		expr n = numerator(a);
		expr d = numerator(a);

		return n.value().doubleValue() / d.value().doubleValue();
	}

	return std::numeric_limits<double>::quiet_NaN();
}

} // namespace alg
