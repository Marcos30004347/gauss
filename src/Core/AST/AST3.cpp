#include "AST3.hpp"

#include "Core/AST/Integer.hpp"
#include "Core/Algebra/Set.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <math.h>
#include <random>
#include <string>
#include <vector>

namespace alg {

expr::expr(expr &&other) {
  //	*this = std::move(other);
  kind_of = other.kind_of;

  switch (kind_of) {
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

  case kind::FUNC: {
    expr_sym = other.expr_sym;
    other.expr_sym = 0;
    expr_childs = std::move(other.expr_childs);

    return;
  }

  case kind::INF:
    return;
  case kind::UNDEF:
    return;
  case kind::FAIL:
    return;

  default: {
    expr_childs = other.expr_childs;
    return;
  }
  }
}

expr::expr(const expr &other) {
  kind_of = other.kind_of;

  switch (kind_of) {
  case kind::SYM: {
    expr_sym = strdup(other.expr_sym);
    return;
  }
  case kind::INT: {
    expr_int = new Int(*other.expr_int);
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

  case kind::INF:
    return;
  case kind::UNDEF:
    return;
  case kind::FAIL:
    return;

  default: {
    expr_childs = other.expr_childs;
    return;
  }
  }
}
expr &expr::operator=(const expr &other) {
  kind_of = other.kind_of;

  switch (kind_of) {

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

  case kind::INF:
    return *this;
  case kind::UNDEF:
    return *this;
  case kind::FAIL:
    return *this;

  default: {
    expr_childs = other.expr_childs;
    return *this;
  }
  }

  return *this;
}

expr &expr::operator=(expr &&other) {
  kind_of = other.kind_of;

  switch (kind_of) {

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

  case kind::INF:
    return *this;
  case kind::UNDEF:
    return *this;
  case kind::FAIL:
    return *this;

  default: {
    expr_childs = std::move(other.expr_childs);
    return *this;
  }
  }

  return *this;
}

expr::expr(enum kind k) { kind_of = k; }

enum kind expr::kind() const { return kind_of; }

expr::expr(list &s) {
  kind_of = kind::LIST;
  expr_list = new list(s);
}

expr::expr(list &&s) {
  kind_of = kind::LIST;
  expr_list = new list(s);
}

expr::expr(set &s) {
  kind_of = kind::SET;
  expr_set = new set(s);
}

expr::expr(set &&s) {
  kind_of = kind::SET;
  expr_set = new set(s);
}

expr::expr(enum kind k, std::initializer_list<expr> &&a) {
  if (k == kind::LIST) {
    this->expr_list = new list(a);
    return;
  }

  if (k == kind::SET) {
    this->expr_set = new set(a);
    return;
  }

  kind_of = k;

  expr_childs = std::move(a);
}

size_t expr::size() { return size_of(this); }

expr::expr(Int v) {
  kind_of = kind::INT;
  this->expr_int = new Int(v);
}

expr::expr(int v) {
  kind_of = kind::INT;
  this->expr_int = new Int(v);
}

expr::expr(long int v) {
  kind_of = kind::INT;
  this->expr_int = new Int(v);
}

expr::expr(long long v) {
  kind_of = kind::INT;
  this->expr_int = new Int(v);
}

expr::expr(std::string v) {
  kind_of = kind::SYM;
  this->expr_sym = strdup(v.c_str());
}

expr::expr() { kind_of = kind::UNDEF; }

expr::~expr() {
  switch (kind_of) {

  case kind::INT: {
    if (expr_int)
      delete expr_int;
    return;
  }

  case kind::SYM: {
    if (expr_sym)
      delete expr_sym;
    return;
  }

  case kind::LIST: {
    if (expr_list)
      delete expr_list;
    return;
  }

  case kind::SET: {
    if (expr_set)
      delete expr_set;
    return;
  }

  case kind::FUNC: {
    if (expr_sym)
      delete expr_sym;
  }

  default:
    return;
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

void expr_set_kind(expr *a, kind kind) {
  if (is(a, kind::SYM)) {
    delete a->expr_sym;
  }

  if (is(a, kind::INT)) {
    delete a->expr_sym;
  }

  if (is(a, kind::LIST)) {
    delete a->expr_list;
  }

  if (is(a, kind::SET)) {
    delete a->expr_set;
  }

  a->kind_of = kind;
}

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
  this->expr_childs.insert(this->expr_childs.begin() + idx, b);
}

void expr::insert(expr &&b, size_t idx) {
  this->expr_childs.insert(this->expr_childs.begin() + idx, std::move(b));
}

void expr::insert(const expr &b) { this->expr_childs.push_back(b); }
void expr::insert(expr &&b) { this->expr_childs.push_back(std::move(b)); }

void expr::remove(list &l) {
  assert(is(this, kind::LIST));
  expr_list->remove(l);
}

void expr::remove(list &&l) {
  assert(is(this, kind::LIST));
  expr_list->remove(std::move(l));
}

void expr::remove(size_t idx) {
  if (is(this, kind::LIST)) {
    this->expr_list->members.erase(this->expr_childs.begin() + idx);
    return;
  }

  if (is(this, kind::SET)) {
    this->expr_set->members.erase(this->expr_childs.begin() + idx);
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

int compare_consts(expr *a, expr *b) {
  assert(is(a, kind::CONST) && is(b, kind::CONST));

  if (is(a, kind::INT) && is(b, kind::INT)) {
		if (get_val(a) == get_val(b)) {
      return 0;
    }

    return get_val(a) > get_val(b) ? 1 : -1;
  }

  if (is(a, kind::FRAC) && is(b, kind::FRAC)) {
    Int na = get_val(operand(a, 0));
    Int da = get_val(operand(a, 1));
    Int nb = get_val(operand(b, 0));
    Int db = get_val(operand(b, 1));

    if (na * db == nb * da) {
      return 0;
    }

    return na * db - nb * da > 0 ? 1 : -1;
  }

  if (is(a, kind::INT) && is(b, kind::FRAC)) {
    Int na = get_val(a);
    Int nb = get_val(operand(b, 0));
    Int db = get_val(operand(b, 1));

    Int ct = na * db;

    if (ct == nb) {
      return 0;
    }

    return ct > nb ? 1 : -1;
  }

  Int nb = get_val(b);
  Int na = get_val(operand(a, 0));
  Int da = get_val(operand(a, 1));

  Int ct = nb * da;

  if (ct == na) {
    return 0;
  }

  return na > ct ? 1 : -1;
}

inline bool should_revert_idx(kind ctx) { return ctx & (kind::ADD); }

inline bool expr_is_zero(expr *a) {
  return (a == nullptr) || (is(a, kind::INT) && get_val(a) == 0);
}

inline int expr_op_cmp(expr *a, expr *b, kind ctx) {
  long m = size_of(a);
  long n = size_of(b);

  long l = std::min(size_of(a), size_of(b));

  m = m - 1;
  n = n - 1;

  if (ctx == kind::ADD) {
    if (is(a, kind::MUL) && is(b, kind::MUL)) {
			for (long i = 0; i < l; i++) {

				int order = compare(operand(a, m - i), operand(b, n - i), ctx);
					// printf("    %s -- %s ", to_string(operand(a, i + j)).c_str(),  to_string(operand(b, i + k)).c_str());
					// printf("result in %i\n", order);

        if (order) {
					return order;
				}
      }

		// 	short j = is(operand(a, 0), kind::CONST) ? 1 : 0;
		// 	short k = is(operand(b, 0), kind::CONST) ? 1 : 0;

		// 	short r = j & k;

		// 	if ((size_of(a) - j) - (size_of(b) - k) != 0) {
    //     return n - m;
    //   }

		// 	printf("cmp %s -- %s\n", to_string(a).c_str(),  to_string(b).c_str());

		// 	for (long i = 0; i < l - r; i++) {

		// 		int order = compare(operand(a, i + j), operand(b, i + k), ctx);
		// 			printf("    %s -- %s ", to_string(operand(a, i + j)).c_str(),  to_string(operand(b, i + k)).c_str());
		// 			printf("result in %i\n", order);

    //     if (order) {

		// 			return order;
		// 		}
    //   }
    // }
		}
    return n - m;
  }

  if (ctx == kind::MUL) {
    for (long i = 0; i < l; i++) {
      int order = kind_of(operand(b, n - i)) - kind_of(operand(a, m - i));

      if (order)
        return order;
    }

    for (long i = 0; i < l; i++) {
      int order = compare(operand(a, m - i), operand(b, n - i), ctx);

      if (order)
        return order;
    }

    return n - m;
  }

  if (size_of(a) != size_of(b)) {
    return (long)size_of(b) - (long)size_of(a);
  }

  for (long i = 0; i < l; i++) {
    int order = compare(operand(a, l - 1 - i), operand(b, l - 1 - i), ctx);

    if (order) {
      return order;
    }
  }

  return 0;
}

inline int compare_idents(expr *a, expr *b) {
  return strcmp(get_id(a), get_id(b));
}

std::string to_string(expr *tree) {
  if (!tree)
    return "null";

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

  if (is(tree, kind::FRAC)) {
    return "sqrt(" + to_string(operand(tree, 0)) + ")";
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
      if (operand(tree, i) &&
          is(operand(tree, i), kind::SUB | kind::ADD )) {

        r += "(";
      }

      r += to_string(operand(tree, i));

      if (operand(tree, i) &&
          is(operand(tree, i), kind::SUB | kind::ADD)) {
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
      if (operand(tree, i) &&
          is(operand(tree, i), kind::SUB | kind::ADD)) {

        r += "(";
      }

      r += to_string(operand(tree, i));

      if (operand(tree, i) &&
          is(operand(tree, i), kind::SUB | kind::ADD)) {
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

      if (is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL)) {
        r += "(";
      }

      r += to_string(operand(tree, i));

      if (is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL)) {
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

  if (is(tree, kind::LIST))
    return to_string(tree->expr_list);
  if (is(tree, kind::SET))
    return to_string(tree->expr_set);

  return "to_string_not_implemented";
}

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
  case kind::SQRT: {
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
  default:
    return "";
  }
}

void expr_print(expr *a, int tabs) {
  printf("%*c<expr ", tabs, ' ');
  printf("address=\"%p\" ", a);
  printf("kind=\"%s\"", kind_of_id(a).c_str());

  if (kind_of(a) == kind::INT) {
    printf(" value=\"%s\"", get_val(a).to_string().c_str());
  }

  if (kind_of(a) == kind::SYM || kind_of(a) == kind::FUNC) {
    printf(" id=\"%s\"", get_id(a));
  }

  if (size_of(a)) {
    printf(" size=\"%li\" >\n", size_of(a));

    for (size_t i = 0; i < size_of(a); i++) {
      expr_print(operand(a, i), tabs + 3);
    }

    printf("%*c</expr>\n", tabs, ' ');
  } else {
    printf(">\n");
  }
}

int compare(expr *const a, expr *const b, kind ctx) {
  if (a == b) {
    return 0;
	}

  if (ctx & kind::MUL) {
    if (is(a, kind::CONST) && is(b, kind::CONST)) {
      return compare_consts(a, b);
    }

    if (is(a, kind::CONST)) {
      return -1;
    }

    if (is(b, kind::CONST)) {
      return +1;
    }

    if (is(a, kind::POW) && is(b, kind::POW)) {
      return compare(operand(a, 0), operand(b, 0), ctx);
    }

    if (is(a, kind::SYM | kind::ADD) && is(b, kind::POW)) {
      int order = compare(a, operand(b, 0), kind::MUL);

      if (order == 0) {
        return kind_of(a) - kind_of(b);
      }

      return order;
    }

    if (is(b, kind::SYM | kind::ADD) && is(a, kind::POW)) {
      int order = compare(operand(a, 0), b, kind::MUL);

      if (order == 0) {
        return kind_of(a) - kind_of(b);
      }

      return order;
    }

    if (is(a, kind::FUNC) && is(b, kind::FUNC)) {
      return strcmp(get_func_id(a), get_func_id(b));
    }

    if (is(a, kind::POW) && is(b, kind::FUNC)) {
      return compare(operand(a, 0), b, kind::MUL);
    }

    if (is(b, kind::POW) && is(a, kind::FUNC)) {
      return compare(b, operand(a, 0), kind::MUL);
    }

    if (is(a, kind::MUL) && is(b, kind::POW | kind::SYM | kind::FUNC)) {
      long k = size_of(a);

      for (long i = k - 1; i >= 0; i--) {
        int order = compare(operand(a, i), b, kind::MUL);

        if (order == 0) {
          return 0;
        }

        if (order < 0)
          break;
      }

      return +k - 1;
    }

    if (is(b, kind::MUL) && is(a, kind::POW | kind::SYM | kind::FUNC)) {
      long k = size_of(b);

      for (long i = k - 1; i >= 0; i--) {
        int order = compare(operand(b, i), a, kind::MUL);

        if (order == 0) {
          return 0;
        }

        if (order < 0)
          break;
      }

      return -k + 1;
    }

    if (is(a, kind::ADD) && is(b, kind::ADD)) {
      return expr_op_cmp(a, b, ctx);
    }

    if (is(a, kind::MUL) && is(b, kind::MUL)) {
      return expr_op_cmp(a, b, ctx);
    }

    if (is(a, kind::LIST) && is(b, kind::LIST)) {
      return a->expr_list->match(b->expr_list);
    }

    if (is(a, kind::SET) && is(b, kind::SET)) {
      return a->expr_set->match(b->expr_set);
    }
  }

  if (ctx & kind::ADD) {
    // printf("--> %s  cmp  %s = ", to_string(a).c_str(), to_string(b).c_str());

    if (is(a, kind::CONST) && is(b, kind::CONST)) {
      // printf("%i\n", compare_consts(b, a));
      return compare_consts(b, a);
    }

    if (is(a, kind::CONST)) {
      // printf("%i\n", +1);
      return +1;
    }

    if (is(b, kind::CONST)) {
      // printf("%i\n", -1);
      return -1;
    }

    if (is(a, kind::POW) && is(b, kind::POW)) {
      int i = compare(operand(a, 1), operand(b, 1), ctx);

      if (i != 0) {
        // printf("%i\n", i);

        return i;
      }

      // printf("%i\n", compare(operand(a, 0), operand(b, 0), ctx));
      return compare(operand(a, 0), operand(b, 0), ctx);
    }

    if(is(a, kind::POW) && is(b, kind::MUL)) {
			// printf("compare = %s %s\n", to_string(a).c_str(), to_string(b).c_str());
			if(size_of(b) == 2) {
    		if(compare(a, operand(b, 1), ctx) == 0) {
    			return 0;
    		}
    	}

			return  (int)kind_of(a) - (int)kind_of(b);
    }

    if (is(b, kind::POW) && is(a, kind::MUL)) {
			// printf("compare = %s %s\n", to_string(a).c_str(), to_string(b).c_str());
      if (size_of(a) == 2) {
        if (compare(b, operand(a, 1), ctx) == 0) {
          return 0;
        }
      }
			return  (int)kind_of(a) - (int)kind_of(b);
    }

    if (is(a, kind::FUNC) && is(b, kind::FUNC)) {
      // printf("%i\n", strcmp(get_func_id(a), get_func_id(b)));
      return strcmp(get_func_id(a), get_func_id(b));
    }

    if (is(a, kind::ADD) && is(b, kind::SYM)) {
      // printf("%i\n", 1);
      return +1;
    }

    if (is(a, kind::SYM) && is(b, kind::ADD)) {
      // printf("%i\n", -1);
      return -1;
    }

    if (is(a, kind::MUL) && is(b, kind::SYM)) {
      long k = size_of(a);

      if (k > 2) {
        // printf("%li\n", -k);
        return -k;
      }

      int order = compare(operand(a, size_of(a) - 1), b, ctx);

      if (order == 0) {
        // printf("%i\n", -1);
        return -1;
      }

      // printf("%i\n", k > 2 ? -1 : order);

      return k > 2 ? -1 : order;
    }

    if (is(b, kind::MUL) && is(a, kind::SYM)) {
      long k = size_of(b);

      if (k > 2) {
        return +k;
      }

      int order = compare(a, operand(b, size_of(b) - 1), ctx);

      if (order == 0) {
        return +1;
      }

      return k > 2 ? +1 : order;
    }

    if (is(a, kind::MUL) && is(b, kind::POW)) {
      long k = size_of(a);

      if (k > 2) {
        // printf("%li\n", -k);
        return -k;
      }

      int order = compare(operand(a, 0), b, ctx);

      if (order == 0) {
        // printf("%li\n", -k);
        return -1;
      }

      // printf("%i\n", k > 2 ? -1 : order);
      return k > 2 ? -1 : order;
    }

    if (is(b, kind::MUL) && is(a, kind::POW)) {
      long k = size_of(b);

      if (k > 2) {
        // printf("%li\n", +k);
        return +k;
      }

      int order = compare(a, operand(b, 0), ctx);

      if (order == 0) {
        // printf("%i\n", +1);
        return +1;
      }

      // printf("%i\n", k > 2 ? -1 : order);
      return k > 2 ? +1 : order;
    }

    if (is(a, kind::ADD) && is(b, kind::ADD)) {
      return expr_op_cmp(a, b, ctx);
    }

    if (is(a, kind::MUL) && is(b, kind::MUL)) {
      return expr_op_cmp(a, b, ctx);
    }

    if (is(a, kind::LIST) && is(b, kind::LIST)) {
      return a->expr_list->match(b->expr_list);
    }

    if (is(a, kind::SET) && is(b, kind::SET)) {
      return a->expr_set->match(b->expr_set);
    }
  }

  if (is(a, kind::FUNC) && is(b, kind::FUNC)) {
    return strcmp(get_func_id(a), get_func_id(b));
  }

  if (is(a, kind::CONST) && is(b, kind::CONST)) {
    return compare_consts(a, b);
  }

  if (is(a, kind::SYM) && is(b, kind::SYM)) {
    return compare_idents(a, b);
  }

  if (is(a, kind::ADD) && is(b, kind::ADD)) {
    return expr_op_cmp(a, b, ctx);
  }

  if (is(a, kind::MUL) && is(b, kind::MUL)) {
    return expr_op_cmp(a, b, ctx);
  }

  if (is(a, kind::POW) && is(b, kind::POW)) {
    int cmp = compare(operand(a, 0), operand(b, 0), ctx);

    if (cmp != 0) {
      return cmp;
    }

    return compare(operand(a, 1), operand(b, 1), ctx);
  }

  if (is(a, kind::DIV) && is(b, kind::DIV)) {
    int cmp = compare(operand(a, 0), operand(b, 0), ctx);

    if (cmp != 0) {
      return cmp;
    }

    return compare(operand(a, 1), operand(b, 1), ctx);
  }

  if (is(a, kind::SQRT) && is(b, kind::SQRT)) {
    return compare(operand(a, 0), operand(b, 0), ctx);
  }

  if (is(a, kind::FACT) && is(b, kind::FACT)) {
    return compare(operand(a, 0), operand(b, 0), ctx);
  }

  if (is(a, kind::LIST) && is(b, kind::LIST)) {
    if (size_of(a) != size_of(b)) {
      return (long)size_of(a) - (long)size_of(b);
    }

    for (size_t i = 0; i < size_of(a); i++) {
      int order =
          compare(&a->expr_list->members[i], &b->expr_list->members[i], ctx);
      if (order != 0) {
        return order;
      }
    }

    return 0;
  }

  if (is(a, kind::SET) && is(b, kind::SET)) {
    if (size_of(a) != size_of(b)) {
      return (long)size_of(a) - (long)size_of(b);
    }

    for (size_t i = 0; i < size_of(a); i++) {
      int order =
          compare(&a->expr_set->members[i], &b->expr_set->members[i], ctx);
      if (order != 0) {
        return order;
      }
    }

    return 0;
  }

  return should_revert_idx(ctx) ? kind_of(a) - kind_of(b)
                                : kind_of(b) - kind_of(a);
}
// int compare(expr *const a, expr *const b, kind ctx) {
// 	int t = compare_rec(a, b, ctx);
// 	if(t == 0) return t;
// 	return t <  0 ? -1 : 1;
// }
long int sort_split(expr *a, long l, long r) {
  long int i = l - 1;

  expr *p = operand(a, r);

  for (long int j = l; j < r; j++) {
    if (compare(operand(a, j), p, kind_of(a)) < 0) {
      std::swap(a->expr_childs[++i], a->expr_childs[j]);
    }
  }

  std::swap(a->expr_childs[++i], a->expr_childs[r]);

  return i;
}

void sort_childs(expr *a, long int l, long int r) {
  if (l < r) {
    long int m = sort_split(a, l, r);

    sort_childs(a, l, m - 1);
    sort_childs(a, m + 1, r);
  }
}

void sort(expr *a) {
  if (is(a, kind::TERMINAL)) {
    return;
  }

  if (is(a, kind::LIST)) {
    for (size_t i = 0; i < size_of(a); i++) {
      sort(&a->expr_list->members[i]);
    }

    return;
  }

  if (is(a, kind::SET)) {
    return a->expr_set->sort();
  }

  size_t i = is(a, kind::SUB) ? 1 : 0;

  for (; i < size_of(a); i++) {
    sort(operand(a, i));
  }

  if (is(a, kind::ORDERED)) {
    return;
  }

  sort_childs(a, 0, size_of(a) - 1);
}

long int sort_split(expr *a, kind k, long l, long r) {
  long int i = l - 1;

  expr *p = operand(a, r);

  for (long int j = l; j < r; j++) {
    if (compare(operand(a, j), p, k) < 0) {
      std::swap(a->expr_childs[++i], a->expr_childs[j]);
    }
  }

  std::swap(a->expr_childs[++i], a->expr_childs[r]);

  return i;
}

void sort_childs(expr *a, kind k, long int l, long int r) {
  if (l < r) {
    long int m = sort_split(a, k, l, r);

    sort_childs(a, k, l, m - 1);
    sort_childs(a, k, m + 1, r);
  }
}

void sort(expr *a, kind k) {
  if (is(a, kind::TERMINAL)) {
    return;
  }

  if (is(a, kind::LIST)) {
    for (size_t i = 0; i < size_of(a); i++) {
      sort(&a->expr_list->members[i], k);
    }

    return;
  }

  if (is(a, kind::SET)) {
    return a->expr_set->sort(k);
  }

  size_t i = is(a, kind::SUB) ? 1 : 0;

  for (; i < size_of(a); i++) {
    sort(operand(a, i), k);
  }

  if (is(a, kind::ORDERED)) {
    return;
  }

  sort_childs(a, k, 0, size_of(a) - 1);
}

inline bool is_negative(expr *a) {
  assert(is(a, kind::CONST));

  if (is(a, kind::INT) && get_val(a) < 0)
    return true;

  if (is(a, kind::FRAC) &&
      ((get_val(operand(a, 0)) < 0 && get_val(operand(a, 1)) > 0) ||
       (get_val(operand(a, 0)) > 0 && get_val(operand(a, 1)) < 0)))
    return true;

  return false;
}

inline void expr_set_to_undefined(expr *a) { expr_set_kind(a, kind::UNDEF); }

inline void expr_set_to_fail(expr *a) { expr_set_kind(a, kind::FAIL); }

inline void expr_set_to_int(expr *a, Int v) {

  expr_set_kind(a, kind::INT);

	a->expr_childs.clear();

  a->expr_int = new Int(v);
}

inline void expr_set_to_inf(expr *a) { expr_set_kind(a, kind::INF); }

inline void expr_set_to_neg_inf(expr *a) {
  expr_set_kind(a, kind::MUL);

  a->insert(integer(-1));
  a->insert(inf());
}

inline bool is_inf(expr *a) { return is(a, kind::INF); }
inline bool is_neg_inf(expr *a) {
  if (!is(a, kind::MUL)) {
    return false;
  }

  bool have_neg = false;
  bool have_inf = false;

  for (size_t i = 0; i < size_of(a); i++) {
    if (is(operand(a, i), kind::CONST) && is_negative(operand(a, i))) {
      have_neg = true;
    }

    if (is_inf(operand(a, i))) {
      have_inf = true;
    }
  }

  return have_neg && have_inf;
}

inline bool is_inv_inf(expr *a) {
  if (!is(a, kind::POW))
    return false;
  return (is_inf(operand(a, 0)) || is_neg_inf(operand(a, 0))) &&
         is(operand(a, 1), kind::CONST) && is_negative(operand(a, 1));
}

inline void expr_set_op_to_int(expr *a, size_t i, Int v) {
  expr_set_to_int(operand(a, i), v);
}

inline void expr_set_to_fra(expr *a, Int u, Int v) {
  a->expr_childs = std::vector<expr>();

  expr_set_kind(a, kind::FRAC);

  a->insert(integer(u));
  a->insert(integer(v));
}

inline void expr_set_op_to_fra(expr *a, size_t i, Int u, Int v) {
  expr_set_to_fra(operand(a, i), u, v);
}

inline void expr_set_to_sym(expr *a, const char *s) {
  expr_set_kind(a, kind::SYM);
  a->expr_sym = strdup(s);
}

inline void expr_set_op_to_sym(expr *a, size_t i, const char *s) {
  expr_set_to_sym(operand(a, i), s);
}

// a = a + b
inline void expr_set_inplace_add_consts(expr *a, expr *b) {
  assert(is(a, kind::CONST));
  assert(is(b, kind::CONST));

  if (is(a, kind::INT) && is(b, kind::INT)) {
    Int x = get_val(a);
    Int y = get_val(b);

    return expr_set_to_int(a, x + y);
  }

  if (is(a, kind::FRAC) && is(b, kind::FRAC)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));
    Int w = get_val(operand(b, 0));
    Int z = get_val(operand(b, 1));

    Int e = x * z + w * y;
    Int k = y * z;

    if (e % k == 0) {
      return expr_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));

    return expr_set_to_fra(a, e / g, k / g);
  }

  if (is(a, kind::FRAC) && is(b, kind::INT)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));
    Int w = get_val(b);

    Int e = x + w * y;
    Int k = y;

    if (e % k == 0) {
      return expr_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));
    return expr_set_to_fra(a, e / g, k / g);
  }

  if (is(b, kind::FRAC) && is(a, kind::INT)) {
    Int x = get_val(operand(b, 0));
    Int y = get_val(operand(b, 1));
    Int w = get_val(a);

    Int e = x + w * y;
    Int k = y;

    if (e % k == 0) {
      return expr_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));

    return expr_set_to_fra(a, e / g, k / g);
  }
}

inline void expr_set_inplace_add_consts(expr *a, Int b) {
  assert(is(a, kind::CONST));

  if (is(a, kind::INT)) {
    Int x = get_val(a);

    return expr_set_to_int(a, x + b);
  }

  if (is(a, kind::FRAC)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));

    Int e = x + b * y;

		if(e == 0) {
			return expr_set_to_int(a, 0);
		}

    if (e % y == 0) {
      return expr_set_to_int(a, e / y);
    }

    Int g = abs(gcd(e, y));
    return expr_set_to_fra(a, e / g, y / g);
  }
}

inline void expr_set_op_inplace_add_consts(expr *a, size_t i, expr *b) {
  expr_set_inplace_add_consts(operand(a, i), b);
}

inline void expr_set_op_inplace_add_consts(expr *a, size_t i, Int b) {
  expr_set_inplace_add_consts(operand(a, i), b);
}

inline void expr_set_inplace_mul_consts(expr *a, expr *b) {
  assert(is(a, kind::CONST));
  assert(is(b, kind::CONST));

  if (is(a, kind::INT) && is(b, kind::INT)) {
    Int x = get_val(a);
    Int y = get_val(b);

    return expr_set_to_int(a, x * y);
  }

  if (is(a, kind::FRAC) && is(b, kind::FRAC)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));
    Int w = get_val(operand(b, 0));
    Int z = get_val(operand(b, 1));

    Int e = x * w;
    Int k = y * z;

		if(e == 0) {
			return expr_set_to_int(a, 0);
		}

    if (e % k == 0) {
      return expr_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));
    return expr_set_to_fra(a, e / g, k / g);
  }

  if (is(a, kind::FRAC) && is(b, kind::INT)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));
    Int w = get_val(b);

    Int e = x * w;
    Int k = y;

		if(e == 0) {
			return expr_set_to_int(a, 0);
		}

		if (e % k == 0) {
      return expr_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));
    return expr_set_to_fra(a, e / g, k / g);
  }

  if (is(b, kind::FRAC) && is(a, kind::INT)) {
    Int x = get_val(operand(b, 0));
    Int y = get_val(operand(b, 1));
    Int w = get_val(a);

    Int e = x * w;
    Int k = y;

		if(e == 0) {
			return expr_set_to_int(a, 0);
		}

    if (e % k == 0) {
      return expr_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));
    return expr_set_to_fra(a, e / g, k / g);
  }
}

inline void expr_set_inplace_mul_consts(expr *a, Int b) {
  assert(is(a, kind::CONST));

  if (is(a, kind::INT)) {
    Int x = get_val(a);

    return expr_set_to_int(a, x * b);
  }

  if (is(a, kind::FRAC)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));

    Int e = x * b;

    if (e % y == 0) {
      return expr_set_to_int(a, e / y);
    }

    Int g = abs(gcd(e, y));
    return expr_set_to_fra(a, e / g, y / g);
  }
}

inline void expr_set_op_inplace_mul_consts(expr *a, size_t i, expr *b) {
  expr_set_inplace_mul_consts(operand(a, i), b);
}

inline void expr_set_op_inplace_mul_consts(expr *a, size_t i, Int b) {
  expr_set_inplace_mul_consts(operand(a, i), b);
}

inline void expr_set_to_mul(Int v, expr *a) {
  if (is(a, kind::MUL)) {
    return a->insert(integer(v), 0);
  }

  a->expr_childs = std::vector<expr>({integer(v), *a});

  expr_set_kind(a, kind::MUL);
}

inline void expr_set_to_mul(expr *a, expr *b) {
  if (is(a, kind::MUL)) {
    return a->insert(*b);
  }

  a->expr_childs = std::vector<expr>({*b, *a});
  expr_set_kind(a, kind::MUL);
}

inline void expr_set_op_to_mul(expr *a, size_t i, Int v) {
  return expr_set_to_mul(v, operand(a, i));
}

inline void expr_set_op_to_mul(expr *a, size_t i, expr *v) {
  return expr_set_to_mul(operand(a, i), v);
}

inline void expr_set_to_pow(expr *a, Int e) {
  a->expr_childs = std::vector<expr>({*a, integer(e)});
  expr_set_kind(a, kind::POW);
}

inline void expr_set_to_add(expr *a, expr *e) {

  if (is(a, kind::ADD)) {
    return a->insert(*e);
  }

  a->expr_childs = std::vector<expr>({*a, *e});

  expr_set_kind(a, kind::ADD);
}

inline void expr_set_op_to_add(expr *a, size_t i, expr *v) {
  expr_set_to_add(operand(a, i), v);
}

inline void expr_set_op_pow_add_to_deg(expr *a, size_t i, expr *e) {
  assert(is(operand(a, i), kind::POW));

  a->expr_childs[i] =
      create(kind::POW, {*operand(operand(a, i), 0),
                         create(kind::ADD, {*operand(operand(a, i), 1), *e})});
}

inline void expr_set_op_pow_add_to_deg(expr *a, size_t i, Int e) {
  assert(is(operand(a, i), kind::POW));

  a->expr_childs[i] = create(
      kind::POW, {*operand(operand(a, i), 0),
                  create(kind::ADD, {*operand(operand(a, i), 1), integer(e)})});
}

inline void expr_set_op_to_pow(expr *a, size_t i, Int v) {
  a->expr_childs[i] = create(kind::POW, {*operand(a, i), integer(v)});
}

inline void expr_set_op_to_pow(expr *a, size_t i, expr *v) {
  a->expr_childs[i] = create(kind::POW, {*operand(a, i), *v});
}

inline bool expr_replace_with(expr *a, expr *t) {
  if (is(t, kind::INT)) {
    expr_set_to_int(a, get_val(t));
    return true;
  }

  if (is(t, kind::SYM)) {
    expr_set_to_sym(a, get_id(t));
    return true;
  }

  if (is(t, kind::FUNC)) {
    a->expr_sym = strdup(get_id(t));
    expr_set_kind(a, kind::FUNC);

    a->expr_childs = std::vector<expr>();

    // TODO: maybe std::move works here
    for (size_t i = 0; i < size_of(t); i++) {
      a->insert(t->expr_childs[i]);
    }

    return a;
  }

  expr_set_kind(a, kind_of(t));

  a->expr_childs = std::vector<expr>();

  // TODO: maybe std::move works here
  for (size_t i = 0; i < size_of(t); i++) {
    a->insert(t->expr_childs[i]);
  }

  return a;
}

inline bool eval_add_consts(expr *u, size_t i, expr *v, size_t j) {
  expr_set_op_inplace_add_consts(u, i, operand(v, j));
  return true;
}

inline bool eval_mul_consts(expr *u, size_t i, expr *v, size_t j) {
  expr_set_op_inplace_mul_consts(u, i, operand(v, j));
  return true;
}

inline bool eval_add_int(expr *a, Int b) {
  if (is(a, kind::INT)) {
    expr_set_to_int(a, get_val(a) + b);

    return true;
  }

  assert(is(a, kind::FRAC));

  Int num = get_val(operand(a, 0));
  Int den = get_val(operand(a, 1));

  num = b * den + num;

  Int cff = abs(gcd(num, den));

  if (den / cff == 1) {
    expr_set_to_int(a, num / cff);

    return true;
  }

  expr_set_to_fra(a, num / cff, den / cff);

  return true;
}

inline bool eval_add_nconst(expr *u, size_t i, expr *v, size_t j) {
  assert(is(u, kind::ADD) && is(v, kind::ADD));

  assert(is(operand(u, i), kind::SUMMABLE));
  assert(is(operand(v, j), kind::SUMMABLE));

  expr *a = operand(u, i);
  expr *b = operand(v, j);

  if (a == 0 || b == 0) {
    return false;
  }

  if (is_inf(a)) {

    if (is_neg_inf(b)) {
      expr_set_to_undefined(a);
      return true;
    }

    expr_set_to_inf(a);

    return true;
  }

  if (is_inf(b)) {

    if (is_neg_inf(a)) {
      expr_set_to_undefined(a);
      return true;
    }

    expr_set_to_inf(a);

    return true;
  }

  if (is(a, kind::UNDEF) || is(b, kind::UNDEF)) {
    expr_set_to_undefined(a);
    return true;
  }

  if (is(a, kind::FAIL) || is(b, kind::FAIL)) {
    expr_set_to_fail(a);
    return true;
  }

  if (is(a, kind::POW) && is(b, kind::SYM)) {
    return false;
  }

  if (is(a, kind::SYM) && is(b, kind::POW)) {
    return false;
  }

  int kind = kind_of(a) & kind_of(b);

  if (kind & (kind::SYM | kind::POW)) {
    if (compare(a, b, kind::ADD) == 0) {
      expr_set_op_to_mul(u, i, 2);
      return true;
    }

    return false;
  }

  long size_a = size_of(a);

  if (kind & kind::MUL) {
    long size_b = size_of(b);

    expr *c;

    long size_c = -1;

    expr *f = operand(a, 0);
    expr *g = operand(b, 0);

    if (is(f, kind::CONST) && is(g, kind::CONST)) {
      if (size_of(a) != size_of(b)) {
        return false;
      }
    }

    else if (is(f, kind::CONST)) {
      if (size_of(a) <= size_of(b)) {
        return false;
      }
    }

    else if (is(g, kind::CONST)) {
      if (size_of(b) <= size_of(a)) {
        return false;
      }
    }

    if (size_b > size_a) {
      c = a;
      a = b;
      b = c;

      size_c = size_b;
      size_b = size_a;
      size_a = size_c;
    }

    if (!(size_a == size_b || size_a == size_b + 1)) {
      return false;
    }

    long size = size_b - (is(operand(b, 0), kind::CONST) ? 1 : 0);

    for (long x = 0; x < size; x++) {
      if (compare(operand(a, size_a - x - 1), operand(b, size_b - x - 1),
                  kind::ADD) != 0) {
        return false;
      }
    }

    int ka = kind_of(operand(a, 0));
    int kb = kind_of(operand(b, 0));

    if ((ka & kind::CONST) && (kb & kind::CONST)) {
      expr_set_op_inplace_add_consts(a, 0, operand(b, 0));

      if (is(operand(a, 0), kind::INT) && get_val(operand(a, 0)) == 0) {
        expr_set_to_int(a, 0);
      }
    } else if (ka & kind::CONST) {
      expr_set_op_inplace_add_consts(a, 0, 1);

      if (is(operand(a, 0), kind::INT) && get_val(operand(a, 0)) == 0) {
        expr_set_to_int(a, 0);
      }
    } else {
      expr_set_to_mul(2, a);
    }

    // expr_set_operand(u, a, i);

    if (size_c != -1) {
      u->expr_childs[i] = *a;
      v->expr_childs[j] = *b;
    }

    return true;
  }

  // if (is(b, kind::MUL)) {
  //   return false;
  // }

  if (is(a, kind::MUL)) {
    if (!is(b, kind::SYM | kind::POW)) {
      return false;
    }

    if (size_of(a) > 2) {
      return false;
    }

    expr *a0 = operand(a, 0);
    expr *a1 = operand(a, 1);

    long ki = kind_of(a1) & kind_of(b);

    if (!is(a0, kind::CONST) || !ki) {
      return false;
    }

    if (compare(b, a1, kind::ADD) == 0) {
      expr_set_op_inplace_add_consts(a, 0, 1);

      if (is(operand(a, 0), kind::INT) && get_val(operand(a, 0)) == 0) {
        expr_set_to_int(a, 0);
      }

      return true;
    }
  }

  if (is(b, kind::MUL)) {
    if (!is(a, kind::SYM | kind::POW)) {
      return false;
    }

    if (size_of(b) > 2) {
      return false;
    }

    expr *b0 = operand(b, 0);
    expr *b1 = operand(b, 1);

    long ki = kind_of(b1) & kind_of(a);

    if (!is(b0, kind::CONST) || !ki) {
      return false;
    }

    if (compare(a, b1, kind::ADD) == 0) {
      expr_replace_with(a, b);

      expr_set_op_inplace_add_consts(a, 0, 1);

      if (is(operand(a, 0), kind::INT) && get_val(operand(a, 0)) == 0) {
        expr_set_to_int(a, 0);
      }

      return true;
    }
  }

  return false;
}

inline bool eval_mul_nconst(expr *u, size_t i, expr *v, size_t j) {
  assert(is(u, kind::MUL) && is(v, kind::MUL));
  assert(is(operand(u, i), kind::MULTIPLICABLE));
  assert(is(operand(v, j), kind::MULTIPLICABLE));

  expr *a = operand(u, i);
  expr *b = operand(v, j);

  if (a == 0 || b == 0) {
    return false;
  }

  if (is_neg_inf(a)) {
    if (is(b, kind::POW)) {
      if (is_inv_inf(b)) {
        expr_set_to_undefined(a);
      }

      expr_set_to_inf(a);
      expr_set_to_neg_inf(a);

      return true;
    }

    if (is_neg_inf(b)) {
      expr_set_to_inf(a);
      return true;
    }

    if (is(b, kind::CONST) && is_negative(b)) {
      expr_set_to_inf(a);
      return true;
    }

    expr_set_to_neg_inf(a);

    return true;
  }

  if (is_neg_inf(b)) {
    if (is(a, kind::POW)) {
      if (is_inv_inf(a)) {
        expr_set_to_undefined(a);
      }

      expr_set_to_neg_inf(a);

      return true;
    }

    if (is_neg_inf(a)) {
      expr_set_to_inf(a);
      return true;
    }

    if (is(a, kind::CONST) && is_negative(a)) {
      expr_set_to_inf(a);
      return true;
    }

    expr_set_to_neg_inf(a);

    return true;
  }

  if (is_inf(a)) {
    if (is(b, kind::POW)) {
      if (is_inv_inf(b)) {
        expr_set_to_undefined(a);
      }

      expr_set_to_inf(a);

      return true;
    }

    if (is_neg_inf(b)) {
      expr_set_to_undefined(a);
      return true;
    }

    if (is(b, kind::CONST) && is_negative(b)) {
      expr_set_to_neg_inf(a);
      return true;
    }

    expr_set_to_inf(a);

    return true;
  }

  if (is_inf(b)) {
    if (is(a, kind::POW)) {
      if (is_inv_inf(a)) {
        expr_set_to_undefined(a);
      }

      expr_set_to_inf(a);

      return true;
    }

    if (is_neg_inf(a)) {
      expr_set_to_undefined(a);
      return true;
    }

    if (is(a, kind::CONST) && is_negative(a)) {
      expr_set_to_neg_inf(a);
      return true;
    }

    expr_set_to_inf(a);

    return true;
  }

  if (is(a, kind::UNDEF) || is(b, kind::UNDEF)) {
    expr_set_to_undefined(a);
    return true;
  }

  if (is(a, kind::FAIL) || is(b, kind::FAIL)) {
    expr_set_to_fail(a);
    return true;
  }

  if (is(a, kind::ADD) && is(b, kind::ADD)) {

    if (compare(a, b, kind::MUL) == 0) {
      expr_set_op_to_pow(u, i, 2);
      return true;
    }

    return false;
  }

  if (is(a, kind::ADD) && is(b, kind::POW)) {
    if (compare(a, operand(b, 0), kind::MUL) == 0) {
      u->expr_childs[i] = create(kind::ADD, {integer(1), *operand(b, 1)});

      reduce(&u->expr_childs[i]);

      return true;
    }

    return false;
  }

  if (is(a, kind::POW) && is(b, kind::ADD)) {
    if (compare(b, operand(a, 0), kind::MUL) == 0) {
      expr_set_op_pow_add_to_deg(u, i, 1);
      reduce(&u->expr_childs[i]);

      return true;
    }

    return false;
  }

  if (is(a, kind::POW) && is(b, kind::POW)) {
    if (compare(operand(a, 0), operand(b, 0), kind::MUL) == 0) {
      expr_set_op_pow_add_to_deg(u, i, operand(b, 1));
      reduce(operand(u, i));
      return true;
    }

    return 0;
  }

  if (is(a, kind::SYM | kind::FUNC) && is(b, kind::SYM | kind::FUNC)) {
    if (compare(a, b, kind::MUL) == 0) {
      expr_set_op_to_pow(u, i, 2);
      return true;
    }

    return false;
  }

  if (is(a, kind::POW) && is(b, kind::SYM | kind::FUNC)) {
    if (compare(operand(a, 0), b, kind::MUL) == 0) {

      expr_set_op_pow_add_to_deg(u, i, 1);
      reduce(operand(u, i));

      return true;
    }

    return false;
  }

  return false;
}

inline bool eval_add_add(expr *a, expr *b) {
  assert(is(a, kind::ADD));
  assert(is(b, kind::ADD));

  size_t j = 0;
  size_t i = 0;

  expr *u = 0;
  expr *v = 0;
	// printf("\n=====> %s\n", to_string(a).c_str());
	// printf("=====> %s\n", to_string(b).c_str());
  while (i < size_of(a) && j < size_of(b)) {
    assert(!is(operand(b, j), kind::ADD));

    u = operand(a, i);
    v = operand(b, j);

    if ((is(u, kind::INT) && get_val(u) == 0)) {
      i++;
      continue;
    }

    if ((is(v, kind::INT) && get_val(v) == 0)) {
      j++;
      continue;
    }

    bool reduced = false;

    if (is(u, kind::CONST) && is(v, kind::CONST)) {
      reduced = eval_add_consts(a, i, b, j);
    } else if (is(u, kind::SUMMABLE) && is(v, kind::SUMMABLE)) {
			// printf("reduced %s and  %s ?", to_string(u).c_str(), to_string(v).c_str());
			reduced = eval_add_nconst(a, i, b, j);
		  // printf("%i\n", reduced);
    }

    if (reduced) {
      if (is(operand(a, i), kind::INT) && get_val(operand(a, i)) == 0) {
        a->remove(i);
      } else {
        i = i + 1;
      }

      j = j + 1;
    } else {
      int order = compare(u, v, kind::ADD);
      // if (!is(u, kind::ADD) && !is(v, kind::ADD)) {
      //   printf("order %s and  %s ? = %i\n", to_string(u).c_str(), to_string(v).c_str(), order);
      // }
      if (order < 0) {
        i = i + 1;
      } else {
        a->insert(*v, i++);
        j = j + 1;
      }
    }
  }
	// printf("======> %s\n", to_string(a).c_str());
	// printf("======> ");
	// for(size_t t = j; t < size_of(b); t++) {
	// 	// printf("%s", to_string(operand(b, t)).c_str());
	// 	if(t < size_of(b) - 1) {
	// 		printf(" + ");
	// 	}
	// }
	// printf("\n");

  while (j < size_of(b)) {
    if (i >= size_of(a)) {
      v = operand(b, j++);
      a->insert(*v, size_of(a));
    } else {
      u = operand(a, i);
      v = operand(b, j);

      if ((is(u, kind::INT) && get_val(u) == 0)) {
        i++;
        continue;
      }

      if ((is(v, kind::INT) && get_val(v) == 0)) {
        j++;
        continue;
      }

      if (v == 0) {
        j++;
      } else if (u == 0 || compare(u, v, kind::ADD) < 0) {
        i++;
      } else {
        a->insert(*v, i++);
        j = j + 1;
      }
    }
  }

  // printf("--===> %s\n", to_string(a).c_str());
  return true;
}

inline bool eval_mul_mul(expr *a, expr *b) {
  assert(is(a, kind::MUL));
  assert(is(b, kind::MUL));

  size_t j = 0;
  size_t i = 0;

  expr *u = 0;
  expr *v = 0;

  while (i < size_of(a) && j < size_of(b)) {
    assert(!is(operand(b, j), kind::MUL));

    u = operand(a, i);
    v = operand(b, j);

    if (u == 0) {
      i++;
      continue;
    }

    if (v == 0) {
      j++;
      continue;
    }

    bool reduced = false;

    if (is(u, kind::CONST) && is(v, kind::CONST)) {
      reduced = eval_mul_consts(a, i, b, j);
    } else if (is(u, kind::MULTIPLICABLE) && is(v, kind::MULTIPLICABLE)) {
      reduced = eval_mul_nconst(a, i, b, j);
    }

    if (reduced) {
      i = i + 1;
      j = j + 1;
    } else {
      int order = compare(u, v, kind::MUL);

      if (order < 0) {
        i = i + 1;
      } else {
        a->insert(*v, i++);
        j = j + 1;
      }
    }
  }

  while (j < size_of(b)) {
    if (i >= size_of(a)) {
      v = operand(b, j++);
      a->insert(*v, size_of(a));
    } else {
      u = operand(a, i);
      v = operand(b, j);

      if (v == 0) {
        j++;
      } else if (u == 0 || compare(u, v, kind::MUL) < 0) {
        i++;
      } else {
        a->insert(*v, i++);
        j = j + 1;
      }
    }
  }

  return a;
}

inline bool eval_mul_int(expr *u, size_t i, Int v) {
  expr *a = operand(u, i);

  if (is(a, kind::INT)) {
    expr_set_op_inplace_mul_consts(u, i, v);
    return true;
  }

  if (is(a, kind::ADD | kind::SUB)) {
    for (size_t j = 0; j < size_of(a); j++) {
      eval_mul_int(a, j, v);
    }

    return true;
  }

  if (is(a, kind::MUL)) {
    expr_set_op_to_mul(u, i, v);
    return true;
  }

  if (is(a, kind::DIV)) {
    eval_mul_int(a, 0, v);
    return true;
  }

  if (is(a, kind::SQRT | kind::POW | kind::FACT | kind::FUNC | kind::SYM)) {
    expr_set_op_to_mul(u, i, v);
    return true;
  }

  if (is(a, kind::CONST)) {
    expr_set_op_inplace_mul_consts(u, i, v);
    return true;
  }

  return true;
}

inline bool expr_raise_to_first_op(expr *a) {
  if (is(operand(a, 0), kind::INT)) {
    expr_set_to_int(a, get_val(operand(a, 0)));
    return true;
  }

  if (is(operand(a, 0), kind::SYM)) {
    expr_set_to_sym(a, get_id(operand(a, 0)));
    return true;
  }

  if (is(operand(a, 0), kind::FUNC)) {
    a->expr_sym = strdup(get_id(operand(a, 0)));
    expr_set_kind(a, kind::FUNC);
    a->expr_childs = expr(a->expr_childs[0]).expr_childs;
    return true;
  }

  *a = expr(a->expr_childs[0]);

  return true;
}

void reduce_add(expr *a) {
  for (size_t i = 0; i < size_of(a); i++) {
    reduce(operand(a, i));
  }

	// printf("--> %s\n", to_string(a).c_str());

	sort_childs(a, 0, size_of(a) - 1);

	// printf("--> %s\n", to_string(a).c_str());

  if (is(operand(a, 0), kind::ADD)) {
    expr t = a->expr_childs[0];

    a->remove(0);

    eval_add_add(a, &t);
  }

  size_t j = 0;

  for (long i = 1; i < (long)size_of(a); i++) {
    if (j > size_of(a)) {
      break;
    }

    expr *aj = operand(a, j);
    expr *ai = operand(a, i);

    bool reduced = 0;

    if (is(ai, kind::FAIL) || is(aj, kind::FAIL)) {
      return expr_set_to_fail(a);
    }

    if (is(ai, kind::UNDEF) || is(aj, kind::UNDEF)) {
      return expr_set_to_undefined(a);
    }

    if (is(ai, kind::INT) && get_val(ai) == 0) {
      a->remove(i--);
      continue;
    }

    if (is(aj, kind::INT) && get_val(aj) == 0) {
      a->remove(j);
      continue;
    }

    if (is(ai, kind::ADD)) {
      expr t = a->expr_childs[i];
      a->remove(i--);
      eval_add_add(a, &t);
    }

    else if (is_inf(ai) && !is_neg_inf(aj)) {
      expr_set_to_inf(aj);
      a->remove(i--);
      reduced = true;
    }

    else if (is_neg_inf(aj) && !is_inf(ai)) {
      a->remove(i--);
      reduced = true;
    }

    else if (is_neg_inf(ai) && !is_inf(aj)) {
      expr_set_to_neg_inf(aj);
      a->remove(i--);
      reduced = true;
    }

    else if (is_inf(aj) && !is_neg_inf(ai)) {
      a->remove(i--);
      reduced = true;
    }

    else if (is(aj, kind::CONST) && is(ai, kind::CONST)) {
      reduced = eval_add_consts(a, j, a, i);

      if (reduced) {
        a->remove(i--);
      }
    }

    else if (is(aj, kind::SUMMABLE) && is(ai, kind::SUMMABLE)) {
      reduced = eval_add_nconst(a, j, a, i);

      if (reduced) {
        a->remove(i--);
      }
    }

    if (reduced && is(operand(a, j), kind::INT) &&
        get_val(operand(a, j)) == 0) {
      a->remove(j);
    }

    if (reduced == false) {
      j = i;
    }
  }

  if (size_of(a) == 0) {
    expr_set_to_int(a, 0);
  } else if (size_of(a) == 1) {
    expr_raise_to_first_op(a);
  }
	// printf("<-- %s\n", to_string(a).c_str());
}

void reduce_mul(expr *a) {
  for (size_t i = 0; i < size_of(a); i++) {
    reduce(operand(a, i));
  }

  sort_childs(a, 0, size_of(a) - 1);

  size_t j = 0;

  if (is(operand(a, 0), kind::MUL) && !is_neg_inf(operand(a, 0))) {
    expr t = a->expr_childs[0];

    a->remove(0);

    eval_mul_mul(a, &t);
  }

  for (long i = 1; i < (long)size_of(a); i++) {
    expr *aj = operand(a, j);
    expr *ai = operand(a, i);

    bool reduced = 0;

    if (is(ai, kind::INT) && get_val(ai) == 0) {
      return expr_set_to_int(a, 0);
    }

    if (is(aj, kind::INT) && get_val(aj) == 0) {
      return expr_set_to_int(a, 0);
    }

    if (is_inf(ai) && is(aj, kind::CONST)) {
      if (is_negative(aj)) {
        expr_set_to_neg_inf(ai);

        a->remove(j);
      } else {
        a->remove(j);
      }
    }

    else if (is_inf(aj) && is(ai, kind::CONST)) {
      if (is_negative(ai)) {
        expr_set_to_neg_inf(ai);

        a->remove(j);
      } else {
        a->remove(j);
      }
    }

    else if (is_neg_inf(ai) && is(aj, kind::CONST)) {
      if (is_negative(aj)) {
        expr_set_to_inf(ai);

        a->remove(j);
      } else {
        a->remove(j);
      }
    }

    else if (is_neg_inf(aj) && is(ai, kind::CONST)) {
      if (is_negative(ai)) {
        expr_set_to_inf(ai);

        a->remove(j);
      } else {
        a->remove(j);
      }
    }

    else if (is(ai, kind::INF) && is(aj, kind::INF)) {
      a->remove(j);
    }

    else if (is(ai, kind::FAIL) || is(aj, kind::FAIL)) {
      return expr_set_to_fail(a);
    }

    else if (is(ai, kind::UNDEF) || is(aj, kind::UNDEF)) {
      return expr_set_to_undefined(a);
    }

    else if (expr_is_zero(operand(a, i)) || expr_is_zero(operand(a, j))) {
      return expr_set_to_int(a, 0);
    }

    else if (is(ai, kind::MUL)) {
      expr t = a->expr_childs[i];
      a->remove(i--);
      eval_mul_mul(a, &t);
    }

    else if (is(aj, kind::CONST) && is(ai, kind::CONST)) {

      reduced = eval_mul_consts(a, j, a, i);

      if (reduced) {
        a->remove(i--);
      }

    }

    else if (is(aj, kind::INT) && get_val(aj) == 1) {
      a->remove(j);
      i = i - 1;
    }

    else if (is(ai, kind::INT) && get_val(ai) == 1) {
      a->remove(i--);
    }

    else if (is(aj, kind::MULTIPLICABLE) && is(ai, kind::MULTIPLICABLE)) {
      reduced = eval_mul_nconst(a, j, a, i);

      if (reduced) {
        a->remove(i--);
      }
    }

    if (reduced == false) {
      j = i;
    }
  }

  if (is(a, kind::MUL) && size_of(a) == 1) {
    expr_raise_to_first_op(a);
  }
}

void reduce_sub(expr *a) {
  for (size_t i = 1; i < size_of(a); i++) {
    eval_mul_int(a, i, -1);
  }

  expr_set_kind(a, kind::ADD);

  reduce(a);
}

void reduce_pow(expr *a) {
  reduce(operand(a, 1));
  reduce(operand(a, 0));

  if (is(operand(a, 0), kind::UNDEF) || is(operand(a, 1), kind::UNDEF)) {
    return expr_set_to_undefined(a);
  }

  if (is(operand(a, 0), kind::FAIL) || is(operand(a, 1), kind::FAIL)) {
    return expr_set_to_undefined(a);
  }

  if (is(operand(a, 1), kind::INT) && get_val(operand(a, 1)) == 0) {
    if (is(operand(a, 0), kind::INT) && get_val(operand(a, 0)) == 0) {
      return expr_set_to_undefined(a);
    }

    return expr_set_to_int(a, 1);
  }

  if (is(operand(a, 0), kind::POW)) {
    expr_set_op_to_mul(operand(a, 0), 1, operand(a, 1));

    expr_raise_to_first_op(a);

    return reduce(a);
  }

  if ((is_inf(operand(a, 0)) || is_neg_inf(operand(a, 0))) &&
      (is_inf(operand(a, 1)) || is_neg_inf(operand(a, 1)))) {
    return expr_set_to_undefined(a);
  }

  if (is(operand(a, 0), kind::MUL)) {
    for (size_t i = 0; i < size_of(operand(a, 0)); i++) {
      expr_set_op_to_pow(operand(a, 0), i, operand(a, 1));
    }

    expr_raise_to_first_op(a);

    return reduce(a);
  }

  if (!is(operand(a, 1), kind::INT)) {
    return;
  }

  if (get_val(operand(a, 1)) == 1) {
    expr_raise_to_first_op(a);
    return;
  }

  if (is(operand(a, 0), kind::INT)) {
    Int b = get_val(operand(a, 0));
    Int c = get_val(operand(a, 1));

    bool n = c < 0;

    c = abs(c);

    Int d = pow(b, c);

    a->remove(1);

    if (!n || d == 1) {
      expr_set_to_int(a, d);
    } else {
      expr_set_to_fra(a, 1, d);
    }

    return;
  }

  if (is(operand(a, 0), kind::FRAC)) {
    Int b = get_val(operand(operand(a, 0), 0));
    Int c = get_val(operand(operand(a, 0), 1));

    Int d = get_val(operand(a, 1));

    bool n = d < 0;

    d = abs(d);

    b = pow(b, d);
    c = pow(c, d);

    Int g = gcd(b, c);

    b = b / g;
    c = c / g;

    if (!n) {
      expr_set_to_fra(a, b, c);
    } else {
      expr_set_to_fra(a, c, b);
    }

    return;
  }

  if (is(operand(a, 0), kind::TERMINAL)) {
    return;
  }

  if (is(operand(a, 0), kind::MUL)) {
    long long y = get_val(operand(a, 1)).longValue();

    expr b = a->expr_childs[0];

    expr_set_to_int(a, 1);

    while (y) {
      if (y % 2 == 1) {
        expr_set_to_mul(a, &b);
        reduce(a);
      }

      y = y >> 1;

      expr t = create(kind::MUL);

      for (size_t i = 0; i < size_of(&b); i++) {
        expr v = create(kind::MUL, {
                                       *operand(&b, i),
                                       *operand(&b, i),
                                   });

        t.insert(v);
      }

      expr_replace_with(&b, &t);
    }

    return;
  }
}

void reduce_div(expr *a) {
  if ((is_inf(operand(a, 0)) || is_neg_inf(operand(a, 0))) &&
      (is_inf(operand(a, 1)) || is_neg_inf(operand(a, 1)))) {
    return expr_set_to_undefined(a);
  }

  if (is_inf(operand(a, 0))) {
    return expr_set_to_inf(a);
  }

  if (is_neg_inf(operand(a, 0))) {
    return expr_set_to_neg_inf(a);
  }

  if (is(operand(a, 1), kind::INT) && get_val(operand(a, 1)) == 0) {
    return expr_set_to_undefined(a);
  }

  if (is(operand(a, 0), kind::INT) && get_val(operand(a, 0)) == 0) {
    return expr_set_to_int(a, 0);
  }

  expr_set_kind(a, kind::MUL);

  expr_set_op_to_pow(a, 1, -1);

  return reduce(a);
}

void reduce_sqr(expr *a) {
  expr_set_kind(a, kind::POW);

  a->insert(fraction(1, 2), 1);

  return reduce(a);
}

void reduce_fac(expr *a) {
  reduce(operand(a, 0));

  if (is(operand(a, 0), kind::FAIL)) {
    return expr_set_to_fail(a);
  }

  if (is(operand(a, 0), kind::UNDEF)) {
    return expr_set_to_undefined(a);
  }

  if (is_inf(operand(a, 0))) {
    return expr_set_to_inf(a);
  }

  if (is_neg_inf(operand(a, 0))) {
    return expr_set_to_neg_inf(a);
  }

  if (is(operand(a, 0), kind::INT)) {
    return expr_set_to_int(a, fact(get_val(operand(a, 0))));
  }

  if (is(operand(a, 0), kind::FRAC)) {
    Int c = fact(get_val(operand(operand(a, 0), 0)));
    Int d = fact(get_val(operand(operand(a, 0), 1)));

    Int g = abs(gcd(c, d));

    return expr_set_to_fra(a, c / g, d / g);
  }
}

expr *reduce_fra(expr *a) {
  Int b = get_val(operand(a, 0));
  Int c = get_val(operand(a, 1));

	if (c == 0) {
    expr_set_to_undefined(a);
    return a;
  }

	if (b % c == 0) {
		expr_set_to_int(a, b / c);
    return a;
	}

  Int d = abs(gcd(b, c));

	if (c / d == 1) {
    expr_set_to_int(a, b / d);
  } else {
    expr_set_to_fra(a, b / d, c / d);
  }
  return a;
}

void reduce(expr *a) {
  if (is(a, kind::LIST)) {
    for (size_t i = 0; i < size_of(a); i++) {
      reduce(&a->expr_list->members[i]);
    }
  } else if (is(a, kind::SET)) {
    for (size_t i = 0; i < size_of(a); i++) {
      reduce(&a->expr_set->members[i]);
    }

    trim(a->expr_set);
  } else if (is(a, kind::FRAC)) {
    reduce_fra(a);
  } else if (is(a, kind::ADD)) {
    reduce_add(a);
  } else if (is(a, kind::MUL)) {
    reduce_mul(a);
  } else if (is(a, kind::SUB)) {
    reduce_sub(a);
  } else if (is(a, kind::DIV)) {
    reduce_div(a);
  } else if (is(a, kind::POW)) {
    reduce_pow(a);
  } else if (is(a, kind::SQRT)) {
    reduce_sqr(a);
  } else if (is(a, kind::FACT)) {
    reduce_fac(a);
  }
}

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

bool expand_pow(expr *u, Int n, expr *a) {
  if (is(u, kind::TERMINAL)) {
    *a = create(kind::POW, {*u, integer(n)});
    return true;
  }

  if (n == 1) {
    *a = *u;
    return true;
  }

  if (n == 0) {
    *a = integer(1);
    return true;
  }

  if (is(u, kind::ADD)) {
    Int c = fact(n);

    expr o = *u;

    expr f = o[0];

    o.remove(0);

    if (size_of(&o) == 0) {
      expr_set_to_int(&o, 0);
    }

    if (size_of(&o) == 1) {
      expr_raise_to_first_op(&o);
    }

    expr s = create(kind::ADD);

    for (Int k = 0; k <= n; k++) {
      expr z = create(kind::MUL, {integer(c / (fact(k) * fact(n - k))),
                                  create(kind::POW, {f, integer(n - k)})});

      expr t;

      if (expand_pow(&o, k, &t)) {
        t = create(kind::MUL, {z, t});
      } else {
        t = z;
      }

      s.insert(t);
    }

    *a = s;

    reduce(a);

    return true;
  }

  return false;
}

int tabs = 0;

void expand(expr *a) {
  if (is(a, kind::TERMINAL)) {
    return;
  }

  if (is(a, kind::SUB | kind::DIV | kind::FACT | kind::FRAC | kind::SQRT)) {
    reduce(a);
  }

  if (is(a, kind::POW)) {
    expand(operand(a, 0));
    expand(operand(a, 1));

    if (is(operand(a, 1), kind::INT)) {
      expand_pow(operand(a, 0), get_val(operand(a, 1)), a);
    }

    reduce(a);
  }

  if (is(a, kind::MUL)) {
    tabs += 3;
    // printf("%*cfrom ----> %s\n",tabs,' ',  to_string(a).c_str());
    while (size_of(a) > 1) {
      expand(operand(a, 0));
      expand(operand(a, 1));

      expr t = expand_mul(a, 0, a, 1);

      a->insert(t, 0);

      a->remove(1);
      a->remove(1);
    }

    expr_raise_to_first_op(a);
    // printf("%*cto ----> %s\n",tabs,' ',  to_string(a).c_str());
    reduce(a);
    // printf("%*cto ----> %s\n",tabs,' ',  to_string(a).c_str());
    tabs -= 3;
  }

  if (is(a, kind::ADD)) {
    for (size_t i = 0; i < size_of(a); i++) {
      expand(operand(a, i));
    }

    reduce(a);
  }
  // printf("from %s\n", to_string(a).c_str());
  // printf("to %s\n", to_string(a).c_str());
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

  expr a = other;
  expr b = *this;

  // printf("a = %s\n", to_string(a).c_str());
  // printf("b = %s\n", to_string(b).c_str());

  sort(&a, kind::UNDEF);

  sort(&b, kind::UNDEF);
  // printf("a = %s\n", to_string(a).c_str());
  // printf("b = %s\n", to_string(b).c_str());

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

  expr b = *this;

  // printf("a = %s\n", to_string(a).c_str());

  sort(&a, kind::UNDEF);

  // printf("a = %s\n", to_string(a).c_str());
  // printf("b = %s\n", to_string(b).c_str());

  sort(&b, kind::UNDEF);
  // printf("b = %s\n", to_string(b).c_str());

  // printf("cmp = %i\n", compare(&a, &b, kind::UNDEF));
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

expr sqrt(const expr &a) { return create(kind::SQRT, {a}); }

expr sqrt(expr &&a) { return create(kind::SQRT, {a}); }

expr fact(const expr &a) { return create(kind::FACT, {a}); }

expr fact(expr &&a) { return create(kind::FACT, {a}); }

expr undefined() { return create(kind::UNDEF); }

expr fail() { return create(kind::FAIL); }

expr inf() { return create(kind::INF); }

expr reduce(expr &a) {
  expr b = a;

  reduce(&b);

  return b;
}

expr reduce(expr &&a) {
  expr b = a;

  reduce(&b);

  return b;
}

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

bool exists(const expr &a, expr &&b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.expr_set;
  set M = *b.expr_set;

  return exists(L, M);
}

bool exists(const expr &a, const expr &b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.expr_set;
  set M = *b.expr_set;

  return exists(L, M);
}

void replace_rec(expr *a, expr *b, expr *c) {
  if (a->match(b)) {
    expr_replace_with(a, c);
    return;
  }

  if (is(a, kind::TERMINAL))
    return;

  for (size_t i = 0; i < size_of(a); i++) {
    replace_rec(operand(a, i), b, c);
  }
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

void list::append(expr &&a) { members.push_back(a); }

void list::append(const expr &a) { members.push_back(a); }
void list::insert(const expr &a, size_t idx) {
  members.insert(members.begin() + idx, a);
}
void list::insert(expr &&a, size_t idx) {
  members.insert(members.begin() + idx, a);
}

list append(list &a, const expr &b) {
  list L = a;
  L.append(b);
  return L;
}

list append(list &a, expr &&b) {
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

  return true;
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

    if (cmp > 0) {
      t.members.push_back(M[j++]);
    }

    if (cmp < 0) {
      t.members.push_back(L[i++]);
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

bool exists(set &L, expr &e) {
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

expr sinh(expr x) { return func_call("sinh", {x}); }

expr cosh(expr x) { return func_call("cosh", {x}); }

expr tanh(expr x) { return func_call("tanh", {x}); }

expr exp(expr x) { return func_call("exp", {x}); }

expr cos(expr x) { return func_call("cos", {x}); }

expr sin(expr x) { return func_call("sin", {x}); }

expr tan(expr x) { return func_call("tan", {x}); }

expr csc(expr x) { return func_call("csc", {x}); }

expr cot(expr x) { return func_call("cot", {x}); }

expr log(expr x) { return func_call("log", {x}); }

expr ln(expr x) { return func_call("ln", {x}); }

expr sec(expr x) { return func_call("sec", {x}); }

expr coth(expr x) { return func_call("coth", {x}); }

expr sech(expr x) { return func_call("sech", {x}); }

expr csch(expr x) { return func_call("csch", {x}); }

expr arccos(expr x) { return func_call("arccos", {x}); }

expr arcsin(expr x) { return func_call("arcsin", {x}); }

expr arctan(expr x) { return func_call("arctan", {x}); }

expr arccot(expr x) { return func_call("arccot", {x}); }

expr arcsec(expr x) { return func_call("arcsec", {x}); }

expr arccsc(expr x) { return func_call("arccsc", {x}); }

expr arccosh(expr x) { return func_call("arccosh", {x}); }

expr arctanh(expr x) { return func_call("arctanh", {x}); }

} // namespace alg
