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

namespace ast_teste {

ast::ast(ast &&other) {
  //	*this = std::move(other);
  kind_of = other.kind_of;

  switch (kind_of) {
  case kind::SYM: {
    ast_sym = other.ast_sym;
    other.ast_sym = 0;
    return;
  }
  case kind::INT: {
    ast_int = other.ast_int;
    other.ast_int = 0;
    return;
  }

  case kind::INF:
    return;
  case kind::NEG_INF:
    return;
  case kind::UNDEF:
    return;
  case kind::FAIL:
    return;

  default: {
    ast_childs = std::move(other.ast_childs);
    return;
  }
  }
}

ast::ast(const ast &other) {
  kind_of = other.kind_of;

  switch (kind_of) {
  case kind::SYM: {
    ast_sym = strdup(other.ast_sym);
    return;
  }
  case kind::INT: {
    ast_int = new Int(*other.ast_int);
    return;
  }
  case kind::INF:
    return;
  case kind::NEG_INF:
    return;
  case kind::UNDEF:
    return;
  case kind::FAIL:
    return;

  default: {
    ast_childs = other.ast_childs;
    return;
  }
  }
}
ast &ast::operator=(const ast &other) {
  kind_of = other.kind_of;

  switch (kind_of) {

  case kind::SYM: {
    ast_sym = strdup(other.ast_sym);
    ast_childs.clear();
    return *this;
  }

  case kind::INT: {
    ast_int = new Int(*other.ast_int);
    ast_childs.clear();
    return *this;
  }

  case kind::INF:
    return *this;
  case kind::NEG_INF:
    return *this;
  case kind::UNDEF:
    return *this;
  case kind::FAIL:
    return *this;

  default: {
    ast_childs = other.ast_childs;
    return *this;
  }
  }

  return *this;
}

ast &ast::operator=(ast &&other) {
  kind_of = other.kind_of;

  switch (kind_of) {

  case kind::SYM: {
    ast_sym = other.ast_sym;
    other.ast_sym = 0;
    return *this;
  }

  case kind::INT: {
    ast_int = other.ast_int;
    other.ast_int = 0;
    return *this;
  }

  case kind::INF:
    return *this;
  case kind::NEG_INF:
    return *this;
  case kind::UNDEF:
    return *this;
  case kind::FAIL:
    return *this;

  default: {
    ast_childs = std::move(other.ast_childs);
    return *this;
  }
  }

  return *this;
}

ast::ast(enum kind k) { kind_of = k; }

ast::ast(list &s) {
  kind_of = kind::LIST;
  ast_list = new list(s);
}

ast::ast(list &&s) {
  kind_of = kind::LIST;
  ast_list = new list(s);
}

ast::ast(set &s) {
  kind_of = kind::SET;
  ast_set = new set(s);
}

ast::ast(set &&s) {
  kind_of = kind::SET;
  ast_set = new set(s);
}

ast::ast(enum kind k, std::initializer_list<ast> &&a) {
  if (k == kind::LIST) {
    this->ast_list = new list(a);
    return;
  }

  if (k == kind::SET) {
    this->ast_set = new set(a);
    return;
  }

  kind_of = k;

  ast_childs = std::move(a);
}

ast::ast(Int v) {
  kind_of = kind::INT;
  this->ast_int = new Int(v);
}

ast::ast(int v) {
  kind_of = kind::INT;
  this->ast_int = new Int(v);
}

ast::ast(long int v) {
  kind_of = kind::INT;
  this->ast_int = new Int(v);
}

ast::ast(long long v) {
  kind_of = kind::INT;
  this->ast_int = new Int(v);
}

ast::ast(std::string v) {
  kind_of = kind::SYM;
  this->ast_sym = strdup(v.c_str());
}

ast::ast() { kind_of = kind::UNDEF; }

ast::~ast() {
  switch (kind_of) {
  case kind::INT: {
    delete ast_int;
    return;
  }
  case kind::SYM: {
    delete ast_sym;
    return;
  }
  default:
    return;
  }
}

ast &ast::operator[](size_t idx) { return ast_childs[idx]; }

ast create(kind kind) { return ast(kind); }

ast create(kind kind, std::initializer_list<ast> &&l) {
  ast u(kind);

  u.ast_childs = std::vector<ast>(std::move(l));

  return u;
}

void ast_set_kind(ast *a, kind kind) { a->kind_of = kind; }

ast symbol(const char *id) {
  ast a = create(kind::SYM);

  a.ast_sym = strdup(id);

  return a;
}

ast integer(Int value) {
  ast a = create(kind::INT);

  a.ast_int = new Int(value);

  return a;
}

ast fraction(Int num, Int den) {
  return create(kind::FRAC, {integer(num), integer(den)});
}

void ast::insert(const ast &b, size_t idx) {
  this->ast_childs.insert(this->ast_childs.begin() + idx, b);
}

void ast::insert(ast &&b, size_t idx) {
  this->ast_childs.insert(this->ast_childs.begin() + idx, std::move(b));
}

void ast::insert(const ast &b) { this->ast_childs.push_back(b); }
void ast::insert(ast &&b) { this->ast_childs.push_back(std::move(b)); }

void ast::remove(list &l) {
  assert(is(this, kind::LIST));
  ast_list->remove(l);
}

void ast::remove(list &&l) {
  assert(is(this, kind::LIST));
  ast_list->remove(std::move(l));
}

void ast::remove(size_t idx) {
  if (is(this, kind::LIST)) {
    this->ast_list->members.erase(this->ast_childs.begin() + idx);
    return;
  }

  if (is(this, kind::SET)) {
    this->ast_set->members.erase(this->ast_childs.begin() + idx);
    return;
  }

  this->ast_childs.erase(this->ast_childs.begin() + idx);
}

void ast::remove() {
  if (is(this, kind::LIST)) {
    this->ast_list->members.pop_back();
    return;
  }

  if (is(this, kind::SET)) {
    this->ast_set->members.pop_back();
    return;
  }

  this->ast_childs.pop_back();
}

int compare_consts(ast *a, ast *b) {
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

inline bool ast_is_zero(ast *a) {
  return (a == nullptr) || (is(a, kind::INT) && get_val(a) == 0);
}

inline int ast_op_cmp(ast *a, ast *b, kind ctx) {
  long m = size_of(a);
  long n = size_of(b);

  long l = std::min(size_of(a), size_of(b));

  m = m - 1;
  n = n - 1;

  if (is(a, kind::CONST) && is(b, kind::CONST)) {
    return compare_consts(a, b);
  }

  if (ctx == kind::ADD) {

    if (is(a, kind::MUL) && is(b, kind::MUL)) {

      if (std::abs(m - n) > 1) {
        return n - m;
      }

      for (long i = 0; i < l; i++) {
        int order = kind_of(operand(a, m - i)) - kind_of(operand(b, n - i));

        if (order)
          return order;
      }

      for (long i = 0; i < l; i++) {
        int order = compare(operand(a, m - i), operand(b, n - i), ctx);

        if (order)
          return order;
      }
    }
  }

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

  return (ctx & kind::ADD) ? m - n : n - m;
}

inline int compare_idents(ast *a, ast *b) {
  return strcmp(get_id(a), get_id(b));
}

std::string to_string(ast *tree) {
  if (!tree)
    return "null";

  if (is(tree, kind::INT)) {
    return tree->ast_int->to_string();
  }

  if (is(tree, kind::SYM)) {
    return std::string(tree->ast_sym);
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

  if (is(tree, kind::NEG_INF)) {
    return "-inf";
  }

  if (is(tree, kind::FRAC)) {
    return to_string(operand(tree, 0)) + "/" + to_string(operand(tree, 1));
  }

  if (is(tree, kind::FRAC)) {
    return "sqrt(" + to_string(operand(tree, 0)) + ")";
  }

  if (is(tree, kind::FUNC)) {
    std::string r = std::string(get_func_id(tree)) + "(";

    if (size_of(tree) > 0) {
      for (size_t i = 0; i < size_of(tree) - 1; i++) {
        r += to_string(operand(tree, i));
        r += ",";
      }

      r += to_string(operand(tree, size_of(tree) - 1));
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
          is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL)) {

        r += "(";
      }

      r += to_string(operand(tree, i));

      if (operand(tree, i) &&
          is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL)) {
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

    for (size_t i = 0; i < size_of(tree) - 1; i++) {
      if (operand(tree, i) &&
          is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL)) {
        r += "(";
      }

      r += to_string(operand(tree, i));

      if (operand(tree, i) &&
          is(operand(tree, i), kind::SUB | kind::ADD | kind::MUL)) {
        r += ")";
      }

      if (i != size_of(tree) - 1) {
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

	if(is(tree, kind::LIST)) return to_string(tree->ast_list);
	if(is(tree, kind::SET)) return to_string(tree->ast_set);

  return "to_string_not_implemented";
}

std::string kind_of_id(ast *a) {
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
  case kind::NEG_INF: {
    return "negative infinity";
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

  default:
    return "";
  }
}

void ast_print(ast *a, int tabs) {
  printf("%*c<ast ", tabs, ' ');
  printf("address=\"%p\" ", a);
  printf("kind=\"%s\"", kind_of_id(a).c_str());

  if (kind_of(a) == kind::INT) {
    printf(" value=\"%s\"", get_val(a).to_string().c_str());
  }

  if (kind_of(a) == kind::SYM) {
    printf(" id=\"%s\"", get_id(a));
  }

  if (size_of(a)) {
    printf(">\n");
    // printf("\n%*c  childs: [\n", tabs, ' ');

    for (size_t i = 0; i < size_of(a); i++) {
      ast_print(operand(a, i), tabs + 3);
    }
    // printf("%*c  ];", tabs, ' ');
    printf("%*c</ast>\n", tabs, ' ');
  } else {
    printf(">\n");
  }
}

int compare(ast *const a, ast *const b, kind ctx) {
  if (a == b)
    return 0;

  if (ctx & kind::MUL) {
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

		if(is(a, kind::LIST) && is(b, kind::LIST)) {
			return a->ast_list->match(b->ast_list);
		}

		if(is(a, kind::SET) && is(b, kind::SET)) {
			return a->ast_set->match(b->ast_set);
		}

  }

  if (ctx & kind::ADD) {

    if (is(a, kind::CONST) && is(b, kind::CONST)) {
      return compare_consts(b, a);
    }

    if (is(a, kind::CONST)) {
      return +1;
    }

    if (is(b, kind::CONST)) {
      return -1;
    }

    if (is(a, kind::POW) && is(b, kind::POW)) {
      int i = compare(operand(a, 1), operand(b, 1), ctx);

      if (i != 0) {
        return i;
      }

      return compare(operand(a, 0), operand(b, 0), ctx);
    }

    if (is(a, kind::FUNC) && is(b, kind::FUNC)) {
      return strcmp(get_func_id(a), get_func_id(b));
    }

    if (is(a, kind::ADD) && is(b, kind::SYM)) {
      return +1;
    }

    if (is(a, kind::SYM) && is(b, kind::ADD)) {
      return -1;
    }

    if (is(a, kind::MUL) && is(b, kind::SYM)) {
      long k = size_of(a);

      if (k > 2)
        return -k;

      int order = compare(operand(a, size_of(a) - 1), b, ctx);

      if (order == 0) {
        return -1;
      }

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

      if (k > 2)
        return -k;

      int order = compare(operand(a, 0), b, ctx);

      if (order == 0)
        return -1;

      return k > 2 ? -1 : order;
    }

    if (is(b, kind::MUL) && is(a, kind::POW)) {
      long k = size_of(b);

      if (k > 2)
        return +k;

      int order = compare(a, operand(b, 0), ctx);

      if (order == 0)
        return +1;

      return k > 2 ? +1 : order;
    }

		if(is(a, kind::LIST) && is(b, kind::LIST)) {
			return a->ast_list->match(b->ast_list);
		}

		if(is(a, kind::SET) && is(b, kind::SET)) {
			return a->ast_set->match(b->ast_set);
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
    return ast_op_cmp(a, b, ctx);
  }

  if (is(a, kind::MUL) && is(b, kind::MUL)) {
    return ast_op_cmp(a, b, ctx);
  }

  if (is(a, kind::POW) && is(b, kind::POW)) {
    return compare(operand(a, 0), operand(b, 0), ctx) ||
           compare(operand(a, 1), operand(b, 1), ctx);
  }

  if (is(a, kind::DIV) && is(b, kind::DIV)) {
    return compare(operand(a, 0), operand(b, 0), ctx) ||
           compare(operand(a, 1), operand(b, 1), ctx);
  }

  return should_revert_idx(ctx) ? kind_of(a) - kind_of(b)
                                : kind_of(b) - kind_of(a);
}

long int sort_split(ast *a, long l, long r) {
  long int i = l - 1;

  ast *p = operand(a, r);

  for (long int j = l; j < r; j++) {
    if (compare(operand(a, j), p, kind_of(a)) < 0) {
      std::swap(a->ast_childs[++i], a->ast_childs[j]);
    }
  }

  std::swap(a->ast_childs[++i], a->ast_childs[r]);

  return i;
}

void sort_childs(ast *a, long int l, long int r) {
  if (l < r) {
    long int m = sort_split(a, l, r);

    sort_childs(a, l, m - 1);
    sort_childs(a, m + 1, r);
  }
}

void sort(ast *a) {
  if (is(a, kind::TERMINAL)) {
    return;
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

inline void ast_set_to_undefined(ast *a) { ast_set_kind(a, kind::UNDEF); }

inline void ast_set_to_fail(ast *a) { ast_set_kind(a, kind::FAIL); }

inline void ast_set_to_int(ast *a, Int v) {
  ast_set_kind(a, kind::INT);

  a->ast_int = new Int(v);
}

inline void ast_set_op_to_int(ast *a, size_t i, Int v) {
  ast_set_to_int(operand(a, i), v);
}

inline void ast_set_to_fra(ast *a, Int u, Int v) {
  a->ast_childs = std::vector<ast>();

  ast_set_kind(a, kind::FRAC);

  a->insert(integer(u));
  a->insert(integer(v));
}

inline void ast_set_op_to_fra(ast *a, size_t i, Int u, Int v) {
  ast_set_to_fra(operand(a, i), u, v);
}

inline void ast_set_to_sym(ast *a, const char *s) {
  ast_set_kind(a, kind::SYM);
  a->ast_sym = strdup(s);
}

inline void ast_set_op_to_sym(ast *a, size_t i, const char *s) {
  ast_set_to_sym(operand(a, i), s);
}

// a = a + b
inline void ast_set_inplace_add_consts(ast *a, ast *b) {
  assert(is(a, kind::CONST));
  assert(is(b, kind::CONST));

  if (is(a, kind::INT) && is(b, kind::INT)) {
    Int x = get_val(a);
    Int y = get_val(b);

    return ast_set_to_int(a, x + y);
  }

  if (is(a, kind::FRAC) && is(b, kind::FRAC)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));
    Int w = get_val(operand(b, 0));
    Int z = get_val(operand(b, 1));

    Int e = x * z + w * y;
    Int k = y * z;

    if (e % k == 0) {
      return ast_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));

    return ast_set_to_fra(a, e / g, k / g);
  }

  if (is(a, kind::FRAC) && is(b, kind::INT)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));
    Int w = get_val(b);

    Int e = x + w * y;
    Int k = y;

    if (e % k == 0) {
      return ast_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));
    return ast_set_to_fra(a, e / g, k / g);
  }

  if (is(b, kind::FRAC) && is(a, kind::INT)) {
    Int x = get_val(operand(b, 0));
    Int y = get_val(operand(b, 1));
    Int w = get_val(a);

    Int e = x + w * y;
    Int k = y;

    if (e % k == 0) {
      return ast_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));

    return ast_set_to_fra(a, e / g, k / g);
  }
}

inline void ast_set_inplace_add_consts(ast *a, Int b) {
  assert(is(a, kind::CONST));

  if (is(a, kind::INT)) {
    Int x = get_val(a);

    return ast_set_to_int(a, x + b);
  }

  if (is(a, kind::FRAC)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));

    Int e = x + b * y;

    if (e % y == 0) {
      return ast_set_to_int(a, e / y);
    }

    Int g = abs(gcd(e, y));
    return ast_set_to_fra(a, e / g, y / g);
  }
}

inline void ast_set_op_inplace_add_consts(ast *a, size_t i, ast *b) {
  ast_set_inplace_add_consts(operand(a, i), b);
}

inline void ast_set_op_inplace_add_consts(ast *a, size_t i, Int b) {
  ast_set_inplace_add_consts(operand(a, i), b);
}

inline void ast_set_inplace_mul_consts(ast *a, ast *b) {
  assert(is(a, kind::CONST));
  assert(is(b, kind::CONST));

  if (is(a, kind::INT) && is(b, kind::INT)) {
    Int x = get_val(a);
    Int y = get_val(b);

    return ast_set_to_int(a, x * y);
  }

  if (is(a, kind::FRAC) && is(b, kind::FRAC)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));
    Int w = get_val(operand(b, 0));
    Int z = get_val(operand(b, 1));

    Int e = x * w;
    Int k = y * z;

    if (e % k == 0) {
      return ast_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));
    return ast_set_to_fra(a, e / g, k / g);
  }

  if (is(a, kind::FRAC) && is(b, kind::INT)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));
    Int w = get_val(b);

    Int e = x * w;
    Int k = y;

    if (e % k == 0) {
      return ast_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));
    return ast_set_to_fra(a, e / g, k / g);
  }

  if (is(b, kind::FRAC) && is(a, kind::INT)) {
    Int x = get_val(operand(b, 0));
    Int y = get_val(operand(b, 1));
    Int w = get_val(a);

    Int e = x * w;
    Int k = y;

    if (e % k == 0) {
      return ast_set_to_int(a, e / k);
    }

    Int g = abs(gcd(e, k));
    return ast_set_to_fra(a, e / g, k / g);
  }
}

inline void ast_set_inplace_mul_consts(ast *a, Int b) {
  assert(is(a, kind::CONST));

  if (is(a, kind::INT)) {
    Int x = get_val(a);

    return ast_set_to_int(a, x * b);
  }

  if (is(a, kind::FRAC)) {
    Int x = get_val(operand(a, 0));
    Int y = get_val(operand(a, 1));

    Int e = x * b;

    if (e % y == 0) {
      return ast_set_to_int(a, e / y);
    }

    Int g = abs(gcd(e, y));
    return ast_set_to_fra(a, e / g, y / g);
  }
}

inline void ast_set_op_inplace_mul_consts(ast *a, size_t i, ast *b) {
  ast_set_inplace_mul_consts(operand(a, i), b);
}

inline void ast_set_op_inplace_mul_consts(ast *a, size_t i, Int b) {
  ast_set_inplace_mul_consts(operand(a, i), b);
}

inline void ast_set_to_mul(Int v, ast *a) {
  if (is(a, kind::MUL)) {
    return a->insert(integer(v), 0);
  }

  a->ast_childs = std::vector<ast>({integer(v), *a});

  ast_set_kind(a, kind::MUL);
}

inline void ast_set_to_mul(ast *a, ast *b) {
  if (is(a, kind::MUL)) {
    return a->insert(*b);
  }

  a->ast_childs = std::vector<ast>({*b, *a});
  ast_set_kind(a, kind::MUL);
}

inline void ast_set_op_to_mul(ast *a, size_t i, Int v) {
  return ast_set_to_mul(v, operand(a, i));
}

inline void ast_set_op_to_mul(ast *a, size_t i, ast *v) {
  return ast_set_to_mul(operand(a, i), v);
}

inline void ast_set_to_pow(ast *a, Int e) {
  a->ast_childs = std::vector<ast>({*a, integer(e)});
  ast_set_kind(a, kind::POW);
}

inline void ast_set_to_add(ast *a, ast *e) {

  if (is(a, kind::ADD)) {
    return a->insert(*e);
  }

  a->ast_childs = std::vector<ast>({*a, *e});

  ast_set_kind(a, kind::ADD);
}

inline void ast_set_op_to_add(ast *a, size_t i, ast *v) {
  ast_set_to_add(operand(a, i), v);
}

inline void ast_set_op_pow_add_to_deg(ast *a, size_t i, ast *e) {
  assert(is(operand(a, i), kind::POW));

  a->ast_childs[i] =
      create(kind::POW, {*operand(operand(a, i), 0),
                         create(kind::ADD, {*operand(operand(a, i), 1), *e})});
}

inline void ast_set_op_pow_add_to_deg(ast *a, size_t i, Int e) {
  assert(is(operand(a, i), kind::POW));

  a->ast_childs[i] = create(
      kind::POW, {*operand(operand(a, i), 0),
                  create(kind::ADD, {*operand(operand(a, i), 1), integer(e)})});
}

inline void ast_set_op_to_pow(ast *a, size_t i, Int v) {
  a->ast_childs[i] = create(kind::POW, {*operand(a, i), integer(v)});
}

inline void ast_set_op_to_pow(ast *a, size_t i, ast *v) {
  a->ast_childs[i] = create(kind::POW, {*operand(a, i), *v});
}

inline bool eval_add_consts(ast *u, size_t i, ast *v, size_t j) {
  ast_set_op_inplace_add_consts(u, i, operand(v, j));
  return true;
}

inline bool eval_mul_consts(ast *u, size_t i, ast *v, size_t j) {
  ast_set_op_inplace_mul_consts(u, i, operand(v, j));
  return true;
}

inline bool eval_add_int(ast *a, Int b) {
  if (is(a, kind::INT)) {
    ast_set_to_int(a, get_val(a) + b);

    return true;
  }

  assert(is(a, kind::FRAC));

  Int num = get_val(operand(a, 0));
  Int den = get_val(operand(a, 1));

  num = b * den + num;

  Int cff = abs(gcd(num, den));

  if (den / cff == 1) {
    ast_set_to_int(a, num / cff);

    return true;
  }

  ast_set_to_fra(a, num / cff, den / cff);

  return true;
}

inline bool eval_add_nconst(ast *u, size_t i, ast *v, size_t j) {
  assert(is(u, kind::ADD) && is(v, kind::ADD));

  assert(is(operand(u, i), kind::SUMMABLE));
  assert(is(operand(v, j), kind::SUMMABLE));

  ast *a = operand(u, i);
  ast *b = operand(v, j);

  if (a == 0 || b == 0) {
    return false;
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
      ast_set_op_to_mul(u, i, 2);
      return true;
    }

    return false;
  }

  long size_a = size_of(a);

  if (kind & kind::MUL) {
    long size_b = size_of(b);

    ast *c;

    long size_c = -1;

    if (is(operand(a, 0), kind::INT) && is(operand(b, 0), kind::INT) &&
        std::abs(size_a - size_b) != 0) {
      return false;
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
        return 0;
      }
    }

    int ka = kind_of(operand(a, 0));
    int kb = kind_of(operand(b, 0));

    if ((ka & kind::CONST) && (kb & kind::CONST)) {
      ast_set_op_inplace_add_consts(a, 0, operand(b, 0));
    } else if (ka & kind::CONST) {
      ast_set_op_inplace_add_consts(a, 0, 1);
    } else {
      ast_set_to_mul(2, a);
    }

    // ast_set_operand(u, a, i);

    if (size_c != -1) {
      u->ast_childs[i] = *a;
      v->ast_childs[j] = *b;
    }

    return true;
  }

  if (is(b, kind::MUL)) {
    return false;
  }

  // a is a mul and b is a sym or a pow
  assert(is(a, kind::MUL));
  assert(is(b, kind::SYM | kind::POW));

  if (size_a > 2) {
    return false;
  }

  ast *a0 = operand(a, 0);
  ast *a1 = operand(a, 1);

  long ki = kind_of(a1) & kind_of(b);

  if (!is(a0, kind::CONST) || !ki) {
    return 0;
  }

  if (compare(b, a1, kind::MUL) == 0) {
    ast_set_op_inplace_add_consts(a, 0, 1);

    return true;
  }

  return false;
}

inline bool eval_mul_nconst(ast *u, size_t i, ast *v, size_t j) {
  assert(is(u, kind::MUL) && is(v, kind::MUL));
  assert(is(operand(u, i), kind::MULTIPLICABLE));
  assert(is(operand(v, j), kind::MULTIPLICABLE));

  ast *a = operand(u, i);
  ast *b = operand(v, j);

  if (a == 0 || b == 0) {
    return false;
  }

  if (is(a, kind::ADD) && is(b, kind::ADD)) {

    if (compare(a, b, kind::MUL) == 0) {
      ast_set_op_to_pow(u, i, 2);
      return true;
    }

    return false;
  }

  if (is(a, kind::ADD) && is(b, kind::POW)) {
    if (compare(a, operand(b, 0), kind::MUL) == 0) {
      u->ast_childs[i] = create(kind::ADD, {integer(1), *operand(b, 1)});

      reduce(&u->ast_childs[i]);

      return true;
    }

    return false;
  }

  if (is(a, kind::POW) && is(b, kind::ADD)) {
    if (compare(b, operand(a, 0), kind::MUL) == 0) {
      ast_set_op_pow_add_to_deg(u, i, 1);
      reduce(&u->ast_childs[i]);

      return true;
    }

    return false;
  }

  if (is(a, kind::POW) && is(b, kind::POW)) {
    if (compare(operand(a, 0), operand(b, 0), kind::MUL) == 0) {
      ast_set_op_pow_add_to_deg(u, i, operand(b, 1));
      reduce(operand(u, i));
      return true;
    }

    return 0;
  }

  if (is(a, kind::SYM | kind::FUNC) && is(b, kind::SYM | kind::FUNC)) {
    if (compare(a, b, kind::MUL) == 0) {
      ast_set_op_to_pow(u, i, 2);
      return true;
    }

    return false;
  }

  if (is(a, kind::POW) && is(b, kind::SYM | kind::FUNC)) {
    if (compare(operand(a, 0), b, kind::MUL) == 0) {

      ast_set_op_pow_add_to_deg(u, i, 1);
      reduce(operand(u, i));

      return true;
    }

    return false;
  }

  return false;
}

inline bool eval_add_add(ast *a, ast *b) {
  assert(is(a, kind::ADD));
  assert(is(b, kind::ADD));

  size_t j = 0;
  size_t i = 0;

  ast *u = 0;
  ast *v = 0;

  while (i < size_of(a) && j < size_of(b)) {
    assert(!is(operand(b, j), kind::ADD));

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
      reduced = eval_add_consts(a, i, b, j);
    } else if (is(u, kind::SUMMABLE) && is(v, kind::SUMMABLE)) {
      reduced = eval_add_nconst(a, i, b, j);
    }

    if (reduced) {
      i = i + 1;
      j = j + 1;
    } else {
      int order = compare(u, v, kind::ADD);

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
      } else if (u == 0 || compare(u, v, kind::ADD) < 0) {
        i++;
      } else {
        a->insert(*v, i++);
        j = j + 1;
      }
    }
  }

  return true;
}

inline bool eval_mul_mul(ast *a, ast *b) {
  assert(is(a, kind::MUL));
  assert(is(b, kind::MUL));

  size_t j = 0;
  size_t i = 0;

  ast *u = 0;
  ast *v = 0;

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

inline bool eval_mul_int(ast *u, size_t i, Int v) {
  ast *a = operand(u, i);

  if (is(a, kind::INT)) {
    ast_set_op_inplace_mul_consts(u, i, v);
    return true;
  }

  if (is(a, kind::ADD | kind::SUB)) {
    for (size_t j = 0; j < size_of(a); j++) {
      ast_set_op_inplace_mul_consts(a, j, v);
    }

    return true;
  }

  if (is(a, kind::MUL)) {
    ast_set_op_to_mul(u, i, v);
    return true;
  }

  if (is(a, kind::DIV)) {
    eval_mul_int(a, 0, v);
    return true;
  }

  if (is(a, kind::SQRT | kind::POW | kind::FACT | kind::FUNC | kind::SYM)) {
    ast_set_op_inplace_mul_consts(u, i, v);
    return true;
  }

  return true;
}

inline bool ast_replace_with(ast *a, ast *t) {
  if (is(t, kind::INT)) {
    ast_set_to_int(a, get_val(t));
    return true;
  }

  if (is(t, kind::SYM)) {
    ast_set_to_sym(a, get_id(t));
    return true;
  }

  if (is(t, kind::FUNC)) {
    a->ast_sym = strdup(get_id(t));
    ast_set_kind(a, kind::FUNC);

    a->ast_childs = std::vector<ast>();

    // TODO: maybe std::move works here
    for (size_t i = 0; i < size_of(t); i++) {
      a->insert(t->ast_childs[i]);
    }

    return a;
  }

  ast_set_kind(a, kind_of(t));

  a->ast_childs = std::vector<ast>();

  // TODO: maybe std::move works here
  for (size_t i = 0; i < size_of(t); i++) {
    a->insert(t->ast_childs[i]);
  }

  return a;
}

inline bool ast_raise_to_first_op(ast *a) {
  if (is(operand(a, 0), kind::INT)) {
    ast_set_to_int(a, get_val(operand(a, 0)));
    return true;
  }

  if (is(operand(a, 0), kind::SYM)) {
    ast_set_to_sym(a, get_id(operand(a, 0)));
    return true;
  }

  if (is(operand(a, 0), kind::FUNC)) {
    a->ast_sym = strdup(get_id(operand(a, 0)));
    ast_set_kind(a, kind::FUNC);
    a->ast_childs = ast(a->ast_childs[0]).ast_childs;
    return true;
  }

  *a = ast(a->ast_childs[0]);

  return true;
}

void reduce_add(ast *a) {
  for (size_t i = 0; i < size_of(a); i++) {
    reduce(operand(a, i));
  }

  sort_childs(a, 0, size_of(a) - 1);

  if (is(operand(a, 0), kind::ADD)) {
    ast t = a->ast_childs[0];

    a->remove(0);

    eval_add_add(a, &t);
  }

  size_t j = 0;

  for (long i = 1; i < (long)size_of(a); i++) {
    ast *aj = operand(a, j);
    ast *ai = operand(a, i);

    bool reduced = 0;

    if (is(ai, kind::FAIL) || is(aj, kind::FAIL)) {
      return ast_set_to_fail(a);
    } else if (is(ai, kind::UNDEF) || is(aj, kind::UNDEF)) {
      return ast_set_to_undefined(a);
    } else if (ast_is_zero(operand(a, j))) {
      a->remove(j);
    } else if (is(ai, kind::ADD)) {
      ast t = a->ast_childs[i];
      a->remove(i--);
      eval_add_add(a, &t);
    } else if (is(aj, kind::CONST) && is(ai, kind::CONST)) {
      reduced = eval_add_consts(a, j, a, i);

      if (reduced) {
        a->remove(i--);
      }
    } else if (is(aj, kind::INT) && get_val(aj) == 0) {
      a->remove(j);
      i = i - 1;
    } else if (is(ai, kind::INT) && get_val(ai) == 0) {
      a->remove(i--);
    } else if (is(aj, kind::SUMMABLE) && is(ai, kind::SUMMABLE)) {
      reduced = eval_add_nconst(a, j, a, i);

      if (reduced) {
        a->remove(i--);
      }
    }

    if (reduced == false) {
      j = i;
    }
  }

  if (size_of(a) == 0) {
    ast_set_to_int(a, 0);
  } else if (size_of(a) == 1) {
    ast_raise_to_first_op(a);
  }
}

void reduce_mul(ast *a) {
  for (size_t i = 0; i < size_of(a); i++) {
    reduce(operand(a, i));
  }

  sort_childs(a, 0, size_of(a) - 1);

  size_t j = 0;

  if (is(operand(a, 0), kind::MUL)) {
    ast t = a->ast_childs[0];
    a->remove(0);
    eval_mul_mul(a, &t);
  }

  for (long i = 1; i < (long)size_of(a); i++) {
    ast *aj = operand(a, j);
    ast *ai = operand(a, i);

    bool reduced = 0;

    if (is(ai, kind::INT) && get_val(ai) == 0) {
      return ast_set_to_int(a, 0);
    }

    if (is(aj, kind::INT) && get_val(aj) == 0) {
      return ast_set_to_int(a, 0);
    }

    if (is(ai, kind::FAIL) || is(aj, kind::FAIL)) {
      return ast_set_to_fail(a);
    } else if (is(ai, kind::UNDEF) || is(aj, kind::UNDEF)) {
      return ast_set_to_undefined(a);
    } else if (ast_is_zero(operand(a, i)) || ast_is_zero(operand(a, j))) {
      return ast_set_to_int(a, 0);
    } else if (is(ai, kind::MUL)) {
      ast t = a->ast_childs[i];
      a->remove(i--);
      eval_mul_mul(a, &t);
    } else if (is(aj, kind::CONST) && is(ai, kind::CONST)) {

      reduced = eval_mul_consts(a, j, a, i);

      if (reduced) {
        a->remove(i--);
      }

    } else if (is(aj, kind::INT) && get_val(aj) == 1) {
      a->remove(j);
      i = i - 1;
    } else if (is(ai, kind::INT) && get_val(ai) == 1) {
      a->remove(i--);
    } else if (is(aj, kind::MULTIPLICABLE) && is(ai, kind::MULTIPLICABLE)) {
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
    ast_raise_to_first_op(a);
  }
}

void reduce_sub(ast *a) {
  for (size_t i = 1; i < size_of(a); i++) {
    eval_mul_int(a, i, -1);
  }

  ast_set_kind(a, kind::ADD);

  reduce(a);
}

void reduce_pow(ast *a) {
  reduce(operand(a, 1));

  if (is(operand(a, 1), kind::INT) && get_val(operand(a, 1)) == 0) {
    return ast_set_to_int(a, 1);
  }

  reduce(operand(a, 0));

  if (is(operand(a, 0), kind::POW)) {
    ast_set_op_to_mul(operand(a, 0), 1, operand(a, 1));
    ast_raise_to_first_op(a);
    return reduce(a);
  }

  if (is(operand(a, 0), kind::MUL)) {
    for (size_t i = 0; i < size_of(operand(a, 0)); i++) {
      ast_set_op_to_pow(operand(a, 0), i, operand(a, 1));
    }
    ast_raise_to_first_op(a);

    return reduce(a);
  }

  if (!is(operand(a, 1), kind::INT)) {
    return;
  }

  if (get_val(operand(a, 1)) == 1) {
    ast_raise_to_first_op(a);
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
      ast_set_to_int(a, d);
    } else {
      ast_set_to_fra(a, 1, d);
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
      ast_set_to_fra(a, b, c);
    } else {
      ast_set_to_fra(a, c, b);
    }

    return;
  }

  if (is(operand(a, 0), kind::TERMINAL)) {
    return;
  }

  if (is(operand(a, 0), kind::MUL)) {
    long long y = get_val(operand(a, 1)).longValue();

    ast b = a->ast_childs[0];

    ast_set_to_int(a, 1);

    while (y) {
      if (y % 2 == 1) {
        ast_set_to_mul(a, &b);
        reduce(a);
      }

      y = y >> 1;

      ast t = create(kind::MUL);

      for (size_t i = 0; i < size_of(&b); i++) {
        ast v = create(kind::MUL, {
                                      *operand(&b, i),
                                      *operand(&b, i),
                                  });

        t.insert(v);
      }

      ast_replace_with(&b, &t);
    }

    return;
  }
}

void reduce_div(ast *a) {
  ast_set_kind(a, kind::MUL);

  ast_set_op_to_pow(a, 1, -1);

  return reduce(a);
}

void reduce_sqr(ast *a) {
  ast_set_kind(a, kind::POW);

  a->insert(fraction(1, 2), 1);

  return reduce(a);
}

ast *reduce_fac(ast *a) {
  if (is(operand(a, 0), kind::INT)) {
    ast_set_to_int(a, fact(get_val(operand(a, 0))));

    return a;
  }

  if (is(operand(a, 0), kind::FRAC)) {
    Int c = fact(get_val(operand(operand(a, 0), 0)));
    Int d = fact(get_val(operand(operand(a, 0), 1)));

    Int g = abs(gcd(c, d));

    ast_set_to_fra(a, c / g, d / g);

    return a;
  }

  return a;
}

ast *reduce_fra(ast *a) {
  Int b = get_val(operand(a, 0));
  Int c = get_val(operand(a, 1));

  Int d = abs(gcd(b, c));

  ast_set_to_fra(a, b / d, c / d);
  return a;
}

void reduce(ast *a) {
  if (is(a, kind::LIST)) {
    for (size_t i = 0; i < size_of(a); i++) {
      reduce(&a->ast_list->members[i]);
    }
  } else if (is(a, kind::SET)) {
    for (size_t i = 0; i < size_of(a); i++) {
      reduce(&a->ast_set->members[i]);
    }

    trim(a->ast_set);
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

ast expand_mul(ast *a, size_t i, ast *b, size_t j) {
  ast *r = operand(a, i);
  ast *s = operand(b, j);

  if (is(r, kind::ADD) && is(s, kind::ADD)) {
    ast u = create(kind::ADD);

    for (size_t k = 0; k < size_of(r); k++) {
      for (size_t t = 0; t < size_of(s); t++) {
        u.insert(
            create(kind::MUL, {ast(r->ast_childs[k]), ast(s->ast_childs[t])}));
      }
    }

    return u;
  }

  if (is(r, kind::ADD)) {
    ast u = create(kind::ADD);

    for (size_t k = 0; k < size_of(r); k++) {
      u.insert(
          create(kind::MUL, {ast(r->ast_childs[k]), ast(b->ast_childs[j])}));
    }

    return u;
  }

  if (is(s, kind::ADD)) {
    ast u = create(kind::ADD);

    for (size_t k = 0; k < size_of(s); k++) {
      u.insert(
          create(kind::MUL, {ast(a->ast_childs[i]), ast(s->ast_childs[k])}));
    }

    return u;
  }

  return create(kind::MUL, {ast(a->ast_childs[i]), ast(b->ast_childs[j])});
}

bool expand_pow(ast *u, Int n, ast *a) {
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

    ast o = *u;

    ast f = o[0];

    o.remove(0);

    if (size_of(&o) == 0) {
      ast_set_to_int(&o, 0);
    }

    if (size_of(&o) == 1) {
      ast_raise_to_first_op(&o);
    }

    ast s = create(kind::ADD);

    for (Int k = 0; k <= n; k++) {
      ast z = create(kind::MUL, {integer(c / (fact(k) * fact(n - k))),
                                 create(kind::POW, {f, integer(n - k)})});

      ast t;

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

void expand(ast *a) {
  if (is(a, kind::TERMINAL)) {
    return;
  }

  if (is(a, kind::SUB | kind::DIV | kind::FACT)) {
    reduce(a);
  }

  if (is(a, kind::POW)) {
    expand(operand(a, 0));
    expand(operand(a, 1));

    if (is(operand(a, 1), kind::INT)) {
      expand_pow(operand(a, 0), get_val(operand(a, 1)), a);
    }
  }

  if (is(a, kind::MUL)) {
    // printf("from ----> %s\n", to_string(a).c_str());
    while (size_of(a) > 1) {
      expand(operand(a, 0));
      expand(operand(a, 1));

      ast t = expand_mul(a, 0, a, 1);

      a->insert(t, 0);

      a->remove(1);
      a->remove(1);
    }

    ast_raise_to_first_op(a);
    // printf("to ----> %s\n", to_string(a).c_str());
  }

  if (is(a, kind::ADD)) {
    for (size_t i = 0; i < size_of(a); i++) {
      expand(operand(a, i));
    }
  }
  // printf("from %s\n", to_string(a).c_str());
  reduce(a);
  // printf("to %s\n", to_string(a).c_str());
}

ast &ast::operator+=(const ast &a) {
  if (is(this, kind::ADD)) {
    this->insert(a);
  } else {
    *this = create(kind::ADD, {*this, a});
  }

  return *this;
}

ast &ast::operator+=(ast &&a) {
  if (is(this, kind::ADD)) {
    this->insert(a);
  } else {
    *this = create(kind::ADD, {*this, a});
  }

  return *this;
}

ast &ast::operator-=(const ast &a) {
  if (is(this, kind::ADD)) {
    this->insert(create(kind::MUL, {integer(-1), a}));
  } else {
    *this = create(kind::ADD, {*this, create(kind::MUL, {integer(-1), a})});
  }

  return *this;
}

ast &ast::operator-=(ast &&a) {
  if (is(this, kind::ADD)) {
    this->insert(create(kind::MUL, {integer(-1), a}));
  } else {
    *this = create(kind::ADD, {*this, create(kind::MUL, {integer(-1), a})});
  }

  return *this;
}

ast ast::operator+(const ast &a) {
  if (is(this, kind::ADD)) {
    ast t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::ADD, {*this, a});
}

ast ast::operator+(ast &&a) {
  if (is(this, kind::ADD)) {
    ast t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::ADD, {*this, a});
}

ast ast::operator-(const ast &a) {
  if (is(this, kind::SUB)) {
    ast t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::SUB, {*this, a});
}

ast ast::operator-(ast &&a) {
  if (is(this, kind::SUB)) {
    ast t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::SUB, {*this, a});
}

ast ast::operator*(const ast &a) {
  if (is(this, kind::MUL)) {
    ast t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::MUL, {*this, a});
}

ast ast::operator*(ast &&a) {
  if (is(this, kind::MUL)) {
    ast t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::MUL, {*this, a});
}

ast ast::operator/(const ast &a) {
  if (is(this, kind::INT) && is(&a, kind::INT)) {
    return create(kind::FRAC, {*this, a});
  }

  if (is(this, kind::DIV)) {
    ast t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::DIV, {*this, a});
}

ast ast::operator/(ast &&a) {
  if (is(this, kind::INT) && is(&a, kind::INT)) {
    return create(kind::FRAC, {*this, a});
  }

  if (is(this, kind::DIV)) {
    ast t = *this;
    t.insert(a);
    return t;
  }

  return create(kind::DIV, {*this, a});
}

bool ast::operator==(const ast &other) {
  ast a = other;
  ast b = *this;

  sort(&a);
  sort(&b);

  enum kind k = kind_of == kind::ADD ? kind::ADD : kind::MUL;

  return ast_op_cmp(&a, &b, k) == 0;
}

bool ast::operator==(ast &&a) {
  ast b = *this;

  sort(&a);
  sort(&b);

  enum kind k = kind_of == kind::ADD ? kind::ADD : kind::MUL;

  return ast_op_cmp(&a, &b, k) == 0;
}

bool ast::operator!=(const ast &other) {
  ast a = other;
  ast b = *this;

  sort(&a);
  sort(&b);

  enum kind k = kind_of == kind::ADD ? kind::ADD : kind::MUL;

  return ast_op_cmp(&a, &b, k) != 0;
}

bool ast::operator!=(ast &&a) {
  ast b = *this;

  sort(&a);
  sort(&b);

  enum kind k = kind_of == kind::ADD ? kind::ADD : kind::MUL;

  return ast_op_cmp(&a, &b, k) != 0;
}

ast ast::operator+() { return *this; }

ast ast::operator-() {
  if (is(this, kind::INF)) {
    return create(kind::NEG_INF);
  }

  if (is(this, kind::NEG_INF)) {
    return create(kind::INF);
  }

  return create(kind::MUL, {integer(-1), *this});
}

ast operator*(Int i, ast &&other) { return integer(i) * other; }

ast operator*(Int i, ast &other) { return integer(i) * other; }

ast operator+(Int i, ast &&other) { return integer(i) + other; }

ast operator+(Int i, ast &other) { return integer(i) + other; }

ast operator-(Int i, ast &&other) { return integer(i) - other; }

ast operator-(Int i, ast &other) { return integer(i) - other; }

ast operator/(Int i, ast &&other) { return integer(i) / other; }

ast operator/(Int i, ast &other) { return integer(i) / other; }

ast operator*(int i, ast &&other) { return integer(i) * other; }

ast operator*(int i, ast &other) { return integer(i) * other; }

ast operator+(int i, ast &&other) { return integer(i) + other; }

ast operator+(int i, ast &other) { return integer(i) + other; }

ast operator-(int i, ast &&other) { return integer(i) - other; }

ast operator-(int i, ast &other) { return integer(i) - other; }

ast operator/(int i, ast &&other) { return integer(i) / other; }

ast operator/(int i, ast &other) { return integer(i) / other; }

ast operator*(long i, ast &&other) { return integer(i) * other; }

ast operator*(long i, ast &other) { return integer(i) * other; }

ast operator+(long i, ast &&other) { return integer(i) + other; }

ast operator+(long i, ast &other) { return integer(i) + other; }

ast operator-(long i, ast &&other) { return integer(i) - other; }

ast operator-(long i, ast &other) { return integer(i) - other; }

ast operator/(long i, ast &&other) { return integer(i) / other; }

ast operator/(long i, ast &other) { return integer(i) / other; }

ast operator*(long long i, ast &&other) { return integer(i) * other; }

ast operator*(long long i, ast &other) { return integer(i) * other; }

ast operator+(long long i, ast &&other) { return integer(i) + other; }

ast operator+(long long i, ast &other) { return integer(i) + other; }

ast operator-(long long i, ast &&other) { return integer(i) - other; }

ast operator-(long long i, ast &other) { return integer(i) - other; }

ast operator/(long long i, ast &&other) { return integer(i) / other; }

ast operator/(long long i, ast &other) { return integer(i) / other; }

ast pow(const ast &a, const ast &b) { return create(kind::POW, {a, b}); }

ast pow(ast &&a, ast &&b) { return create(kind::POW, {a, b}); }

ast pow(const ast &a, ast &&b) { return create(kind::POW, {a, b}); }

ast pow(ast &&a, const ast &b) { return create(kind::POW, {a, b}); }

ast sqrt(const ast &a) { return create(kind::SQRT, {a}); }

ast sqrt(ast &&a) { return create(kind::SQRT, {a}); }

ast fact(const ast &a) { return create(kind::FACT, {a}); }

ast fact(ast &&a) { return create(kind::FACT, {a}); }

ast undefined() { return create(kind::UNDEF); }

ast fail() { return create(kind::FAIL); }

ast inf() { return create(kind::INF); }

ast reduce(ast &a) {
  ast b = a;

  reduce(&b);

  return b;
}

ast reduce(ast &&a) {
  ast b = a;

  reduce(&b);

  return b;
}

ast expand(ast &a) {
  ast b = a;

  expand(&b);

  return b;
}

ast expand(ast &&a) {
  ast b = a;

  expand(&b);

  return b;
}

ast first(ast &a) {
  if (is(&a, kind::LIST)) {
    return a.ast_list->members[0];
  }

  if (is(&a, kind::SET)) {
    return a.ast_set->members[0];
  }

  assert(!is(&a, kind::TERMINAL));

  return a[0];
}

ast rest(ast &a) {
  if (is(&a, kind::LIST)) {
    list L = *a.ast_list;
    return ast(rest(L));
  }

  if (is(&a, kind::SET)) {
    set L = *a.ast_set;
    return ast(rest(L));
  }

  assert(!is(&a, kind::TERMINAL));

  return a[0];
}

ast append(const ast &a, const ast &b) {
  assert(is(&a, kind::LIST) && is(&b, kind::LIST));

  list L = *a.ast_list;
  list M = *b.ast_list;

  return ast(append(L, M));
}

ast append(const ast &a, ast &&b) {
  assert(is(&a, kind::LIST) && is(&b, kind::LIST));

  list L = *a.ast_list;
  list M = *b.ast_list;

  return ast(append(L, M));
}

ast join(const ast &a, const ast &b) {
  assert(is(&a, kind::LIST) && is(&b, kind::LIST));

  list L = *a.ast_list;
  list M = *b.ast_list;

  return ast(join(L, M));
}

ast join(const ast &a, ast &&b) {
  assert(is(&a, kind::LIST) && is(&b, kind::LIST));

  list L = *a.ast_list;
  list M = *b.ast_list;

  return ast(join(L, M));
}

ast difference(const ast &a, ast &&b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.ast_set;
  set M = *b.ast_set;

  return ast(difference(L, M));
}

ast difference(const ast &a, const ast &b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.ast_set;
  set M = *b.ast_set;

  return ast(difference(L, M));
}

ast unnification(const ast &a, ast &&b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.ast_set;
  set M = *b.ast_set;

  return ast(unnification(L, M));
}

ast uniffication(const ast &a, const ast &b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.ast_set;
  set M = *b.ast_set;

  return ast(difference(L, M));
}
ast intersection(const ast &a, ast &&b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.ast_set;
  set M = *b.ast_set;

  return ast(intersection(L, M));
}

ast intersection(const ast &a, const ast &b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.ast_set;
  set M = *b.ast_set;

  return ast(intersection(L, M));
}
int exists(const ast &a, ast &&b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.ast_set;
  set M = *b.ast_set;

  return exists(L, M);
}

int exists(const ast &a, const ast &b) {
  assert(is(&a, kind::SET) && is(&b, kind::SET));

  set L = *a.ast_set;
  set M = *b.ast_set;

  return exists(L, M);
}

void replace_rec(ast *a, ast *b, ast *c) {
  if (ast_op_cmp(a, b, a->kind_of == kind::ADD ? kind::ADD : kind::MUL) == 0) {
    *a = *c;
    return;
  }

  if (is(a, kind::TERMINAL))
    return;

  for (size_t i = 0; i < size_of(a); i++) {
    replace_rec(operand(a, i), b, c);
  }
}

ast replace(ast &a, ast &b, ast &c) {
  ast d = a;

  replace_rec(&d, &b, &c);

  return d;
}

ast replace(ast &a, ast &&b, ast &&c) {
  ast d = a;

  replace_rec(&d, &b, &c);

  return d;
}

ast replace(ast &a, ast &&b, ast &c) {
  ast d = a;

  replace_rec(&d, &b, &c);

  return d;
}

ast replace(ast &a, ast &b, ast &&c) {
  ast d = a;

  replace_rec(&d, &b, &c);

  return d;
}

ast map(ast &u, ast &v, ast (*f)(ast &, ast &)) {
  if (is(&u, kind::TERMINAL)) {
    return f(u, v);
  }

  if (size_of(&u) == 0) {
    return f(u, v);
  }

  ast t = create(kind_of(&u));

  if (is(&u, kind::FUNC)) {
    u.ast_sym = strdup(u.ast_sym);
  }

  for (size_t i = 0; i < size_of(&u); i++)
    t.insert(f(u[i], v));

  return t;
}

ast map(ast &u, ast (*f)(ast &)) {
  if (is(&u, kind::TERMINAL)) {
    return f(u);
  }

  if (size_of(&u) == 0) {
    return f(u);
  }

  ast t = create(kind_of(&u));

  if (is(&u, kind::FUNC)) {
    u.ast_sym = strdup(u.ast_sym);
  }

  for (size_t i = 0; i < size_of(&u); i++) {
    t.insert(f(u[i]));
  }

  return t;
}

bool free_of_rec(ast *a, ast *b) {
  if (ast_op_cmp(a, b, a->kind_of == kind::ADD ? kind::ADD : kind::MUL) == 0) {
    return false;
  }

  for (size_t i = 0; i < size_of(a); i++) {
    if (free_of_rec(operand(a, i), b) == false) {
      return false;
    }
  }
  return true;
}

bool ast::freeOf(ast &a) { return free_of_rec(this, &a); }
bool ast::freeOf(ast &&a) { return free_of_rec(this, &a); }

list::list(std::initializer_list<ast> &&a) { members = a; }

list::list(std::vector<ast> &&a) { members = std::move(a); }

list::list(std::vector<ast> &a) { members = a; }

void list::append(ast &&a) { members.push_back(a); }

void list::append(const ast &a) { members.push_back(a); }

list append(list &a, const ast &b) {
  list L = a;
  L.append(b);
  return L;
}

list append(list &a, ast &&b) {
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
  std::vector<ast> L;

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
  std::vector<ast> L;

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
  std::vector<ast> l;

  for (size_t i = from; i < size(); i++) {
    l.push_back(members[i]);
  }

  return l;
}

list rest(list &a, size_t from) {
  std::vector<ast> l;

  for (size_t i = from; i < a.size(); i++) {
    l.push_back(a[i]);
  }

  return l;
}

ast fist(list &l) { return l[0]; }

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

bool list::match(list* a) {
	if(size() != a->size()) {
		return (long)size() - (long)a->size();
	}

	for(size_t i = 0; i < size(); i++) {
		bool found = false;

		for(size_t j = 0; j < size(); j++) {
			long cmp = compare(&members[i], &a->members[j], kind::ADD);

			if(cmp == 0) {
				found = true;
				break;
			}
		}

		if(found == false) {
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

std::string to_string(list& a) {
	std::string s = "{";
	for(size_t i = 0; i < a.size(); i++) {
		s += to_string(&a.members[i]);
		if(i < a.size() - 1) {
			s += ", ";
		}
	}

	s += "}";

	return s;
}

std::string to_string(list* a) {
	std::string s = "{";
	for(size_t i = 0; i < a->size(); i++) {
		s += to_string(&a->members[i]);

		if(i < a->size() - 1) {
			s += ", ";
		}
	}

	s += "}";

	return s;
}
bool list::operator==(list& a) {
	return this->match(&a) == 0;
}

bool list::operator==(list&& a) {
	return this->match(&a) == 0;
}

bool list::operator!=(list& a) {
	return this->match(&a) != 0;
}

bool list::operator!=(list&& a) {
	return this->match(&a) != 0;
}

long int sort_split(std::vector<ast> &a, long l, long r) {
  long int i = l - 1;

  ast &p = a[r];

  for (long int j = l; j < r; j++) {
    if (compare(&a[j], &p, kind::ADD) < 0) {
      std::swap(a[++i], a[j]);
    }
  }

  i = i + 1;

  std::swap(a[i], a[r]);

  return i;
}

void sort(std::vector<ast> &a, long int l, long int r) {
  if (l < r) {
    long int m = sort_split(a, l, r);

    sort(a, l, m - 1);
    sort(a, m + 1, r);
  }
}

set::set(std::initializer_list<ast> &&a) {
  members = std::move(a);
  trim(this);
}

set::set(std::vector<ast> &&a) {
  members = std::move(a);
  trim(this);
}

set::set(std::vector<ast> &a) {
  members = a;
  trim(this);
}

set difference(set &L, set &M) {
  set t = {};

  size_t j = 0;
  size_t i = 0;

	while(true) {
		if (j >= M.size() || i >= L.size()) {
			break;
		}

    int cmp = compare(&L[i], &M[j], kind::ADD);

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

set unnification(set &L, set &M) {
  set t = {};

  size_t j = 0;
  size_t i = 0;

	while(true) {
		if (j >= M.size() || i >= L.size()) {
			break;
		}

    int cmp = compare(&L[i], &M[j], kind::ADD);

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

	while(true) {
		if (j >= M.size() || i >= L.size()) {
			break;
		}

    int cmp = compare(&L[i], &M[j], kind::ADD);

		if(cmp < 0) {
			i = i + 1;
		}

		if(cmp > 0) {
			j = j + 1;
		}

    if (cmp == 0) {
      t.members.push_back(L[i++]);
			j = j + 1;
    }
	}

	return t;
}

int search(std::vector<ast> &a, ast &x, int l, int r) {
  if (r >= l) {
    int m = l + (r - l) / 2;

    int cmp = compare(&a[m], &x, kind::ADD);

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

int exists(set &L, ast &e) { return search(L.members, e, 0, L.size() - 1); }

set rest(set &a, size_t from) {
  std::vector<ast> l;

  for (size_t i = from; i < a.size(); i++) {
    l.push_back(a[i]);
  }

  return l;
}

ast fist(set &l) { return l[0]; }

void trim(set *s) {
	sort(s->members, 0, s->members.size() - 1);

  for (long i = 0; i < ((long)s->size()) - 1; i++) {
    int cmp = compare(&s->members[i], &s->members[i + 1], kind::ADD);

    if (cmp == 0) {
      s->members.erase(s->members.begin() + i);
			i = i - 1;
    }
  }
}

long set::match(set* a) {
	if(size() != a->size()) {
		return (long)size() - (long)a->size();
	}
	for(size_t i = 0; i < size(); i++) {
		long cmp = compare(&members[i], &a->members[i], kind::ADD);
		if(cmp != 0) {
			return cmp;
		}
	}

	return 0;
}

bool set::operator==(set& a) {
	return this->match(&a) == 0;
}

bool set::operator==(set&& a) {
	return this->match(&a) == 0;
}

bool set::operator!=(set& a) {
	return this->match(&a) != 0;
}

bool set::operator!=(set&& a) {
	return this->match(&a) != 0;
}

std::string to_string(set& a) {
	std::string s = "{";
	for(size_t i = 0; i < a.size(); i++) {
		s += to_string(&a.members[i]);
		if(i < a.size() - 1) {
			s += ", ";
		}
	}

	s += "}";

	return s;
}

std::string to_string(set* a) {
	std::string s = "{";
	for(size_t i = 0; i < a->size(); i++) {
		s += to_string(&a->members[i]);

		if(i < a->size() - 1) {
			s += ", ";
		}
	}

	s += "}";

	return s;
}

} // namespace ast_teste
