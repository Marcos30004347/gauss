#include "AST3.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <math.h>
#include <random>
#include <string>



namespace ast_teste {

ast *ast_create(ast::kind kind){
  ast *a = (ast *)malloc(sizeof(ast));

	a->ast_kind = kind;

  a->ast_childs =
      (ast **)(calloc(sizeof(void (ast::*)()), ast::childs_margin + 1));

  a->ast_format = ast::default_format;

  a->ast_size = 0;

  a->ast_reserved_size = ast::childs_margin + 1;

	return a;
}

ast *ast_create(ast::kind kind, std::initializer_list<ast *> l) {
  ast *a = (ast *)malloc(sizeof(ast));

  a->ast_kind = kind;

  a->ast_format = ast::default_format;

  a->ast_size = l.size();

  a->ast_reserved_size = a->ast_size;

  a->ast_childs = (ast **)(malloc(sizeof(void (ast::*)()) * l.size()));

  memcpy(a->ast_childs, l.begin(), (size_t)l.end() - (size_t)l.begin());

  return a;
}

void ast_set_kind(ast *a, ast::kind kind) { a->ast_kind = kind; }

void ast_set_operand(ast *a, ast *v, size_t i) { a->ast_childs[i] = v; }

void ast_set_size(ast *a, size_t s) {
  assert(a->ast_reserved_size >= s);
  a->ast_size = s;
}

void ast_delete_operands(ast *a) {
  for (size_t i = 0; i < ast_size(a); i++) {
    ast *op = ast_operand(a, i);
    if (op)
      ast_delete(op);
  }

  if (a->ast_childs) {
    a->ast_size = 0;
    a->ast_reserved_size = 0;

    free(a->ast_childs);
  }
}

void ast_delete(ast *a) {
  if (a == nullptr)
    return;

  ast_delete_operands(a);

  if (ast_is_kind(a, ast::symbol)) {
    free(a->ast_sym);
  }

  if (ast_is_kind(a, ast::integer)) {
    delete a->ast_int;
  }

  free(a);
}

ast *ast_symbol(const char *id) {
  ast *a = ast_create(ast::symbol);

  a->ast_sym = strdup(id);

  return a;
}

ast *ast_integer(Int value) {
  ast *a = ast_create(ast::integer);

  a->ast_int = new Int(value);

  return a;
}

ast *ast_fraction(Int num, Int den) {
  return ast_create(ast::fraction, {ast_integer(num), ast_integer(den)});
}

ast *ast_operand(ast *a, size_t i) { return a->ast_childs[i]; }

void ast_alocate_operands(ast *a, size_t at_least) {
  if (ast_size(a)) {
    ast_delete_operands(a);
  }

  a->ast_reserved_size = at_least + ast::childs_margin;

  a->ast_childs =
      (ast **)calloc(sizeof(void (ast::*)()), at_least + ast::childs_margin);
}

void ast_insert(ast *a, ast *b, size_t idx) {
  assert(a->ast_size >= idx);

  const size_t ast_size = sizeof(void (ast::*)());

  if (a->ast_reserved_size == a->ast_size) {

    ast **childs = (ast **)calloc(ast_size, (a->ast_size + ast::childs_margin));

    a->ast_reserved_size = a->ast_size + ast::childs_margin;

    memcpy(childs, a->ast_childs, ast_size * idx);

    childs[idx] = b;

    if (idx > a->ast_size) {
      memcpy(childs + idx + 1, a->ast_childs + idx,
             ast_size * (a->ast_size - idx));
    }

    free(a->ast_childs);

    a->ast_childs = childs;
  } else {
    if (idx == a->ast_size) {
      a->ast_childs[idx] = b;
    } else {
      memcpy(a->ast_childs + idx + 1, a->ast_childs + idx,
             ast_size * (a->ast_size - idx));
      a->ast_childs[idx] = b;
    }
  }

  a->ast_size += 1;
}

void ast_remove(ast *a, size_t idx) {
  if (ast_operand(a, idx)) {
    ast_delete(ast_operand(a, idx));
  }

  if (idx != a->ast_size - 1) {
    memcpy(a->ast_childs + idx, a->ast_childs + idx + 1,
           sizeof(void (ast::*)()) * (a->ast_size - idx - 1));
  }

  a->ast_size -= 1;

  if (a->ast_size - a->ast_reserved_size > ast::childs_margin) {
    ast **childs = (ast **)calloc(sizeof(void (ast::*)()),
                                  (a->ast_size + ast::childs_margin));

    memcpy(childs, a->ast_childs, sizeof(void (ast::*)()) * a->ast_size);

    free(a->ast_childs);

    a->ast_childs = childs;
  }
}

int ast_cmp_consts(ast *a, ast *b) {
  assert(ast_is_kind(a, ast::constant) && ast_is_kind(b, ast::constant));

  if (ast_is_kind(a, ast::integer) && ast_is_kind(b, ast::integer)) {
    if (ast_value(a) == ast_value(b)) {
      return 0;
    }

    return ast_value(a) > ast_value(b) ? 1 : -1;
  }

  if (ast_is_kind(a, ast::fraction) && ast_is_kind(b, ast::fraction)) {
    Int na = ast_value(a->ast_childs[0]);
    Int da = ast_value(a->ast_childs[1]);
    Int nb = ast_value(b->ast_childs[0]);
    Int db = ast_value(b->ast_childs[1]);

    if (na * db == nb * da)
      return 0;

    return na * db - nb * da > 0 ? 1 : -1;
  }

  if (ast_is_kind(a, ast::integer) && ast_is_kind(b, ast::fraction)) {
    Int na = ast_value(a);
    Int nb = ast_value(b->ast_childs[0]);
    Int db = ast_value(b->ast_childs[1]);

    Int ct = na * db;

    if (ct == nb) {
      return 0;
    }

    return ct > nb ? 1 : -1;
  }

  Int nb = ast_value(b);
  Int na = ast_value(a->ast_childs[0]);
  Int da = ast_value(a->ast_childs[1]);

  Int ct = nb * da;

  if (ct == na) {
    return 0;
  }

  return na > ct ? 1 : -1;
}

bool should_revert_idx(ast::kind ctx) { return ctx & (ast::add); }

inline int ast_op_cmp(ast *a, ast *b, ast::kind ctx) {
  long m = ast_size(a);
  long n = ast_size(b);

  long l = std::min(ast_size(a), ast_size(b));

  m = m - 1;
  n = n - 1;

  if (ast_is_kind(a, ast::constant) && ast_is_kind(b, ast::constant)) {
    return ast_cmp_consts(a, b);
  }

  if (ctx == ast::add) {

    if (ast_is_kind(a, ast::mul) && ast_is_kind(b, ast::mul)) {

      if (std::abs(m - n) > 1) {
        return n - m;
      }

      for (long i = 0; i < l; i++) {
        int order =
            ast_kind(ast_operand(a, m - i)) - ast_kind(ast_operand(b, n - i));

        if (order)
          return order;
      }

      for (long i = 0; i < l; i++) {
        int order = ast_cmp(ast_operand(a, m - i), ast_operand(b, n - i), ctx);

        if (order)
          return order;
      }
    }
  }

  for (long i = 0; i < l; i++) {
    int order =
        ast_kind(ast_operand(b, n - i)) - ast_kind(ast_operand(a, m - i));

    if (order)
      return order;
  }

  for (long i = 0; i < l; i++) {
    int order = ast_cmp(ast_operand(a, m - i), ast_operand(b, n - i), ctx);

    if (order)
      return order;
  }

  return (ctx & ast::add) ? m - n : n - m;
}

inline int ast_cmp_idents(ast *a, ast *b) {
  return strcmp(ast_id(a), ast_id(b));
}

int ast_cmp(ast *a, ast *b, ast::kind ctx) {

  if (ctx & ast::mul) {
    if (ast_is_kind(a, ast::constant)) {
      return -1;
    }

    if (ast_is_kind(b, ast::constant)) {
      return +1;
    }

    if (ast_is_kind(a, ast::pow) && ast_is_kind(b, ast::pow)) {
      return ast_cmp(ast_operand(a, 0), ast_operand(b, 0), ctx);
    }

    if (ast_is_kind(a, ast::symbol) && ast_is_kind(b, ast::pow)) {
      int order = ast_cmp(a, ast_operand(b, 0), ast::mul);

      if (order == 0) {
        return ast_kind(a) - ast_kind(b);
      }

      return order;
    }

    if (ast_is_kind(b, ast::symbol) && ast_is_kind(a, ast::pow)) {
      int order = ast_cmp(ast_operand(a, 0), b, ast::mul);

      if (order == 0) {
        return ast_kind(a) - ast_kind(b);
      }

      return order;
    }

    if (ast_is_kind(a, ast::funcall) && ast_is_kind(b, ast::funcall)) {
      return strcmp(ast_funname(a), ast_funname(b));
    }

    if (ast_is_kind(a, ast::pow) && ast_is_kind(b, ast::funcall)) {
      return ast_cmp(ast_operand(a, 0), b, ast::mul);
    }

    if (ast_is_kind(b, ast::pow) && ast_is_kind(a, ast::funcall)) {
      return ast_cmp(b, ast_operand(a, 0), ast::mul);
    }

    if (ast_is_kind(a, ast::mul) &&
        ast_is_kind(b, ast::pow | ast::symbol | ast::funcall)) {
      long k = ast_size(a);

      for (long i = k - 1; i >= 0; i--) {
        int order = ast_cmp(ast_operand(a, i), b, ast::mul);

        if (order == 0) {
          return 0;
        }

        if (order < 0)
          break;
      }

      return -k + 1;
    }

    if (ast_is_kind(b, ast::mul) &&
        ast_is_kind(a, ast::pow | ast::symbol | ast::funcall)) {
      long k = ast_size(b);

      for (long i = k - 1; i >= 0; i--) {
        int order = ast_cmp(ast_operand(b, i), a, ast::mul);

        if (order == 0) {
          return 0;
        }

        if (order < 0)
          break;
      }

      return +k - 1;
    }
  }

  if (ctx & ast::add) {

    if (ast_is_kind(a, ast::constant) && ast_is_kind(b, ast::constant)) {
      return ast_cmp_consts(b, a);
    }

    if (ast_is_kind(a, ast::constant)) {
      return +1;
    }

    if (ast_is_kind(b, ast::constant)) {
      return -1;
    }

    if (ast_is_kind(a, ast::pow) && ast_is_kind(b, ast::pow)) {
      int i = ast_cmp(ast_operand(a, 0), ast_operand(b, 0), ctx);

      if (i != 0) {
        return i;
      }

      return ast_cmp(ast_operand(a, 1), ast_operand(b, 1), ctx);
    }

    if (ast_is_kind(a, ast::funcall) && ast_is_kind(b, ast::funcall)) {
      return strcmp(ast_funname(a), ast_funname(b));
    }

    if (ast_is_kind(a, ast::add) && ast_is_kind(b, ast::symbol)) {
      return +1;
    }

    if (ast_is_kind(a, ast::symbol) && ast_is_kind(b, ast::add)) {
      return -1;
    }

    if (ast_is_kind(a, ast::mul) && ast_is_kind(b, ast::symbol)) {
      long k = ast_size(a);

      if (k > 2)
        return -k;

      int order = ast_cmp(ast_operand(a, ast_size(a) - 1), b, ctx);

      if (order == 0) {
        return -1;
      }

      return k > 2 ? -1 : order;
    }

    if (ast_is_kind(b, ast::mul) && ast_is_kind(a, ast::symbol)) {
      long k = ast_size(b);

      if (k > 2) {
        return +k;
      }

      int order = ast_cmp(a, ast_operand(b, ast_size(b) - 1), ctx);

      if (order == 0) {
        return +1;
      }

      return k > 2 ? +1 : order;
    }

    if (ast_is_kind(a, ast::mul) && ast_is_kind(b, ast::pow)) {
      long k = ast_size(a);

      if (k > 2)
        return -k;

      int order = ast_cmp(ast_operand(a, 0), b, ctx);

      if (order == 0)
        return -1;

      return k > 2 ? -1 : order;
    }

    if (ast_is_kind(b, ast::mul) && ast_is_kind(a, ast::pow)) {
      long k = ast_size(b);

      if (k > 2)
        return +k;

      int order = ast_cmp(a, ast_operand(b, 0), ctx);

      if (order == 0)
        return +1;

      return k > 2 ? +1 : order;
    }
  }

  if (ast_is_kind(a, ast::funcall) && ast_is_kind(b, ast::funcall)) {
    return strcmp(ast_funname(a), ast_funname(b));
  }

  if (ast_is_kind(a, ast::constant) && ast_is_kind(b, ast::constant)) {
    return ast_cmp_consts(a, b);
  }

  if (ast_is_kind(a, ast::symbol) && ast_is_kind(b, ast::symbol)) {
    return ast_cmp_idents(a, b);
  }

  if (ast_is_kind(a, ast::add) && ast_is_kind(b, ast::add)) {
    return ast_op_cmp(a, b, ctx);
  }

  if (ast_is_kind(a, ast::mul) && ast_is_kind(b, ast::mul)) {
    return ast_op_cmp(a, b, ctx);
  }

  if (ast_is_kind(a, ast::pow) && ast_is_kind(b, ast::pow)) {
    return ast_cmp(ast_operand(a, 0), ast_operand(b, 0), ctx) ||
           ast_cmp(ast_operand(a, 1), ast_operand(b, 1), ctx);
  }

  if (ast_is_kind(a, ast::div) && ast_is_kind(b, ast::div)) {
    return ast_cmp(ast_operand(a, 0), ast_operand(b, 0), ctx) ||
           ast_cmp(ast_operand(a, 1), ast_operand(b, 1), ctx);
  }

  return should_revert_idx(ctx) ? ast_kind(a) - ast_kind(b)
                                : ast_kind(b) - ast_kind(a);
}

long int ast_sort_split(ast *a, long l, long r) {
  long int i = l - 1;

  ast *p = a->ast_childs[r];

  for (long int j = l; j < r; j++) {
    if (ast_cmp(a->ast_childs[j], p, ast_kind(a)) < 0) {
      i = i + 1;

      // swap i and j
      ast *t = a->ast_childs[i];

      a->ast_childs[i] = a->ast_childs[j];
      a->ast_childs[j] = t;
    }
  }

  ast *t = a->ast_childs[i + 1];

  a->ast_childs[i + 1] = a->ast_childs[r];
  a->ast_childs[r] = t;

  return i + 1;
}

void ast_sort_childs(ast *a, long int l, long int r) {
  if (l < r) {
    long int m = ast_sort_split(a, l, r);

    ast_sort_childs(a, l, m - 1);
    ast_sort_childs(a, m + 1, r);
  }
}

void ast_sort(ast *a) {
  if (ast_is_kind(a, ast::integer | ast::fraction | ast::infinity |
                         ast::negative_infinity | ast::symbol | ast::fail |
                         ast::undefined))
    return;

  size_t i = ast_is_kind(a, ast::sub) ? 1 : 0;

  for (; i < a->ast_size; i++) {
    ast_sort(ast_operand(a, i));
  }

  if (ast_is_kind(a, ast::pow | ast::div))
    return;

  ast_sort_childs(a, 0, ast_size(a) - 1);
}

ast *ast_copy(ast *a) {
  ast *b = (ast *)malloc(sizeof(ast));

  b->ast_kind = a->ast_kind;

  b->ast_size = a->ast_size;
  b->ast_reserved_size = a->ast_reserved_size;
  b->ast_format = a->ast_format;

  if (ast_is_kind(a, ast::integer)) {
    b->ast_int = new Int(*a->ast_int);
    return b;
  }

  if (ast_is_kind(a, ast::symbol)) {
    b->ast_sym = (char *)malloc(strlen(a->ast_sym));
    strcpy(b->ast_sym, a->ast_sym);
    return b;
  }

  b->ast_childs =
      (ast **)malloc(sizeof(void (ast::*)()) * b->ast_reserved_size);

  for (size_t i = 0; i < ast_size(a); i++) {
    b->ast_childs[i] = ast_copy(ast_operand(a, i));
  }

  return b;
}

ast *ast_copy_from(ast *a, size_t from) {
  assert(ast_is_kind(a, ast::add | ast::mul | ast::sub));

  ast *b = (ast *)malloc(sizeof(ast));

  b->ast_kind = a->ast_kind;

  b->ast_size = a->ast_size - 1;
  b->ast_reserved_size = a->ast_size - 1 + ast::childs_margin;
  b->ast_format = a->ast_format;

  b->ast_childs = (ast **)calloc(sizeof(void (ast::*)()), b->ast_reserved_size);

  for (size_t i = from; i < ast_size(a); i++) {
    b->ast_childs[i] = ast_copy(ast_operand(a, i));
  }

  return b;
}

std::string ast_to_string(ast *tree) {
  if (!tree)
    return "null";

  if (ast_is_kind(tree, ast::integer)) {
    return tree->ast_int->to_string();
  }

  if (ast_is_kind(tree, ast::symbol)) {
    return std::string(tree->ast_sym);
  }

  if (ast_is_kind(tree, ast::undefined)) {
    return "undefined";
  }

  if (ast_is_kind(tree, ast::fail)) {
    return "fail";
  }

  if (ast_is_kind(tree, ast::infinity)) {
    return "inf";
  }

  if (ast_is_kind(tree, ast::negative_infinity)) {
    return "-inf";
  }

  if (ast_is_kind(tree, ast::fraction)) {
    return ast_to_string(ast_operand(tree, 0)) + "/" +
           ast_to_string(ast_operand(tree, 1));
  }

  if (ast_is_kind(tree, ast::fraction)) {
    return "sqrt(" + ast_to_string(ast_operand(tree, 0)) + ")";
  }

  if (ast_is_kind(tree, ast::funcall)) {
    std::string r = std::string(ast_funname(tree)) + "(";

    if (ast_size(tree) > 0) {
      for (size_t i = 0; i < ast_size(tree) - 1; i++) {
        r += ast_to_string(ast_operand(tree, i));
        r += ",";
      }

      r += ast_to_string(ast_operand(tree, ast_size(tree) - 1));
    }

    r += ")";

    return r;
  }

  if (ast_is_kind(tree, ast::pow)) {
    std::string r = "";

    if (ast_is_kind(ast_operand(tree, 0),
                    ast::sub | ast::add | ast::mul | ast::div)) {
      r += "(";
    }

    r += ast_to_string(ast_operand(tree, 0));

    if (ast_is_kind(ast_operand(tree, 0),
                    ast::sub | ast::add | ast::mul | ast::div)) {
      r += ")";
    }

    r += "^";

    if (ast_is_kind(ast_operand(tree, 1),
                    ast::sub | ast::add | ast::mul | ast::div)) {
      r += "(";
    }

    r += ast_to_string(ast_operand(tree, 1));

    if (ast_is_kind(ast_operand(tree, 1),
                    ast::sub | ast::add | ast::mul | ast::div)) {
      r += ")";
    }

    return r;
  }

  if (ast_is_kind(tree, ast::div)) {
    std::string r = "";

    if (ast_is_kind(ast_operand(tree, 0),
                    ast::sub | ast::add | ast::mul | ast::div)) {
      r += "(";
    }

    r += ast_to_string(ast_operand(tree, 0));

    if (ast_is_kind(ast_operand(tree, 0),
                    ast::sub | ast::add | ast::mul | ast::div)) {
      r += ")";
    }

    r += " ÷ ";

    if (ast_is_kind(ast_operand(tree, 1),
                    ast::sub | ast::add | ast::mul | ast::div)) {
      r += "(";
    }

    r += ast_to_string(ast_operand(tree, 1));

    if (ast_is_kind(ast_operand(tree, 1),
                    ast::sub | ast::add | ast::mul | ast::div)) {
      r += ")";
    }

    return r;
  }

  if (ast_is_kind(tree, ast::add)) {
    std::string r = "";

    for (size_t i = 0; i < ast_size(tree); i++) {
      if (ast_operand(tree, i) &&
          ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {

        r += "(";
      }

      r += ast_to_string(ast_operand(tree, i));

      if (ast_operand(tree, i) &&
          ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
        r += ")";
      }

      if (i < ast_size(tree) - 1) {
        r += " + ";
      }
    }

    return r;
  }

  if (ast_is_kind(tree, ast::sub)) {
    std::string r = "";

    for (size_t i = 0; i < ast_size(tree) - 1; i++) {
      if (ast_operand(tree, i) &&
          ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
        r += "(";
      }

      r += ast_to_string(ast_operand(tree, i));

      if (ast_operand(tree, i) &&
          ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
        r += ")";
      }

      if (i != ast_size(tree) - 1) {
        r += " - ";
      }
    }

    return r;
  }

  if (ast_is_kind(tree, ast::mul)) {
    std::string r = "";

    for (size_t i = 0; i < ast_size(tree); i++) {

      if (ast_operand(tree, i) == nullptr) {
        continue;
      }

      if (ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
        r += "(";
      }

      r += ast_to_string(ast_operand(tree, i));

      if (ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
        r += ")";
      }

      if (i < ast_size(tree) - 1) {
        r += "⋅";
      }
    }

    return r;
  }

  if (ast_is_kind(tree, ast::fact)) {
    return ast_to_string(ast_operand(tree, 0)) + "!";
  }

  return "to_string_not_implemented";
}

inline ast *ast_stole_operand(ast *u, size_t j, ast *v, size_t i) {
  ast_insert(u, v->ast_childs[i], j);

  v->ast_childs[i] = nullptr;

  return u;
}

inline ast *eval_add_consts(ast *u, size_t i, ast *v, size_t j) {
  ast *a = ast_operand(u, i);
  ast *b = ast_operand(v, j);

  if (a == nullptr || b == nullptr) {
    return u;
  }

  assert(ast_is_kind(a, ast::integer | ast::fraction));
  assert(ast_is_kind(b, ast::integer | ast::fraction));

  if (ast_is_kind(a, ast::integer) && ast_is_kind(b, ast::integer)) {

    Int t = ast_value(a) + ast_value(b);

    if (t == 0) {
      ast_remove(u, i);
    } else {

      delete a->ast_int;

      a->ast_int = new Int(t);
    }

    return u;
  }

  if (ast_is_kind(a, ast::fraction) && ast_is_kind(b, ast::fraction)) {
    Int num_a = ast_value(ast_operand(a, 0));
    Int den_a = ast_value(ast_operand(a, 1));
    Int num_b = ast_value(ast_operand(b, 0));
    Int den_b = ast_value(ast_operand(b, 1));

    Int num = num_a * den_b + num_b * den_a;
    Int den = den_b * den_a;

    if (num == 0) {
      ast_remove(u, i);
      return u;
    }

    Int cff = abs(gcd(num, den));

    if (den / cff == 1) {
      ast_delete(ast_operand(a, 0));
      ast_delete(ast_operand(a, 1));

      ast_set_size(a, 0);
      ast_set_kind(a, ast::integer);

      a->ast_int = new Int(num / cff);

      return u;
    }

    ast_delete(ast_operand(a, 0));
    ast_delete(ast_operand(a, 1));

    ast_set_operand(a, ast_integer(num / cff), 0);
    ast_set_operand(a, ast_integer(den / cff), 1);

    return u;
  }

  ast *c = ast_is_kind(a, ast::integer) ? a : b;
  ast *f = ast_is_kind(a, ast::integer) ? b : a;

  Int val = ast_value(c);

  Int num = ast_value(ast_operand(f, 0));
  Int den = ast_value(ast_operand(f, 1));

  num = val * den + num;

  if (num == 0) {
    ast_remove(u, i);
    return u;
  }

  Int cff = abs(gcd(num, den));

  if (den / cff == 1) {
    ast_delete(ast_operand(a, 0));
    ast_delete(ast_operand(a, 1));

    ast_set_size(a, 0);
    ast_set_kind(a, ast::integer);

    a->ast_int = new Int(num / cff);

    return u;
  }

  ast_delete(ast_operand(a, 0));
  ast_delete(ast_operand(a, 1));

  ast_set_operand(a, ast_integer(num / cff), 0);
  ast_set_operand(a, ast_integer(den / cff), 1);

  return u;
}

inline ast *eval_mul_consts(ast *u, size_t i, ast *v, size_t j) {
  ast *a = ast_operand(u, i);
  ast *b = ast_operand(v, j);

  if (a == nullptr || b == nullptr) {
    return u;
  }

  assert(ast_is_kind(a, ast::integer | ast::fraction));
  assert(ast_is_kind(b, ast::integer | ast::fraction));

  if (ast_is_kind(a, ast::integer) && ast_is_kind(b, ast::integer)) {

    Int t = ast_value(a) * ast_value(b);

    if (t == 0) {
      ast_remove(u, i);
    } else {
      delete a->ast_int;
      a->ast_int = new Int(t);
    }

    return u;
  }

  if (ast_is_kind(a, ast::fraction) && ast_is_kind(b, ast::fraction)) {
    Int num_a = ast_value(ast_operand(a, 0));
    Int den_a = ast_value(ast_operand(a, 1));
    Int num_b = ast_value(ast_operand(b, 0));
    Int den_b = ast_value(ast_operand(b, 1));

    Int num = num_a * num_b;
    Int den = den_b * den_a;

    if (num == 0) {
      ast_remove(u, i);
      return u;
    }

    Int cff = abs(gcd(num, den));

    if (den / cff == 1) {
      ast_delete(ast_operand(a, 0));
      ast_delete(ast_operand(a, 1));

      ast_set_size(a, 0);
      ast_set_kind(a, ast::integer);

      a->ast_int = new Int(num / cff);

      return u;
    }

    ast_delete(ast_operand(a, 0));
    ast_delete(ast_operand(a, 1));

    ast_set_operand(a, ast_integer(num / cff), 0);
    ast_set_operand(a, ast_integer(den / cff), 1);

    return u;
  }

  ast *c = ast_is_kind(a, ast::integer) ? a : b;
  ast *f = ast_is_kind(a, ast::integer) ? b : a;

  Int val = ast_value(c);

  Int num = ast_value(ast_operand(f, 0));
  Int den = ast_value(ast_operand(f, 1));

  num = val * num;

  if (num == 0) {
    ast_remove(u, i);
    return u;
  }

  Int cff = abs(gcd(num, den));

  if (den / cff == 1) {
    ast_delete(ast_operand(a, 0));
    ast_delete(ast_operand(a, 1));

    ast_set_size(a, 0);
    ast_set_kind(a, ast::integer);

    a->ast_int = new Int(num / cff);

    return u;
  }

  ast_delete(ast_operand(a, 0));
  ast_delete(ast_operand(a, 1));

  ast_set_operand(a, ast_integer(num / cff), 0);
  ast_set_operand(a, ast_integer(den / cff), 1);

  return u;
}

inline ast *eval_add_int(ast *a, Int b) {
  if (ast_is_kind(a, ast::integer)) {
    Int *c = new Int(Int(*a->ast_int) + b);

    delete a->ast_int;

    a->ast_int = c;

    return a;
  }

  assert(ast_is_kind(a, ast::fraction));

  Int num = ast_value(ast_operand(a, 0));
  Int den = ast_value(ast_operand(a, 1));

  num = b * den + num;

  Int cff = abs(gcd(num, den));

  if (den / cff == 1) {
    ast_delete_operands(a);

    a->ast_kind = ast::integer;
    a->ast_int = new Int(num / cff);

    return a;
  }

  delete a->ast_childs[0];
  delete a->ast_childs[1];

  a->ast_childs[0] = ast_integer(num / cff);
  a->ast_childs[1] = ast_integer(den / cff);

  return a;
}

inline bool eval_add_nconst(ast *u, size_t i, ast *v, size_t j) {
  assert(ast_is_kind(u, ast::add) && ast_is_kind(v, ast::add));
  assert(ast_is_kind(ast_operand(u, i), ast::summable));
  assert(ast_is_kind(ast_operand(v, j), ast::summable));

  ast *a = ast_operand(u, i);
  ast *b = ast_operand(v, j);

  if (a == nullptr || b == nullptr) {
    return false;
  }

  if (ast_is_kind(a, ast::pow) && ast_is_kind(b, ast::symbol)) {
    return false;
  }

  if (ast_is_kind(a, ast::symbol) && ast_is_kind(b, ast::pow)) {
    return false;
  }

  int kind = ast_kind(a) & ast_kind(b);

  // both are a sym or both are a pow
  if (kind & (ast::symbol | ast::pow)) {
    if (ast_cmp(a, b, ast::add) == 0) {
      ast *c = ast_integer(2);

      ast_set_kind(a, ast::mul);

      ast_set_size(a, 2);

      ast_set_operand(a, c, 0);
      ast_set_operand(a, b, 1);

      ast_set_operand(v, nullptr, j);

      return true;
    }

    return false;
  }

  long size_a = ast_size(a);

  // both are mul
  if (kind & ast::mul) {
    long size_b = ast_size(b);

    ast *c;

    long size_c = -1;

    if (std::abs(size_a - size_b) > 1) {
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

    if (!ast_is_kind(ast_operand(a, 0), ast::constant)) {
      return false;
    }

    long start_a = ast_is_kind(ast_operand(a, 0), ast::integer) ? 1 : 0;
    long start_b = ast_is_kind(ast_operand(b, 0), ast::integer) ? 1 : 0;

    start_a = start_a + (size_b > size_a ? 1 : 0);

    for (long x = 0; x < size_b - start_b; x++) {
      if (ast_cmp(ast_operand(a, x + start_a), ast_operand(b, x + start_b),
                  ast::mul) != 0) {
        return false;
      }
    }

    if (ast_is_kind(ast_operand(a, 0), ast::constant) &&
        ast_is_kind(ast_operand(b, 0), ast::constant)) {
      eval_add_consts(a, 0, b, 0);
    } else if (ast_is_kind(ast_operand(a, 0), ast::constant)) {
      eval_add_int(ast_operand(a, 0), 1);
    } else {
      c = ast_integer(2);

      ast_set_size(a, 2);
      ast_set_kind(a, ast::mul);

      ast_set_operand(a, c, 0);
      ast_set_operand(a, b, 1);

      ast_set_operand(v, nullptr, j);
    }

    if (size_c != -1) {
      ast_set_operand(u, a, i);
      ast_set_operand(v, b, j);
    }

    return true;
  }

  if (ast_is_kind(b, ast::mul)) {
    return false;
  }

  // a is a mul and b is a sym or a pow
  assert(ast_is_kind(a, ast::mul));
  assert(ast_is_kind(b, ast::symbol | ast::pow));

  if (size_a > 2) {
    return false;
  }

  ast *a0 = ast_operand(a, 0);
  ast *a1 = ast_operand(a, 1);

  long ki = ast_kind(a1) & ast_kind(b);

  if (!ast_is_kind(a0, ast::constant) || !ki) {
    return false;
  }

  if (ast_cmp(b, a1, ast::mul) == 0) {
    eval_add_int(a0, 1);
  } else {
    return false;
  }

  return true;
}

inline bool eval_mul_nconst(ast *u, size_t i, ast *v, size_t j) {
  assert(ast_is_kind(u, ast::mul) && ast_is_kind(v, ast::mul));

  assert(ast_is_kind(ast_operand(u, i), ast::multiplicable));
  assert(ast_is_kind(ast_operand(v, j), ast::multiplicable));

  ast *a = ast_operand(u, i);
  ast *b = ast_operand(v, j);

  if (a == nullptr || b == nullptr) {
    return false;
  }

  if (ast_is_kind(a, ast::constant) || ast_is_kind(b, ast::constant)) {
    return false;
  }

  if (ast_is_kind(a, ast::add) || ast_is_kind(b, ast::add)) {
    return false;
  }

  if (ast_is_kind(a, ast::pow) && ast_is_kind(b, ast::pow)) {
    if (ast_cmp(ast_operand(a, 0), ast_operand(b, 0), ast::mul) == 0) {

      ast_set_operand(
          a, ast_create(ast::add, {ast_operand(a, 1), ast_operand(b, 1)}), 1);
      ast_set_operand(a, ast_eval(ast_operand(a, 1)), 1);

      if (ast_is_kind(ast_operand(a, 1), ast::integer) &&
          ast_value(ast_operand(a, 1)) == 0) {
        ast_remove(u, i);
      }

      return true;
    }

    return false;
  }

  if (ast_is_kind(a, ast::symbol) && ast_is_kind(b, ast::symbol)) {
    if (ast_cmp(a, b, ast::mul) == 0) {
      ast_set_operand(u, ast_create(ast::pow, {a, ast_integer(2)}), i);

      return true;
    }

    return false;
  }

  if (ast_is_kind(a, ast::funcall) && ast_is_kind(b, ast::funcall)) {
    if (ast_cmp(a, b, ast::mul) == 0) {
      ast_set_operand(u, ast_create(ast::pow, {a, ast_integer(2)}), i);

      return true;
    }

    return false;
  }

  if (ast_is_kind(a, ast::pow) && ast_is_kind(b, ast::symbol | ast::funcall)) {
    if (ast_cmp(ast_operand(a, 0), b, ast::mul) == 0) {

      ast_set_operand(
          a, ast_create(ast::add, {ast_operand(a, 1), ast_integer(1)}), 1);
      ast_set_operand(a, ast_eval(ast_operand(a, 1)), 1);

      if (ast_is_kind(ast_operand(a, 1), ast::integer) &&
          ast_value(ast_operand(a, 1)) == 0) {
        ast_remove(u, i);
      }

      return true;
    }

    return false;
  }

  return false;
}

inline bool eval_add_add(ast *a, ast *b) {
  assert(ast_is_kind(a, ast::add));
  assert(ast_is_kind(b, ast::add));

  size_t j = 0;
  size_t i = 0;

  ast *u = 0;
  ast *v = 0;

  while (i < ast_size(a) && j < ast_size(b)) {
    assert(!ast_is_kind(ast_operand(b, j), ast::add));

    u = ast_operand(a, i);
    v = ast_operand(b, j);

    if (u == 0) {
      i++;
      continue;
    }

    if (v == 0) {
      j++;
      continue;
    }

    bool cond = false;

    if (ast_is_kind(u, ast::constant) && ast_is_kind(v, ast::constant)) {
      cond = eval_add_consts(a, i, b, j);

      if (cond) {
        ast_set_operand(b, 0, j);
      }
    } else if (ast_is_kind(u, ast::summable) && ast_is_kind(v, ast::summable)) {
      cond = eval_add_nconst(a, i, b, j);

      if (cond) {
        ast_set_operand(b, 0, j);
      }
    }

    if (cond) {
      i = i + 1;
      j = j + 1;
    } else {
      int order = ast_cmp(u, v, ast::add);

      if (order < 0) {
        i = i + 1;
      } else {
        ast_insert(a, v, i++);
        ast_set_operand(b, 0, j++);
      }
    }
  }

  while (j < ast_size(b)) {
    if (i >= ast_size(a)) {
      ast_insert(a, ast_operand(b, j++), ast_size(a));
    } else {
      u = ast_operand(a, i);
      v = ast_operand(b, j);

      if (v == 0) {
        j++;
      } else if (u == 0) {
        i++;
      } else if (ast_cmp(u, v, ast::add) < 0) {
        i++;
      } else {
        ast_insert(a, v, i++);
        ast_set_operand(b, 0, j++);
      }
    }
  }

  return true;
}

inline bool eval_mul_mul(ast *a, ast *b) {
  assert(ast_is_kind(a, ast::mul));
  assert(ast_is_kind(b, ast::mul));

  size_t j = 0;
  size_t i = 0;

  ast *u = 0;
  ast *v = 0;

  while (i < ast_size(a) && j < ast_size(b)) {
    assert(!ast_is_kind(ast_operand(b, j), ast::mul));

    u = ast_operand(a, i);
    v = ast_operand(b, j);

    if (u == 0) {
      i++;
      continue;
    }

    if (v == 0) {
      j++;
      continue;
    }

    bool cond = false;

    if (ast_is_kind(u, ast::constant) && ast_is_kind(v, ast::constant)) {
      cond = eval_mul_consts(a, i, b, j);

      if (cond) {
        ast_set_operand(b, 0, j);
      }
    } else if (ast_is_kind(u, ast::summable) && ast_is_kind(v, ast::summable)) {
      cond = eval_mul_nconst(a, i, b, j);

      if (cond) {
        ast_set_operand(b, 0, j);
      }
    }

    if (cond) {
      i = i + 1;
      j = j + 1;
    } else {
      int order = ast_cmp(u, v, ast::mul);

      if (order < 0) {
        i = i + 1;
      } else {
        ast_insert(a, v, i++);
        ast_set_operand(b, 0, j++);
      }
    }
  }

  while (j < ast_size(b)) {
    if (i >= ast_size(a)) {
      ast_insert(a, ast_operand(b, j++), ast_size(a));
    } else {
      u = ast_operand(a, i);
      v = ast_operand(b, j);

      if (v == 0) {
        j++;
      } else if (u == 0) {
        i++;
      } else if (ast_cmp(u, v, ast::mul) < 0) {
        i++;
      } else {
        ast_insert(a, v, i++);
        ast_set_operand(b, 0, j++);
      }
    }
  }

  return true;
}

bool eval_mul_int(ast *u, size_t i, Int v) {
  ast *a = ast_operand(u, i);

  if (ast_is_kind(a, ast::integer)) {
    Int c = ast_value(a);

    delete a->ast_int;

    a->ast_int = new Int(c * v);

    return true;
  }

  if (ast_is_kind(a, ast::integer)) {
    return eval_mul_int(a, 0, v);
  }

  if (ast_is_kind(a, ast::add | ast::sub)) {
    for (size_t i = 0; i < ast_size(a); i++) {
      eval_mul_int(a, i, v);
    }

    return true;
  }

  if (ast_is_kind(a, ast::mul)) {
    if (ast_is_kind(ast_operand(a, 0), ast::constant)) {
      return eval_mul_int(a, 0, v);
    }

    ast_insert(a, ast_integer(v), 0);

    return true;
  }

  if (ast_is_kind(a, ast::div)) {
    return eval_mul_int(a, 0, v);
  }

  if (ast_is_kind(a, ast::sqrt | ast::pow | ast::fact | ast::funcall |
                         ast::symbol)) {
    ast *p = ast_create(ast::mul, {ast_integer(v), a});
    ast_set_operand(u, p, i);
    return true;
  }
  assert(ast_is_kind(a, ast::undefined | ast::fail));

  return true;
}

ast *ast_raise_to_first_op(ast *a) {
  if (ast_is_kind(ast_operand(a, 0), ast::integer)) {
    a->ast_int = new Int(ast_value(ast_operand(a, 0)));

    ast_set_kind(a, ast::integer);

    ast_delete(ast_operand(a, 0));

    ast_set_operand(a, 0, 0);

    ast_set_size(a, 0);

    return a;
  }

  if (ast_is_kind(ast_operand(a, 0), ast::symbol)) {
    a->ast_sym = strdup(ast_id(a));

    ast_set_kind(a, ast::symbol);

    ast_delete(ast_operand(a, 0));

    ast_set_operand(a, 0, 0);

    ast_set_size(a, 0);

    return a;
  }

  if (ast_is_kind(ast_operand(a, 0), ast::funcall)) {
    a->ast_sym = strdup(ast_id(a));

    ast_set_kind(a, ast::funcall);

    ast *t = ast_operand(a, 0);

    ast **t_childs = a->ast_childs;
    size_t t_size = a->ast_size;
    size_t t_rsize = a->ast_reserved_size;

    a->ast_childs = t->ast_childs;
    a->ast_size = t->ast_size;
    a->ast_reserved_size = t->ast_reserved_size;

    t->ast_childs = t_childs;
    t->ast_size = t_size;
    t->ast_reserved_size = t_rsize;

    ast_delete(t);

    return a;
  }

  ast *t = ast_operand(a, 0);

  ast_set_kind(a, ast_kind(t));

  ast **t_childs = a->ast_childs;
  size_t t_size = a->ast_size;
  size_t t_rsize = a->ast_reserved_size;

  a->ast_childs = t->ast_childs;
  a->ast_size = t->ast_size;
  a->ast_reserved_size = t->ast_reserved_size;

  t->ast_childs = t_childs;
  t->ast_size = t_size;
  t->ast_reserved_size = t_rsize;

  ast_delete(t);

  return a;
}

	ast *ast_eval(ast *a, bool print, ast *parent) {
  if (a == nullptr)
    return nullptr;

  if (ast_is_kind(a, ast::add)) {

    for (size_t i = 0; i < ast_size(a); i++) {
      a->ast_childs[i] = ast_eval(a->ast_childs[i], print, parent ? parent : a);
    }

    ast_sort_childs(a, 0, ast_size(a) - 1);

    size_t j = 0;

    for (long i = 1; i < (long)ast_size(a); i++) {
      ast *aj = ast_operand(a, j);
      ast *ai = ast_operand(a, i);

      bool cond = false;

      if (ast_is_kind(ai, ast::add)) {
        ast_set_operand(a, 0, i);
        ast_remove(a, i);
        cond = eval_add_add(a, ai);
        i--;
      } else if (ast_is_kind(aj, ast::constant) &&
                 ast_is_kind(ai, ast::constant)) {

        cond = eval_add_consts(a, j, a, i);

        if (cond) {
          ast_remove(a, i--);
        }
      } else if (ast_is_kind(aj, ast::summable) &&
                 ast_is_kind(ai, ast::summable)) {
        cond = eval_add_nconst(a, j, a, i);

        if (cond) {
          ast_remove(a, i--);
        }
      }

      if (cond == false) {
        j = i;
      }

      // TODO: code like this could be used to return the history of operations
      if (cond && print) {
        printf("%s\n", ast_to_string(parent ? parent : a).c_str());
      }
    }

    if (ast_size(a) == 1) {
      a = ast_raise_to_first_op(a);
    }

    return a;
  }

  if (ast_is_kind(a, ast::mul)) {

    for (size_t i = 0; i < ast_size(a); i++) {
      a->ast_childs[i] = ast_eval(a->ast_childs[i], print, parent ? parent : a);
    }

    ast_sort_childs(a, 0, ast_size(a) - 1);

    size_t j = 0;

    for (long i = 1; i < (long)ast_size(a); i++) {
      ast *aj = ast_operand(a, j);
      ast *ai = ast_operand(a, i);

      bool cond = false;

      if (ast_is_kind(ai, ast::mul)) {
        ast_set_operand(a, 0, i);
        ast_remove(a, i);
        cond = eval_mul_mul(a, ai);
        i--;
      } else if (ast_is_kind(aj, ast::constant) &&
                 ast_is_kind(ai, ast::constant)) {
        cond = eval_mul_consts(a, j, a, i);
        if (cond)
          ast_remove(a, i--);
      } else if (ast_is_kind(aj, ast::multiplicable) &&
                 ast_is_kind(ai, ast::multiplicable)) {
        cond = eval_mul_nconst(a, j, a, i);
        if (cond)
          ast_remove(a, i--);
      }

      if (cond == false) {
        j = i;
      }

      // TODO: code like this could be used to return the history of operations
      if (cond && print) {
        printf("%s\n", ast_to_string(parent ? parent : a).c_str());
      }
    }

    if (ast_size(a) == 1) {
      ast_raise_to_first_op(a);
    }

    return a;
  }

  if (ast_is_kind(a, ast::sub)) {
    for (size_t i = 1; i < ast_size(a); i++) {
      eval_mul_int(a, i, -1);
    }

    ast_set_kind(a, ast::add);

    return ast_eval(a, print, parent ? parent : a);
  }

  return a;
}

} // namespace ast_teste