#include "AST3.hpp"

#include <initializer_list>
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <random>
#include <string>
#include <cmath>

namespace ast_teste {

ast *ast_create(ast::kind kind) {
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

void ast_set_operand(ast *a, ast *v, size_t i) {
  // if (a->ast_childs[i])
  //   ast_delete(a->ast_childs[i]);
  a->ast_childs[i] = v;
}

void ast_set_size(ast* a, size_t s) {
	assert(a->ast_reserved_size >= s);
	a->ast_size = s;
}

void ast_delete_operands(ast *a) {
  for (size_t i = 0; i < a->ast_reserved_size; i++) {
    ast *op = ast_operand(a, i);
    if (op) ast_delete(op);
  }

  if (a->ast_childs) {
    a->ast_size = 0;
    a->ast_reserved_size = 0;

    free(a->ast_childs);
  }
}

void ast_delete(ast *a) {
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
  a->ast_sym = (char *)malloc(strlen(id));
  strcpy(a->ast_sym, id);

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

    ast **childs =
			(ast **)calloc(ast_size, (a->ast_size + ast::childs_margin));

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
  ast *b = a->ast_childs[idx];

  ast_delete(b);

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
      int order = ast_cmp(a, ast_operand(b, 0), ctx);

      if (order == 0) {
        return -1;
      }

      return order;
    }

    if (ast_is_kind(b, ast::symbol) && ast_is_kind(a, ast::pow)) {
      int order = ast_cmp(ast_operand(a, 0), b, ctx);

      if (order == 0) {
        return +1;
      }

      return order;
    }

    if (ast_is_kind(a, ast::mul) && ast_is_kind(b, ast::pow | ast::symbol)) {
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

    if (ast_is_kind(b, ast::mul) && ast_is_kind(a, ast::pow | ast::symbol)) {
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
    if (ast_is_kind(a, ast::constant)) {
      return +1;
    }

    if (ast_is_kind(b, ast::constant)) {
      return -1;
    }

    if (ast_is_kind(a, ast::pow) && ast_is_kind(b, ast::pow)) {
      int i = ast_cmp(ast_operand(a, 0), ast_operand(b, 0), ctx);

      if (i != 0)
        return i;

      return ast_cmp(ast_operand(a, 1), ast_operand(b, 1), ctx);
    }

    if (ast_is_kind(a, ast::add) && ast_is_kind(b, ast::symbol)) {
      return +1;
    }

    if (ast_is_kind(a, ast::symbol) && ast_is_kind(b, ast::add)) {
      return -1;
    }

    if (ast_is_kind(a, ast::mul) && ast_is_kind(b, ast::symbol | ast::pow)) {
      long k = ast_size(a);

			if (k > 2) return -k;

			int order = ast_cmp(ast_operand(a, ast_size(a) - 1), b, ctx);

      if (order == 0)
        return -1;

      return k > 2 ? -1 : order;
    }

    if (ast_is_kind(b, ast::mul) && ast_is_kind(a, ast::symbol | ast::pow)) {
      long k = ast_size(b);

			if (k > 2) return +k;

      int order = ast_cmp(a, ast_operand(b, ast_size(b) - 1), ctx);

      if (order == 0)
        return +1;

      return k > 2 ? +1 : order;
    }
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
    return ast_cmp(ast_operand(a, 1), ast_operand(b, 1), ctx) ||
           ast_cmp(ast_operand(a, 0), ast_operand(b, 0), ctx);
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

  if (ast_is_kind(a, ast::pow | ast::div | ast::integral | ast::derivative))
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

  b->ast_childs =
      (ast **)malloc(sizeof(void (ast::*)()) * b->ast_reserved_size);

  for (size_t i = from; i < ast_size(a); i++) {
    b->ast_childs[i] = ast_copy(ast_operand(a, i));
  }

  return b;
}

std::string ast_to_string(ast *tree) {
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
      if (ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
        r += "(";
      }

      r += ast_to_string(ast_operand(tree, i));

      if (ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
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
      if (ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
        r += "(";
      }

      r += ast_to_string(ast_operand(tree, i));

      if (ast_is_kind(ast_operand(tree, i), ast::sub | ast::add | ast::mul)) {
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

// ast *eval_add_consts(ast *a, ast *b) {
//   if (ast_is_kind(a, ast::integer) && ast_is_kind(b, ast::integer)) {
//     return ast_integer(ast_value(a) + ast_value(b));
//   }

//   if (ast_is_kind(a, ast::fraction) && ast_is_kind(b, ast::fraction)) {
//     Int num_a = ast_value(ast_operand(a, 0));
//     Int den_a = ast_value(ast_operand(a, 1));
//     Int num_b = ast_value(ast_operand(b, 0));
//     Int den_b = ast_value(ast_operand(b, 1));

//     Int num = num_a * den_b + num_b * den_a;
//     Int den = den_b * den_a;

//     Int cff = abs(gcd(num, den));

//     if (den / cff == 1) {
//       return ast_integer(num / cff);
//     }

//     return ast_create(ast::fraction,
//                       {ast_integer(num / cff), ast_integer(den / cff)});
//   }

//   ast *i = ast_is_kind(a, ast::integer) ? a : b;
//   ast *f = ast_is_kind(a, ast::integer) ? b : a;

//   Int val = ast_value(i);
//   Int num = ast_value(ast_operand(f, 0));
//   Int den = ast_value(ast_operand(f, 1));

//   num = val * den + num;

//   Int cff = abs(gcd(num, den));

//   if (den / cff == 1) {
//     return ast_integer(num / cff);
//   }

//   return ast_create(ast::fraction,
//                     {ast_integer(num / cff), ast_integer(den / cff)});
// }

inline void ast_stole_operand(ast *u, size_t j, ast *v, size_t i) {
  ast_insert(u, v->ast_childs[i], j);
  v->ast_childs[i] = nullptr;
}

ast *eval_add_consts(ast *u, size_t i, ast *v, size_t j) {
  ast *a = ast_operand(u, i);
  ast *b = ast_operand(v, j);

  assert(ast_is_kind(a, ast::integer | ast::fraction));
  assert(ast_is_kind(b, ast::integer | ast::fraction));

  if (ast_is_kind(a, ast::integer) && ast_is_kind(b, ast::integer)) {
    Int *c = new Int(ast_value(a) + ast_value(b));

    delete a->ast_int;

    a->ast_int = c;

    return a;
  }

  if (ast_is_kind(a, ast::fraction) && ast_is_kind(b, ast::fraction)) {
    Int num_a = ast_value(ast_operand(a, 0));
    Int den_a = ast_value(ast_operand(a, 1));
    Int num_b = ast_value(ast_operand(b, 0));
    Int den_b = ast_value(ast_operand(b, 1));

    Int num = num_a * den_b + num_b * den_a;
    Int den = den_b * den_a;

    Int cff = abs(gcd(num, den));

    if (den / cff == 1) {
      ast_delete_operands(a);

      a->ast_kind = ast::integer;
      a->ast_int = new Int(num / cff);

      return a;
    }

    ast_delete(ast_operand(a, 0));
    ast_delete(ast_operand(a, 1));

    a->ast_childs[0] = ast_integer(num / cff);
    a->ast_childs[1] = ast_integer(den / cff);

    return a;
  }

  ast *c = ast_is_kind(a, ast::integer) ? a : b;
  ast *f = ast_is_kind(a, ast::integer) ? b : a;

  Int val = ast_value(c);

  Int num = ast_value(ast_operand(f, 0));
  Int den = ast_value(ast_operand(f, 1));

  num = val * den + num;

  Int cff = abs(gcd(num, den));

  if (den / cff == 1) {
    ast_delete_operands(a);

    a->ast_kind = ast::integer;
    a->ast_int = new Int(num / cff);

    return a;
  }

  ast_delete(ast_operand(a, 0));
  ast_delete(ast_operand(a, 1));

  a->ast_childs[0] = ast_integer(num / cff);
  a->ast_childs[1] = ast_integer(den / cff);

  return a;
}

ast *eval_add_int(ast *a, Int b) {
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

ast *eval_add_nconst(ast *u, size_t i, ast *v, size_t j, bool stole) {
  assert(ast_is_kind(u, ast::add) && ast_is_kind(v, ast::add));

  assert(ast_is_kind(ast_operand(u, i), ast::mul | ast::symbol));
  assert(ast_is_kind(ast_operand(v, i), ast::mul | ast::symbol));

  ast *a = ast_operand(u, i);
  ast *b = ast_operand(v, j);

  int kind = ast_kind(a) & ast_kind(b);

  // both are a sym or both are a pow
  if (kind & (ast::symbol | ast::pow)) {
    if (ast_cmp(a, b, ast::undefined) == 0) {
      ast *c = ast_integer(2);

      ast_set_kind(a, ast::mul);

      ast_set_size(a, 2);

      ast_set_operand(a, c, 0);
      ast_set_operand(a, b, 1);

      ast_set_operand(v, nullptr, j);

      return u;
    }

    if (stole) {
      ast_stole_operand(u, i + 1, v, j);
      return u;
    }

    return nullptr;
  }

  long size_a = ast_size(a);

  // both are mul
  if (kind & ast::mul) {
    long size_b = ast_size(b);

    ast *c;

    long size_c = -1;

    if (std::abs(size_a - size_b) > 1) {
      if (stole) {
        ast_stole_operand(u, i + 1, v, j);
        return u;
      }

      return nullptr;
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
      if (stole) {
        ast_stole_operand(u, i + 1, v, j);
        return u;
      }

      return nullptr;
    }

    long start_a = ast_is_kind(ast_operand(a, 0), ast::integer) ? 1 : 0;
    long start_b = ast_is_kind(ast_operand(b, 0), ast::integer) ? 1 : 0;

    start_a = start_a + (size_b > size_a ? 1 : 0);

    for (long x = 0; x < size_b - start_b; x++) {
      if (ast_cmp(ast_operand(a, x + start_a), ast_operand(b, x + start_b),
                  ast::mul) != 0) {
        if (stole) {
          ast_stole_operand(u, i + 1, v, j);
          return u;
        }

        return nullptr;
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

    return u;
  }

  // a is a mul and b is a sym or a pow
  assert(ast_is_kind(a, ast::mul));
  assert(ast_is_kind(b, ast::symbol | ast::pow));

  if (size_a > 2) {
    if (stole) {
      ast_stole_operand(u, i + 1, v, j);
      return u;
    }

    return nullptr;
  }

  ast *a0 = ast_operand(a, 0);
  ast *a1 = ast_operand(a, 1);

  long ki = ast_kind(a1) & ast_kind(b);

  if (!ast_is_kind(a0, ast::constant) || !ki) {
    if (stole) {
      ast_stole_operand(u, i + 1, v, j);
      return u;
    }

    return nullptr;
  }

  if (ast_cmp(b, a1, ast::mul) == 0) {
    eval_add_int(a0, 1);
  }

  return u;
}

ast *eval_add_add(ast *a, ast *b) {
  assert(ast_is_kind(a, ast::add));
  assert(ast_is_kind(b, ast::add));

  size_t j = 0, i = 0;

  while (i < ast_size(a) && j < ast_size(b)) {
    assert(!ast_is_kind(ast_operand(b, j), ast::add));

    ast *u = ast_operand(a, i);
    ast *v = ast_operand(b, j);

    bool cond = false;

    if (ast_is_kind(u, ast::constant) && ast_is_kind(v, ast::constant)) {
      cond = eval_add_consts(a, i, b, j);
    } else if (ast_is_kind(u, ast::summable) && ast_is_kind(v, ast::summable)) {
      cond = eval_add_nconst(a, i, b, j, false);
    }

    if (cond) {
      i = i + 1;
      j = j + 1;
    } else {
      if (ast_cmp(u, v, ast::add) > 0) {
        i = i + 1;
      } else {
				ast_stole_operand(a, i, b, j++);
      }
    }
  }

	while(j < ast_size(b)) {
		if(i >= ast_size(a)) {
			ast_stole_operand(a, ast_size(a) - 1, b, j++);
		} else {
			if(ast_cmp(ast_operand(a, i), ast_operand(b, j), ast::add) > 0) {
				i++;
			} else {
				ast_stole_operand(a, i, b, j++);
			}
		}
	}

  return a;
}

ast *ast_eval(ast *a) {
  if (a == nullptr)
    return nullptr;

  if (ast_is_kind(a, ast::add)) {
    long k = ast_size(a);

    for (long i = 0; i < k; i++) {
      a->ast_childs[i] = ast_eval(a->ast_childs[i]);
    }

    ast_sort(a);

    ast *u = ast_create(ast::add);


		for (long i = 1; i < k; i++) {
      if (ast_size(u) == 0) {
        ast_stole_operand(u, 0, a, i);
      }

      ast *ai = a->ast_childs[i];
      ast *it = u->ast_childs[ast_size(u) - 1];

      if (ast_is_kind(ai, ast::add)) {
        eval_add_add(u, ai);
      } else if (ast_is_kind(it, ast::constant) &&
                 ast_is_kind(ai, ast::constant)) {
        eval_add_consts(u, ast_size(u) - 1, a, i);
      } else if (ast_is_kind(it, ast::summable) &&
								 ast_is_kind(ai, ast::summable)) {
        eval_add_nconst(u, ast_size(u) - 1, a, i, true);
      } else {
				ast_stole_operand(u, ast_size(u) - 1, a, i);
			}
    }
  }

  return a;
}

} // namespace ast_teste
