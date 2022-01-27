#include "AST4.hpp"
#include "Core/AST/AST.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <iterator>
#include <math.h>
#include <random>
#include <string>
#include <vector>

namespace expression {

expr::expr(expr &&other) {
  type = other.type;

  switch (type) {

	case kind::symbol: {
    if(other.id) {
			id = strdup(other.id);
		} else {
			id = 0;
		}
		return;
  }
  case kind::integer: {
    if(other.val) {
			val = new Int(*other.val);
		} else {
			val = 0;
		}
		return;
  }

	case kind::infinity: return;
	case kind::negative_infinity: return;
	case kind::undefined: return;
	case kind::fail: return;

  default: {
    ops = std::move(other.ops);
		return;
  }
  }
}

expr::expr(const expr &other) {
  type = other.type;

  switch (type) {
	case kind::symbol: {
    if(other.id) {
			id = strdup(other.id);
		} else {
			id = 0;
		}
		return;
  }
  case kind::integer: {
    if(other.val) {
			val = new Int(*other.val);
		} else {
			val = 0;
		}
		return;
  }
	case kind::infinity: return;
	case kind::negative_infinity: return;
	case kind::undefined: return;
	case kind::fail: return;

  default: {
    ops = other.ops;
		return;
  }
  }
}

expr::expr(kind k) {
 	type = k;
}

expr::expr() {
	type = kind::undefined;
}


expr& expr::operator[](size_t idx) {
	return this->ops[idx];
}

// expr expr_create(kind k, size_t reserve) {
//   expr a(k);

//   a.ops = std::vector<expr>(reserve);

//   return a;
// }

expr expr_create(kind kind) {
	expr a(kind);

	a.ops = std::vector<expr>();

  return a;
}

expr expr_create(kind kind, std::initializer_list<expr>&& l) {
  expr a;

  a.type = kind;

	a.ops = std::vector<expr>(l);

  return a;
}

expr expr::operator=(const expr &other) {
	return expr(other);
}

expr expr::operator=(expr &&other) {
	return expr(other);
}

void expr_set_kind(expr& a, kind k) { a.type = k; }

// void expr_set_operand(expr& a, expr& v, size_t i) {
//   a.ops[i] = v;
// }

expr::~expr() {
	switch (type) {
	case kind::symbol: {
		if(id) free(id);
		return;
	}
	case kind::integer: {
		if(val) delete val;
		return;
	}

	default: return;
	}
}

expr expr_symbol(const char *id) {
  expr a = expr_create(kind::symbol);

  a.id = strdup(id);

  return a;
}

expr expr_integer(Int value) {
  expr a = expr_create(kind::integer);

  a.val = new Int(value);

  return a;
}

expr expr_fraction(Int num, Int den) {
  return expr_create(kind::fraction, { expr_integer(num), expr_integer(den) });
}

void expr_insert(expr& a, expr& b, size_t idx) {
	a.ops.insert(a.ops.begin() + idx, b);
}

void expr_insert(expr& a, expr&& b, size_t idx) {
	a.ops.insert(a.ops.begin() + idx, b);
}

void expr_insert(expr& a, expr& b) {
	a.ops.push_back(b);
}

void expr_insert(expr& a, expr&& b) {
	a.ops.push_back(b);
}


void expr_remove(expr& a, size_t idx) {
	a.ops.erase(a.ops.begin() + idx);
}

int expr_cmp_consts(expr& a, expr& b) {
  assert(expr_is_kind(a, kind::constant) && expr_is_kind(b, kind::constant));

  if (expr_is_kind(a, kind::integer) && expr_is_kind(b, kind::integer)) {
    if (expr_value(a) == expr_value(b)) {
      return 0;
    }

    return expr_value(a) > expr_value(b) ? 1 : -1;
  }

  if (expr_is_kind(a, kind::fraction) && expr_is_kind(b, kind::fraction)) {
    Int na = expr_value(a[0]);
    Int da = expr_value(a[1]);
    Int nb = expr_value(b[0]);
    Int db = expr_value(b[1]);

    if (na * db == nb * da)
      return 0;

    return na * db - nb * da > 0 ? 1 : -1;
  }

  if (expr_is_kind(a, kind::integer) && expr_is_kind(b, kind::fraction)) {
    Int na = expr_value(a);
    Int nb = expr_value(b[0]);
    Int db = expr_value(b[1]);

    Int ct = na * db;

    if (ct == nb) {
      return 0;
    }

    return ct > nb ? 1 : -1;
  }

  Int nb = expr_value(b);
  Int na = expr_value(a[0]);
  Int da = expr_value(a[1]);

  Int ct = nb * da;

  if (ct == na) {
    return 0;
  }

  return na > ct ? 1 : -1;
}

inline int expr_op_cmp(expr& a, expr& b, kind ctx) {
  long m = expr_size(a);
  long n = expr_size(b);

  long l = std::min(expr_size(a), expr_size(b));

  m = m - 1;
  n = n - 1;

  if (expr_is_kind(a, kind::constant) && expr_is_kind(b, kind::constant)) {
    return expr_cmp_consts(a, b);
  }

  if (ctx == kind::add) {

    if (expr_is_kind(a, kind::mul) && expr_is_kind(b, kind::mul)) {

      if (std::abs(m - n) > 1) {
        return n - m;
      }

      for (long i = 0; i < l; i++) {
        int order =
            expr_kind(a[m - i]) - expr_kind(b[n - i]);

        if (order)
          return order;
      }

      for (long i = 0; i < l; i++) {
        int order = expr_cmp(a[m - i], b[n - i], ctx);

        if (order)
          return order;
      }
    }
  }

  for (long i = 0; i < l; i++) {
    int order =
        expr_kind(b[n - i]) - expr_kind(a[m - i]);

    if (order)
      return order;
  }

  for (long i = 0; i < l; i++) {
    int order = expr_cmp(a[m - i], b[n - i], ctx);

    if (order)
      return order;
  }

  return (ctx & kind::add) ? m - n : n - m;
}

inline int expr_cmp_idents(expr& a, expr& b) {
  return strcmp(expr_id(a), expr_id(b));
}

std::string expr_to_string(expr& tree) {
  if (expr_is_kind(tree, kind::integer)) {
    return tree.val->to_string();
  }

  if (expr_is_kind(tree, kind::symbol)) {
    return std::string(tree.id);
  }

  if (expr_is_kind(tree, kind::undefined)) {
    return "undefined";
  }

  if (expr_is_kind(tree, kind::fail)) {
    return "fail";
  }

  if (expr_is_kind(tree, kind::infinity)) {
    return "inf";
  }

  if (expr_is_kind(tree, kind::negative_infinity)) {
    return "-inf";
  }

  if (expr_is_kind(tree, kind::fraction)) {
    return expr_to_string(tree[0]) + "/" +
           expr_to_string(tree[1]);
  }

  if (expr_is_kind(tree, kind::fraction)) {
    return "sqrt(" + expr_to_string(tree[0]) + ")";
  }

  // if (expr_is_kind(tree, kind::funcall)) {
  //   std::string r = std::string(expr_funname(tree)) + "(";

  //   if (expr_size(tree) > 0) {
  //     for (size_t i = 0; i < expr_size(tree) - 1; i++) {
  //       r += expr_to_string(expr_operand(tree, i));
  //       r += ",";
  //     }

  //     r += expr_to_string(expr_operand(tree, expr_size(tree) - 1));
  //   }

  //   r += ")";

  //   return r;
  // }

  if (expr_is_kind(tree, kind::power)) {
    std::string r = "";

    if (expr_is_kind(tree[0], kind::sub | kind::add | kind::mul | kind::div)) {
      r += "(";
    }

    r += expr_to_string(tree[0]);

    if (expr_is_kind(tree[0], kind::sub | kind::add | kind::mul | kind::div)) {
      r += ")";
    }

    r += "^";

    if (expr_is_kind(tree[1], kind::sub | kind::add | kind::mul | kind::div)) {
      r += "(";
    }

    r += expr_to_string(tree[1]);

    if (expr_is_kind(tree[1], kind::sub | kind::add | kind::mul | kind::div)) {
      r += ")";
    }

    return r;
  }

  if (expr_is_kind(tree, kind::div)) {
    std::string r = "";

    if (expr_is_kind(tree[0], kind::sub | kind::add | kind::mul | kind::div)) {
      r += "(";
    }

    r += expr_to_string(tree[0]);

    if (expr_is_kind(tree[0], kind::sub | kind::add | kind::mul | kind::div)) {
      r += ")";
    }

    r += " ÷ ";

    if (expr_is_kind(tree[1], kind::sub | kind::add | kind::mul | kind::div)) {
      r += "(";
    }

    r += expr_to_string(tree[1]);

    if (expr_is_kind(tree[1], kind::sub | kind::add | kind::mul | kind::div)) {
      r += ")";
    }

    return r;
  }

  if (expr_is_kind(tree, kind::add)) {
    std::string r = "";

    for (size_t i = 0; i < expr_size(tree); i++) {
      if (expr_is_kind(tree[i], kind::sub | kind::add | kind::mul)) {

        r += "(";
      }

      r += expr_to_string(tree[i]);

      if (expr_is_kind(tree[i], kind::sub | kind::add | kind::mul)) {
        r += ")";
      }

      if (i < expr_size(tree) - 1) {
        r += " + ";
      }
    }

    return r;
  }

  if (expr_is_kind(tree, kind::sub)) {
    std::string r = "";

    for (size_t i = 0; i < expr_size(tree) - 1; i++) {
      if (expr_is_kind(tree[i], kind::sub | kind::add | kind::mul)) {
        r += "(";
      }

      r += expr_to_string(tree[i]);

      if (expr_is_kind(tree[i], kind::sub | kind::add | kind::mul)) {
        r += ")";
      }

      if (i != expr_size(tree) - 1) {
        r += " - ";
      }
    }

    return r;
  }

  if (expr_is_kind(tree, kind::mul)) {
    std::string r = "";

    for (size_t i = 0; i < expr_size(tree); i++) {

      if (expr_is_kind(tree[i], kind::sub | kind::add | kind::mul)) {
        r += "(";
      }

      r += expr_to_string(tree[i]);

      if (expr_is_kind(tree[i], kind::sub | kind::add | kind::mul)) {
        r += ")";
      }

      if (i < expr_size(tree) - 1) {
        r += "⋅";
      }
    }

    return r;
  }

  if (expr_is_kind(tree, kind::factorial)) {
    return expr_to_string(tree[0]) + "!";
  }

	assert(false);

	return "";
}

std::string expr_kind_id(expr a) {
  switch (expr_kind(a)) {

  case kind::integer: {
    return "integer";
  }
  case kind::symbol: {
    return "symbol";
  }
  // case kind::funcall: {
  //   return "funcall";
  // }
  case kind::factorial: {
    return "fact";
  }
  case kind::power: {
    return "pow";
  }
  case kind::mul: {
    return "mul";
  }
  case kind::add: {
    return "add";
  }
  case kind::sub: {
    return "div";
  }
  case kind::sqrt: {
    return "sqrt";
  }
  case kind::infinity: {
    return "infinity";
  }
  case kind::negative_infinity: {
    return "negative infinity";
  }

  case kind::undefined: {
    return "undefined";
  }

  case kind::fail: {
    return "fail";
  }

  case kind::fraction: {
    return "fraction";
  }

  case kind::div: {
    return "div";
  }

  default:
    return "";
  }
}

void expr_print(expr& a, int tabs) {
  printf("%*c<ast ", tabs, ' ');
  printf("kind=\"%s\"", expr_kind_id(a).c_str());

  if (expr_kind(a) == kind::integer) {
    printf(" value=\"%s\"", expr_value(a).to_string().c_str());
  }

  if (expr_kind(a) == kind::symbol) {
    printf(" id=\"%s\"", expr_id(a));
  }

  if (expr_size(a)) {
    printf(">\n");
    for (size_t i = 0; i < expr_size(a); i++) {
      expr_print(a[i], tabs + 3);
    }
    printf("%*c</ast>\n", tabs, ' ');
  } else {
    printf(">\n");
  }
}

int expr_cmp(expr& a, expr& b, kind ctx) {

  if (ctx & kind::mul) {
    if (expr_is_kind(a, kind::constant)) {
      return -1;
    }

    if (expr_is_kind(b, kind::constant)) {
      return +1;
    }

    if (expr_is_kind(a, kind::power) && expr_is_kind(b, kind::power)) {
      return expr_cmp(a[0], b[0], ctx);
    }
    if (expr_is_kind(a, kind::symbol | kind::add) && expr_is_kind(b, kind::power)) {

      int order = expr_cmp(a, b[0], kind::mul);

      if (order == 0) {
        return expr_kind(a) - expr_kind(b);
      }

      return order;
    }

    if (expr_is_kind(b, kind::symbol | kind::add) && expr_is_kind(a, kind::power)) {
      int order = expr_cmp(a[0], b, kind::mul);

      if (order == 0) {
        return expr_kind(a) - expr_kind(b);
      }

      return order;
    }

    // if (expr_is_kind(a, kind::funcall) && expr_is_kind(b, kind::funcall)) {
    //   return strcmp(expr_funname(a), expr_funname(b));
    // }

    // if (expr_is_kind(a, kind::power) && expr_is_kind(b, kind::funcall)) {
    //   return expr_cmp(a[0], b, kind::mul);
    // }

    // if (expr_is_kind(b, kind::power) && expr_is_kind(a, kind::funcall)) {
    //   return expr_cmp(b, a[0], kind::mul);
    // }

    if (expr_is_kind(a, kind::mul) &&
        expr_is_kind(b, kind::power | kind::symbol /*| kind::funcall*/)) {
      long k = expr_size(a);

      for (long i = k - 1; i >= 0; i--) {
        int order = expr_cmp(a[i], b, kind::mul);

        if (order == 0) {
          return 0;
        }

        if (order < 0)
          break;
      }

      return +k - 1;
    }

    if (expr_is_kind(b, kind::mul) &&
        expr_is_kind(a, kind::power | kind::symbol /* | kind::funcall */)) {
      long k = expr_size(b);

      for (long i = k - 1; i >= 0; i--) {
        int order = expr_cmp(b[i], a, kind::mul);

        if (order == 0) {
          return 0;
        }

        if (order < 0)
          break;
      }

      return -k + 1;
    }
  }

  if (ctx & kind::add) {

    if (expr_is_kind(a, kind::constant) && expr_is_kind(b, kind::constant)) {
      return expr_cmp_consts(b, a);
    }

    if (expr_is_kind(a, kind::constant)) {
      return +1;
    }

    if (expr_is_kind(b, kind::constant)) {
      return -1;
    }

    if (expr_is_kind(a, kind::power) && expr_is_kind(b, kind::power)) {
      int i = expr_cmp(a[1], b[1], ctx);

      if (i != 0) {
        return i;
      }

      return expr_cmp(a[0], b[0], ctx);
    }

    // if (expr_is_kind(a, kind::funcall) && expr_is_kind(b, kind::funcall)) {
    //   return strcmp(expr_funname(a), expr_funname(b));
    // }

    if (expr_is_kind(a, kind::add) && expr_is_kind(b, kind::symbol)) {
      return +1;
    }

    if (expr_is_kind(a, kind::symbol) && expr_is_kind(b, kind::add)) {
      return -1;
    }

    if (expr_is_kind(a, kind::mul) && expr_is_kind(b, kind::symbol)) {
      long k = expr_size(a);

      if (k > 2)
        return -k;

      int order = expr_cmp(a[expr_size(a) - 1], b, ctx);

      if (order == 0) {
        return -1;
      }

      return k > 2 ? -1 : order;
    }

    if (expr_is_kind(b, kind::mul) && expr_is_kind(a, kind::symbol)) {
      long k = expr_size(b);

      if (k > 2) {
        return +k;
      }

      int order = expr_cmp(a, b[expr_size(b) - 1], ctx);

      if (order == 0) {
        return +1;
      }

      return k > 2 ? +1 : order;
    }

    if (expr_is_kind(a, kind::mul) && expr_is_kind(b, kind::power)) {
      long k = expr_size(a);

      if (k > 2)
        return -k;

      int order = expr_cmp(a[0], b, ctx);

      if (order == 0)
        return -1;

      return k > 2 ? -1 : order;
    }

    if (expr_is_kind(b, kind::mul) && expr_is_kind(a, kind::power)) {
      long k = expr_size(b);

      if (k > 2)
        return +k;

      int order = expr_cmp(a, b[0], ctx);

      if (order == 0)
        return +1;

      return k > 2 ? +1 : order;
    }
  }

  // if (expr_is_kind(a, kind::funcall) && expr_is_kind(b, kind::funcall)) {
  //   return strcmp(expr_funname(a), expr_funname(b));
  // }

  if (expr_is_kind(a, kind::constant) && expr_is_kind(b, kind::constant)) {
    return expr_cmp_consts(a, b);
  }

  if (expr_is_kind(a, kind::symbol) && expr_is_kind(b, kind::symbol)) {
    return expr_cmp_idents(a, b);
  }

  if (expr_is_kind(a, kind::add) && expr_is_kind(b, kind::add)) {
    return expr_op_cmp(a, b, ctx);
  }

  if (expr_is_kind(a, kind::mul) && expr_is_kind(b, kind::mul)) {
    return expr_op_cmp(a, b, ctx);
  }

  if (expr_is_kind(a, kind::power) && expr_is_kind(b, kind::power)) {
    return expr_cmp(a[0], b[0], ctx) ||
           expr_cmp(a[1], b[1], ctx);
  }

  if (expr_is_kind(a, kind::div) && expr_is_kind(b, kind::div)) {
    return expr_cmp(a[0], b[0], ctx) ||
           expr_cmp(a[1], b[1], ctx);
  }

  return (ctx & kind::add) ? expr_kind(a) - expr_kind(b)
                                : expr_kind(b) - expr_kind(a);
}

long int expr_sort_split(expr& a, long l, long r) {
  long int i = l - 1;

  expr p = a[r];

  for (long int j = l; j < r; j++) {
    if (expr_cmp(a[j], p, expr_kind(a)) < 0) {
      i = i + 1;

      // swap i and j
      expr t = a[i];

      a[i] = a[j];
      a[j] = t;
    }
  }

  expr t = a[i + 1];

  a[i + 1] = a[r];
  a[r] = t;

  return i + 1;
}

void expr_sort_childs(expr& a, long int l, long int r) {
  if (l < r) {
    long int m = expr_sort_split(a, l, r);

    expr_sort_childs(a, l, m - 1);
    expr_sort_childs(a, m + 1, r);
  }

  return;
}

void expr_sort(expr& a) {
  if (expr_is_kind(a, kind::terminal)) {
    return;
  }

  size_t i = expr_is_kind(a, kind::sub) ? 1 : 0;

  for (; i < expr_size(a); i++) {
    expr_sort(a[i]);
  }

  if (expr_is_kind(a, kind::ordered)) {
    return;
  }

  expr_sort_childs(a, 0, expr_size(a) - 1);
}

inline expr& expr_set_to_undefined(expr& a) {
  expr_set_kind(a, kind::undefined);

	a.ops = std::vector<expr>();

	return a;
}

inline expr& expr_set_to_fail(expr& a) {
  expr_set_kind(a, kind::fail);

	a.ops = std::vector<expr>();

	return a;
}

inline expr& expr_set_to_int(expr& a, Int v) {
  expr_set_kind(a, kind::integer);

	a.ops = std::vector<expr>();

  a.val = new Int(v);

	return a;
}

inline expr& expr_set_op_to_int(expr& a, size_t i, Int v) {
	expr_set_to_int(a[i], v);
	return a;
}

inline expr& expr_set_to_fra(expr& a, Int u, Int v) {
  expr_set_kind(a, kind::fraction);

	a.ops = std::vector<expr>({expr_integer(u), expr_integer(v)});
	return a;
}

inline expr expr_set_op_to_fra(expr& a, size_t i, Int u, Int v) {
	expr_set_to_fra(a[i], u, v);
	return a;
}

inline expr& expr_set_to_sym(expr& a, const char *s) {
  expr_set_kind(a, kind::symbol);

	a.ops = std::vector<expr>();

	a.id = strdup(s);

	return a;
}

inline expr& expr_set_op_to_sym(expr& a, size_t i, const char *s) {
  expr_set_to_sym(a[i], s);
	return a;
}

// a = a + b
inline expr& expr_set_inplace_add_consts(expr& a, expr& b) {
  assert(expr_is_kind(a, kind::constant));
  assert(expr_is_kind(b, kind::constant));

  if (expr_is_kind(a, kind::integer) && expr_is_kind(b, kind::integer)) {
    Int x = expr_value(a);
    Int y = expr_value(b);

    return expr_set_to_int(a, x + y);
  }

  if (expr_is_kind(a, kind::fraction) && expr_is_kind(b, kind::fraction)) {
    Int x = expr_value(a[0]);
    Int y = expr_value(a[1]);
    Int w = expr_value(b[0]);
    Int z = expr_value(b[1]);

    Int e = x * z + w * y;
    Int k = y * z;

    if (e % k == 0) {
      a = expr_set_to_int(a, e / k);
    } else {

      Int g = abs(gcd(e, k));

      a = expr_set_to_fra(a, e / g, k / g);
    }

    return a;
  }

  if (expr_is_kind(a, kind::fraction) && expr_is_kind(b, kind::integer)) {
    Int x = expr_value(a[0]);
    Int y = expr_value(a[1]);
    Int w = expr_value(b);

    Int e = x + w * y;
    Int k = y;

    if (e % k == 0) {
      a = expr_set_to_int(a, e / k);
    } else {

      Int g = abs(gcd(e, k));

      a = expr_set_to_fra(a, e / g, k / g);
    }

    return a;
  }

  if (expr_is_kind(b, kind::fraction) && expr_is_kind(a, kind::integer)) {
    Int x = expr_value(b[0]);
    Int y = expr_value(b[1]);
    Int w = expr_value(a);

    Int e = x + w * y;
    Int k = y;

    if (e % k == 0) {
      a = expr_set_to_int(a, e / k);
    } else {

      Int g = abs(gcd(e, k));

      a = expr_set_to_fra(a, e / g, k / g);
    }

    return a;
  }

  return a;
}

inline expr& expr_set_inplace_add_consts(expr& a, Int b) {
  assert(expr_is_kind(a, kind::constant));

  if (expr_is_kind(a, kind::integer)) {
    Int x = expr_value(a);

    return expr_set_to_int(a, x + b);
  }

  if (expr_is_kind(a, kind::fraction)) {
    Int x = expr_value(a[0]);
    Int y = expr_value(a[1]);

    Int e = x + b * y;

    if (e % y == 0) {
      a = expr_set_to_int(a, e / y);
    } else {

      Int g = abs(gcd(e, y));

      a = expr_set_to_fra(a, e / g, y / g);
    }

    return a;
  }

  return a;
}

inline expr& expr_set_op_inplace_add_consts(expr& a, size_t i, expr b) {
  expr_set_inplace_add_consts(a[i], b);
  return a;
}

inline expr& expr_set_op_inplace_add_consts(expr &a, size_t i, Int b) {
  expr_set_inplace_add_consts(a[i], b);
  return a;
}

inline expr& expr_set_inplace_mul_consts(expr& a, expr& b) {
  assert(expr_is_kind(a, kind::constant));
  assert(expr_is_kind(b, kind::constant));

  if (expr_is_kind(a, kind::integer) && expr_is_kind(b, kind::integer)) {
    Int x = expr_value(a);
    Int y = expr_value(b);

    return expr_set_to_int(a, x * y);
  }

  if (expr_is_kind(a, kind::fraction) && expr_is_kind(b, kind::fraction)) {
    Int x = expr_value(a[0]);
    Int y = expr_value(a[1]);
    Int w = expr_value(b[0]);
    Int z = expr_value(b[1]);

    Int e = x * w;
    Int k = y * z;

    if (e % k == 0) {
      a = expr_set_to_int(a, e / k);
    } else {

      Int g = abs(gcd(e, k));

      a = expr_set_to_fra(a, e / g, k / g);
    }

    return a;
  }

  if (expr_is_kind(a, kind::fraction) && expr_is_kind(b, kind::integer)) {
    Int x = expr_value(a[0]);
    Int y = expr_value(a[1]);
    Int w = expr_value(b);

    Int e = x * w;
    Int k = y;

    if (e % k == 0) {
      a = expr_set_to_int(a, e / k);
    } else {

      Int g = abs(gcd(e, k));

      a = expr_set_to_fra(a, e / g, k / g);
    }

    return a;
  }

  if (expr_is_kind(b, kind::fraction) && expr_is_kind(a, kind::integer)) {
    Int x = expr_value(b[0]);
    Int y = expr_value(b[1]);
    Int w = expr_value(a);

    Int e = x * w;
    Int k = y;

    if (e % k == 0) {
      a = expr_set_to_int(a, e / k);
    } else {

      Int g = abs(gcd(e, k));

      a = expr_set_to_fra(a, e / g, k / g);
    }

    return a;
  }

  return a;
}

inline expr& expr_set_inplace_mul_consts(expr& a, Int b) {
  assert(expr_is_kind(a, kind::constant));

  if (expr_is_kind(a, kind::integer)) {
    Int x = expr_value(a);

    return expr_set_to_int(a, x * b);
  }

  if (expr_is_kind(a, kind::fraction)) {
    Int x = expr_value(a[0]);
    Int y = expr_value(a[1]);

    Int e = x * b;

    if (e % y == 0) {
      a = expr_set_to_int(a, e / y);
    } else {

      Int g = abs(gcd(e, y));

      a = expr_set_to_fra(a, e / g, y / g);
    }

    return a;
  }
  return a;
}

inline expr& expr_set_op_inplace_mul_consts(expr& a, size_t i, expr b) {
  expr_set_inplace_mul_consts(a[i], b);
  return a;
}

inline expr& expr_set_op_inplace_mul_consts(expr& a, size_t i, Int b) {
  expr_set_inplace_mul_consts(a[i], b);
	return a;
}

inline expr& expr_set_to_mul(Int v, expr& a) {
  if (expr_is_kind(a, kind::mul)) {
    expr_insert(a, expr_integer(v), 0);
  } else {
		a = expr_create(kind::mul, {expr_integer(v), a});
  }

  return a;
}

inline expr& expr_set_to_mul(expr& a, expr& b) {
  if (expr_is_kind(a, kind::mul)) {
    expr_insert(a, b, 0);
  } else {
    a = expr_create(kind::mul, {a, b});
  }

  return a;
}

inline expr& expr_set_op_to_mul(expr &a, size_t i, Int v) {
  expr_set_to_mul(v, a[i]);
  return a;
}

inline expr& expr_set_op_to_mul(expr& a, size_t i, expr v) {
  expr_set_to_mul(v, a[i]);
  return a;
}

inline expr& expr_set_to_pow(expr& a, Int e) {
	a = expr_create(kind::power, {a, expr_integer(e)});
  return a;
}

inline expr& expr_set_to_add(expr& a, expr& e) {
  if (expr_is_kind(a, kind::add)) {
    expr_insert(a, e);

    return a;
  }

	a = expr_create(kind::add, {a, e});

	return a;
}

inline expr& expr_set_op_to_add(expr& a, size_t i, expr v) {
  expr_set_to_add(a[i], v);
  return a;
}

inline expr& expr_set_op_pow_add_to_deg(expr& a, size_t i, expr& e) {
  assert(expr_is_kind(a[i], kind::power));

  a[i] =
      expr_create(kind::power, {a[i][0], expr_create(kind::add, {a[i][1], e})});

  return a;
}

inline expr& expr_set_op_to_pow(expr& a, size_t i, Int v) {
  a[i] = expr_create(kind::power, {a[i], expr_integer(v)});
  return a;
}

inline expr& expr_set_op_to_pow(expr& a, size_t i, expr& v) {
	a[i] = expr_create(kind::power, {a[i], v});
  return a;
}

inline expr expr_detatch_operand(expr& a, size_t i) {
  expr b = a[i];

	expr_remove(a, i);

  return b;
}

inline expr& eval_add_consts(expr& u, size_t i, expr& v, size_t j) {
  return expr_set_op_inplace_add_consts(u, i, v[j]);
}

inline expr& eval_mul_consts(expr& u, size_t i, expr& v, size_t j) {
  return expr_set_op_inplace_mul_consts(u, i, v[j]);
}

inline expr& eval_add_int(expr& a, Int b) {
  if (expr_is_kind(a, kind::integer)) {
    return expr_set_to_int(a, expr_value(a) + b);
  }

  assert(expr_is_kind(a, kind::fraction));

  Int num = expr_value(a[0]);
  Int den = expr_value(a[1]);

  num = b * den + num;

  Int cff = abs(gcd(num, den));

  if (den / cff == 1) {
    return expr_set_to_int(a, num / cff);
  }

  return expr_set_to_fra(a, num / cff, den / cff);
}

expr Fail = expr_create(kind::fail);

inline expr& eval_add_nconst(expr& u, size_t i, expr& v, size_t j) {
  assert(expr_is_kind(u, kind::add) && expr_is_kind(v, kind::add));

  assert(expr_is_kind(u[i], kind::summable));
  assert(expr_is_kind(v[j], kind::summable));

  // expr a = u[i];
  // expr b = v[j];

	if (expr_is_kind(u[i], kind::power) && expr_is_kind(v[j], kind::symbol)) {
    return Fail;
  }

  if (expr_is_kind(u[i], kind::symbol) && expr_is_kind(v[j], kind::power)) {
    return Fail;
  }

  int kind = expr_kind(u[i]) & expr_kind(v[j]);

  if (kind & (kind::symbol | kind::power)) {
    if (expr_cmp(u[i], v[j], kind::add) == 0) {
      return expr_set_op_to_mul(u, i, 2);
    }

    return Fail;
  }

  long size_a = expr_size(u[i]);

  if (kind & kind::mul) {
    long size_b = expr_size(v[j]);

    long size_c = -1;

    if (expr_is_kind(u[i][0], kind::integer) &&
				expr_is_kind(v[j][0], kind::integer) &&
        std::abs(size_a - size_b) != 0) {
      return Fail;
    }

    if (size_b > size_a) {
      expr c = u[i];

			u[i] = v[j];
      v[j] = c;

      size_c = size_b;
      size_b = size_a;
      size_a = size_c;
    }

    assert(size_a == size_b || size_a == size_b + 1);

    long size =
        size_b - (expr_is_kind(v[j][0], kind::constant) ? 1 : 0);

    for (long x = 0; x < size; x++) {
      if (expr_cmp(u[i][size_a - x - 1], v[j][size_b - x - 1], kind::add) != 0) {
        return Fail;
      }
    }

    int ka = expr_kind(u[i][0]);
    int kb = expr_kind(v[j][0]);

    if ((ka & kind::constant) && (kb & kind::constant)) {
      u[i] = expr_set_op_inplace_add_consts(u[i], 0, v[j][0]);
    } else if (ka & kind::constant) {
      u[i] = expr_set_op_inplace_add_consts(u[i], 0, 1);
    } else {
      u[i] = expr_set_to_mul(2, u[i]);
    }

    if (size_c != -1) {
			expr c = u[i];

			u[i] = v[j];
			v[j] = c;
    }

    return u;
  }

  if (expr_is_kind(v[j], kind::mul)) {
    return Fail;
  }

  assert(expr_is_kind(u[i], kind::mul));
  assert(expr_is_kind(v[j], kind::symbol | kind::power));

  if (size_a > 2) {
    return Fail;
  }

  long ki = expr_kind(u[i][1]) & expr_kind(v[j]);

  if (!expr_is_kind(u[i][0], kind::constant) || !ki) {
    return Fail;
  }

  if (expr_cmp(v[j], u[i][1], kind::mul) == 0) {
    u[i] = expr_set_op_inplace_add_consts(u[i], 0, 1);
    return u;
  }

  return Fail;
}

inline expr& eval_mul_nconst(expr& u, size_t i, expr& v, size_t j) {
  assert(expr_is_kind(u, kind::mul) && expr_is_kind(v, kind::mul));

  assert(expr_is_kind(u[i], kind::multiplicable));
  assert(expr_is_kind(v[j], kind::multiplicable));

  if (expr_is_kind(u[i], kind::add) && expr_is_kind(v[j], kind::add)) {

    if (expr_cmp(u[i], v[j], kind::mul) == 0) {
      return expr_set_op_to_pow(u, i, 2);
    }

    return Fail;
  }

  if (expr_is_kind(u[i], kind::add) && expr_is_kind(v[j], kind::power)) {
    if (expr_cmp(u[i], v[j][0], kind::mul) == 0) {
      expr e = expr_create(kind::add, { expr_integer(1), v[j][1] });

      u = expr_set_op_to_pow(u, i, e);

			expr_eval(u[i]);

			return u;
    }

    return Fail;
  }

  if (expr_is_kind(u[i], kind::power) && expr_is_kind(v[j], kind::add)) {
    if (expr_cmp(v[j], u[i][0], kind::mul) == 0) {
			expr one = expr_integer(1);

			u = expr_set_op_pow_add_to_deg(u, i, one);

			expr_eval(u[i]);

			return u;
    }

    return Fail;
  }

  if (expr_is_kind(u[i], kind::power) && expr_is_kind(v[j], kind::power)) {
    if (expr_cmp(u[i][0], v[j][0], kind::mul) == 0) {
			u = expr_set_op_pow_add_to_deg(u, i, v[j][1]);

			expr_eval(u[i]);

      return u;
    }

    return Fail;
  }

  if (expr_is_kind(u[i], kind::symbol /* | kind::funcall */) &&
      expr_is_kind(v[j], kind::symbol /* | kind::funcall */)) {
    if (expr_cmp(u[i], v[j], kind::mul) == 0) {
      return expr_set_op_to_pow(u, i, 2);
    }

    return Fail;
  }

  if (expr_is_kind(u[i], kind::power) && expr_is_kind(v[j], kind::symbol /*| kind::funcall */)) {
    if (expr_cmp(u[i][0], v[j], kind::mul) == 0) {
      expr one = expr_integer(1);

      u = expr_set_op_pow_add_to_deg(u, i, one);

      expr_eval(u[i]);

			return u;
    }

    return Fail;
  }

  return Fail;
}

expr& eval_add_add(expr& a, expr& b) {
  assert(expr_is_kind(a, kind::add));
  assert(expr_is_kind(b, kind::add));

  size_t j = 0;
  size_t i = 0;

  while (i < expr_size(a) && j < expr_size(b)) {
    assert(!expr_is_kind(b[j], kind::add));

    expr& tmp = Fail;

    if (expr_is_kind(a[i], kind::constant) && expr_is_kind(b[j], kind::constant)) {
      tmp = eval_add_consts(a, i, b, j);
    } else if (expr_is_kind(a[i], kind::summable) && expr_is_kind(b[j], kind::summable)) {
      tmp = eval_add_nconst(a, i, b, j);
    }

    if (!expr_is_kind(tmp, kind::fail)) {
      a = tmp;

      i = i + 1;
      j = j + 1;

		} else {
      int order = expr_cmp(a[i], b[j], kind::add);

      if (order < 0) {
        i = i + 1;
      } else {
        expr_insert(a, b[j++], i++);
      }
    }
  }

  while (j < expr_size(b)) {
    if (i >= expr_size(a)) {
      expr_insert(a, b[j++], expr_size(a));
    } else {
      if (expr_cmp(a[i], b[j], kind::add) < 0) {
        i++;
      } else {
        expr_insert(a, b[j++], i++);
      }
    }
  }

  return a;
}

inline expr& eval_mul_mul(expr& a, expr& b) {
  assert(expr_is_kind(a, kind::mul));
  assert(expr_is_kind(b, kind::mul));

  size_t j = 0;
  size_t i = 0;

  while (i < expr_size(a) && j < expr_size(b)) {
    assert(!expr_is_kind(b[j], kind::mul));

    expr& tmp = Fail;

    if (expr_is_kind(a[i], kind::constant) && expr_is_kind(b[j], kind::constant)) {
      tmp = eval_mul_consts(a, i, b, j);
    } else if (expr_is_kind(a[i], kind::multiplicable) &&
               expr_is_kind(b[j], kind::multiplicable)) {
      tmp = eval_mul_nconst(a, i, b, j);
    }

    if (!expr_is_kind(tmp, kind::fail)) {
			a = tmp;
			i = i + 1;
      j = j + 1;
    } else {
      int order = expr_cmp(a[i], b[j], kind::mul);

      if (order < 0) {
        i = i + 1;
      } else {
        expr_insert(a, b[j++], i++);
      }
    }
  }

  while (j < expr_size(b)) {
    if (i >= expr_size(a)) {
      expr_insert(a, b[j++], expr_size(a));
    } else {
      if (expr_cmp(a[i], b[j], kind::mul) < 0) {
        i++;
      } else {
        expr_insert(a, b[j++], i++);
      }
    }
  }

  return a;
}

expr& eval_mul_int(expr& u, size_t i, Int v) {
  if (expr_is_kind(u[i], kind::integer)) {
    return expr_set_op_inplace_mul_consts(u, i, v);
  }

  if (expr_is_kind(u[i], kind::add | kind::sub)) {
    for (size_t j = 0; j < expr_size(u[i]); j++) {
      u[i] = expr_set_op_inplace_mul_consts(u[i], j, v);
    }

    return u;
  }

  if (expr_is_kind(u[i], kind::mul)) {
    return expr_set_op_to_mul(u, i, v);
  }

  if (expr_is_kind(u[i], kind::div)) {
    return eval_mul_int(u[i], 0, v);
  }

  if (expr_is_kind(u[i], kind::sqrt | kind::power | kind::factorial /* | kind::funcall */ |
                         kind::symbol)) {

    return expr_set_op_inplace_mul_consts(u, i, v);
  }

  assert(expr_is_kind(u[i], kind::undefined | kind::fail));

  return u;
}

expr& expr_replace_with(expr& a, expr& t) {
  if (expr_is_kind(t, kind::integer)) {
    return expr_set_to_int(a, expr_value(t));
  }

  if (expr_is_kind(t, kind::symbol)) {
    return expr_set_to_sym(a, expr_id(t));
  }

  // if (expr_is_kind(t, kind::funcall)) {
  //   a.id = strdup(expr_id(t));

  //   expr_set_kind(a, kind::funcall);

	// 	expr_allocate_operands(a, expr_size(t));

	// 	for(size_t i = 0; i < expr_size(t); i++) {
	// 		a[i] = expr_inc_ref(t[i]);
	// 	}

  //   return a;
  // }

  expr_set_kind(a, expr_kind(t));

	a.ops = t.ops;

	return a;
}


expr& expr_raise_to_first_op(expr& a) {

  if (expr_is_kind(a[0], kind::integer)) {
    return expr_set_to_int(a, expr_value(a[0]));
  }

  if (expr_is_kind(a[0], kind::symbol)) {
    return expr_set_to_sym(a, expr_id(a[0]));
  }

  // if (expr_is_kind(a[0], kind::funcall)) {
  //   a.id = strdup(expr_id(a[0]));

  //   expr_set_kind(a, kind::funcall);

  //   expr t = a[0];

  //   expr *t_childs = a;
  //   size_t t_size = a->expr_size;
  //   size_t t_rsize = a->expr_reserved_size;

  //   a = t;
  //   a->expr_size = t->expr_size;
  //   a->expr_reserved_size = t->expr_reserved_size;

  //   t = t_childs;
  //   t->expr_size = t_size;
  //   t->expr_reserved_size = t_rsize;

  //   expr_delete(t);

  //   return a;
  // }

  expr_set_kind(a, expr_kind(a[0]));
	a.ops = a[0].ops;

  return a;
}

expr& expr_eval_add(expr& a) {
	for (size_t i = 0; i < expr_size(a); i++) {
    expr_eval(a[i]);
  }

  expr_sort_childs(a, 0, expr_size(a) - 1);

  if (expr_is_kind(a[0], kind::add)) {
    expr t = expr_detatch_operand(a, 0);

    eval_add_add(a, t);
  }

  size_t j = 0;

  for (long i = 1; i < (long)expr_size(a); i++) {
    expr& t = Fail;

    if (expr_is_kind(a[i], kind::fail) || expr_is_kind(a[j], kind::fail)) {
      return expr_set_to_fail(a);
    } else if (expr_is_kind(a[i], kind::undefined) ||
               expr_is_kind(a[j], kind::undefined)) {
      return expr_set_to_undefined(a);
    } else if (expr_is_kind(a[j], kind::integer) && expr_value(a[j]) == 0) {
			expr_remove(a, j);
    } else if (expr_is_kind(a[i], kind::add)) {
      expr g = expr_detatch_operand(a, i--);
      a = eval_add_add(a, g);
    } else if (expr_is_kind(a[j], kind::constant) &&
               expr_is_kind(a[i], kind::constant)) {

      t = eval_add_consts(a, j, a, i);

      if (!expr_is_kind(t, kind::fail)) {
        a = t;
        expr_remove(a, i--);
      }

    } else if (expr_is_kind(a[j], kind::summable) &&
               expr_is_kind(a[i], kind::summable)) {
      t = eval_add_nconst(a, j, a, i);

      if (!expr_is_kind(t, kind::fail)) {
        a = t;
        expr_remove(a, i--);
      }
    }

    if (expr_is_kind(t, kind::fail))
      j = i;
  }

  if (expr_size(a) == 0) {
		return expr_set_to_int(a, 0);
  }

	if (expr_size(a) == 1) {
    return expr_raise_to_first_op(a);
  }

  return a;
}

expr& expr_eval_mul(expr& a) {
  for (size_t i = 0; i < expr_size(a); i++) {
		expr_eval(a[i]);
  }


	expr_sort_childs(a, 0, expr_size(a) - 1);

  size_t j = 0;

  if (expr_is_kind(a[0], kind::mul)) {
    expr t = expr_detatch_operand(a, 0);
    eval_mul_mul(a, t);
  }

  for (long i = 1; i < (long)expr_size(a); i++) {
    expr &t = Fail;

    if (expr_is_kind(a[i], kind::fail) || expr_is_kind(a[j], kind::fail)) {
      return expr_set_to_fail(a);
    } else if (expr_is_kind(a[i], kind::undefined) ||
               expr_is_kind(a[j], kind::undefined)) {
      return expr_set_to_undefined(a);
    } else if ((expr_is_kind(a[i], kind::integer) && expr_value(a[i]) == 0) ||
               (expr_is_kind(a[j], kind::integer) && expr_value(a[j]) == 0)) {
      return expr_set_to_int(a, 0);
    } else if (expr_is_kind(a[i], kind::mul)) {
      expr g = expr_detatch_operand(a, i--);
      a = eval_mul_mul(a, g);
    } else if (expr_is_kind(a[j], kind::constant) &&
               expr_is_kind(a[i], kind::constant)) {
      t = eval_mul_consts(a, j, a, i);

      if (!expr_is_kind(t, kind::fail)) {
        a = t;
        expr_remove(a, i--);
      }
    } else if (expr_is_kind(a[j], kind::multiplicable) &&
               expr_is_kind(a[i], kind::multiplicable)) {
      t = eval_mul_nconst(a, j, a, i);

      if (!expr_is_kind(t, kind::fail)) {
        a = t;
        expr_remove(a, i--);
      }
    }

    if (expr_is_kind(t, kind::fail)) {
      j = i;
    }
  }

        if(expr_size(a) == 0) {
		return expr_set_to_int(a, 1);
	}

	if (expr_is_kind(a, kind::mul) && expr_size(a) == 1) {
		return expr_raise_to_first_op(a);
  }

  return a;
}

expr& expr_eval_sub(expr& a) {
  for (size_t i = 1; i < expr_size(a); i++) {
    eval_mul_int(a, i, -1);
  }

  expr_set_kind(a, kind::add);

  return expr_eval(a);
}

expr& expr_eval_pow(expr& a) {
	a[1] = expr_eval(a[1]);

  // TODO: if expoent is zero return 1, if expoent is 1 return base

	a[0] = expr_eval(a[0]);

	if (!expr_is_kind(a[1], kind::integer)) {
    return a;
  }

  if (expr_value(a[1]) == 1) {
		return expr_raise_to_first_op(a);
  }

  if (expr_value(a[1]) == 0) {
    return expr_set_to_int(a, 1);
  }

  if (expr_is_kind(a[0], kind::integer)) {
    Int b = expr_value(a[0]);
    Int c = expr_value(a[1]);

    bool n = c < 0;

    c = abs(c);

    Int d = pow(b, c);

    if (!n || d == 1) {
      a = expr_set_to_int(a, d);
    } else {
      a = expr_set_to_fra(a, 1, d);
    }

    return a;
  }
  if (expr_is_kind(a[0], kind::fraction)) {
    Int b = expr_value(a[0][0]);
    Int c = expr_value(a[0][1]);

    Int d = expr_value(a[1]);

    bool n = d < 0;

    d = abs(d);

    b = pow(b, d);
    c = pow(c, d);

    Int g = gcd(b, c);

    b = b / g;
    c = c / g;

    a = !n ? expr_set_to_fra(a, b, c) : expr_set_to_fra(a, c, b);

    return a;
  }

  if (expr_is_kind(a[0], kind::terminal)) {
    return a;
  }

  if (expr_is_kind(a[0], kind::mul)) {

    long long y = expr_value(a[1]).longValue();

    expr b = expr_detatch_operand(a, 0);

    a = expr_set_to_int(a, 1);

    while (y) {
      if (y % 2 == 1) {
        a = expr_create(kind::mul, { a, b });
        a = expr_eval(a);
      }

      y = y >> 1;

      expr t = expr_create(kind::mul);

			// expr_set_size(t, expr_size(b));
      for (size_t i = 0; i < expr_size(b); i++) {
        expr_insert(t, expr_create(kind::mul, {
							b[i],
							b[i],
						}));
      }

			expr_replace_with(b, t);
    }

    return a;
  }

  return a;
}

expr& expr_eval_div(expr& a) {
  expr_set_kind(a, kind::mul);

  a = expr_set_op_to_pow(a, 1, -1);

  return expr_eval(a);
}

expr& expr_eval_sqr(expr& a) {
  expr_set_kind(a, kind::power);

  expr_insert(a, expr_fraction(1, 2), 1);

  return expr_eval(a);
}

expr& expr_eval_fac(expr& a) {
  if (expr_is_kind(a[0], kind::integer)) {
    return expr_set_to_int(a, fact(expr_value(a[0])));
  }

  if (expr_is_kind(a[0], kind::fraction)) {
    Int c = fact(expr_value(a[0][0]));
    Int d = fact(expr_value(a[0][1]));

    Int g = abs(gcd(c, d));

    return expr_set_to_fra(a, c / g, d / g);
  }

  return a;
}

expr& expr_eval_fra(expr& a) {
  Int b = expr_value(a[0]);
  Int c = expr_value(a[1]);

  Int d = abs(gcd(b, c));

  return expr_set_to_fra(a, b / d, c / d);
}

expr& expr_eval(expr& a) {
  if (expr_is_kind(a, kind::fraction)) {
    return expr_eval_fra(a);
  }

	if (expr_is_kind(a, kind::add)) {
    return expr_eval_add(a);
  }

	if (expr_is_kind(a, kind::mul)) {
    return expr_eval_mul(a);
  }

	if (expr_is_kind(a, kind::sub)) {
    return expr_eval_sub(a);
  }

	if (expr_is_kind(a, kind::div)) {
    return expr_eval_div(a);
  }

	if (expr_is_kind(a, kind::power)) {
    return expr_eval_pow(a);
  }

	if (expr_is_kind(a, kind::sqrt)) {
    return expr_eval_sqr(a);
  }

	if (expr_is_kind(a, kind::factorial)) {
    return expr_eval_fac(a);
  }

  return a;
}

// expr expr_expand_mul(expr a, size_t i, expr b, size_t j) {
//   expr r = a[i];
//   expr s = b[j];

//   if (expr_is_kind(r, kind::add) && expr_is_kind(s, kind::add)) {
//     expr u = expr_create(kind::add, expr_size(r) * expr_size(s));

//     expr_set_size(u, expr_size(r) * expr_size(s));

//     for (size_t k = 0; k < expr_size(r); k++) {
//       for (size_t t = 0; t < expr_size(s); t++) {
//         expr v = expr_create(kind::mul, {expr_inc_ref(r[k]),
//                                        expr_inc_ref(s[t])});

//         expr_set_operand(u, v, k * expr_size(s) + t);
//       }
//     }

//     return u;
//   }

//   if (expr_is_kind(r, kind::add)) {
//     expr u = expr_create(kind::add, expr_size(r));

//     expr_set_size(u, expr_size(r));

//     for (size_t k = 0; k < expr_size(r); k++) {
//       expr v = expr_create(kind::mul,
//                           {expr_inc_ref(r[k]), expr_inc_ref(s)});

//       expr_set_operand(u, v, k);
//     }

//     return u;
//   }

//   if (expr_is_kind(s, kind::add)) {
//     expr u = expr_create(kind::add, expr_size(s));

//     expr_set_size(u, expr_size(s));

//     for (size_t k = 0; k < expr_size(s); k++) {
//       expr v = expr_create(kind::mul,
//                           {expr_inc_ref(r), expr_inc_ref(s[k])});

//       expr_set_operand(u, v, k);
//     }

//     return u;
//   }

//   return expr_create(kind::mul, {expr_inc_ref(r), expr_inc_ref(s)});
// }

// expr expr_expand_pow(expr u, Int n) {
//   if (n == 1)
//     return expr_inc_ref(u);
//   if (n == 0)
//     return expr_integer(1);

//   if (expr_is_kind(u, kind::add)) {
//     Int c = fact(n);

//     expr o = expr_copy(u);

//     expr f = expr_detatch_operand(o, 0);

//     if (expr_size(o) == 0)
//       o = expr_set_to_int(o, 0);
//     if (expr_size(o) == 1)
//       o = expr_raise_to_first_op(o);

//     expr s = expr_create(kind::add, n.longValue() + 1);

//     expr_set_size(s, n.longValue() + 1);

//     for (Int k = 0; k <= n; k++) {
//       expr z = expr_create(
//           kind::mul,
//           {expr_integer(c / (fact(k) * fact(n - k))),
//            expr_create(kind::power, {expr_inc_ref(f), expr_integer(n - k)})});

//       expr t = expr_expand_pow(o, k);

//       expr q = t ? expr_create(kind::mul, {z, t}) : z;

//       expr_set_operand(s, q, k.longValue());
//     }

//     expr_delete(f);
//     expr_delete(o);

//     return s;
//   }

//   if (expr_is_kind(u, kind::terminal)) {
//     return 0;
//   }

//   return expr_create(kind::power, {
//                                   expr_inc_ref(u),
//                                   expr_integer(n),
//                               });

//   // return expr_eval(t);
// }

// expr expr_expand(expr a) {
//   expr_assign_modifiable_ref(a, &a);

//   if (expr_is_kind(a, kind::terminal)) {
//     return a;
//   }

//   if (expr_is_kind(a, kind::sub | kind::div | kind::factorial)) {
//     a = expr_eval(a);
//   }

//   if (expr_is_kind(a, kind::power)) {
// 		expr_set_operand(a, expr_expand(a[0]), 0);
//     expr_set_operand(a, expr_expand(a[1]), 1);

//     if (expr_is_kind(a[1], kind::integer)) {
//       expr t = expr_expand_pow(a[0], expr_value(a[1]));

//       if (t) {
//         expr_delete(a);
//         a = t;
//       }
//     }
// 	}

//   if (expr_is_kind(a, kind::mul)) {
//     while (expr_size(a) > 1) {
//       expr_set_operand(a, expr_expand(a[0]), 0);
//       expr_set_operand(a, expr_expand(a[1]), 1);

//       expr_insert(a, expr_expand_mul(a, 0, a, 1), 0);

//       expr_remove(a, 1);
//       expr_remove(a, 1);
//     }

//     a = expr_raise_to_first_op(a);
//   }

//   if (expr_is_kind(a, kind::add)) {
//     for (size_t i = 0; i < expr_size(a); i++) {
//       expr_set_operand(a, expr_expand(a[i]), i);
//     }
//   }

//   return expr_eval(a);
// }
} // namespace expr_teste
