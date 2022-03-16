#include "Reduction.hpp"

#include "Utils.hpp"
#include "Sorting.hpp"
#include "Expression.hpp"


namespace alg {

using namespace utils;

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

    if (e == 0) {
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

  // TODO: this may be unnecessary, remove it if thats the case
  set_to_unreduced(a);
}

inline void expr_set_op_inplace_add_consts(expr *a, size_t i, Int b) {
  expr_set_inplace_add_consts(operand(a, i), b);

  // TODO: this may be unnecessary, remove it if thats the case
  set_to_unreduced(a);
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

    if (e == 0) {
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

    if (e == 0) {
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

    if (e == 0) {
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

  // TODO: this may be unnecessary, remove it if thats the case
  set_to_unreduced(a);
}

inline void expr_set_op_inplace_mul_consts(expr *a, size_t i, Int b) {
  expr_set_inplace_mul_consts(operand(a, i), b);

  // TODO: this may be unnecessary, remove it if thats the case
  set_to_unreduced(a);
}

inline void expr_set_to_mul(Int v, expr *a) {
  set_to_unreduced(a);

  if (is(a, kind::MUL)) {
    return a->insert(integer(v), 0);
  }

  a->expr_childs = std::vector<expr>({integer(v), *a});

  expr_set_kind(a, kind::MUL);
}

inline void expr_set_to_mul(expr *a, expr *b) {
  set_to_unreduced(a);

  if (is(a, kind::MUL)) {
    return a->insert(*b);
  }

  a->expr_childs = std::vector<expr>({*b, *a});

  expr_set_kind(a, kind::MUL);
}

inline void expr_set_op_to_mul(expr *a, size_t i, Int v) {
  set_to_unreduced(a);

  return expr_set_to_mul(v, operand(a, i));
}

inline void expr_set_op_to_mul(expr *a, size_t i, expr *v) {
  set_to_unreduced(a);

  return expr_set_to_mul(operand(a, i), v);
}

inline void expr_set_to_pow(expr *a, Int e) {
  set_to_unreduced(a);

  a->expr_childs = std::vector<expr>({*a, integer(e)});

  expr_set_kind(a, kind::POW);
}

inline void expr_set_to_add(expr *a, expr *e) {
  set_to_unreduced(a);

  if (is(a, kind::ADD)) {
    return a->insert(*e);
  }

  a->expr_childs = std::vector<expr>({*a, *e});

  expr_set_kind(a, kind::ADD);
}

inline void expr_set_op_to_add(expr *a, size_t i, expr *v) {
  set_to_unreduced(a);

  expr_set_to_add(operand(a, i), v);
}

inline void expr_set_op_pow_add_to_deg(expr *a, size_t i, expr *e) {
  assert(is(operand(a, i), kind::POW));

  set_to_unreduced(a);

  a->expr_childs[i] =
      create(kind::POW, {*operand(operand(a, i), 0),
                         create(kind::ADD, {*operand(operand(a, i), 1), *e})});
}

inline void expr_set_op_pow_add_to_deg(expr *a, size_t i, Int e) {
  assert(is(operand(a, i), kind::POW));

  set_to_unreduced(a);

  a->expr_childs[i] = create(
      kind::POW, {*operand(operand(a, i), 0),
                  create(kind::ADD, {*operand(operand(a, i), 1), integer(e)})});
}

inline void expr_set_op_to_pow(expr *a, size_t i, Int v) {
  set_to_unreduced(a);

  a->expr_childs[i] = create(kind::POW, {*operand(a, i), integer(v)});
}

inline void expr_set_op_to_pow(expr *a, size_t i, expr *v) {
  set_to_unreduced(a);

  a->expr_childs[i] = create(kind::POW, {*operand(a, i), *v});
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

      u->expr_childs[i] = pow(*a, degree(*b) + 1);

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
      // printf("reduced %s and  %s ?", to_string(u).c_str(),
      // to_string(v).c_str());
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
      //   printf("order %s and  %s ? = %i\n", to_string(u).c_str(),
      //   to_string(v).c_str(), order);
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

  if (is(a, kind::MAT)) {
    matrix *A = u->expr_mat;

    u->expr_mat = matrix::mul_ptr(A, v.longValue());

    delete A;

    return true;
  }

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

void reduce_add(expr *a) {
  for (size_t i = 0; i < size_of(a); i++) {
    reduce(operand(a, i));
  }

  sort_childs(a, 0, size_of(a) - 1);

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

    else if (is(aj, kind::MAT) && is(ai, kind::MAT)) {
      if (columns(ai) == columns(aj) && rows(ai) == rows(aj)) {
        matrix *A = aj->expr_mat;
        matrix *B = ai->expr_mat;

        matrix *C = matrix::add_ptr(A, B);

        reduced = true;

        delete aj->expr_mat;

        aj->expr_mat = C;

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

    if (is(ai, kind::ADD) || is(aj, kind::ADD)) {
      set_to_unexpanded(a);
    }

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

    else if (is(aj, kind::MAT) && is(ai, kind::MAT)) {
      if (is(aj, kind::MAT) && is(ai, kind::MAT)) {
        if (columns(aj) == rows(ai)) {
          matrix *t = matrix::mul_ptr(aj->expr_mat, ai->expr_mat);

          delete aj->expr_mat;

          aj->expr_mat = t;

          reduced = true;

          a->remove(i--);
        }
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



  if (is(operand(a, 1), kind::INT)) {
    set_to_unexpanded(a);
  }

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

    set_to_unexpanded(a);

    return reduce(a);
  }

  if ((is_inf(operand(a, 0)) || is_neg_inf(operand(a, 0))) &&
      (is_inf(operand(a, 1)) || is_neg_inf(operand(a, 1)))) {
    return expr_set_to_undefined(a);
  }

	if(is(operand(a, 0), kind::SYM) && strcmp(operand(a, 0)->expr_sym, "i") == 0) {
		if(is(operand(a, 1), kind::INT)) {

			Int d = get_val(operand(a, 1));

			if(d == 0) {
				return expr_set_to_int(a, 1);
			}

			if(d == 1) {
				return expr_set_to_sym(a, "i");
			}

			if(d % 2 == 0) {
				return expr_set_to_int(a, pow(-1, d/2));
			}

			expr i = symbol("i");

			if(d % 2 == 1) {
				expr t = pow(-1, (d - 1)/2) * i;
				expr_replace_with(a, &t);
				return;
			}
		}
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

  if (is(operand(a, 0), kind::MAT)) {
    if (is(operand(a, 1), kind::CONST)) {
      if (is(operand(a, 1), kind::INT)) {
        expr_set_to_mat(
            a, matrix::div_ptr(operand(a, 0)->expr_mat,
                               operand(a, 1)->expr_int->doubleValue()));
      } else if (is(operand(a, 1), kind::FRAC)) {
        double n = operand(operand(a, 1), 0)->expr_int->doubleValue() /
                   operand(operand(a, 1), 1)->expr_int->doubleValue();
        matrix *t = matrix::div_ptr(operand(a, 0)->expr_mat, n);
        expr_set_to_mat(a, t);
      }

      return;
    }

    if (is(operand(a, 1), kind::MAT)) {
      double eps = std::numeric_limits<double>::epsilon();

      if (rows(operand(a, 1)) == columns(operand(a, 1)) &&
          matrix::det_ptr(operand(a, 1)->expr_mat) > eps) {
        // TODO: if condition above are not true a error should be
        //  emited, maybe? ...
        matrix *t = matrix::inv_ptr(operand(a, 1)->expr_mat);

        a->kind_of = kind::MUL;

        expr_set_to_mat(operand(a, 1), t);

        return reduce(a);
      }
    }

    return;
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

	a->expr_childs[1] = expr(kind::FRAC, { 1, a->expr_childs[1] });

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
  if (is_reduced(a)) {
    return;
  }

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

  set_to_reduced(a);
}


expr reduce(expr a) {
  reduce(&a);

  return a;
}

}
