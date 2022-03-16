#include "Utils.hpp"

namespace alg {

namespace utils {

Int rows(expr *a) {
  assert(is(a, kind::MAT));

  return a->expr_mat->lines();
}

Int columns(expr *a) {
  assert(is(a, kind::MAT));

  return a->expr_mat->columns();
}

bool is_reduced(expr *a) { return a->expr_info & info::REDUCED; }
bool is_sorted(expr *a, enum kind k) {
  return a->expr_info & info::SORTED && a->sort_kind == k;
}
bool expr_is_zero(expr *a) {
  return (a == nullptr) || (is(a, kind::INT) && get_val(a) == 0);
}

void set_to_reduced(expr *a) { a->expr_info |= info::REDUCED; }

void set_to_unsorted(expr *a) {
  a->expr_info &= ~(1UL << SORTED_BIT);
  a->sort_kind = kind::UNDEF;
}

void set_to_unreduced(expr *a) {
  a->expr_info &= ~(1UL << REDUCED_BIT);
  set_to_unsorted(a);
}

void set_to_unexpanded(expr *a) {
  a->expr_info &= ~(1UL << EXPANDED_BIT);
}


void set_to_sorted(expr *a, enum kind k) {
  a->expr_info |= info::SORTED;
  a->sort_kind = k & (kind::ADD | kind::MUL) ? k : kind::UNDEF;
}

bool is_negative(expr *a) {
  assert(is(a, kind::CONST));

  if (is(a, kind::INT) && get_val(a) < 0)
    return true;

  if (is(a, kind::FRAC) &&
      ((get_val(operand(a, 0)) < 0 && get_val(operand(a, 1)) > 0) ||
       (get_val(operand(a, 0)) > 0 && get_val(operand(a, 1)) < 0)))
    return true;

  return false;
}

void set_to_expanded(expr *a) { a->expr_info |= info::EXPANDED; }

void expr_set_kind(expr *a, kind kind) {
  if (is(a, kind::SYM | kind::FUNC | kind::ERROR)) {
    free(a->expr_sym);
    a->expr_sym = 0;
  }

  if (is(a, kind::INT)) {
    delete a->expr_int;
    a->expr_int = 0;
  }

  if (is(a, kind::LIST)) {
    delete a->expr_list;
    a->expr_list = 0;
  }

  if (is(a, kind::SET)) {
    delete a->expr_set;
    a->expr_set = 0;
  }

  if (is(a, kind::MAT)) {
    delete a->expr_mat;
    a->expr_set = 0;
  }

  a->kind_of = kind;
}

void expr_set_to_int(expr *a, Int v) {
  expr_set_kind(a, kind::INT);

  a->expr_childs.clear();

  a->expr_int = new Int(v);
}

void expr_set_to_mat(expr *a, matrix *v) {
  expr_set_kind(a, kind::MAT);

  a->expr_childs.clear();

  a->expr_mat = v;
}


void expr_set_op_to_int(expr *a, size_t i, Int v) {
  expr_set_to_int(operand(a, i), v);

  set_to_unreduced(a);
}

void expr_set_to_fra(expr *a, Int u, Int v) {
  expr_set_kind(a, kind::FRAC);

  a->expr_childs = std::vector<expr>();

  a->insert(integer(u));
  a->insert(integer(v));

  set_to_unreduced(a);
}

void expr_set_op_to_fra(expr *a, size_t i, Int u, Int v) {
  expr_set_to_fra(operand(a, i), u, v);

  set_to_unreduced(a);
}

void expr_set_to_sym(expr *a, const char *s) {
  expr_set_kind(a, kind::SYM);

  a->expr_sym = strdup(s);

  set_to_reduced(a);
}

void expr_set_op_to_sym(expr *a, size_t i, const char *s) {
  expr_set_to_sym(operand(a, i), s);

  set_to_unreduced(a);
}

void expr_set_to_undefined(expr *a) { expr_set_kind(a, kind::UNDEF); }

void expr_set_to_fail(expr *a) { expr_set_kind(a, kind::FAIL); }

void expr_set_to_inf(expr *a) { expr_set_kind(a, kind::INF); }

void expr_set_to_neg_inf(expr *a) {
  expr_set_kind(a, kind::MUL);

  a->insert(integer(-1));
  a->insert(inf());

  set_to_reduced(a);
  set_to_expanded(a);
}

bool is_inf(expr *a) { return is(a, kind::INF); }

bool is_neg_inf(expr *a) {
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

bool is_inv_inf(expr *a) {
  if (!is(a, kind::POW))
    return false;
  return (is_inf(operand(a, 0)) || is_neg_inf(operand(a, 0))) &&
         is(operand(a, 1), kind::CONST) && is_negative(operand(a, 1));
}


bool expr_raise_to_first_op(expr *a) {
  if (is(operand(a, 0), kind::INT)) {
    expr_set_to_int(a, get_val(operand(a, 0)));
    return true;
  }

  if (is(operand(a, 0), kind::SYM)) {
    expr_set_to_sym(a, get_id(operand(a, 0)));
    return true;
  }

  if (is(operand(a, 0), kind::FUNC)) {
    a->expr_info = operand(a, 0)->expr_info;

    a->expr_sym = strdup(get_id(operand(a, 0)));
    expr_set_kind(a, kind::FUNC);
    a->expr_childs = expr(a->expr_childs[0]).expr_childs;
    return true;
  }

  expr k = a->expr_childs[0];

  *a = k;

  return true;
}

bool expr_replace_with(expr *a, expr *t) {
  a->expr_info = t->expr_info;

  // if(is(a, kind::INT)) {
  // 	delete a->expr_int;
  // }

  // if(is(a, kind::LIST)) {
  // 	delete a->expr_list;
  // }

  // if(is(a, kind::SET)) {
  // 	delete a->expr_set;
  // }

  // if(is(a, kind::MAT)) {
  // 	delete a->expr_mat;
  // }

  // if(is(a, kind::SYM)) {
  // 	delete a->expr_sym;
  // }

  if (is(t, kind::INT)) {
    expr_set_to_int(a, get_val(t));
    return true;
  }

  if (is(t, kind::SYM)) {
    expr_set_to_sym(a, get_id(t));
    return true;
  }

  if (is(t, kind::MAT)) {
    expr_set_kind(a, kind::MAT);
    a->expr_mat = new matrix(*t->expr_mat);
  }

  if (is(t, kind::FUNC)) {
    expr_set_kind(a, kind::FUNC);

    a->expr_sym = strdup(get_id(t));

    a->expr_childs = std::vector<expr>();

    for (size_t i = 0; i < size_of(t); i++) {
      a->insert(t->expr_childs[i]);
    }

    return true;
  }

  expr_set_kind(a, kind_of(t));

  a->expr_childs = std::vector<expr>();

  for (size_t i = 0; i < size_of(t); i++) {
    a->insert(t->expr_childs[i]);
  }

  return true;
}


}
}
