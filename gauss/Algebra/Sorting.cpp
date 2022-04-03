#include "Sorting.hpp"
#include "Expression.hpp"
#include "Utils.hpp"
#include <cstddef>

namespace alg {

using namespace utils;

inline bool should_revert_idx(kind ctx) { return ctx & (kind::ADD); }

inline int compare_idents(expr *a, expr *b) {
  return strcmp(get_id(a), get_id(b));
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
        // printf("    %s -- %s ", to_string(operand(a, i + j)).c_str(),
        // to_string(operand(b, i + k)).c_str()); printf("result in %i\n",
        // order);

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

      // 	printf("cmp %s -- %s\n", to_string(a).c_str(),
      // to_string(b).c_str());

      // 	for (long i = 0; i < l - r; i++) {

      // 		int order = compare(operand(a, i + j), operand(b, i +
      // k),
      // ctx); 			printf("    %s -- %s ", to_string(operand(a, i +
      // j)).c_str(), to_string(operand(b, i + k)).c_str());
      // printf("result in %i\n", order);

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
			int order = strcmp(get_func_id(a), get_func_id(b));

			if(order) return order;

			order = size_of(a) - size_of(b);

			if(order) return order;

			for(size_t i = 0; i < size_of(a); i++) {
				order = compare(operand(a, i), operand(b, i), kind::MUL);
				if(order) return order;
			}

			return 0;
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

    if (is(a, kind::MAT) && is(b, kind::MAT)) {
      return 0;
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

    if (is(a, kind::POW) && is(b, kind::MUL)) {
      if (size_of(b) == 2) {
        if (compare(a, operand(b, 1), ctx) == 0) {
          return 0;
        }
      }

      return (int)kind_of(a) - (int)kind_of(b);
    }

    if (is(b, kind::POW) && is(a, kind::MUL)) {
      if (size_of(a) == 2) {
        if (compare(b, operand(a, 1), ctx) == 0) {
          return 0;
        }
      }
      return (int)kind_of(a) - (int)kind_of(b);
    }

		if (is(a, kind::FUNC) && is(b, kind::FUNC)) {
			int order = strcmp(get_func_id(a), get_func_id(b));

			if(order) return order;

			order = size_of(a) - size_of(b);

			if(order) return order;

			for(size_t i = 0; i < size_of(a); i++) {
				order = compare(operand(a, i), operand(b, i), kind::MUL);
				if(order) return order;
			}

			return 0;
    }

    if (is(a, kind::ADD) && is(b, kind::SYM)) {
      return +1;
    }

    if (is(a, kind::SYM) && is(b, kind::ADD)) {
      return -1;
    }

    if (is(a, kind::MUL) && is(b, kind::SYM)) {
      long k = size_of(a);

      if (k > 2) {
        return -k;
      }

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

      if (k > 2) {
        return -k;
      }

      int order = compare(operand(a, 0), b, ctx);

      if (order == 0) {
        return -1;
      }

      return k > 2 ? -1 : order;
    }

    if (is(b, kind::MUL) && is(a, kind::POW)) {
      long k = size_of(b);

      if (k > 2) {
        return +k;
      }

      int order = compare(a, operand(b, 0), ctx);

      if (order == 0) {
        return +1;
      }

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

    if (is(a, kind::MAT) && is(b, kind::MAT)) {
      if (rows(a) - rows(b) == 0) {
        if (columns(a) - columns(b) == 0) {
          return 0;
        }

        return (columns(a) - columns(b)).longValue();
      }

      return (rows(a) - rows(b)).longValue();
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

  if (is(a, kind::ROOT) && is(b, kind::ROOT)) {
		int r = compare(operand(a, 1), operand(b, 1), ctx);

		if(r != 0) return r;

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

  enum kind k = is(a, kind::ADD | kind::MUL) ? kind_of(a) : kind::UNDEF;

  // NOTE: this may be causing bugs
  if (is_sorted(a, k)) {
    return;
  }

  set_to_sorted(a, k);

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

  if (is_sorted(a, k)) {
    return;
  }

  set_to_sorted(a, k);

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

}
