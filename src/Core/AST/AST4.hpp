#include "Integer.hpp"
#include <cstddef>
#include <cstring>
#include <initializer_list>
#include <string>

namespace expression {

enum kind {
  factorial = (1 << 0),
  power = (1 << 1),
  mul = (1 << 2),
  add = (1 << 3),
  sub = (1 << 4),
  div = (1 << 5),
  sqrt = (1 << 6),
  infinity = (1 << 7),
  negative_infinity = (1 << 8),
  undefined = (1 << 9),
  integer = (1 << 10),
  fraction = (1 << 11),
  symbol = (1 << 12),
  fail = (1 << 13),
  // funcall = (1 << 14),
  constant = integer | fraction,
  summable = mul | power | symbol,
  multiplicable = power | symbol | add,
  nonconstant = symbol /* | funcall */,
  terminal = fail | undefined | fail | infinity | negative_infinity | symbol |
             integer | fraction,
  ordered = power | div | sqrt /* | funcall */
};

struct expr {
  kind type;

  union {
    char* id;
    Int* val;
  };

	std::vector<expr> ops;

	expr();
	expr(kind k);
	expr(expr&& other);
	expr(const expr& other);

	expr& operator=(const expr&);
	expr& operator=(expr&&);

	~expr();

	expr& operator[](size_t idx);
};

void expr_insert(expr& a, expr& b, size_t idx);
void expr_insert(expr& a, expr&& b, size_t idx);
void expr_insert(expr& a, expr& b);
void expr_insert(expr& a, expr&& b);

void expr_remove(expr& a, size_t idx);

void expr_sort(expr& a);

expr expr_create(kind kind);

expr expr_create(kind kind, std::initializer_list<expr>&&);

expr expr_symbol(const char *id);

expr expr_integer(Int value);

expr expr_fraction(Int num, Int den);

expr expr_rational(Int num, Int den);

inline int expr_is_kind(expr& a, int k) { return a.type & k; }

inline size_t expr_size(expr& ast) { return ast.ops.size(); }

inline kind expr_kind(expr& ast) { return ast.type; }

inline char *expr_id(expr& ast) { return ast.id; }

inline Int expr_value(expr& ast) { return *ast.val; }

inline char *expr_funname(expr& ast) { return ast.id; }

std::string expr_to_string(expr& a);

int expr_cmp(expr& a, expr& b, kind ctx);

expr& expr_eval(expr& a);

void expr_expand(expr& a);

void expr_print(expr& a, int tabs = 0);

}
