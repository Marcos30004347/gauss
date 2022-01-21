#include "Integer.hpp"
#include <cstddef>
#include <initializer_list>
#include <string>

namespace ast_teste {

struct ast {
  // *****************
  // Node metadata
  // *****************
	long long ref_count = 0;

  enum kind {
    fact = (1 << 0),
    pow = (1 << 1),
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
    funcall = (1 << 14),
    constant = integer | fraction,
		summable = mul | pow | symbol,
		multiplicable = pow | symbol,
    nonconstant = symbol | funcall,
		terminal = fail | undefined | fail | infinity | negative_infinity | symbol | integer | fraction,
		ordered = pow | div | sqrt | funcall
  };

  enum format {
    default_format = 0,
    poly_expr_format = 1,
  };

  const static size_t childs_margin = 8;

  // **********************
  // Node metadata members
  // **********************

  // The kind of the Node
  kind ast_kind = kind::undefined;

  // format that the ast is formated
  format ast_format = format::default_format;

  // the number of operands on child array
  size_t ast_size = 0;

  // the number of pre allocated ast* in the childs array
  size_t ast_reserved_size = 0;

  // array of childs of this node
  ast **ast_childs;

  // ******************
  // Node data members
  // ******************

  union {
    char *ast_sym;
    Int *ast_int;
  };
};

inline ast *ast_inc_ref(ast *a) {
  a->ref_count = a->ref_count + 1;
  return a;
}

inline ast *ast_dec_ref(ast *a) {
	assert(a->ref_count > 0);
  a->ref_count = a->ref_count - 1;
	return a;
}

void ast_delete(ast *a);

void ast_insert(ast *a, ast *b, size_t idx);

void ast_remove(ast *a, size_t idx);

ast* ast_sort(ast *a);

ast *ast_create(ast::kind kind);

ast *ast_create(ast::kind kind, std::initializer_list<ast *>);

ast *ast_symbol(const char *id);

ast *ast_integer(Int value);

ast *ast_fraction(Int num, Int den);

ast *ast_rational(Int num, Int den);

ast *ast_copy(ast *a);

inline ast *ast_operand(ast *a, size_t i) { return a->ast_childs[i]; }

inline int ast_is_kind(ast *a, int k) { return a->ast_kind & k; }

inline size_t ast_size(ast *ast) { return ast->ast_size; }

inline ast::kind ast_kind(ast *ast) { return ast->ast_kind; }

inline char *ast_id(ast *ast) { return ast->ast_sym; }

inline Int ast_value(ast *ast) { return Int(*ast->ast_int); }

inline char *ast_funname(ast *ast) { return ast->ast_sym; }

std::string ast_to_string(ast *a);

int ast_cmp(ast *a, ast *b, ast::kind ctx);

ast* ast_eval(ast* a, bool print = false, ast* parent = nullptr);

ast* ast_expand(ast* a);

void ast_print(ast* a, int tabs = 0);

} // namespace ast_teste
