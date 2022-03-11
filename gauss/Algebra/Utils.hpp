#include "Expression.hpp"

namespace alg {
namespace utils {
Int rows(expr *a);
Int columns(expr *a);
bool is_sorted(expr *a, enum kind k);
bool is_reduced(expr *a);
bool expr_is_zero(expr *a);
void set_to_reduced(expr *a);
void set_to_unsorted(expr *a);
void set_to_unreduced(expr *a) ;
void set_to_unexpanded(expr *a);
void set_to_sorted(expr *a, enum kind k);
bool is_negative(expr *a);
void set_to_expanded(expr *a);
void expr_set_kind(expr *a, kind kind);
void expr_set_to_int(expr *a, Int v);
void expr_set_to_mat(expr *a, matrix *v);
void expr_set_op_to_int(expr *a, size_t i, Int v);
void expr_set_to_fra(expr *a, Int u, Int v);
void expr_set_op_to_fra(expr *a, size_t i, Int u, Int v);
void expr_set_to_sym(expr *a, const char *s);
void expr_set_op_to_sym(expr *a, size_t i, const char *s);
void expr_set_to_undefined(expr *a);
void expr_set_to_fail(expr *a);
void expr_set_to_inf(expr *a);
void expr_set_to_neg_inf(expr *a);
bool is_inf(expr *a);
bool is_neg_inf(expr *a);
bool is_inv_inf(expr *a);
bool expr_raise_to_first_op(expr *a);
bool expr_replace_with(expr *a, expr *t);

}
}
