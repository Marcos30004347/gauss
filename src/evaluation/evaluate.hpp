#include "core/expression.hpp"

namespace evaluation {

core::expr* evaluate_product(core::expr* u, core::expr* v);
core::expr* evaluate_summation(core::expr* u, core::expr* v);
core::expr* evaluate_difference(core::expr* u, core::expr* v);
core::expr* evaluate_quotient(core::expr* v, core::expr* w);
core::expr* evaluate_power(core::expr* v, core::expr* w);


} // simplification