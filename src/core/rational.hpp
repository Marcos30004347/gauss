#ifndef CORE_RATIONAL_H
#define CORE_RATIONAL_H

#include "expression.hpp"

namespace core {

expr* numerator(const expr* u);
expr* denominator(const expr* u);
expr* gcd(const expr* a, const expr* b);
bool is_rational_number_expression(const expr* u);

}
#endif