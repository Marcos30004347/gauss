#ifndef ALGEBRA_RATIONAL_H
#define ALGEBRA_RATIONAL_H

#include "expression.hpp"

namespace algebra {

expression* numerator(const expression* u);
expression* denominator(const expression* u);
expression* gcd(const expression* a, const expression* b);

bool is_rational_number_expression(const expression* u);

}
#endif