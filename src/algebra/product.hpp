#ifndef ALGEBRA_PRODUCT_H
#define ALGEBRA_PRODUCT_H

#include "expression.hpp"

namespace algebra {

expression* constant_coefficient(expression* u);
expression* nonconstant_coefficient(expression* u);

}

#endif