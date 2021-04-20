#ifndef SIMPLIFICATION_POWER_H
#define SIMPLIFICATION_POWER_H

#include "algebra/expression.hpp"

namespace simplification {

// core::expression* simplify_integer_power(const core::expression* v, const core::expression* n);
algebra::expression* simplify_power(const algebra::expression* u);

}

#endif