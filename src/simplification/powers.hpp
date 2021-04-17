#ifndef SIMPLIFICATION_POWER_H
#define SIMPLIFICATION_POWER_H

#include "core/expression.hpp"

namespace simplification {

// core::expr* simplify_integer_power(const core::expr* v, const core::expr* n);
core::expr* simplify_power(const core::expr* u);

}

#endif