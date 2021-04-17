#include "core/expression.hpp"

namespace ordering {

/**
 *  Both u and v need to be and ASAE(automatically simplified algebraic expression)
 */
bool order_relation(const core::expr* u, const core::expr* v);

}