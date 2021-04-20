#ifndef ORDERING_H
#define ORDERING_H

#include "algebra/expression.hpp"

namespace ordering {

/**
 *  Both u and v need to be and ASAE(automatically simplified algebraic expression)
 */
bool order_relation(const algebra::expression* u, const algebra::expression* v);

}

#endif
