#ifndef EVALUATE_H
#define EVALUATE_H

#include "algebra/expression.hpp"

namespace evaluation {

algebra::expression* evaluate_product(algebra::expression* u, algebra::expression* v);
algebra::expression* evaluate_summation(algebra::expression* u, algebra::expression* v);
algebra::expression* evaluate_difference(algebra::expression* u, algebra::expression* v);
algebra::expression* evaluate_quotient(algebra::expression* v, algebra::expression* w);
algebra::expression* evaluate_power(algebra::expression* v, algebra::expression* w);


} // simplification

#endif