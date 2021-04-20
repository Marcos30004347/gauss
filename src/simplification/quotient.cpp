#include "quotient.hpp"
#include "products.hpp"
#include "powers.hpp"

using namespace algebra;

namespace simplification {

// a/b = a * b^-1
expression* simplify_quotient(const expression* u) {
    expression* u1 = operand(u,0);
    expression* u2 = operand(u,1);
    return simplify_product(product(u1, simplify_power(power(u2, integer(-1)))));
}

}