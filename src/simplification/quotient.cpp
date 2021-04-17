#include "quotient.hpp"
#include "products.hpp"
#include "powers.hpp"

using namespace core;

namespace simplification {

// a/b = a * b^-1
expr* simplify_quotient(const expr* u) {
    expr* u1 = operand(u,0);
    expr* u2 = operand(u,1);
    return simplify_product(product(u1, simplify_power(power(u2, integer(-1)))));
}

}