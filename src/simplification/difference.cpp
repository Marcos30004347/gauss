#include "difference.hpp"
#include "products.hpp"
#include "summations.hpp"

using namespace core;

namespace simplification {

expr* simplify_difference(const expr* u) {
    expr* sum = construct(expr::ALG_OP_SUMMATION);
    include_operand(sum, operand(u,0));
    for(int i=1; i<number_of_operands(u); i++) {
        include_operand(sum, simplify_product(product(integer(-1), operand(u,i))));
    }

    return simplify_summation(sum);
}

}