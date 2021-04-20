#include "powers.hpp"
#include "rationals.hpp"
#include "products.hpp"
#include "algebra/power.hpp"

using namespace algebra;

namespace simplification {

expression* simplify_integer_power(const expression* v, const expression* n) {
    //(a/b)^n
    if(kind(v) == expression::FRACTION || kind(v) == expression::INTEGER) {
        expression* r = construct(expression::ALG_OP_POWER);
        include_operand(r, copy(v));
        include_operand(r, copy(n));
        return simplify_rational_number_expression(r);
    }

    // (v^0) = 1
    if(kind(n) == expression::INTEGER && integer_value(n) == 0)
        return integer(1);

    // (n^1) = n
    if(kind(n) == expression::INTEGER && integer_value(n) == 1)
        return copy(v);
    
    // (n^i)^k = n^(i*k)
    if(kind(v) == expression::ALG_OP_POWER) {
        expression* p = construct(expression::ALG_OP_PRODUCT);

        expression* r = operand(v, 0);
        expression* s = operand(v, 1);

        include_operand(p, copy(n));
        include_operand(p, copy(s));

        p = simplify_product(p);

        expression* v = construct(expression::ALG_OP_POWER);
        include_operand(v, r);
        include_operand(v, p);
    
        if(kind(p) == expression::INTEGER)
            return simplify_integer_power(operand(v,0), operand(v, 1));

        return v;
    }

    // (a*b*c)^n = a^n * b^n * c^n
    if(kind(v) == expression::ALG_OP_PRODUCT) {
        expression* r = binary_map(v, n, simplify_integer_power);
        return simplify_product(r);
    }

    expression* r = construct(expression::ALG_OP_POWER);
    include_operand(r, copy(v));
    include_operand(r, copy(n));
    return r;
}

expression* simplify_power(const expression* u) {
    const expression* v = base(u);
    const expression* n = expoent(u);
    
    if(kind(v) == expression::UNDEFINED)
        return undefined();

    if(kind(n) == expression::UNDEFINED)
        return undefined();
    
    if(kind(v) == expression::INTEGER && integer_value(v) == 0) {
        if(kind(n) == expression::INTEGER && integer_value(n) > 0)
            return integer(0); // 0^0 = 0
        else
            return undefined();  // 0^0 = 0
    }

    if(kind(v) == expression::INTEGER && integer_value(v) == 1)
        return integer(1);

    if(kind(n) == expression::INTEGER)
        return simplify_integer_power(v, n);

    return copy(u);
}


}