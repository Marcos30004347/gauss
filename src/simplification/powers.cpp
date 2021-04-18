#include "powers.hpp"
#include "rationals.hpp"
#include "products.hpp"
#include "core/power.hpp"

using namespace core;

namespace simplification {

expr* simplify_integer_power(const expr* v, const expr* n) {
    //(a/b)^n
    if(kind(v) == expr::FRACTION || kind(v) == expr::INTEGER) {
        expr* r = construct(expr::ALG_OP_POWER);
        include_operand(r, copy(v));
        include_operand(r, copy(n));
        return simplify_rational_number_expression(r);
    }

    // (v^0) = 1
    if(kind(n) == expr::INTEGER && integer_value(n) == 0)
        return integer(1);

    // (n^1) = n
    if(kind(n) == expr::INTEGER && integer_value(n) == 1)
        return copy(v);
    
    // (n^i)^k = n^(i*k)
    if(kind(v) == expr::ALG_OP_POWER) {
        expr* p = construct(expr::ALG_OP_PRODUCT);

        expr* r = operand(v, 0);
        expr* s = operand(v, 1);

        include_operand(p, copy(n));
        include_operand(p, copy(s));

        p = simplify_product(p);

        expr* v = construct(expr::ALG_OP_POWER);
        include_operand(v, r);
        include_operand(v, p);
    
        if(kind(p) == expr::INTEGER)
            return simplify_integer_power(operand(v,0), operand(v, 1));

        return v;
    }

    // (a*b*c)^n = a^n * b^n * c^n
    if(kind(v) == expr::ALG_OP_PRODUCT) {
        expr* r = binary_map(v, n, simplify_integer_power);
        return simplify_product(r);
    }

    expr* r = construct(expr::ALG_OP_POWER);
    include_operand(r, copy(v));
    include_operand(r, copy(n));
    return r;
}

expr* simplify_power(const expr* u) {
    const expr* v = base(u);
    const expr* n = expoent(u);
    
    if(kind(v) == expr::UNDEFINED)
        return undefined();

    if(kind(n) == expr::UNDEFINED)
        return undefined();
    
    if(kind(v) == expr::INTEGER && integer_value(v) == 0) {
        if(kind(n) == expr::INTEGER && integer_value(n) > 0)
            return integer(0); // 0^0 = 0
        else
            return undefined();  // 0^0 = 0
    }

    if(kind(v) == expr::INTEGER && integer_value(v) == 1)
        return integer(1);

    if(kind(n) == expr::INTEGER)
        return simplify_integer_power(v, n);

    return copy(u);
}


}