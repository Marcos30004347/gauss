#include "rational.hpp"
#include "power.hpp"

namespace core {
int gcd_rec(int a, int b) {
    if (b == 0)
        return a;

    return gcd_rec(b, a % b);
}

expr* gcd(const expr* a, const expr* b) {
    return integer(gcd_rec(integer_value(a), integer_value(b)));
}


expr* numerator(const expr* u) {
    if(kind(u) != expr::FRACTION && kind(u) != expr::ALG_OP_QUOTIENT)
        return copy(u);
    
    return operand(u, 0);
}

expr* denominator(const expr* u) {
    if(kind(u) != expr::FRACTION && kind(u) != expr::ALG_OP_QUOTIENT)
        return integer(1);
    
    return operand(u, 1);
}

bool is_rational_number_expression(const expr* u) {
    if(kind(u) == expr::INTEGER)
        return true;

    if(kind(u) == expr::FRACTION)
        return is_constant(operand(u,0)) && is_constant(operand(u,1));

    if(kind(u) == expr::ALG_OP_SUMMATION && number_of_operands(u) <= 2) {
        bool result = true;
        for(int i=0; i<number_of_operands(u); i++)
            result = result && is_rational_number_expression(operand(u,i));
        return result;
    }

    if(kind(u) == expr::ALG_OP_DIFFERENCE && number_of_operands(u) <= 2) {
        bool result = true;
        for(int i=0; i<number_of_operands(u); i++)
            result = result && is_rational_number_expression(operand(u,i));
        return result;
    }

    if(kind(u) == expr::ALG_OP_PRODUCT && number_of_operands(u) == 2) {
        bool result = true;
        for(int i=0; i<number_of_operands(u); i++)
            result = result && is_rational_number_expression(operand(u,i));
        return result;
    }

    if(kind(u) == expr::ALG_OP_QUOTIENT) {
        bool result = true;
        for(int i=0; i<number_of_operands(u); i++)
            result = result && is_rational_number_expression(operand(u,i));
        return result;
    }

    if(kind(u) == expr::ALG_OP_POWER) {
        return is_rational_number_expression(base(u)) && kind(expoent(u)) == expr::INTEGER;
    }

    return false;
}

}