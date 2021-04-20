#include "rational.hpp"
#include "power.hpp"

namespace algebra {
int gcd_rec(int a, int b) {
    if (b == 0)
        return a;

    return gcd_rec(b, a % b);
}

expression* gcd(const expression* a, const expression* b) {
    return integer(gcd_rec(integer_value(a), integer_value(b)));
}


expression* numerator(const expression* u) {
    if(kind(u) != expression::FRACTION && kind(u) != expression::ALG_OP_QUOTIENT)
        return copy(u);
    
    return operand(u, 0);
}

expression* denominator(const expression* u) {
    if(kind(u) != expression::FRACTION && kind(u) != expression::ALG_OP_QUOTIENT)
        return integer(1);
    
    return operand(u, 1);
}

bool is_rational_number_expression(const expression* u) {
    if(kind(u) == expression::INTEGER)
        return true;

    if(kind(u) == expression::FRACTION)
        return is_constant(operand(u,0)) && is_constant(operand(u,1));

    if(kind(u) == expression::ALG_OP_SUMMATION && number_of_operands(u) <= 2) {
        bool result = true;
        for(int i=0; i<number_of_operands(u); i++)
            result = result && is_rational_number_expression(operand(u,i));
        return result;
    }

    if(kind(u) == expression::ALG_OP_DIFFERENCE && number_of_operands(u) <= 2) {
        bool result = true;
        for(int i=0; i<number_of_operands(u); i++)
            result = result && is_rational_number_expression(operand(u,i));
        return result;
    }

    if(kind(u) == expression::ALG_OP_PRODUCT && number_of_operands(u) == 2) {
        bool result = true;
        for(int i=0; i<number_of_operands(u); i++)
            result = result && is_rational_number_expression(operand(u,i));
        return result;
    }

    if(kind(u) == expression::ALG_OP_QUOTIENT) {
        bool result = true;
        for(int i=0; i<number_of_operands(u); i++)
            result = result && is_rational_number_expression(operand(u,i));
        return result;
    }

    if(kind(u) == expression::ALG_OP_POWER) {
        return is_rational_number_expression(base(u)) && kind(expoent(u)) == expression::INTEGER;
    }

    return false;
}

}