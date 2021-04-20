#include <assert.h>

#include "rationals.hpp"
#include "evaluation/evaluate.hpp"
#include "algebra/rational.hpp"

using namespace algebra;
using namespace evaluation;

namespace simplification {


expression* simplyfy_quotient(const expression* u, const expression* v) {
    return integer(integer_value(u) / integer_value(v));
}


expression* simplify_rational_number(const expression* u) {
    if(kind(u) == expression::INTEGER)
        return copy(u);

    if(kind(u) == expression::FRACTION) {
        expression* n = operand(u, 0);
        expression* d = operand(u, 1);
    
        assert(is_constant(u));
        assert(is_constant(d));

        if(equals(d, integer(1)))
            return copy(n);

        if(integer_value(n) % integer_value(d) == 0)
            return simplyfy_quotient(n, d);
        else {
            expression* g = gcd(n, d);
            if(integer_value(d) > 0) {
                if(integer_value(g) > 0)
                    return fraction(
                        simplyfy_quotient(n, g),
                        simplyfy_quotient(d, g)
                    );
                else 
                    return fraction(
                        simplyfy_quotient(n, evaluate_product(integer(-1), g)),
                        simplyfy_quotient(d, evaluate_product(integer(-1), g))
                    );
            }
            else if(integer_value(d) < 0) {
                // maybe avaluate
                return fraction(
                    simplyfy_quotient(evaluate_product(n, integer(-1)), g),
                    simplyfy_quotient(evaluate_product(d, integer(-1)), g)
                );
            }
        }
    }

    return copy(u);    
}

expression* simplify_rational_number_expression_rec(const expression* u) {
    assert(is_rational_number_expression(u) == true);

    if(kind(u) == expression::INTEGER)
        return copy(u);

    if(kind(u) == expression::FRACTION) {
        if(equals(denominator(u), integer(0))) {
            return undefined();
        } else {
            return copy(u);
        }
    }

    if(number_of_operands(u) == 1) {
        expression* v = simplify_rational_number_expression_rec(operand(u,0));
        if(kind(v) == expression::UNDEFINED)
            return undefined();
        else if(kind(u) == expression::ALG_OP_SUMMATION)
            return v;
        else if(kind(u) == expression::ALG_OP_DIFFERENCE)
            return evaluate_product(integer(-1), v);
        else
            return copy(u);
    }

    if(number_of_operands(u) < 2)
        return copy(operand(u,0));

    if(
        kind(u) == expression::ALG_OP_SUMMATION ||
        kind(u) == expression::ALG_OP_PRODUCT ||
        kind(u) == expression::ALG_OP_DIFFERENCE ||
        kind(u) == expression::ALG_OP_QUOTIENT
    ) {
        expression* v = simplify_rational_number_expression_rec(operand(u,0));
        expression* w = simplify_rational_number_expression_rec(operand(u,1));

        if(kind(v) == expression::UNDEFINED || kind(w) == expression::UNDEFINED)
            return undefined();
        
        if(kind(u) == expression::ALG_OP_SUMMATION)
            return evaluate_summation(v,w);

        if(kind(u) == expression::ALG_OP_DIFFERENCE)
            return evaluate_difference(v,w);

        if(kind(u) == expression::ALG_OP_PRODUCT)
            return evaluate_product(v,w);

        if(kind(u) == expression::ALG_OP_QUOTIENT)
            return evaluate_quotient(v,w);

    } else if(kind(u) == expression::ALG_OP_POWER) {
        expression* v = simplify_rational_number_expression_rec(operand(u, 0));

        if(kind(v) == expression::UNDEFINED)
            return undefined();
        else if(kind(operand(u,1)) == expression::INTEGER)
            return evaluate_power(v, operand(u,1));
        else
            return copy(u);
    }

    return copy(u);
}

expression* simplify_rational_number_expression(const expression* u) {
    if(!is_rational_number_expression(u))
        return copy(u);

    expression* v = simplify_rational_number_expression_rec(u);

    if(kind(v) == expression::UNDEFINED)
        return undefined();

    return simplify_rational_number(v);
}

}
