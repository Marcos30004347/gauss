#include "evaluate.hpp"
#include "algebra/rational.hpp"
#include "simplification/rationals.hpp"

using namespace algebra;
using namespace simplification;

namespace evaluation {


expression* remainder(const expression* u, const expression* v) {
    return integer(integer_value(u) % integer_value(v));
}

expression* evaluate_product(expression* u, expression* v) {
    if(is_constant(u) && is_constant(v)) {
        if(kind(u) == expression::INTEGER && kind(v) == expression::INTEGER)
            return integer(integer_value(u) * integer_value(v));

        
        return simplify_rational_number_expression(fraction(
            evaluate_product(numerator(u), numerator(v)),
            evaluate_product(denominator(u), denominator(v))  
        ));
    }

    return product(u, v);
}

expression* evaluate_summation(expression* u, expression* v) {
    if(kind(u) == expression::INTEGER && kind(v) == expression::INTEGER)
        return integer(integer_value(u) + integer_value(v));

    if(kind(u) == expression::FRACTION && kind(v) == expression::INTEGER)
        return simplify_rational_number_expression(
            fraction(
                evaluate_summation(
                    numerator(u),
                    evaluate_product(v, denominator(u))
                ),
                denominator(u)
            )
        );

    if(kind(v) == expression::FRACTION && kind(u) == expression::INTEGER)
        return simplify_rational_number_expression(
            fraction(
                evaluate_summation(
                    numerator(v),
                    evaluate_product(u, denominator(v))
                ),
                denominator(v)
            )
        );

    if(kind(v) == expression::FRACTION && kind(u) == expression::FRACTION)
        return simplify_rational_number_expression(fraction(
            evaluate_summation(
                evaluate_product(numerator(u), denominator(v)),
                evaluate_product(denominator(u), numerator(v))
            ),
            evaluate_product(denominator(u), denominator(v))  
        ));

    return summation(u, v);
}

expression* evaluate_difference(expression* u, expression* v) {
    if(kind(u) == expression::INTEGER && kind(v) == expression::INTEGER)
        return integer(integer_value(u) - integer_value(v));

    if(kind(u) == expression::FRACTION && kind(v) == expression::INTEGER)
        return simplify_rational_number_expression(
            fraction(
                evaluate_difference(
                    numerator(u),
                    evaluate_product(v, denominator(u))
                ),
                denominator(u)
            )
        );

    if(kind(v) == expression::FRACTION && kind(u) == expression::INTEGER)
        return simplify_rational_number_expression(
            fraction(
                evaluate_difference(
                    numerator(v),
                    evaluate_product(u, denominator(v))
                ),
                denominator(v)
            )
        );

    if(kind(v) == expression::FRACTION && kind(u) == expression::FRACTION)
        return simplify_rational_number_expression(fraction(
            evaluate_difference(
                evaluate_product(numerator(u), denominator(v)),
                evaluate_product(denominator(u), numerator(v))
            ),
            evaluate_product(denominator(u), denominator(v))  
        ));

    return difference(u, v);
}

expression* evaluate_quotient(expression* v, expression* w) {
    if(!is_constant(v) || !is_constant(w))
        return quotient(v, w);

    if(equals(numerator(w), integer(0)))
        return undefined();

    if(equals(numerator(w), integer(1)))
        return copy(v);

    return simplify_rational_number_expression(
        fraction(
            evaluate_product(numerator(v), denominator(w)),
            evaluate_product(numerator(w), denominator(v))
        )
    );
}

expression* evaluate_power(expression* v, expression* w) {
    if(kind(w) != expression::INTEGER)
        return power(v,w);
    
    long long n = integer_value(w);

    if(!is_constant(v))
        return power(v, integer(n));

    if(!equals(numerator(v), integer(0))) {
        if(n > 0) {
            expression* s = evaluate_power(v, integer(n - 1));
            return evaluate_product(s, v);
        } else if(n == 0) {
            return (integer(1));
        } else if(n == -1) {
            return simplify_rational_number_expression(
                fraction(denominator(v), numerator(v))
            );
        } else {
            expression* s = simplify_rational_number_expression(
                fraction(denominator(v), numerator(v))
            );
            return evaluate_power(s, integer(-1 * n));
        }
    }

    if(n >= 1) 
        return integer(0);
    else if(n <= 0) 
        return undefined();
    
    return power(v, w);
}


}