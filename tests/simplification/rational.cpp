#include <assert.h>
#include "simplification/rationals.hpp"

using namespace algebra;
using namespace simplification;

void should_simplify_rational_number_expression() {
    expression* e0 = simplify_rational_number_expression(
        fraction(integer(2), integer(4))
    );
    assert(equals(e0, fraction(integer(1), integer(2))));

    expression* e1 = simplify_rational_number_expression(
        fraction(integer(4), integer(2))
    );
    assert(equals(e1, integer(2)));

    expression* e2 = simplify_rational_number_expression(
        quotient(summation(integer(4), integer(4)), integer(2))
    );
    assert(equals(e2, integer(4)));

    expression* e3 = simplify_rational_number_expression(
        quotient(summation(integer(4), integer(4)), summation(integer(4), integer(4)))
    );
    assert(equals(e3, integer(1)));

    expression* e4 = simplify_rational_number_expression(
        quotient(product(integer(4), integer(4)), summation(integer(4), integer(4)))
    );
    assert(equals(e4, integer(2)));

    expression* e5 = simplify_rational_number_expression(
        quotient(product(integer(4), integer(4)), fraction(integer(4), integer(4)))
    );
    assert(equals(e5, integer(16)));
}

int main() {
    should_simplify_rational_number_expression();
    return 0;
}