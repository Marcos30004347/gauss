#include <assert.h>
#include "evaluation/evaluate.hpp"

using namespace algebra;

using namespace evaluation;

void should_evaluate_power() {
    expression* b0 = integer(2);
    expression* e0 = integer(2);
    expression* r0 = evaluate_power(b0, e0);
    assert(equals(r0, integer(4)) == true);

    expression* b1 = fraction(integer(2), integer(2));
    expression* e1 = integer(2);
    expression* r1 = evaluate_power(b1, e1);
    assert(equals(r1, integer(1)) == true);

    expression* b2 = fraction(integer(4), integer(2));
    expression* e2 = integer(2);
    expression* r2 = evaluate_power(b2, e2);
    assert(equals(r2, integer(4)) == true);

    expression* b3 = fraction(integer(4), integer(2));
    expression* e3 = integer(-1);
    expression* r3 = evaluate_power(b3, e3);
    assert(equals(r3, fraction(integer(1), integer(2))) == true);
}

void should_evaluate_summations() {
    expression* b0 = integer(2);
    expression* e0 = integer(2);
    expression* r0 = evaluate_summation(b0, e0);
    assert(equals(r0, integer(4)) == true);

    expression* b1 = fraction(integer(2), integer(2));
    expression* e1 = integer(2);
    expression* r1 = evaluate_summation(b1, e1);
    assert(equals(r1, integer(3)) == true);

    expression* b2 = fraction(integer(1), integer(2));
    expression* e2 = fraction(integer(1), integer(2));
    expression* r2 = evaluate_summation(b2, e2);
    assert(equals(r2, integer(1)) == true);


    expression* b3 = fraction(integer(2), integer(4));
    expression* e3 = fraction(integer(3), integer(2));
    expression* r3 = evaluate_summation(b3, e3);
    assert(equals(r3, integer(2)) == true);

    expression* b4 = fraction(integer(1), integer(2));
    expression* e4 = fraction(integer(3), integer(5));
    expression* r4 = evaluate_summation(b4, e4);
    assert(equals(r4, fraction(integer(11), integer(10))) == true);
}

void should_evaluate_differences() {
    expression* b0 = integer(2);
    expression* e0 = integer(2);
    expression* r0 = evaluate_difference(b0, e0);
    assert(equals(r0, integer(0)) == true);

    expression* b1 = fraction(integer(2), integer(2));
    expression* e1 = integer(2);
    expression* r1 = evaluate_difference(b1, e1);
    assert(equals(r1, integer(-1)) == true);

    expression* b2 = fraction(integer(1), integer(2));
    expression* e2 = fraction(integer(1), integer(2));
    expression* r2 = evaluate_difference(b2, e2);
    assert(equals(r2, integer(0)) == true);

    expression* b3 = fraction(integer(2), integer(4));
    expression* e3 = fraction(integer(3), integer(2));
    expression* r3 = evaluate_difference(b3, e3);
    assert(equals(r3, integer(-1)) == true);

    expression* b4 = fraction(integer(1), integer(2));
    expression* e4 = fraction(integer(3), integer(5));
    expression* r4 = evaluate_difference(b4, e4);
    assert(equals(r4, fraction(integer(-1), integer(10))) == true);
}

void should_evaluate_products() {
    expression* b0 = integer(2);
    expression* e0 = integer(2);
    expression* r0 = evaluate_product(b0, e0);
    assert(equals(r0, integer(4)) == true);

    expression* b1 = integer(2);
    expression* e1 = fraction(integer(2), integer(4));
    expression* r1 = evaluate_product(b1, e1);
    assert(equals(r1, integer(1)) == true);

    expression* b2 = fraction(integer(1), integer(7));
    expression* e2 = fraction(integer(5), integer(3));
    expression* r2 = evaluate_product(b2, e2);
    assert(equals(r2, fraction(integer(5), integer(21))) == true);
}

void should_evaluate_quotients() {
    expression* b0 = integer(2);
    expression* e0 = integer(2);
    expression* r0 = evaluate_quotient(b0, e0);
    assert(equals(r0, integer(1)) == true);

    expression* b1 = integer(2);
    expression* e1 = fraction(integer(2), integer(4));
    expression* r1 = evaluate_quotient(b1, e1);
    assert(equals(r1, integer(4)) == true);

    expression* b2 = fraction(integer(2), integer(4));
    expression* e2 = integer(2);
    expression* r2 = evaluate_quotient(b2, e2);
    assert(equals(r2, fraction(integer(1), integer(4))) == true);

    expression* b3 = fraction(integer(1), integer(7));
    expression* e3 = fraction(integer(5), integer(3));
    expression* r3 = evaluate_quotient(b3, e3);
    assert(equals(r3, fraction(integer(3), integer(35))) == true);
}


int main() {
    should_evaluate_power();
    should_evaluate_summations();
    should_evaluate_differences();
    should_evaluate_products();
    should_evaluate_quotients();
    return 0;
}