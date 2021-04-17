#include <assert.h>
#include "evaluation/evaluate.hpp"

using namespace core;
using namespace evaluation;

void should_evaluate_power() {
    expr* b0 = integer(2);
    expr* e0 = integer(2);
    expr* r0 = evaluate_power(b0, e0);
    assert(equals(r0, integer(4)) == true);

    expr* b1 = fraction(integer(2), integer(2));
    expr* e1 = integer(2);
    expr* r1 = evaluate_power(b1, e1);
    assert(equals(r1, integer(1)) == true);

    expr* b2 = fraction(integer(4), integer(2));
    expr* e2 = integer(2);
    expr* r2 = evaluate_power(b2, e2);
    assert(equals(r2, integer(4)) == true);

    expr* b3 = fraction(integer(4), integer(2));
    expr* e3 = integer(-1);
    expr* r3 = evaluate_power(b3, e3);
    assert(equals(r3, fraction(integer(1), integer(2))) == true);
}

void should_evaluate_summations() {
    expr* b0 = integer(2);
    expr* e0 = integer(2);
    expr* r0 = evaluate_summation(b0, e0);
    assert(equals(r0, integer(4)) == true);

    expr* b1 = fraction(integer(2), integer(2));
    expr* e1 = integer(2);
    expr* r1 = evaluate_summation(b1, e1);
    assert(equals(r1, integer(3)) == true);

    expr* b2 = fraction(integer(1), integer(2));
    expr* e2 = fraction(integer(1), integer(2));
    expr* r2 = evaluate_summation(b2, e2);
    assert(equals(r2, integer(1)) == true);


    expr* b3 = fraction(integer(2), integer(4));
    expr* e3 = fraction(integer(3), integer(2));
    expr* r3 = evaluate_summation(b3, e3);
    assert(equals(r3, integer(2)) == true);

    expr* b4 = fraction(integer(1), integer(2));
    expr* e4 = fraction(integer(3), integer(5));
    expr* r4 = evaluate_summation(b4, e4);
    assert(equals(r4, fraction(integer(11), integer(10))) == true);
}

void should_evaluate_differences() {
    expr* b0 = integer(2);
    expr* e0 = integer(2);
    expr* r0 = evaluate_difference(b0, e0);
    assert(equals(r0, integer(0)) == true);

    expr* b1 = fraction(integer(2), integer(2));
    expr* e1 = integer(2);
    expr* r1 = evaluate_difference(b1, e1);
    assert(equals(r1, integer(-1)) == true);

    expr* b2 = fraction(integer(1), integer(2));
    expr* e2 = fraction(integer(1), integer(2));
    expr* r2 = evaluate_difference(b2, e2);
    assert(equals(r2, integer(0)) == true);

    expr* b3 = fraction(integer(2), integer(4));
    expr* e3 = fraction(integer(3), integer(2));
    expr* r3 = evaluate_difference(b3, e3);
    assert(equals(r3, integer(-1)) == true);

    expr* b4 = fraction(integer(1), integer(2));
    expr* e4 = fraction(integer(3), integer(5));
    expr* r4 = evaluate_difference(b4, e4);
    assert(equals(r4, fraction(integer(-1), integer(10))) == true);
}

void should_evaluate_products() {
    expr* b0 = integer(2);
    expr* e0 = integer(2);
    expr* r0 = evaluate_product(b0, e0);
    assert(equals(r0, integer(4)) == true);

    expr* b1 = integer(2);
    expr* e1 = fraction(integer(2), integer(4));
    expr* r1 = evaluate_product(b1, e1);
    assert(equals(r1, integer(1)) == true);

    expr* b2 = fraction(integer(1), integer(7));
    expr* e2 = fraction(integer(5), integer(3));
    expr* r2 = evaluate_product(b2, e2);
    assert(equals(r2, fraction(integer(5), integer(21))) == true);
}

void should_evaluate_quotients() {
    expr* b0 = integer(2);
    expr* e0 = integer(2);
    expr* r0 = evaluate_quotient(b0, e0);
    assert(equals(r0, integer(1)) == true);

    expr* b1 = integer(2);
    expr* e1 = fraction(integer(2), integer(4));
    expr* r1 = evaluate_quotient(b1, e1);
    assert(equals(r1, integer(4)) == true);

    expr* b2 = fraction(integer(2), integer(4));
    expr* e2 = integer(2);
    expr* r2 = evaluate_quotient(b2, e2);
    assert(equals(r2, fraction(integer(1), integer(4))) == true);

    expr* b3 = fraction(integer(1), integer(7));
    expr* e3 = fraction(integer(5), integer(3));
    expr* r3 = evaluate_quotient(b3, e3);
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