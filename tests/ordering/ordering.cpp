#include <assert.h>
#include "ordering/ordering.hpp"

using namespace algebra;
using namespace ordering;

void should_order_relate_integers() {
    assert(order_relation(integer(0), integer(1)) == true);
}

void should_order_relate_symbols() {
    assert(order_relation(symbol("x"), symbol("y")) == true);
    assert(order_relation(symbol("abc"), symbol("def")) == true);
}

void should_order_relate_integers_and_symbols() {
    assert(order_relation(integer(1), symbol("x")) == true);
}

void should_order_relate_powers() {
    assert(order_relation(power(integer(2),integer(2)), power(integer(2),integer(3))) == true);
    assert(order_relation(power(integer(3),integer(2)), power(integer(4),integer(2))) == true);
}

void should_order_relate_summations() {
    assert(order_relation(summation(integer(1),integer(2)),  summation(integer(1),integer(3))) == true);
    assert(order_relation(summation(integer(1),integer(4)), summation(integer(2),integer(4))) == true);
}

void should_order_relate_products() {
    assert(order_relation(product(integer(1),integer(2)), product(integer(1),integer(3))) == true);
    assert(order_relation(product(integer(1),integer(4)), product(integer(2),integer(4))) == true);
}

int main() {
    should_order_relate_integers();
    should_order_relate_symbols();
    should_order_relate_integers_and_symbols();
    should_order_relate_powers();
    should_order_relate_summations();
    should_order_relate_products();
    return 0;
}
