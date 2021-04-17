#include <assert.h>
#include "simplification/summations.hpp"

using namespace core;
using namespace simplification;

void should_simplify_summation() {
    expr* e0 = simplify_summation(
        summation(integer(2), integer(4))
    );
    assert(equals(e0, integer(6)));

    expr* e1 = simplify_summation(
        summation(fraction(integer(1),integer(2)), fraction(integer(1),integer(2)))
    );
    assert(equals(e1, integer(1)) == true);

    expr* e2 = simplify_summation(
        summation(summation(summation(integer(1), integer(2)), summation(integer(3), integer(4))), summation(summation(integer(5), integer(6)), summation(integer(7), integer(8))))
    );
    assert(equals(e2, integer(36)) == true);

    expr* e3 = simplify_summation(
        summation(summation(summation(symbol("x"), integer(2)), summation(integer(3), integer(4))), summation(summation(integer(5), integer(6)), summation(integer(7), integer(8))))
    );
    assert(equals(e3, summation(integer(35), symbol("x")) ) == true);

    expr* e4 = simplify_summation(summation(symbol("x"), symbol("x")));
    print(e4);
    assert(equals(e4,product(integer(2), symbol("x"))) == true);

    expr* e5 = simplify_summation(summation(symbol("x"), summation(symbol("x"), symbol("x"))));
    print(e5);
    assert(equals(e5,product(integer(3), symbol("x"))) == true);

}

int main() {
    should_simplify_summation();
    return 0;
}