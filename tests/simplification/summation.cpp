#include <assert.h>
#include <cstdio>
#include "simplification/summations.hpp"

using namespace core;
using namespace simplification;

void should_simplify_summation() {
    expr* e0 = simplify_summation(
        summation(integer(2), integer(4))
    );
    print(e0);
    printf("\n");
    assert(equals(e0, integer(6)));

    expr* e1 = simplify_summation(
        summation(fraction(integer(1),integer(2)), fraction(integer(1),integer(2)))
    );
    print(e1);
    printf("\n");
    assert(equals(e1, integer(1)) == true);

    expr* e2 = simplify_summation(
        summation(
            summation(
                summation(integer(1), integer(2)),
                summation(integer(3), integer(4))
            ),
            summation(
                summation(integer(5), integer(6)),
                summation(integer(7), integer(8))
            )
        )
    );
    print(e2);
    printf("\n");
    assert(equals(e2, integer(36)) == true);

    expr* e4 = simplify_summation(
        summation(
            summation(
                summation(symbol("x"), integer(2)),
                summation(integer(3), integer(4))
            ),
            summation(
                summation(integer(5), symbol("x")),
                summation(integer(7), symbol("x"))
            )
        )
    );

    print(e4);
    printf("\n");
    assert(equals(e4, summation(integer(21), product(integer(3), symbol("x")))) == true);



    expr* e5 = simplify_summation(
        summation(
            summation(
                summation(symbol("x"), symbol("y")),
                summation(integer(3), integer(4))
            ),
            summation(
                summation(symbol("y"), symbol("x")),
                summation(integer(7), symbol("x"))
            )
        )
    );

    print(e5);
    printf("\n");

    assert(equals(e5, summation({integer(14), product(integer(3), symbol("x")), product(integer(2), symbol("y"))})) == true);

    expr* e6 = simplify_summation(
        summation(
            summation(
                summation(symbol("x"), symbol("y")),
                summation(integer(1), integer(2))
            ),
            summation(
                summation(symbol("y"), symbol("x")),
                summation(
                    summation(
                        summation(integer(3), integer(4)),
                        integer(5)
                    ),
                    summation(
                        summation(integer(6), integer(7)),
                        integer(8)
                    )
                )
            )
        )
    );

    print(e6);
    printf("\n");




    // expr* e3 = simplify_summation(
    //     summation(summation(summation(symbol("x"), integer(2)), summation(integer(3), integer(4))), summation(summation(integer(5), integer(6)), summation(integer(7), integer(8))))
    // );
    // assert(equals(e3, summation(integer(35), symbol("x")) ) == true);

    // expr* e4 = simplify_summation(summation(symbol("x"), symbol("x")));
    // print(e4);
    // assert(equals(e4,product(integer(2), symbol("x"))) == true);

    // expr* e5 = simplify_summation(summation(symbol("x"), summation(symbol("x"), symbol("x"))));
    // print(e5);
    // assert(equals(e5,product(integer(3), symbol("x"))) == true);

}

int main() {
    should_simplify_summation();
    return 0;
}