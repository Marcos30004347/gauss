#include <assert.h>
#include <cstdio>

#include "simplification/products.hpp"
#include "simplification/summations.hpp"

using namespace core;
using namespace simplification;

void should_simplify_products() {
    expr* e0 = simplify_product(product(integer(3), integer(2)));
    print(e0);
    printf("\n");

    expr* e1 = simplify_product(product(product(integer(1), integer(2)), product(integer(3), integer(4))));
    print(e1);
    printf("\n");


    expr* e2 = simplify_product(
        product(
            product(
                product(integer(1), integer(2)),
                product(integer(3), integer(4))
            ),
            product(
                product(integer(5), integer(6)),
                product(integer(7), integer(8))
            )
        )
    );

    print(e2);
    printf("\n");


    expr* e3 = simplify_product(
        product(
            summation(integer(2), integer(3)),
            summation(integer(4), integer(5))
        )
    );

    print(e3);
    printf("\n");

    expr* e4 = simplify_product(
        product(
            summation(symbol("x"), integer(3)),
            summation(symbol("x"), integer(5))
        )
    );

    print(e4);
    printf("\n");



    expr* e5 = simplify_product(
        product(
            product(
                symbol("x"),
                symbol("y")
            ),
            product(
                summation(symbol("x"), integer(3)),
                summation(symbol("x"), integer(5))
            )
        )
    );

    print(e5);
    printf("\n");


    // expr* e6 =  simplify_product(
    //     product(
    //         product(
    //             summation(symbol("x"), symbol("x")),
    //             summation(symbol("x"), symbol("y"))
    //         ),
    //         product(
    //             summation(symbol("x"), integer(3)),
    //             summation(symbol("x"), integer(5))
    //         )
    //     )
    // );

    // print(e6);
    // printf("\n");

}


int main() {
    should_simplify_products();
    return 0;
}