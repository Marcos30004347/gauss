
#include <vector>
#include <assert.h>
#include "products.hpp"
#include "powers.hpp"
#include "summations.hpp"
#include "rationals.hpp"
#include "ordering/ordering.hpp"
#include "algebra/product.hpp"

#include <cstdio>

using namespace algebra;
using namespace ordering;

namespace simplification {

std::vector<expression*> simplify_summation_rec(std::vector<expression*> L);
std::vector<expression*> adjoin_summations(expression* p, std::vector<expression*> q);
std::vector<expression*> merge_summations(std::vector<expression*> p, std::vector<expression*> q);

std::vector<expression*> transform_summation(expression* a, expression* b);
std::vector<expression*> transform_summation_associative(std::vector<expression*> L);
std::vector<expression*> transform_summation_expand(std::vector<expression*> L);

// bool is_unifyable(expression* u, expression* v);

// bool is_summation_unifyable(expression* u, expression* v) {
//     return (
//         (is_constant(u) && is_constant(v)) ||
//         equals(
//             nonconstant_coefficient(u),
//             nonconstant_coefficient(v)
//         )
//     );
// }

std::vector<expression*> adjoin_summations(expression* p, std::vector<expression*> q) {
    std::vector<expression*> tmp = std::vector<expression*>(0);

    int i = 0;
    bool inc = false;
    
    while(i != q.size()) {
        // if(!inc && is_summation_unifyable(p, q[i])) {
        //     std::vector<core::expression *> res = transform_summation(p, q[i]);
        //     inc = true;
        //     tmp.push_back(res[0]);
        //     i++;
        // } else 
        if(!inc && order_relation(p, q[i])) {
            inc = true;
            tmp.push_back(p);
        }
        else {
            tmp.push_back(q[i++]);
        }
    }

    if(!q.size() || !inc) {
        tmp.push_back(p);
    }
    
    return tmp;
}

std::vector<expression*> merge_summations(std::vector<expression*> p, std::vector<expression*> q) {
    if(p.size() == 0)
        return q;

    if(q.size() == 0)
        return p;

    std::vector<expression*> H = simplify_summation_rec({p[0], q[0]});
    std::vector<expression*> R = merge_summations(rest(p), rest(q));

    for(int i=0; i<H.size(); i++) {
        // printf("adjoining ");
        // print(H[i]);
        // printf(" with [");
        // for(int j=0; j<R.size(); j++) {
        //     print(R[j]);
        //     if(j != R.size() - 1)
        //         printf(", ");
        // }
        // printf("]\n");
        R = adjoin_summations(H[i], R);
    }

    return (R);
}

std::vector<expression*> transform_summation(expression* a, expression* b) {
    // printf("a: ");
    // print(a);
    // printf("\nb: ");
    // print(b);
    // printf("\n");
    // printf("\n");

    if(is_constant(a) && is_constant(b))
        return { simplify_rational_number_expression(summation(a, b)) };
   
    // printf("\n ncst a: ");
    // // print(constant_coefficient(a));
    // // printf("\n");
    // print(nonconstant_coefficient(a));
    // printf("\n ncst b: ");
    // // print(constant_coefficient(b));
    // // printf("\n");
    // print(nonconstant_coefficient(b));
    // printf("\n");
    // printf("%i %i\n", kind(nonconstant_coefficient(a)), kind(nonconstant_coefficient(b)));
    // printf("%i\n", equals(
    //     nonconstant_coefficient(a),
    //     nonconstant_coefficient(b)
    // ));

    if(equals(
        nonconstant_coefficient(a),
        nonconstant_coefficient(b)
    )) {
        return {
            simplify_product(
                product(
                    simplify_rational_number_expression(
                        summation(
                            constant_coefficient(a),
                            constant_coefficient(b)
                        )
                    ),
                    nonconstant_coefficient(a)
                )
            )
        };
    }

    if(order_relation(a, b))
        return { a, b };

    return { b, a };
}

// simplify((a+b)+c) = simplify(a+b+c)
std::vector<expression*> transform_summation_associative(std::vector<expression*> L) {
    std::vector<expression*> _L;
    // printf("L[0]: ");
    // print(L[0]);
    // printf("\nL[1]: ");
    // if(L.size() > 1)
    // print(L[1]);
    // printf("\n");

    for(int j=0; j<number_of_operands(L[0]); j++) 
        _L.push_back(operand(L[0], j));

    return simplify_summation_rec(
        merge_summations(_L, simplify_summation_rec(rest(L)))
    );        
}

std::vector<expression*> transform_summation_expand(std::vector<expression*> L) {
    // printf("transform_expand\n");
    if(L.size() == 2 && kind(L[1]) != expression::ALG_OP_SUMMATION)
        return transform_summation(L[0], L[1]);

    return merge_summations({L[0]}, simplify_summation_rec(rest(L))); 
}








std::vector<expression*> simplify_summation_rec(std::vector<expression*> L) {
    if(L.size() == 0)
        return L;

    if(L.size() == 1 && kind(L[0]) != expression::ALG_OP_SUMMATION)
        return L;

    // 0+a = a
    if(equals(L[0], integer(0)))
        return simplify_summation_rec(rest(L));

    // (a+b)+(c+d) = a+b+c+d
    if(kind(L[0]) == expression::ALG_OP_SUMMATION)
        return transform_summation_associative(L);

    return transform_summation_expand(L);       
}



expression* simplify_summation(const expression* u) {
    // assert(kind(u) == expression::ALG_OP_SUMMATION);

    if(kind(u) == expression::UNDEFINED)
        return undefined();
    
    for(int i=0; i<number_of_operands(u); i++) {
        expression* o = operand(u,i);
        if(kind(o) == expression::INTEGER && integer_value(o) == 0)
            return integer(0);
    }

    if(number_of_operands(u) == 1)
        return operand(u, 0);

    std::vector<expression*> L;
    for(int i=0; i<number_of_operands(u); i++)
        L.push_back(operand(u,i));

    std::vector<expression*> R = simplify_summation_rec(L);
    
    if(R.size() == 1)
        return R[0];

    if(R.size() == 0)
        return integer(0);
    
    return summation(R);
}

}