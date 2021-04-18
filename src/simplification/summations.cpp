#include "summations.hpp"
#include "products.hpp"
#include "rationals.hpp"
#include <assert.h>
#include <cstdio>
#include "ordering/ordering.hpp"
#include "core/power.hpp"

using namespace core;
using namespace ordering;

namespace simplification {

std::vector<expr*> simplify_summation_rec(std::vector<expr*> L);

std::vector<expr*> adjoin_summations(expr* p, std::vector<expr*> q) {
    std::vector<expr*> tmp = std::vector<expr*>(0);

    int i = 0;
    bool inc = false;

    while(i != q.size()) {
        if(!inc && order_relation(p, q[i])) {
            tmp.push_back(p);
            inc = true;
        } else {
            tmp.push_back(q[i++]);
        }
    }

    if(!q.size() || !inc) {
        tmp.push_back(p);
    }

    for(int i=0; i<tmp.size(); i++) {
        for(int j= i + 1; j<tmp.size(); j++) {
    
            if(equals(tmp[i], tmp[j])) {
                expr* a = tmp[i];
                tmp.erase(tmp.begin() + i);
                tmp.erase(tmp.begin() + j - 1);

                tmp.push_back(simplify_product(product(integer(2), copy(a))));
                continue;
            }

            if(kind(tmp[i]) == expr::ALG_OP_PRODUCT && kind(tmp[j]) == expr::ALG_OP_PRODUCT) {
                if(equals(operand(tmp[i], 1), operand(tmp[j], 1))) {
                    expr* a = tmp[i];
                    expr* b = tmp[j];

                    tmp.erase(tmp.begin() + i);
                    tmp.erase(tmp.begin() + j - 1);

                    tmp.push_back(simplify_product(
                        product(
                            simplify_summation(
                                summation(
                                    operand(a, 0),
                                    operand(b, 0)
                                )
                            ),
                            copy(operand(b,1))
                        )
                    ));
                    continue;
                }
            }

            if(kind(tmp[i]) == expr::ALG_OP_PRODUCT && kind(tmp[j]) != expr::ALG_OP_PRODUCT) {
                if(equals(operand(tmp[i], 1), tmp[j])) {
                    expr* a = tmp[i];
                    expr* b = tmp[j];

                    tmp.erase(tmp.begin() + i);
                    tmp.erase(tmp.begin() + j - 1);

                    tmp.push_back(simplify_product(product(simplify_summation(summation(operand(a, 0), integer(1))), copy(b))));
                    continue;
                }
            }

            if(kind(tmp[i]) != expr::ALG_OP_PRODUCT && kind(tmp[j]) == expr::ALG_OP_PRODUCT) {
                if(equals(operand(tmp[j], 1), tmp[i])) {
                    expr* a = tmp[j];
                    expr* b = tmp[i];

                    tmp.erase(tmp.begin() + i);
                    tmp.erase(tmp.begin() + j - 1);

                    tmp.push_back(simplify_product(product(simplify_summation(summation(operand(a, 0), integer(1))), copy(b))));
                    continue;

                }
            }
        }
    }

    return tmp;
}

std::vector<expr*> merge_summations(std::vector<expr*> p, std::vector<expr*> q) {
    if(p.size() == 0)
        return q;
    if(q.size() == 0)
        return p;

    std::vector<expr*> L = std::vector<expr*>(0);

    L.push_back(p[0]);
    L.push_back(q[0]);

    // printf("L's: \n");

    // print(L[0]);
    // printf("\n");

    // print(L[1]);
    // printf("\n");
    // printf("\n");

    std::vector<expr*> H = simplify_summation_rec(L);
    std::vector<expr*> R = merge_summations(rest(p), rest(q));

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
    
    return simplify_summation_rec(R);
}


std::vector<expr*> simplify_summation_rec(std::vector<expr*> L) {
    // printf("simplify_summation_rec\n");
    // for(int i=0; i<L.size(); i++) {
    //     print(L[i]);
    //     printf(" - ");
    // }
    // printf("\n");

    if(L.size() == 0)
        return L;

    if(L.size() == 1 && kind(L[0]) != expr::ALG_OP_SUMMATION)
        return L;


    std::vector<expr*> _L;

    // compare equality
    // printf("asdasdasdasd\n");
    // if(equals(L[0], L[1])) {
    //     printf("asdasdasdasd\n");
    //     expr* prod = simplify_product(product(integer(2), L[0]));
    //     _L.push_back(prod);

    //     return merge_summations(_L, rest(L, 2));
    // }
    // a + 0 = a

    //

    if(equals(L[0], integer(0))) {
        return simplify_summation_rec(rest(L));
    }

    // if(is_constant(L[0]) && is_constant(L[1])) {
    //     expr* e = simplify_rational_number_expression(summation(L[0], L[1]));
    //     _L.push_back(e);
    //     return merge_summations(_L, rest(L, 2));
    // }
    
    if(is_constant(L[0]) && is_constant(L[1])) {
        expr* e = simplify_rational_number_expression(summation(L[0], L[1]));
        _L.push_back(e);
        return merge_summations(_L, rest(L, 2));
    }

    // if(L.size() == 2 && kind(L[0]) != expr::ALG_OP_SUMMATION && kind(L[1]) != expr::ALG_OP_SUMMATION) {

    //     // if(order_relation(L[0], L[1])) {
    //     //     _L.push_back(L[0]);
    //     //     _L.push_back(L[1]);
    //     // } else {
    //     //     _L.push_back(L[1]);
    //     //     _L.push_back(L[0]);
    //     // }

    //     return _L;
    // }

    if(kind(L[0]) == expr::ALG_OP_SUMMATION) {
    
        for(int j=0; j<number_of_operands(L[0]); j++) 
            _L.push_back(operand(L[0], j));

        std::vector<expr*> rst = rest(L);

        return merge_summations(_L, simplify_summation_rec(rest(L)));        
    }

    // _L.push_back(L[0]);

    return adjoin_summations(L[0], simplify_summation_rec(rest(L)));        

    // return L;
}


expr* simplify_summation(const expr* u) {
    assert(kind(u) == expr::ALG_OP_SUMMATION);

    if(kind(u) == expr::UNDEFINED)
        return undefined();

    if(number_of_operands(u) == 1)
        return operand(u, 0);

    std::vector<expr*> L;
    for(int i=0; i<number_of_operands(u); i++)
        L.push_back(operand(u,i));

    std::vector<expr*> R = simplify_summation_rec(L);

    if(R.size() == 1)
        return R[0];

    if(R.size() == 0)
        return integer(0);
    
    return summation(R);
}

}