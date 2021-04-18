#include "summations.hpp"
#include "products.hpp"
#include "rationals.hpp"

#include "ordering/ordering.hpp"
#include "core/power.hpp"

using namespace core;
using namespace ordering;

namespace simplification {

std::vector<expr*> simplify_summation_rec(std::vector<expr*> L);

std::vector<expr*> adjoin_summations(expr* p, std::vector<expr*> q) {
    
    std::vector<expr*> res = std::vector<expr*>(0);
    std::vector<expr*> tmp = std::vector<expr*>(0);

    for(int i=0; i<number_of_operands(p); i++) {
        tmp.push_back(operand(p, i));
    }

    for(int i=0; i<q.size(); i++) {
        tmp.push_back(q[i]);
    }

    std::vector<bool> included = std::vector<bool>(tmp.size(), false);

    for(int i=0; i< tmp.size() && !included[i]; i++) {
        int c = 1;

        for(int j=i + 1; j<tmp.size(); j++) {
            if(equals(tmp[i], tmp[j])) {
                included[j] = true;
                c++;
            }
        }
        res.push_back(simplify_product(product(integer(c), tmp[i])));
    }

    // printf("\n");
    // print(p);
    // printf("\n");
    // for(int i=0; i < q.size(); i++) {
    //     print(q[i]);
    //     printf(" ");
    // }
    // printf("\n");

    // for(int i=0; i < q.size(); i++) {
    //     included = false;

    //     for(int j=0; j<number_of_operands(p); j++) {
    //         if(equals(operand(p, j), q[i])) {
    //             included = true;
    //             counts[j]++;
    //         }
    //     }

    //     if(!included)
    //         res.push_back(copy(q[i]));
    // }
    // for(int i=0; i < res.size(); i++) {
    //     print(res[i]);
    //     printf(" ");
    // }

    // printf("\n");
    // for(int i=0; i<number_of_operands(p); i++) {
    //     res.push_back(simplify_product(product(integer(counts[i]), operand(p,i))));
    // }

    // for(int i=0; i < res.size(); i++) {
    //     print(res[i]);
    //     printf(" ");
    // }

    // printf("\n");

    return res;
}

std::vector<expr*> merge_summations(std::vector<expr*> p, std::vector<expr*> q) {
    if(p.size() == 0)
        return q;
    if(q.size() == 0)
        return p;

    std::vector<expr*> L = std::vector<expr*>(0);

    L.push_back(p[0]);
    L.push_back(q[0]);

    std::vector<expr*> H = simplify_summation_rec(L);

    if(H.size() == 0) {
        return merge_summations(rest(p), rest(q));
    } else if(H.size() == 1) {
        return adjoin_summations(H[0], merge_summations(rest(p), rest(q)));
    }

    return merge_summations(H, merge_summations(rest(p), rest(q)));
}


std::vector<expr*> simplify_summation_rec(std::vector<expr*> L) {
    if(L.size() == 2) {

        expr* u1 = L[0];
        expr* u2 = L[1];
    
        if(kind(u1) != expr::ALG_OP_SUMMATION && kind(u2) != expr::ALG_OP_SUMMATION) {
            if(is_constant(u1) && is_constant(u2)) {
                expr* P = simplify_rational_number_expression(summation(u1, u2));
                std::vector<expr*> res = std::vector<expr*>(0);
        
                if(kind(P) != expr::INTEGER  || integer_value(P) != 1)
                    res.push_back(P);
        
                return res;
            }


            if(kind(u1) == expr::INTEGER && integer_value(u1) == 0) {
                std::vector<expr*> res = std::vector<expr*>(0);
                res.push_back(u2);
                return res;
            }

            if(kind(u2) == expr::INTEGER && integer_value(u2) == 0) {
                std::vector<expr*> res = std::vector<expr*>(0);
                res.push_back(u1);
                return res;
            }
            
            if(equals(u1, u2)) {
                expr* prod = simplify_product(product(integer(2), u1));
                print(u1);
                printf("  - \n");
                print(u2);
                printf("  - \n");

                std::vector<expr*> res = std::vector<expr*>(0);

                if(kind(prod) != expr::INTEGER || integer_value(u1) != 0)
                    res.push_back(prod);

                return res;
            }

            if(order_relation(u2, u1)) {
                std::vector<expr*> res = std::vector<expr*>(0);
                res.push_back(u2);
                res.push_back(u1);
                return res;
            }

            return L;
        }

        if(kind(u1) == expr::ALG_OP_SUMMATION || kind(u2) == expr::ALG_OP_SUMMATION) {
            printf("ASASA\n");
            print(u1);
            printf(" + ");
            print(u2);
            printf("\n");

            if(kind(u1) == expr::ALG_OP_SUMMATION && kind(u2) == expr::ALG_OP_SUMMATION) {
                std::vector<expr*> L0 = std::vector<expr*>(0);
                std::vector<expr*> L1 = std::vector<expr*>(0);

                expr* U1 = simplify_summation(u1);
                expr* U2 = simplify_summation(u2);
        
                // print(U1);
                // printf(" + ");
                // print(U2);
                // printf("\n");
        
                // printf("L0 = ");
                for(int i=0; i<number_of_operands(U1); i++) {
                    // print(operand(U1, i));
                    // printf(" ");
                    L0.push_back(operand(U1, i));
                }
    
                // printf("L1 = ");
                for(int i=0; i<number_of_operands(U2); i++) {
                    // print(operand(U2, i));
                    // printf(" ");
                    L1.push_back(operand(U2, i));
                }
                // printf("\n");

                std::vector<expr*> L = merge_summations(L0, L1);
                
                return L;
            }

            // printf("asdsad\n");

            if(kind(u1) != expr::ALG_OP_SUMMATION && kind(u2) == expr::ALG_OP_SUMMATION) {
                std::vector<expr*> L0 = std::vector<expr*>(0);
                std::vector<expr*> L1 = std::vector<expr*>(0);
                L0.push_back(u1);

                for(int i=0; i<number_of_operands(u2); i++)
                    L1.push_back(operand(u2, i));

                return merge_summations(L0, L1);
            }

            if(kind(u1) == expr::ALG_OP_SUMMATION && kind(u2) != expr::ALG_OP_SUMMATION) {
                std::vector<expr*> L0 = std::vector<expr*>(0);
                std::vector<expr*> L1 = std::vector<expr*>(0);
                for(int i=0; i<number_of_operands(u1); i++)
                    L0.push_back(operand(u1, i));
                L1.push_back(u2);
                return merge_summations(L0, L1);
            }
        }
    } else if(L.size() > 2) {
        std::vector<expr*> w = simplify_summation_rec(rest(L));
        expr* u1 = L[0];
        if(kind(u1) == expr::ALG_OP_SUMMATION) {
            std::vector<expr*> L0 = std::vector<expr*>(0);
            for(int i=0; i<number_of_operands(u1); i++)
                L0.push_back(operand(u1, i));
            return merge_summations(L0, w);
        }

        std::vector<expr*> L0 = std::vector<expr*>(0);
        L0.push_back(u1);
        return merge_summations(L0, w);
    }

    return L;
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

    // printf("ASDAS\n");
    // printf("ASDAS\n");
    // printf("ASDAS\n");
    // printf("ASDAS\n");

    std::vector<expr*> R = simplify_summation_rec(L);

    // for(int i=0; i<R.size(); i++) {
    //     print(R[i]);
    //     printf(" ");
    // }
    // printf("\n");

    if(R.size() == 1)
        return R[0];

    if(R.size() == 0)
        return integer(1);
    
    return summation(R);
}

}