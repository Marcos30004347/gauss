
#include <vector>
#include <assert.h>
#include "products.hpp"
#include "powers.hpp"
#include "summations.hpp"
#include "rationals.hpp"
#include "ordering/ordering.hpp"
#include "core/power.hpp"

using namespace core;
using namespace ordering;

namespace simplification {

std::vector<expr*> simplify_product_rec(std::vector<expr*> L);

std::vector<expr*> adjoin_products(expr* p, std::vector<expr*> q) {
    bool included = false;
    std::vector<expr*> res = std::vector<expr*>(0);

    for(int i=0; i < q.size(); i++) {
        expr* o = q[i];

        // if(equals(o, p)) {
        //     res.push_back(simplify_power(power(copy(o), integer(2))));
        //     included = true;
        // } else
            res.push_back(copy(o));
    }

    // if(!included)
        res.push_back(p);

    return res;
}

std::vector<expr*> merge_products(std::vector<expr*> p, std::vector<expr*> q) {

    if(p.size() == 0)
        return q;
    if(q.size() == 0)
        return p;

    std::vector<expr*> L = std::vector<expr*>(0);

    L.push_back(p[0]);
    L.push_back(q[0]);

    std::vector<expr*> H = simplify_product_rec(L);

    if(H.size() == 0) {
        return merge_products(rest(p), rest(q));
    }
    
    if(H.size() == 1)
        return adjoin_products(H[0], merge_products(rest(p), rest(q)));

    if(H.size() == 2) {
        if(equals(H[0], p[0]) && equals(H[1], q[0]))
            return adjoin_products(p[0], merge_products(rest(p), q));

        if(equals(H[0], p[0]) && equals(H[1], q[0]))
            return adjoin_products(q[0], merge_products(p, rest(q)));
    }

    std::vector<expr*> r = std::vector<expr*>(0);
    
    for(int i=0; i < p.size(); i++)
        adjoin_products(p[i], r);

    for(int i=0; i < q.size(); i++)
        adjoin_products(q[i], r);
    
    return r;
}

std::vector<expr*> simplify_product_rec(std::vector<expr*> L) {
    if(L.size() == 2) {

        expr* u1 = L[0];
        expr* u2 = L[1];
        if(kind(u1) != expr::ALG_OP_PRODUCT && kind(u2) != expr::ALG_OP_PRODUCT) {
    
            if(is_constant(u1) && is_constant(u2)) {
                expr* P = simplify_rational_number_expression(product(u1, u2));
                std::vector<expr*> res = std::vector<expr*>(0);

                if(kind(P) != expr::INTEGER  || (kind(P) == expr::INTEGER && integer_value(P) != 1))
                    res.push_back(P);
        
                return res;
            }

            if(kind(u1) == expr::INTEGER && integer_value(u1) == 1) {
                std::vector<expr*> res = std::vector<expr*>(0);
                res.push_back(u2);
                return res;
            }

            if(kind(u2) == expr::INTEGER && integer_value(u2) == 1) {
                std::vector<expr*> res = std::vector<expr*>(0);
                res.push_back(u1);
                return res;
            }

            if(equals(base(u1), base(u2))) {
                expr* S = simplify_summation(summation(expoent(u1), expoent(u2)));

                expr* P = simplify_power(power(base(u1), S));
                std::vector<expr*> res = std::vector<expr*>(0);
                
                if(kind(P) == expr::INTEGER && integer_value(P) == 1) {
                    return res;
                }

                res.push_back(P);
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

        if(kind(u1) == expr::ALG_OP_PRODUCT || kind(u2) == expr::ALG_OP_PRODUCT) {
            if(kind(u1) == expr::ALG_OP_PRODUCT && kind(u2) == expr::ALG_OP_PRODUCT) {
                std::vector<expr*> L0 = std::vector<expr*>(0);
                std::vector<expr*> L1 = std::vector<expr*>(0);
                for(int i=0; i<number_of_operands(u1); i++)
                    L0.push_back(operand(u1, i));
                for(int i=0; i<number_of_operands(u2); i++)
                    L1.push_back(operand(u2, i));
                return merge_products(L0, L1);
            }
            if(kind(u1) != expr::ALG_OP_PRODUCT && kind(u2) == expr::ALG_OP_PRODUCT) {
                std::vector<expr*> L0 = std::vector<expr*>(0);
                std::vector<expr*> L1 = std::vector<expr*>(0);
                L0.push_back(u1);
                for(int i=0; i<number_of_operands(u2); i++)
                    L1.push_back(operand(u2, i));
                return merge_products(L0, L1);
            }
            if(kind(u1) == expr::ALG_OP_PRODUCT && kind(u2) != expr::ALG_OP_PRODUCT) {
                std::vector<expr*> L0 = std::vector<expr*>(0);
                std::vector<expr*> L1 = std::vector<expr*>(0);
                for(int i=0; i<number_of_operands(u1); i++)
                    L0.push_back(operand(u1, i));
                L1.push_back(u2);
                return merge_products(L0, L1);
            }
        }
    } else if(L.size() > 2) {
        std::vector<expr*> w = simplify_product_rec(rest(L));
        expr* u1 = L[0];
        if(kind(u1) == expr::ALG_OP_PRODUCT) {
            std::vector<expr*> L0 = std::vector<expr*>(0);
            for(int i=0; i<number_of_operands(u1); i++)
                L0.push_back(operand(u1, i));
            return merge_products(L0, w);
        }

        std::vector<expr*> L0 = std::vector<expr*>(0);
        L0.push_back(u1);
        return merge_products(L0, w);
    }

    return L;
}

expr* simplify_product(const expr* u) {
    assert(kind(u) == expr::ALG_OP_PRODUCT);

    if(kind(u) == expr::UNDEFINED)
        return undefined();
    
    for(int i=0; i<number_of_operands(u); i++) {
        expr* o = operand(u,i);
        if(kind(o) == expr::INTEGER && integer_value(o) == 0)
            return integer(0);
    }

    if(number_of_operands(u) == 1)
        return operand(u, 0);

    std::vector<expr*> L;
    for(int i=0; i<number_of_operands(u); i++)
        L.push_back(operand(u,i));

    std::vector<expr*> R = simplify_product_rec(L);
    
    if(R.size() == 1)
        return R[0];
    if(R.size() == 0)
        return integer(1);
    
    return product(R);
}

}