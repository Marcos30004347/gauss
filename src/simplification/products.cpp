
#include <vector>
#include <assert.h>
#include "products.hpp"
#include "powers.hpp"
#include "summations.hpp"
#include "rationals.hpp"
#include "ordering/ordering.hpp"
#include "core/power.hpp"
#include <cstdio>

using namespace core;
using namespace ordering;

namespace simplification {

std::vector<expr*> simplify_product_rec(std::vector<expr*> L);
std::vector<expr*> adjoin_products(expr* p, std::vector<expr*> q);
std::vector<expr*> merge_products(std::vector<expr*> p, std::vector<expr*> q);

std::vector<expr*> transform_product_distributive(expr* a, expr* b);
std::vector<expr*> transform_product(expr* a, expr* b);
std::vector<expr*> transform_product_associative(std::vector<expr*> L);
std::vector<expr*> transform_expand(std::vector<expr*> L);

std::vector<expr*> adjoin_products(expr* p, std::vector<expr*> q) {
    std::vector<expr*> tmp = std::vector<expr*>(0);

    int i = 0;
    bool inc = false;

    while(i != q.size()) {
        if(!inc && order_relation(p, q[i])) {
            tmp.push_back(p);
            inc = true;
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

std::vector<expr*> merge_products(std::vector<expr*> p, std::vector<expr*> q) {
    if(p.size() == 0)
        return q;

    if(q.size() == 0)
        return p;

    std::vector<expr*> H = simplify_product_rec({p[0], q[0]});
    std::vector<expr*> R = merge_products(rest(p), rest(q));

    for(int i=0; i<H.size(); i++) {
        R = adjoin_products(H[i], R);
    }

    return (R);
}

std::vector<expr*> transform_product_distributive(expr* a, expr* b) {
    if(kind(a) == expr::ALG_OP_SUMMATION && kind(b) == expr::ALG_OP_SUMMATION) {
        expr* u = construct(expr::ALG_OP_SUMMATION);
        for(int i=0; i<number_of_operands(a); i++) 
            for(int j=0; j<number_of_operands(b); j++) 
                include_operand(u, simplify_product(
                    product(
                        merge_products(
                            {operand(a, i)},
                            {operand(b, j)}
                        )
                    )
                ));

        return { simplify_summation(u) };
    }

    if(kind(a) != expr::ALG_OP_SUMMATION && kind(b) == expr::ALG_OP_SUMMATION) {
        expr* u = construct(expr::ALG_OP_SUMMATION);
        
        for(int j=0; j<number_of_operands(b); j++) 
            include_operand(u, simplify_product(
                product(
                    merge_products(
                        {operand(b, j)},
                        {a}
                    )
                )
            ));

        return { simplify_summation(u) };
    }

    if(kind(a) == expr::ALG_OP_SUMMATION && kind(b) != expr::ALG_OP_SUMMATION) {
        expr* u = construct(expr::ALG_OP_SUMMATION);
        
        for(int j=0; j<number_of_operands(a); j++)
            include_operand(u, simplify_product(
                product(
                    merge_products(
                        {operand(a, j)},
                        {b}
                    )
                )
            ));

        return { simplify_summation(u) };
    }

    return { a, b };
}

// simplify(a*a) = a^2, simplify(a*b) = a*b
std::vector<expr*> transform_product(expr* a, expr* b) {
    if(equals(base(a), base(b)))
        return {simplify_power(
            power(
                base(a),
                simplify_summation(
                    summation(
                        expoent(a),
                        expoent(b)
                    )
                )
            )
        )};

    // printf("a: ");
    // print(a);
    // printf("\nb: ");
    // print(b);
    // printf("\n");

    // printf("%i %i\n", kind(a), kind(b));

    if(is_constant(a) && is_constant(b))
        return { simplify_rational_number_expression(product(a, b)) };

    if(order_relation(a, b))
        return { a, b };

    return { b, a };
}

// simplify((a*b)*c) = simplify(a*b*c)
std::vector<expr*> transform_product_associative(std::vector<expr*> L) {
    std::vector<expr*> _L;

    for(int j=0; j<number_of_operands(L[0]); j++) 
        _L.push_back(operand(L[0], j));

    std::vector<expr*> rst = rest(L);

    return simplify_product_rec(
        merge_products(_L, simplify_product_rec(rest(L)))
    );        
}

std::vector<expr*> transform_expand(std::vector<expr*> L) {
    // printf("transform_expand\n");
    if(L.size() == 2) {
        if(
            kind(L[0]) == expr::ALG_OP_SUMMATION ||
            kind(L[1]) == expr::ALG_OP_SUMMATION
        ) return transform_product_distributive(L[0], L[1]);
    
        if(
            kind(L[0]) != expr::ALG_OP_PRODUCT &&
            kind(L[1]) != expr::ALG_OP_PRODUCT
        ) return transform_product(L[0], L[1]);
    }

    return merge_products({L[0]}, simplify_product_rec(rest(L))); 
}

std::vector<expr*> simplify_product_rec(std::vector<expr*> L) {
    // printf("simplify_product_rec\n");
    // for(int i=0; i<L.size(); i++) {
    //     print(L[i]);
    //     printf(" - ");
    // }
    // printf("\n");
    if(L.size() == 0)
        return L;

    if(L.size() == 1 && kind(L[0]) != expr::ALG_OP_PRODUCT)
        return L;

    // 0*a = 0
    if(equals(L[0], integer(0)))
        return {};

    // 1*a = a
    if(equals(L[0], integer(1)))
        return simplify_product_rec(rest(L));

    // if(is_constant(L[0]) && is_constant(L[1])) {
        // return simplify_rational_number_expression(product(L[0], L[1])), rest(L,2));
    // }

    // (a*b)*(c*d) = a*b*c*d
    if(kind(L[0]) == expr::ALG_OP_PRODUCT)
        return transform_product_associative(L);

    return transform_expand(L);       
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
        return integer(0);
    
    return product(R);
}

}