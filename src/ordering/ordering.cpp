#include "ordering.hpp"
#include "core/power.hpp"
#include "core/rational.hpp"
#include "evaluation/evaluate.hpp"
#include <string.h>

using namespace core;
using namespace evaluation;

namespace ordering {

bool compare_symbols(const char* a, const char* b) {
    return std::lexicographical_compare(a, a+strlen(a), b, b+strlen(b));
}

bool compare_constants(const expr* u, const expr* v) {
    if(kind(u) == expr::INTEGER && kind(v) == expr::INTEGER)
        return integer_value(u) < integer_value(v);

    expr* d = gcd(u, v);
    if(
        kind(d) == expr::INTEGER &&
        kind(numerator(u)) == expr::INTEGER &&
        kind(numerator(v)) == expr::INTEGER
    ) {
        return  integer_value(evaluate_product(d, numerator(u))) < integer_value(evaluate_product(d, numerator(v)));
    }

    return false;
}

bool compare_products_and_summations(const expr* u, const expr* v) {
    long m = number_of_operands(u) - 1;
    long n = number_of_operands(v) - 1;

    for(int k=0; k <= std::min(m, n); k++) {
        if(!equals(operand(u, m - k), operand(v, n - k))) {
            return order_relation(operand(u,m - k), operand(v,n - k));
        }
    }
    
    return m < n;
}

bool compare_powers(const expr* u, const expr* v) {
    if(!equals(base(u), base(v)))
        return order_relation(base(u), base(v));

    return order_relation(expoent(u), expoent(v));
}

bool compare_factorials(const expr* u, const expr* v) {
    return order_relation(operand(u,0), operand(v,0));
}

bool compare_functions(const expr* u, const expr* v) {
    if(strcmp(function_name(u), function_name(v)) != 0)
        return order_relation(operand(u,0), operand(v, 0));

    // operand(function, 1) is the arguments list
    expr* argsu = operand(u, 1);
    expr* argsv = operand(v, 1);

    if(number_of_operands(argsu) >= 1 && number_of_operands(argsv) >= 1) {
        long m = number_of_operands(argsu) - 1;
        long n = number_of_operands(argsv) - 1;

        for(int k=0; k <= std::min(m, n); k++)
            if(!equals(operand(argsu, m - k), operand(argsv, n - k))) 
                return order_relation(operand(argsu, m - k), operand(argsv, n - k));
        
        return m < n;
    }

    return false;
}

bool order_relation(const expr* u, const expr* v) {
    if(is_constant(u) && is_constant(v))
        return compare_constants(u, v);

    if(kind(u) == expr::SYMBOL && kind(v) == expr::SYMBOL)
        return compare_symbols(symbol_value(u), symbol_value(v));

    if(kind(u) == expr::ALG_OP_SUMMATION && kind(v) == expr::ALG_OP_SUMMATION)
        return compare_products_and_summations(u, v);

    if(kind(u) == expr::ALG_OP_PRODUCT && kind(v) == expr::ALG_OP_PRODUCT)
        return compare_products_and_summations(u, v);
    
    if(kind(u) == expr::ALG_OP_POWER && kind(v) == expr::ALG_OP_POWER)
        return compare_powers(u, v);

    if(kind(u) == expr::ALG_OP_FACTORIAL && kind(v) == expr::ALG_OP_FACTORIAL)
        return compare_factorials(u, v);

    if(kind(u) == expr::FUNCTION && kind(v) == expr::FUNCTION)
        return compare_functions(u, v);

    if(is_constant(u))
        return true;

    if(kind(u) == expr::ALG_OP_PRODUCT && (
        kind(v) == expr::ALG_OP_POWER ||
        kind(v) == expr::ALG_OP_SUMMATION ||
        kind(v) == expr::ALG_OP_FACTORIAL ||
        kind(v) == expr::FUNCTION ||
        kind(v) == expr::SYMBOL
    )) return order_relation(u, product(v));

    if(kind(u) == expr::ALG_OP_POWER && (
        kind(v) == expr::ALG_OP_SUMMATION ||
        kind(v) == expr::ALG_OP_FACTORIAL ||
        kind(v) == expr::FUNCTION ||
        kind(v) == expr::SYMBOL
    )) return order_relation(u, power(v, integer(1)));

    if(kind(u) == expr::ALG_OP_SUMMATION && (
        kind(v) == expr::ALG_OP_FACTORIAL ||
        kind(v) == expr::FUNCTION ||
        kind(v) == expr::SYMBOL
    )) return order_relation(u, summation(v));


    if(kind(u) == expr::ALG_OP_FACTORIAL && (
        kind(v) == expr::FUNCTION ||
        kind(v) == expr::SYMBOL
    )) {
        if(equals(operand(u,0), v)) {
            return false;
        } else {
            return order_relation(u, factorial(v));
        }
    }

    if(kind(u) == expr::FUNCTION && kind(v) == expr::SYMBOL) {
        if(strcmp(function_name(u), symbol_value(v)) == 0) {
            return false;
        } else {
            return order_relation(operand(u,0), v);
        }
    }

    return false;
}

}