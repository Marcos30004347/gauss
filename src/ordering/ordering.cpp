#include "ordering.hpp"
#include "algebra/power.hpp"
#include "algebra/rational.hpp"
#include "evaluation/evaluate.hpp"
#include <string.h>
#include <cstdio>

using namespace algebra;
using namespace evaluation;

namespace ordering {

bool compare_symbols(const char* a, const char* b) {
    return std::lexicographical_compare(a, a+strlen(a), b, b+strlen(b));
}

bool compare_constants(const expression* u, const expression* v) {
    if(kind(u) == expression::INTEGER && kind(v) == expression::INTEGER)
        return integer_value(u) < integer_value(v);

    expression* d = gcd(u, v);
    if(
        kind(d) == expression::INTEGER &&
        kind(numerator(u)) == expression::INTEGER &&
        kind(numerator(v)) == expression::INTEGER
    ) {
        return  integer_value(evaluate_product(d, numerator(u))) < integer_value(evaluate_product(d, numerator(v)));
    }

    return false;
}

bool compare_products_and_summations(const expression* u, const expression* v) {
    long m = number_of_operands(u) - 1;
    long n = number_of_operands(v) - 1;

    for(int k=0; k <= std::min(m, n); k++) {
        if(!equals(operand(u, m - k), operand(v, n - k))) {
            return order_relation(operand(u,m - k), operand(v,n - k));
        }
    }
    
    return m < n;
}

bool compare_powers(const expression* u, const expression* v) {
    if(!equals(base(u), base(v)))
        return order_relation(base(u), base(v));

    return order_relation(expoent(u), expoent(v));
}

bool compare_factorials(const expression* u, const expression* v) {
    return order_relation(operand(u,0), operand(v,0));
}

bool compare_functions(const expression* u, const expression* v) {
    if(strcmp(function_name(u), function_name(v)) != 0)
        return order_relation(operand(u,0), operand(v, 0));

    // operand(function, 1) is the arguments list
    expression* argsu = operand(u, 1);
    expression* argsv = operand(v, 1);

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

bool order_relation(const expression* u, const expression* v) {
    if(is_constant(u) && is_constant(v))
        return compare_constants(u, v);

    if(kind(u) == expression::SYMBOL && kind(v) == expression::SYMBOL)
        return compare_symbols(symbol_value(u), symbol_value(v));

    if(kind(u) == expression::ALG_OP_SUMMATION && kind(v) == expression::ALG_OP_SUMMATION)
        return compare_products_and_summations(u, v);

    if(kind(u) == expression::ALG_OP_PRODUCT && kind(v) == expression::ALG_OP_PRODUCT)
        return compare_products_and_summations(u, v);
    
    if(kind(u) == expression::ALG_OP_POWER && kind(v) == expression::ALG_OP_POWER)
        return compare_powers(u, v);

    if(kind(u) == expression::ALG_OP_FACTORIAL && kind(v) == expression::ALG_OP_FACTORIAL)
        return compare_factorials(u, v);

    if(kind(u) == expression::FUNCTION && kind(v) == expression::FUNCTION)
        return compare_functions(u, v);

    if(is_constant(u))
        return true;

    // if(is_constant(v))
    //     return false;

    if(kind(u) == expression::ALG_OP_PRODUCT && (
        kind(v) == expression::ALG_OP_POWER ||
        kind(v) == expression::ALG_OP_SUMMATION ||
        kind(v) == expression::ALG_OP_FACTORIAL ||
        kind(v) == expression::FUNCTION ||
        kind(v) == expression::SYMBOL
    )) return order_relation(u, product(v));

    if(kind(u) == expression::ALG_OP_POWER && (
        kind(v) == expression::ALG_OP_SUMMATION ||
        kind(v) == expression::ALG_OP_FACTORIAL ||
        kind(v) == expression::FUNCTION ||
        kind(v) == expression::SYMBOL
    )) return order_relation(u, power(v, integer(1)));

    if(kind(u) == expression::ALG_OP_SUMMATION && (
        kind(v) == expression::ALG_OP_FACTORIAL ||
        kind(v) == expression::FUNCTION ||
        kind(v) == expression::SYMBOL
    )) return order_relation(u, summation(v));

    if(kind(u) == expression::ALG_OP_FACTORIAL && (
        kind(v) == expression::FUNCTION ||
        kind(v) == expression::SYMBOL
    )) {
        if(equals(operand(u,0), v)) {
            return false;
        } else {
            return order_relation(u, factorial(v));
        }
    }

    if(kind(u) == expression::FUNCTION && kind(v) == expression::SYMBOL) {
        if(strcmp(function_name(u), symbol_value(v)) == 0) {
            return false;
        } else {
            return order_relation(operand(u,0), v);
        }
    }
    // printf("NOT COMPARABLE %i %i\n", kind(u), kind(v));
    return !order_relation(v, u);
}

}