#include <assert.h>
#include <string.h>
#include "core/expression.hpp"
#include "core/rational.hpp"

using namespace core;

void should_create_integer_expression() {
    expr* e = integer(3);
    
    assert(kind(e) == expr::INTEGER);
    assert(integer_value(e) == 3);
    
    destroy(e);
}

void should_create_fraction_expression() {
    expr* e = fraction(integer(3), integer(2));

    assert(kind(e) == expr::FRACTION);
    assert(equals(numerator(e), integer(3)) == true);
    assert(equals(denominator(e), integer(2)) == true);

    destroy(e);
}

void should_create_quotient_expression() {
    expr* e = quotient(integer(3), integer(2));

    assert(kind(e) == expr::ALG_OP_QUOTIENT);
    assert(equals(numerator(e), integer(3)) == true);
    assert(equals(denominator(e), integer(2)) == true);

    destroy(e);
}

void should_create_product_expression() {
    expr* e0 = product(integer(3), integer(4));

    assert(kind(e0) == expr::ALG_OP_PRODUCT);
    assert(equals(operand(e0,0), integer(3)) == true);
    assert(equals(operand(e0,1), integer(4)) == true);

    destroy(e0);

    expr* e1 = product(integer(3));

    assert(kind(e1) == expr::ALG_OP_PRODUCT);
    assert(equals(operand(e1,0), integer(3)) == true);

    destroy(e1);
    
    expr* e2 = product(integer(3));

    assert(kind(e2) == expr::ALG_OP_PRODUCT);
    assert(equals(operand(e2,0), integer(3)) == true);

    destroy(e2);
}

void should_create_power_expression() {
    expr* e1 = power(integer(3), integer(4));

    assert(kind(e1) == expr::ALG_OP_POWER);
    assert(equals(operand(e1,0), integer(3)) == true);
    assert(equals(operand(e1,1), integer(4)) == true);

    destroy(e1);
}

void should_create_difference_expression() {
    expr* e1 = difference(integer(3), integer(4));

    assert(kind(e1) == expr::ALG_OP_DIFFERENCE);
    assert(equals(operand(e1,0), integer(3)) == true);
    assert(equals(operand(e1,1), integer(4)) == true);

    destroy(e1);
}

void should_create_summation_expression() {
    expr* e1 = summation(integer(3), integer(4));

    assert(kind(e1) == expr::ALG_OP_SUMMATION);
    assert(equals(operand(e1,0), integer(3)) == true);
    assert(equals(operand(e1,1), integer(4)) == true);

    destroy(e1);
}

void should_create_symbol_expression() {
    expr* e = symbol("function");

    assert(kind(e) == expr::SYMBOL);
    assert(strcmp("function", (const char*)e->_data) == 0);

    destroy(e);
}

void should_be_free_of() {
    expr* e = summation(quotient(integer(3), integer(4)), product(integer(5), integer(6)));

    assert(free_of(e, quotient(integer(3), integer(4))) == false); // LEAK
    assert(free_of(e, product(integer(5), integer(6))) == false); // LEAK

    // integer(5) is not a complete subtree of 'e'
    assert(free_of(e, integer(3)) == true);
    assert(free_of(operand(e, 0), integer(3)) == false);
    assert(free_of(operand(e, 1), integer(5)) == false);
}

void should_substitute() {
    expr* e0 = summation(quotient(integer(4), integer(2)), product(integer(5), integer(6)));

    expr* e1 = substitute(e0, quotient(integer(4), integer(2)), integer(2)); // LEAK
    expr* e2 = substitute(e0, product(integer(5), integer(6)), integer(30)); // LEAK
    
    assert(equals(operand(e1, 0), integer(2)) == true); // LEAK
    assert(equals(operand(e2, 1), integer(30)) == true); // LEAK
}

expr* map_handler_0(expr* e) {
    return product(copy(e), integer(2));
}

expr* map_mul_by_2(const expr* e) {
    return product(copy(e), integer(2));
}

expr* map_mul_by_e(const expr* e0, const expr* e1) {
    return product(copy(e0), copy(e1));
}

void should_map() {
    expr* e0 = integer(2);
    expr* e1 = integer(3);

    expr* e2 = unary_map(e0, map_mul_by_2);
    expr* e3 = binary_map(e0, e1, map_mul_by_e);

    expr* e4 = product(integer(2), integer(2));
    expr* e5 = product(integer(2), integer(3));

    assert(equals(e2, e4) == true);
    assert(equals(e3, e5) == true);

    destroy(e0);
    destroy(e1);
    destroy(e2);
    destroy(e3);
    destroy(e4);
    destroy(e5);
}


void should_be_equal() {
    expr* e0 = integer(3);
    expr* e1 = integer(3);

    assert(equals(e0, e1) == true);
    assert(equals(e0, e0) == true);
    assert(equals(e1, e1) == true);

    destroy(e0);
    destroy(e1);

    expr* e2 = summation(summation(integer(3), integer(4)), product(integer(5), integer(6)));
    expr* e3 = summation(product(integer(5), integer(6)), summation(integer(3), integer(4)));

    assert(equals(e2, e3));
    
    destroy(e2);
    destroy(e3);

    expr* e4 = power(integer(3), integer(2));
    expr* e5 = power(integer(2), integer(3));

    assert(equals(e4, e5) == false);

    destroy(e4);
    destroy(e5);
}

int main() {
    should_create_integer_expression();
    should_create_product_expression();
    should_create_fraction_expression();
    should_create_quotient_expression();
    should_create_power_expression();
    should_create_difference_expression();
    should_create_summation_expression();
    should_create_symbol_expression();
    should_be_equal();
    should_be_free_of();
    should_substitute();
    should_map();
    return 0;
}