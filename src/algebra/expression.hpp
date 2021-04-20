#ifndef ALGEBRA_EXPRESSION_H
#define ALGEBRA_EXPRESSION_H

#include <vector>

namespace algebra {

struct expression;

class argument {
public:
    argument* _next;
    argument* _prev;
    
    expression* _operation;
    expression* _operand;
};

void include_operand(expression* operation, expression* operand);
void remove_operand(expression* operation, expression* operand);

struct expression {
    enum kind {
        UNDEFINED = 0,
        INTEGER,
        FRACTION,
        SYMBOL,

        // Algebraic operators
        ALG_OP_SUMMATION,
        ALG_OP_DIFFERENCE,
        ALG_OP_PRODUCT,
        ALG_OP_QUOTIENT,
        ALG_OP_POWER,
        ALG_OP_FACTORIAL,

        // Relational operators
        REL_OP_EQUAL,
        REL_OP_NOT_EQUAL,
        REL_OP_LESS,
        REL_OP_LESS_EQUAL,
        REL_OP_GREATER,
        REL_OP_GREATER_EQUAL,

        // Function forms
        FUNCTION,
    };

    // The kind of the epression
    kind _kind;

    // The data associated with the expression.
    // When the kind is symbols or constants, _data
    // will hold the value for those symbols.
    void* _data;

    // std::vector<expression*> _operands;
    argument* _operands;
    unsigned _operands_count;
};

/**
 * Return and undefined expression
 */
expression* undefined();

/**
 * Returns the expression that corresponds to a given integer.
 */
expression* integer(long long i);

/**
 * Returns a fraction expression.
 */
expression* fraction(expression* numerator, expression* denominator);

/**
 * Returns a quotient expression.
 */
expression* quotient(expression* numerator, expression* denominator);

/**
 * Returns a product expression.
 */
expression* product(const std::vector<expression*> operands);
expression* product(const expression* u, const expression* v);
expression* product(const expression* u);

/**
 * Returns a power u^v expression.
 */
expression* power(const expression* u, const expression* v);

/**
 * Returns a difference expression.
 */
expression* difference(const std::vector<expression*> operands);
expression* difference(expression* u, expression* v);

/**
 * Returns a summation expression.
 */
expression* summation(const std::vector<expression*> operands);
expression* summation(const expression* u, const expression* v);
expression* summation(const expression* u);

/**
 * Returns the expression that corresponds to a factorial of u
 */
expression* factorial(const expression* u);

/**
 * Returns the expression that corresponds to a given symbol.
 */
expression* symbol(const char* i);

/**
 * Return a expression with main operator of a given type and operands.
 * 
 * EXAMPLE: construct('+', [a,b,c]) -> a + b + c;
 */
expression* construct(expression::kind kind, std::vector<expression*> operands);

/**
 * Return a expression with main operator of a given type and operands.
 */
expression* construct(expression::kind kind);

/**
 * Destroys a given expression dealocating all its resources and sub trees.
 */
void destroy(expression* u);

/**
 * Return the i'th operand of u.
 * 
 * EXAMPLE: operand(a+ b + c, 2) -> b;
 */
expression* operand(const expression* u, unsigned i);

/**
 * Return the kind of the main operator of u.
 * 
 * EXAMPLE: kind(a + b + 3c) -> +;
 */
expression::kind kind(const expression* u);

/**
 * Return the number of operands of the main operator of u.
 * 
 * EXAMPLE: number_of_operands(a + 2b + 4(a+b)) -> 3;
 */
unsigned number_of_operands(const expression* u);

/**
 * Return true when u is equal to v and false otherwise.
 * 
 * EXAMPLE: equals(a+b, a+b) -> true;
 */
bool equals(const expression* u, const expression* v);

/**
 * Return false when v is a complete sub-expression of u,
 * and true otherwise.
 * 
 * EXAMPLE: free_of((a+b)c, a+b) -> false;
 */
bool free_of(const expression* u, const expression* v);

/**
 * Substitutes the t expression in u by r.
 * 
 * EXAMPLE: substitute(a+2*b+c, 2*b, d) -> a+d+c;
 * 
 * NOTE: The substitute only substitures complete subtrees,
 * so the substitution substitute(a+2*b+c, 2*b+c, d) wont 
 * work because 2*b+c is not a complete subtree of u;
 */
expression* substitute(const expression* u, const expression* t, expression* r);

/**
 * Return expression of kind kind(u) with operands F(operand(u,i))
 * for each i > 0 and i < number_of_operands(u). The interface for
 * the map handler is void (*)(const expression*), it will transform the 
 * given expression to a new one without modifying the original, the
 * original is deleted at the end of each map.
 * 
 * EXAMPLE: map(f(x) = (x+2)^2, x + 2) -> (x+2)^2 + (2+2)^2;
 * EXAMPLE: map(f(x, y) = (x+y)^2, x + 2, 3) -> (x+3)^2 + (2+3)^2;
 */
expression* unary_map(const expression* u, expression* (*f)(const expression*));
expression* binary_map(const expression* u, const expression* v, expression* (*f)(const expression*, const expression*));

/**
 * Return a full copy of a given expression.
 */
expression* copy(const expression* u);

/**
 * Returns the value of a integer expression.
 */
long long integer_value(const expression* u);

/**
 * Returns the value of a symbol expression.
 */
const char* symbol_value(const expression* u);

/**
 * Returns the value of a symbol expression.
 */
const char* function_name(const expression* u);

/**
 * Returns the value of a symbol expression.
 */
const char* symbol_value(const expression* u);

bool is_constant(const expression* u);

void print(const expression* u);

std::vector<expression*> rest(std::vector<expression*> p, int from = 1);

// std::vector<expression*> equal_operands(expression* u, expression* v);
// std::vector<expression*> different_operands(expression* u, expression* v);


}

#endif