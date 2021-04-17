#ifndef CORE_EXPRESSION_H
#define CORE_EXPRESSION_H

#include <vector>

namespace core {

struct expr;

class argument {
public:
    argument* _next;
    argument* _prev;
    
    expr* _operation;
    expr* _operand;
};

void include_operand(expr* operation, expr* operand);
void remove_operand(expr* operation, expr* operand);

struct expr {
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

    // std::vector<expr*> _operands;
    argument* _operands;
    unsigned _operands_count;
};

/**
 * Return and undefined expression
 */
expr* undefined();

/**
 * Returns the expression that corresponds to a given integer.
 */
expr* integer(long long i);

/**
 * Returns a fraction expression.
 */
expr* fraction(expr* numerator, expr* denominator);

/**
 * Returns a quotient expression.
 */
expr* quotient(expr* numerator, expr* denominator);

/**
 * Returns a product expression.
 */
expr* product(const std::vector<expr*> operands);
expr* product(const expr* u, const expr* v);
expr* product(const expr* u);

/**
 * Returns a power u^v expression.
 */
expr* power(const expr* u, const expr* v);

/**
 * Returns a difference expression.
 */
expr* difference(const std::vector<expr*> operands);
expr* difference(expr* u, expr* v);

/**
 * Returns a summation expression.
 */
expr* summation(const std::vector<expr*> operands);
expr* summation(const expr* u, const expr* v);
expr* summation(const expr* u);

/**
 * Returns the expression that corresponds to a factorial of u
 */
expr* factorial(const expr* u);

/**
 * Returns the expression that corresponds to a given symbol.
 */
expr* symbol(const char* i);

/**
 * Return a expression with main operator of a given type and operands.
 * 
 * EXAMPLE: construct('+', [a,b,c]) -> a + b + c;
 */
expr* construct(expr::kind kind, std::vector<expr*> operands);

/**
 * Return a expression with main operator of a given type and operands.
 */
expr* construct(expr::kind kind);

/**
 * Destroys a given expression dealocating all its resources and sub trees.
 */
void destroy(expr* u);

/**
 * Return the i'th operand of u.
 * 
 * EXAMPLE: operand(a+ b + c, 2) -> b;
 */
expr* operand(const expr* u, unsigned i);

/**
 * Return the kind of the main operator of u.
 * 
 * EXAMPLE: kind(a + b + 3c) -> +;
 */
expr::kind kind(const expr* u);

/**
 * Return the number of operands of the main operator of u.
 * 
 * EXAMPLE: number_of_operands(a + 2b + 4(a+b)) -> 3;
 */
unsigned number_of_operands(const expr* u);

/**
 * Return true when u is equal to v and false otherwise.
 * 
 * EXAMPLE: equals(a+b, a+b) -> true;
 */
bool equals(const expr* u, const expr* v);

/**
 * Return false when v is a complete sub-expression of u,
 * and true otherwise.
 * 
 * EXAMPLE: free_of((a+b)c, a+b) -> false;
 */
bool free_of(const expr* u, const expr* v);

/**
 * Substitutes the t expression in u by r.
 * 
 * EXAMPLE: substitute(a+2*b+c, 2*b, d) -> a+d+c;
 * 
 * NOTE: The substitute only substitures complete subtrees,
 * so the substitution substitute(a+2*b+c, 2*b+c, d) wont 
 * work because 2*b+c is not a complete subtree of u;
 */
expr* substitute(const expr* u, const expr* t, expr* r);

/**
 * Return expression of kind kind(u) with operands F(operand(u,i))
 * for each i > 0 and i < number_of_operands(u). The interface for
 * the map handler is void (*)(const expr*), it will transform the 
 * given expression to a new one without modifying the original, the
 * original is deleted at the end of each map.
 * 
 * EXAMPLE: map(f(x) = (x+2)^2, x + 2) -> (x+2)^2 + (2+2)^2;
 * EXAMPLE: map(f(x, y) = (x+y)^2, x + 2, 3) -> (x+3)^2 + (2+3)^2;
 */
expr* unary_map(const expr* u, expr* (*f)(const expr*));
expr* binary_map(const expr* u, const expr* v, expr* (*f)(const expr*, const expr*));

/**
 * Return a full copy of a given expression.
 */
expr* copy(const expr* u);

/**
 * Returns the value of a integer expression.
 */
long long integer_value(const expr* u);

/**
 * Returns the value of a symbol expression.
 */
const char* symbol_value(const expr* u);

/**
 * Returns the value of a symbol expression.
 */
const char* function_name(const expr* u);

/**
 * Returns the value of a symbol expression.
 */
const char* symbol_value(const expr* u);

bool is_constant(const core::expr* u);

void print(const core::expr* u);
std::vector<expr*> rest(std::vector<expr*> p);
}

#endif