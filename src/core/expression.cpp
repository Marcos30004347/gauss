#include "expression.hpp"

#include <string.h>

namespace core {

void include_operand(expr* operation, expr* operand) {
    if(!operation) return;
    if(!operand) return;

    operation->_operands_count++;
    
    if(!operation->_operands) {
        argument* arg = new argument();
        arg->_operand = operand;
        arg->_operation = operation;

        arg->_next = nullptr;
        arg->_prev = nullptr;
    
        operation->_operands = arg;
        return;
    }

    argument* arg = operation->_operands;

    while(arg->_next)
        arg = arg->_next;
    
    arg->_next = new argument();
    arg->_next->_operand = operand;
    arg->_next->_operation = operation;

    arg->_next->_next = nullptr;
    arg->_next->_prev = arg;
}

void remove_operand(expr* operation, expr* operand) {
    argument* arg = operation->_operands;

    while(arg && !equals(arg->_operand, operand))
        arg = arg->_next;
    
    if(arg) {
        operation->_operands_count--;

        if(arg->_next) arg->_next->_prev = arg->_prev;
        if(arg->_prev) arg->_prev->_next = arg->_next;
        
        destroy(arg->_operand);
        free(arg);
    }
}

expr* undefined() {
    return construct(expr::UNDEFINED);
}

expr* integer(long long i) {
    expr* u = construct(expr::INTEGER);
    u->_data = new long long(i);
    
    return u;
}

expr* symbol(const char* i) {
    expr* u = construct(expr::SYMBOL);

    u->_data = new char[strlen(i)];
    strcpy((char*)u->_data, i);

    return u;
}

expr* fraction(expr* numerator, expr* denominator) {
    expr* f = construct(expr::FRACTION);
    include_operand(f, numerator);
    include_operand(f, denominator);
    return f;
}

expr* quotient(expr* numerator, expr* denominator) {
    expr* f = construct(expr::ALG_OP_QUOTIENT);
    include_operand(f, numerator);
    include_operand(f, denominator);
    return f;
}

expr* product(const std::vector<expr*> operands) {
    return construct(expr::ALG_OP_PRODUCT, operands);
}

expr* difference(const std::vector<expr*>& operands) {
    return construct(expr::ALG_OP_DIFFERENCE, operands);

}

expr* difference(expr* u, expr* v) {
    expr* f = construct(expr::ALG_OP_DIFFERENCE);
    include_operand(f, u);
    include_operand(f, v);
    return f;
}

expr* summation(const std::vector<expr*> operands) {
    return construct(expr::ALG_OP_SUMMATION, operands);
}

expr* construct(expr::kind kind, std::vector<expr*> operands) {
    expr* u = construct(kind);

    for(int i=0; i < operands.size(); i++)
        include_operand(u, operands[i]);

    return u;
}

expr* construct(expr::kind kind) {
    expr* u = new expr();

    u->_kind = kind;
    u->_data = nullptr;
    u->_operands_count = 0;
    u->_operands = nullptr;

    return u;
}

expr* product(const expr* u, const expr* v) {
    expr* r = construct(expr::ALG_OP_PRODUCT);
    include_operand(r, copy(u));
    include_operand(r, copy(v));
    return r;
}

expr* product(const expr* u) {
    expr* r = construct(expr::ALG_OP_PRODUCT);
    include_operand(r, copy(u));
    return r;
}


expr* power(const expr* u, const expr* v) {
    expr* r = construct(expr::ALG_OP_POWER);
    include_operand(r, copy(u));
    include_operand(r, copy(v));
    return r;
}

expr* summation(const expr* u, const expr* v) {
    expr* r = construct(expr::ALG_OP_SUMMATION);
    include_operand(r, copy(u));
    include_operand(r, copy(v));
    return r; 
}
expr* summation(const expr* u) {
    expr* r = construct(expr::ALG_OP_SUMMATION);
    include_operand(r, copy(u));
    return r;
}

expr* factorial(const expr* u) {
    expr* r = construct(expr::ALG_OP_FACTORIAL);
    include_operand(r, copy(u));
    return r;
}

void destroy(expr* u) {
    if(!u) return;

    if(kind(u) == expr::INTEGER) {
        delete (long long*)u->_data;
    } else if(kind(u) == expr::SYMBOL) {
        delete (char*)u->_data;
    } else {
        for(int i=0; i < number_of_operands(u); i++)
            remove_operand(u, operand(u, i));
    }

    
    delete u;
}

expr::kind kind(const expr* u) {
    return u->_kind;
}

unsigned number_of_operands(const expr* u) {
    if(
        kind(u) == expr::INTEGER ||
        kind(u) == expr::SYMBOL ||
        kind(u) == expr::FUNCTION
    ) return 1;

    return u->_operands_count;
}

expr* operand(const expr* u,unsigned i) {
    if(
        kind(u) == expr::INTEGER ||
        kind(u) == expr::SYMBOL ||
        kind(u) == expr::FUNCTION
    ) return copy(u);

    argument* arg = u->_operands;

    for(int j=0; j<i; j++) {
        arg = arg->_next;
    }

    return copy(arg->_operand);
}

bool free_of(const expr* u, const expr* v) {
    for(int i=0; i<number_of_operands(u); i++) {
        if(equals(operand(u, i), v)) {
            return false;
        }
    }

    return true;
}

bool equals(const expr* a, const expr* b) {

    if(a->_kind != b->_kind)
        return false;
    
    if(number_of_operands(a) != number_of_operands(b))
        return false;

    expr::kind exp_kind = kind(a);

    // compare expressions that have meaningfull data
    if(exp_kind == expr::INTEGER)
        return integer_value(a) == integer_value(b);
    if(exp_kind == expr::SYMBOL)
        return strcmp(symbol_value(a), symbol_value(b)) == 0;

    // order of the operators dont matter
    if(
        exp_kind == expr::ALG_OP_SUMMATION ||
        exp_kind == expr::ALG_OP_DIFFERENCE ||
        exp_kind == expr::ALG_OP_PRODUCT
    ) {
        long matches = 0;
        long match = false;

        for(int i=0; i < number_of_operands(a); i++) {
            for(int j=0; j < number_of_operands(b); j++) {
                if(equals(operand(a,i), operand(b,j))) {
                    matches++;
                    match = true;
                    break;
                }
            }

            if(match) {
                match = false;
                continue;
            }
        }

        return matches == number_of_operands(a);    
    }


    // order of the operators does matter
    for(int i=0; i < number_of_operands(a); i++) 
        if(!equals(operand(a,i), operand(b,i)))
            return false;
    
    return true;
}

expr* substitute(const expr* u, const expr* t, expr* r) {
    expr* v = construct(kind(u));

    for(int i=0; i < number_of_operands(u); i++) {
        expr* o = operand(u, i);
        if(equals(o, t))
            include_operand(v, copy(r));
        else
            include_operand(v, o);
    }

    return v;
}

expr* unary_map(const expr* u, expr* (*f)(const expr*)) {    
    if(kind(u) == expr::INTEGER || kind(u) == expr::SYMBOL) 
        return f(u);

    if(number_of_operands(u) == 0)
        return f(u);

    expr* v = construct(kind(u));

    for(int i=0; i < number_of_operands(u); i++)
        include_operand(v, f(operand(u,i)));

    return v;
}

expr* binary_map(const expr* u, const expr* v, expr* (*f)(const expr*, const expr*)) {
    if(kind(u) == expr::INTEGER || kind(u) == expr::SYMBOL) 
        return f(u, v);

    if(number_of_operands(u) == 0)
        return f(u, v);

    expr* r = construct(kind(u));

    for(int i=0; i< number_of_operands(u); i++)
        include_operand(r, f(operand(u,i), v));

    return r;
}

expr* copy(const expr* u) {
    expr* v = construct(kind(u));

    switch (kind(u)) {
    case expr::SYMBOL:
        v->_data = new char[strlen((char*)u->_data)];
        strcpy((char*)v->_data, (char*)u->_data);
        break;
    case expr::INTEGER:
        v->_data = new long long(*((long long*)u->_data));
        break;
    default:
        for(int i=0; i < number_of_operands(u); i++)
            include_operand(v, copy(operand(u, i)));    
        break;
    }

    return v;
}

long long integer_value(const expr* u) {
    return *(long long*)u->_data;
}

const char* symbol_value(const expr* u) {
    return (const char*)u->_data;
}


const char* function_name(const expr* u) {
    return symbol_value(operand(u,0));
}

bool is_constant(const expr* u) {
    return kind(u) == expr::INTEGER || kind(u) == expr::FRACTION;
}

void print(const core::expr* u) {

    switch(kind(u)) {
        case expr::ALG_OP_SUMMATION:
        printf("(");
        for(int i=0; i<number_of_operands(u); i++) {
            print(operand(u, i));
            if(i != number_of_operands(u) -1)
                printf(" + ");
        }
        printf(")");
        break;
        
        case expr::ALG_OP_DIFFERENCE:
        printf("(");
        for(int i=0; i<number_of_operands(u); i++) {
            print(operand(u, i));
            if(i != number_of_operands(u) -1)
                printf(" - ");
        }
        printf(")");
        break;
        
        case expr::ALG_OP_POWER:
        printf("(");
        printf("(");
        print(operand(u, 0));
        printf(")^");
        print(operand(u, 1));
        printf(")");
        break;
        
        case expr::ALG_OP_PRODUCT:
        printf("(");
        for(int i=0; i<number_of_operands(u); i++) {
            print(operand(u, i));
            if(i != number_of_operands(u) -1)
                printf(" * ");
        }
        printf(")");
        break;
        
        case expr::ALG_OP_QUOTIENT:
        printf("(");
        printf("(");
        print(operand(u, 0));
        printf(")");
        printf(" / ");
        printf("(");
        print(operand(u, 1));
        printf(")");
        printf(")");
        break;
        
        case expr::ALG_OP_FACTORIAL:
        printf("(");
        printf("!");
        printf("(");
        print(operand(u, 0));
        printf(")");
        printf(")");
        break;
        
        case expr::FRACTION:
        printf("(");
        print(operand(u, 0));
        printf(")");
        printf(" / ");
        printf("(");
        print(operand(u, 1));
        printf(")");
        break;
        
        case expr::INTEGER:
        printf("%lld", *(long long*)u->_data);
        break;
        
        case expr::SYMBOL:
        printf("%s", (const char*)u->_data);
        break;
    
        default:
            break;
    }
}
std::vector<expr*> rest(std::vector<expr*> p) {
    std::vector<expr*> a;
    for(int i=1; i <= p.size() - 1; i++) 
        a.push_back(copy(p[i]));
    return a;
}


} // namespace core