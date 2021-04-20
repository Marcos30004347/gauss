#include "expression.hpp"

#include <string.h>
#include <cstdlib>
#include <cstdio>

#include "ordering/ordering.hpp"

namespace algebra {

void include_operand(expression* operation, expression* operand) {
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

void remove_operand(expression* operation, expression* operand) {
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

expression* undefined() {
    return construct(expression::UNDEFINED);
}

expression* integer(long long i) {
    expression* u = construct(expression::INTEGER);
    u->_data = new long long(i);
    
    return u;
}

expression* symbol(const char* i) {
    expression* u = construct(expression::SYMBOL);

    u->_data = new char[strlen(i)];
    strcpy((char*)u->_data, i);

    return u;
}

expression* fraction(expression* numerator, expression* denominator) {
    expression* f = construct(expression::FRACTION);
    include_operand(f, numerator);
    include_operand(f, denominator);
    return f;
}

expression* quotient(expression* numerator, expression* denominator) {
    expression* f = construct(expression::ALG_OP_QUOTIENT);
    include_operand(f, numerator);
    include_operand(f, denominator);
    return f;
}

expression* product(const std::vector<expression*> operands) {
    return construct(expression::ALG_OP_PRODUCT, operands);
}

expression* difference(const std::vector<expression*>& operands) {
    return construct(expression::ALG_OP_DIFFERENCE, operands);

}

expression* difference(expression* u, expression* v) {
    expression* f = construct(expression::ALG_OP_DIFFERENCE);
    include_operand(f, u);
    include_operand(f, v);
    return f;
}

expression* summation(const std::vector<expression*> operands) {
    return construct(expression::ALG_OP_SUMMATION, operands);
}

expression* construct(expression::kind kind, std::vector<expression*> operands) {
    expression* u = construct(kind);

    for(int i=0; i < operands.size(); i++)
        include_operand(u, operands[i]);

    return u;
}

expression* construct(expression::kind kind) {
    expression* u = new expression();

    u->_kind = kind;
    u->_data = nullptr;
    u->_operands_count = 0;
    u->_operands = nullptr;

    return u;
}

expression* product(const expression* u, const expression* v) {
    expression* r = construct(expression::ALG_OP_PRODUCT);
    include_operand(r, copy(u));
    include_operand(r, copy(v));
    return r;
}

expression* product(const expression* u) {
    expression* r = construct(expression::ALG_OP_PRODUCT);
    include_operand(r, copy(u));
    return r;
}


expression* power(const expression* u, const expression* v) {
    expression* r = construct(expression::ALG_OP_POWER);
    include_operand(r, copy(u));
    include_operand(r, copy(v));
    return r;
}

expression* summation(const expression* u, const expression* v) {
    expression* r = construct(expression::ALG_OP_SUMMATION);
    include_operand(r, copy(u));
    include_operand(r, copy(v));
    return r; 
}
expression* summation(const expression* u) {
    expression* r = construct(expression::ALG_OP_SUMMATION);
    include_operand(r, copy(u));
    return r;
}

expression* factorial(const expression* u) {
    expression* r = construct(expression::ALG_OP_FACTORIAL);
    include_operand(r, copy(u));
    return r;
}

void destroy(expression* u) {
    if(!u) return;

    if(kind(u) == expression::INTEGER) {
        delete (long long*)u->_data;
    } else if(kind(u) == expression::SYMBOL) {
        delete (char*)u->_data;
    } else {
        for(int i=0; i < number_of_operands(u); i++)
            remove_operand(u, operand(u, i));
    }

    
    delete u;
}

expression::kind kind(const expression* u) {
    return u->_kind;
}

unsigned number_of_operands(const expression* u) {
    if(
        kind(u) == expression::INTEGER ||
        kind(u) == expression::SYMBOL ||
        kind(u) == expression::FUNCTION
    ) return 1;

    return u->_operands_count;
}

expression* operand(const expression* u,unsigned i) {
    if(
        kind(u) == expression::INTEGER ||
        kind(u) == expression::SYMBOL ||
        kind(u) == expression::FUNCTION
    ) return copy(u);

    argument* arg = u->_operands;

    for(int j=0; j<i; j++) {
        arg = arg->_next;
    }

    return copy(arg->_operand);
}

bool free_of(const expression* u, const expression* v) {
    for(int i=0; i<number_of_operands(u); i++) {
        if(equals(operand(u, i), v)) {
            return false;
        }
    }

    return true;
}

bool equals(const expression* a, const expression* b) {

    if(a->_kind != b->_kind)
        return false;
    
    if(number_of_operands(a) != number_of_operands(b))
        return false;

    expression::kind exp_kind = kind(a);

    // compare expressions that have meaningfull data
    if(exp_kind == expression::INTEGER)
        return integer_value(a) == integer_value(b);
    if(exp_kind == expression::SYMBOL)
        return strcmp(symbol_value(a), symbol_value(b)) == 0;

    // order of the operators dont matter
    if(
        exp_kind == expression::ALG_OP_SUMMATION ||
        exp_kind == expression::ALG_OP_DIFFERENCE ||
        exp_kind == expression::ALG_OP_PRODUCT
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

expression* substitute(const expression* u, const expression* t, expression* r) {
    expression* v = construct(kind(u));

    for(int i=0; i < number_of_operands(u); i++) {
        expression* o = operand(u, i);
        if(equals(o, t))
            include_operand(v, copy(r));
        else
            include_operand(v, o);
    }

    return v;
}

expression* unary_map(const expression* u, expression* (*f)(const expression*)) {    
    if(kind(u) == expression::INTEGER || kind(u) == expression::SYMBOL) 
        return f(u);

    if(number_of_operands(u) == 0)
        return f(u);

    expression* v = construct(kind(u));

    for(int i=0; i < number_of_operands(u); i++)
        include_operand(v, f(operand(u,i)));

    return v;
}

expression* binary_map(const expression* u, const expression* v, expression* (*f)(const expression*, const expression*)) {
    if(kind(u) == expression::INTEGER || kind(u) == expression::SYMBOL) 
        return f(u, v);

    if(number_of_operands(u) == 0)
        return f(u, v);

    expression* r = construct(kind(u));

    for(int i=0; i< number_of_operands(u); i++)
        include_operand(r, f(operand(u,i), v));

    return r;
}

expression* copy(const expression* u) {
    expression* v = construct(kind(u));

    switch (kind(u)) {
    case expression::SYMBOL:
        v->_data = new char[strlen((char*)u->_data)];
        strcpy((char*)v->_data, (char*)u->_data);
        break;
    case expression::INTEGER:
        v->_data = new long long(*((long long*)u->_data));
        break;
    default:
        for(int i=0; i < number_of_operands(u); i++)
            include_operand(v, copy(operand(u, i)));    
        break;
    }

    return v;
}

long long integer_value(const expression* u) {
    return *(long long*)u->_data;
}

const char* symbol_value(const expression* u) {
    return (const char*)u->_data;
}


const char* function_name(const expression* u) {
    return symbol_value(operand(u,0));
}

bool is_constant(const expression* u) {
    return kind(u) == expression::INTEGER || kind(u) == expression::FRACTION;
}

void print(const expression* u) {

    switch(kind(u)) {
        case expression::ALG_OP_SUMMATION:
        for(int i=0; i<number_of_operands(u); i++) {
            print(operand(u, i));
            if(i != number_of_operands(u) -1)
                printf(" + ");
        }
        break;
        
        case expression::ALG_OP_DIFFERENCE:
        for(int i=0; i<number_of_operands(u); i++) {
            print(operand(u, i));
            if(i != number_of_operands(u) -1)
                printf(" - ");
        }
        break;
        
        case expression::ALG_OP_POWER:
        printf("(");
        printf("(");
        print(operand(u, 0));
        printf(")^");
        print(operand(u, 1));
        printf(")");
        break;
        
        case expression::ALG_OP_PRODUCT:
        for(int i=0; i<number_of_operands(u); i++) {
            if(
                kind(operand(u,i)) == expression::ALG_OP_SUMMATION || 
                kind(operand(u,i)) == expression::ALG_OP_DIFFERENCE 
            ) printf("(");
            print(operand(u, i));
            if(
                kind(operand(u,i)) == expression::ALG_OP_SUMMATION || 
                kind(operand(u,i)) == expression::ALG_OP_DIFFERENCE 
            ) printf(")");

            if(i != number_of_operands(u) -1)
                printf("⋅");
                // ×
        }
        break;
        
        case expression::ALG_OP_QUOTIENT:
        printf("(");
        print(operand(u, 0));
        printf(")");
        printf(" / ");
        printf("(");
        print(operand(u, 1));
        printf(")");
        break;
        
        case expression::ALG_OP_FACTORIAL:
        printf("!");
        printf("(");
        print(operand(u, 0));
        printf(")");
        break;
        
        case expression::FRACTION:
        print(operand(u, 0));
        printf(" / ");
        print(operand(u, 1));
        break;
        
        case expression::INTEGER:
        printf("%lld", *(long long*)u->_data);
        break;
        
        case expression::SYMBOL:
        printf("%s", (const char*)u->_data);
        break;
    
        default:
            break;
    }
}


std::vector<expression*> rest(std::vector<expression*> p, int from) {
    std::vector<expression*> a;
    for(int i=from; i <= p.size() - 1; i++) 
        a.push_back(copy(p[i]));
    return a;
}

std::vector<expression*> equal_operands(expression* u, expression* v) {
    std::vector<expression*> m;

    for(int i=0; i<number_of_operands(u); i++) {
        for(int j=i; j<number_of_operands(v); j++) {
            if(equals(operand(u, i), operand(v, j)))
                m.push_back(operand(u, i));
        }
    }

    return m;
}

// std::vector<expression*> different_operands(expression* u, expression* v) {
//     std::vector<expression*> r; // difference between a and b
//     std::vector<expression*> m = equal_operands(u, v);
    
//     bool e = false;
//     for(int i=0; i<number_of_operands(u); i++) {
//         for(int j=0; j<m.size(); j++) {
//             if(equals(operand(u,i), m[j])) {
//                 e = true;
//                 break;
//             }
//         }
//         if(e) continue;
//         r.push_back(operand(u, i));
//     }

//     return r;
// }

} // namespace core