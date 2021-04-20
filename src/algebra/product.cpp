#include "product.hpp"
#include <assert.h>
namespace algebra {

expression* constant_coefficient(expression* u) {
    if(kind(u) == expression::ALG_OP_PRODUCT) {
        std::vector<expression*> terms;
        for(int i=0; i<number_of_operands(u); i++) {
            if(is_constant(operand(u,i))) {
                terms.push_back(operand(u,i));
            }
        }

        if(terms.size() == 0)
            return integer(1);

        if(terms.size() > 1)
            return product(terms);
        
        return terms[0];
    }

    return integer(1);
}

expression* nonconstant_coefficient(expression* u) {
    // assert(kind(u) == expression::ALG_OP_PRODUCT);
    if(kind(u) == expression::ALG_OP_PRODUCT) {
        std::vector<expression*> terms;
        for(int i=0; i<number_of_operands(u); i++) {
            if(!is_constant(operand(u,i))) {
                terms.push_back(operand(u,i));
            }
        }

        if(terms.size() == 0)
            return undefined();
    
        if(terms.size() > 1)
            return product(terms);
        
        return terms[0];
    }

    return copy(u);
}

}