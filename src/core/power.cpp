#include "power.hpp"

namespace core {

expr* base(const expr* u) {
    if(kind(u) == expr::ALG_OP_POWER) 
        return operand(u, 0);
    
    return (expr*)u;
}

expr* expoent(const expr* u) {
    if(kind(u) == expr::ALG_OP_POWER) 
        return operand(u, 1);

    return integer(1);
}

// bool monomial_non_negative_expoent_property(expr* u) {
//     switch (kind(u)) {
//         case expr::ALG_OP_PRODUCT:
//             for(int i=0; i<number_of_operands(u); i++) {
//                 if(kind(operand(u, i)) != expr::ALG_OP_POWER) {
//                     expr* expoent = operand(operand(u,i), 1);
//                     if(
//                         kind(expoent) != expr::INTEGER ||
//                         (kind(expoent) == expr::INTEGER && *(long long*)expoent->_data <= 0)
//                     ) return false;
//                 }
//             }
//             return true;

//         case expr::ALG_OP_POWER:
//             expr* expoent = operand(u,0);
//             if(
//                 kind(expoent) != expr::INTEGER ||
//                 (kind(expoent) == expr::INTEGER && *(long long*)expoent->_data <= 0)
//             ) return false;
//     }
// }

// bool monomial_indenpendence_property(expr* u) {
//     // TODO: if the operands of u are always sorted
//     // this nested loop can be optmized.

//     // For a polynomial of type c_1*c_2*...c_r*(x_1^n_1)*(x_2^n_2)*...*(x_n^n_n)
//     // where the c_i's are coefficients and the x_j's are the variables, this
//     // function computes if for each possible i, free_of(c_i, x_j) for each possible j.
//     switch (kind(u)) {
//         case expr::ALG_OP_PRODUCT:
//             for(int i=0; i<number_of_operands(u); i++) {
//                 if(kind(operand(u, i)) != expr::ALG_OP_POWER) {
//                     for(int j=0; j<number_of_operands(u); j++) {
//                         if(kind(operand(u, i)) == expr::ALG_OP_POWER) {
//                             if(!free_of(operand(u,i), operand(u,j)))
//                                 return false;
//                         }
//                     }
//                 }
//             }
//         break;
//     }

//     return true;
// }

// bool is_constant(expr* u) {
//     if(kind(u) == expr::FRACTION || kind(u) == expr::INTEGER)
//         return true;

//     return false;
// }

// bool general_monomial_gpe(expr* u, const expr* v) {
//     if(!is_constant(u)) return false;
//     if(!monomial_non_negative_expoent_property(u)) return false;
//     if(!monomial_indenpendence_property(u)) return false;


//     switch (kind(u)) {
//         case expr::ALG_OP_PRODUCT:
//             for(int i=0; i < number_of_operands(u); i++) {
//                 if(general_monomial_gpe(u, v))
//                     return true;
//             }

//         case expr::ALG_OP_POWER:
//             expr* base = operand(u, 0);
//             return equals(base, v);
//     }

//     return equals(u, v);
// }

// bool polynomial_gpe(expr* u, std::vector<const expr*> v) {
//     // A GPE is a GME(general monomial expression) or a sum of GME's
    
//     if (kind(u) == expr::ALG_OP_SUMMATION) {
//         for(int i=0; i<number_of_operands(u); i++) {
//             for(int j=0; j<v.size(); j++) {
//                 if(!general_monomial_gpe(operand(u,i), v[j]))
//                     return false;
//             }
//         }
//         return true;
//     }

//     if (kind(u) == expr::ALG_OP_PRODUCT) {
//         for(int j=0; j<v.size(); j++) {
//             if(!general_monomial_gpe(u, v[j]))
//                 return false;
//         }
//         return true;
//     }

//     return false;
// }
}