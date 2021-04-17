#ifndef CORE_POLYNOMIAL_H
#define CORE_POLYNOMIAL_H

#include "expression.hpp"

namespace core {
expr* base(const expr* u);
expr* expoent(const expr* u);
// /**
//  * Return true when u is a GME(General Polynomial Expression)
//  * in the set v.
//  * 
//  * EXAMPLE: polynomial_gpe(x^2 + y^2, [x,y]) -> true;
//  */
// bool polynomial_gpe(const expr* u, std::vector<const expr*> v);


// /**
//  * If v is and polynomial in [v], degree_gpe returns the maximum degree
//  * of v in u.
//  * 
//  * EXAMPLE:  
// expr* degree_gpe(const expr* u, const expr* v);

}
#endif