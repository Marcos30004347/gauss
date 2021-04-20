#ifndef ALGEBRA_POLYNOMIAL_H
#define ALGEBRA_POLYNOMIAL_H

#include "expression.hpp"

namespace algebra {
expression* base(const expression* u);
expression* expoent(const expression* u);
// /**
//  * Return true when u is a GME(General Polynomial Expression)
//  * in the set v.
//  * 
//  * EXAMPLE: polynomial_gpe(x^2 + y^2, [x,y]) -> true;
//  */
// bool polynomial_gpe(const expression* u, std::vector<const expression*> v);


// /**
//  * If v is and polynomial in [v], degree_gpe returns the maximum degree
//  * of v in u.
//  * 
//  * EXAMPLE:  
// expression* degree_gpe(const expression* u, const expression* v);

}
#endif