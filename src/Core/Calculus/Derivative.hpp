#ifndef CALCULUS_H
#define CALCULUS_H

#include "Core/AST/AST3.hpp"

namespace calculus {

// alg::expr derivative(alg::expr u, alg::expr x);
alg::expr derivate(alg::expr u, alg::expr x);

}

#endif
