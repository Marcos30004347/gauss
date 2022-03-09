#ifndef SVD_HPP
#define SVD_HPP

#include <tuple>
#include "MathSystem/Algebra/Matrix.hpp"

namespace alg {
	std::tuple<matrix, matrix, matrix> svd(matrix a, double tol = 2.22e-15);
};


#endif
