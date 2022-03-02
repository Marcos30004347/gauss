

#ifdef WASM_BUILD

#include "MathSystem/Algebra/Expression.hpp"
#include "MathSystem/Polynomial/Polynomial.hpp"
#include "MathSystem/Factorization/Wang.hpp"

#define API_FUNC __attribute__( ( visibility( "default" ) ))

API_FUNC alg::expr nondivisors_api(Int G, alg::expr F, Int c, alg::expr L, alg::expr K) {
	return factorization::nondivisors(G, F, c, L, K);
}


#endif
