#include "gauss/Gauss.hpp"

#include "test.hpp"

using gauss::polynomial::factorPoly;
using gauss::algebra::pow;

using namespace gauss;
using namespace algebra;

void should_factorize_polynomials() {
	expr x = symbol("x");
	expr y = symbol("y");
	expr z = symbol("z");

	expr f = algebra::pow(x, 2)*(algebra::pow(y, 2)*algebra::pow(z, 2)) + -9;

	assert(factorPoly(9*x) == 9*x);
	assert(factorPoly(f) == (x*y*z + -3)*(x*y*z + 3));
}

int main() {
	TEST(should_factorize_polynomials)
		return 0;
}
