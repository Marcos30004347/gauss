#include "Core/AST/Real.hpp"

#include <assert.h>
#include<iostream>

void should_perform_basic_operations() {
	Real a = 2000;
	Real b = 2000;

	Real c = a + b;

	assert(c == 4000);

	assert(c / Int(2) == 2000);
	assert(c * Int(3) == 12000);

	Real d = 3.14;

	assert(d.numerator() == 157);
	assert(d.denominator() == 50);
}

void should_compute_euler_constant()
{
	Real eps = Real(Int(1), Int("10000000000000000000000000000000"));

	Real e = Real::computeEulerConstant(eps);

	assert(e < 2.7182818284590459);
	assert(e > 2.7182818284590450);

	Real pi = Real::computePiConstant(eps);

	assert(pi < 3.141592653589799);
	assert(pi > 3.141592653589790);
}

void should_compute_natural_log()
{
	printf("aa\n");
	Real l = ln(3);
	printf("%s\n", l.to_string().c_str());
}

int main() {
	should_perform_basic_operations();
	should_compute_euler_constant();
	should_compute_natural_log();

}
