#include "Core/AST/Real.hpp"

#include <assert.h>
#include<iostream>

void should_perform_basic_operations() {
	Real a = 2000;
	Real b = 2000;

	Real c = a + b;

	assert(c == 4000);

	assert(c / 2 == 2000);
	assert(c * 3 == 12000);

	Real d = 3.14;

	assert(d.numerator() == 157);
	assert(d.denominator() == 50);

	std::cout << d.to_string() << std::endl;
}

void should_compute_euler_constant()
{
	Real e = Real::computeEulerConstant(Real(Int(1), Int(10000000000)));
	std::cout << e.to_string() << std::endl;
	Real pi = Real::computePiConstant(Real(Int(1), Int(100)));
	std::cout << pi.to_string() << std::endl;

}

int main() {
	should_perform_basic_operations();
	should_compute_euler_constant();
}
