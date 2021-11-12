#include "Core/AST/Integer.hpp"

#include <assert.h>
#include <iostream>

void should_perform_basic_operations() {
	Int a = 2000;
	Int b = 2000;

	Int c = a + b;

	assert(c == 4000);

	assert(c / 2 == 2000);
	assert(c * 3 == 12000);

}

int main() {
	should_perform_basic_operations();
}
