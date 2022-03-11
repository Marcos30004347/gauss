#include<cstdlib>

#include "test.hpp"

#include "gauss/Algebra/Parser.hpp"
#include "gauss/Algebra/Expression.hpp"


using namespace alg;

void should_parse_exprs() {
	Parser p0("3 + 4.4");

	expr e0 = p0.parse();

	assert(e0 == expr(3) + expr(22)/expr(5));
}

int main() {
	TEST(should_parse_exprs)
	return 0;
}
