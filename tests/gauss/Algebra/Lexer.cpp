#include<cstdlib>

#include "test.hpp"

#include "gauss/Algebra/Lexer.hpp"

using namespace alg;

void should_lex_expressions() {
	Lexer l0("3 + 4.321");

	Token a = l0.getToken();
	Token b = l0.getToken();
	Token c = l0.getToken();

	assert(a.type == Token::TOKEN_INT_LITERAL);
	assert(a.value == "3");

	assert(b.type == Token::TOKEN_PLUS);

	assert(c.type == Token::TOKEN_FLOAT_LITERAL);
	assert(c.value == "4.321");
}

int main() {
	TEST(should_lex_expressions)
	return 0;
}
