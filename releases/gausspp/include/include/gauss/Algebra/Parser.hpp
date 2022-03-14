#ifndef PARSER_HPP
#define PARSER_HPP

#include "Expression.hpp"
#include "Lexer.hpp"

namespace alg {

class Parser {
private:
  Lexer lexer;

	Token current;
public:
  Parser(std::string src);

  void readToken(Token::kind tokenType);

  void readNumeric();

	bool isNumeric();

	void moveToNextToken();

	inline Token currentToken() { return this->current; }

  expr parse();
};

} // namespace alg
#endif
