#ifndef PARSER_LEXER_H
#define PARSER_LEXER_H

#include <string>

namespace alg {

class Token {
public:
  enum kind {
		TOKEN_UNDEFINED,
    TOKEN_DIV,
    TOKEN_TIMES,
    TOKEN_MINUS,
    TOKEN_PLUS,
    TOKEN_POW,
    TOKEN_OPEN_PARENTESIS,
    TOKEN_CLOSE_PARENTESIS,
    TOKEN_OPEN_SQUARE_BRACKETS,
    TOKEN_CLOSE_SQUARE_BRACKETS,
    TOKEN_FLOAT_LITERAL,
    TOKEN_INT_LITERAL,
    TOKEN_STRING_LITERAL,
		TOKEN_COMMA,
    TOKEN_EOF,
  };

  kind type;

	std::string value;

  unsigned int line;
  unsigned int position;

  Token(Token::kind type, std::string value, unsigned line, unsigned pos);
};

const char* tokenToString(Token::kind t);

class Lexer {
private:
	std::string source;

  unsigned int head;
  char character;
  unsigned int line;
  unsigned int source_size;
  unsigned int eof;

  bool skipSpaces();

  void advance();

  Token collectIdentifier();
  Token collectStringLiteral();
  Token collectNumberLiteral();

public:
  Lexer(std::string src);

	Token getToken();
};

} // namespace alg

#endif
