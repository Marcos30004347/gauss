#ifndef PARSER_LEXER_H
#define PARSER_LEXER_H

#include "String.hpp"

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
    TOKEN_EOF,
  };

  kind type;

  string value;

  unsigned int line;
  unsigned int position;

  Token(Token::kind type, string value, unsigned line, unsigned pos);
};

class Lexer {
private:
  string source;

  unsigned int head;
  char character;
  unsigned int line;
  unsigned int source_size;
  unsigned int eof;

	Token current;

  bool skipSpaces();

  void advance();

  Token collectIdentifier();
  Token collectStringLiteral();
  Token collectNumberLiteral();

public:
  Lexer(string src);

	inline Token currentToken() { return current; }

	Token getToken();
};

} // namespace alg

#endif
