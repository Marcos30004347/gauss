
#include "Parser.hpp"

#include "gauss/Algebra/Expression.hpp"
#include "gauss/Algebra/Lexer.hpp"

#include <cstdlib>
#include <ctype.h>
#include <limits>
#include <stdio.h>
#include <stdlib.h>

using namespace alg;

Parser::Parser(std::string src)
    : lexer{src}, current{Token::TOKEN_UNDEFINED, "", 0, 0} {
  moveToNextToken();
}

void Parser::moveToNextToken() { this->current = this->lexer.getToken(); }

bool Parser::isNumeric() {
  return this->currentToken().type == Token::TOKEN_FLOAT_LITERAL ||
         this->currentToken().type == Token::TOKEN_INT_LITERAL;
}

void Parser::readNumeric() {
  if (this->currentToken().type != Token::TOKEN_FLOAT_LITERAL &&
      this->currentToken().type != Token::TOKEN_INT_LITERAL) {
    printf("Expecting Numeric value at line '%u' position '%u'\n",
           this->currentToken().line, this->currentToken().position);
  }

  this->readToken(this->currentToken().type);
}

void Parser::readToken(Token::kind tokenType) {
  Token token = currentToken();

  if (token.type != tokenType) {
    printf("Unexpected token at line '%s' and position '%u', expecting %s!\n",
           tokenToString(token.type), token.position, tokenToString(tokenType));

    exit(-1);
  }

  moveToNextToken();
}

expr parseExpression(Parser *parser);
expr parseTerm(Parser *parser);
expr parseFactor(Parser *parser);
expr parseUnary(Parser *parser);
expr parsePrimary(Parser *parser);
expr parseLiteral(Parser *parser);

expr parseExpression(Parser *parser) { return parseTerm(parser); }

expr parseTerm(Parser *parser) {
  expr root = parseFactor(parser);

  Token curr = parser->currentToken();

  if (curr.type == Token::TOKEN_PLUS) {
    parser->readToken(parser->currentToken().type);
    return root + parseTerm(parser);
  }

  if (curr.type == Token::TOKEN_MINUS) {
    parser->readToken(parser->currentToken().type);
    return root - parseTerm(parser);
  }

  return root;
}

expr parseFactor(Parser *parser) {

  expr root = parseUnary(parser);

  Token curr = parser->currentToken();

  if (curr.type == Token::TOKEN_DIV) {
    parser->readToken(Token::TOKEN_DIV);

    return root / parseFactor(parser);
  }

  if (curr.type == Token::TOKEN_TIMES) {
    parser->readToken(Token::TOKEN_DIV);

    return root * parseFactor(parser);
  }

  return root;
}

expr parseUnary(Parser *parser) {
  Token curr = parser->currentToken();

  if (curr.type == Token::TOKEN_MINUS) {
    parser->readToken(Token::TOKEN_MINUS);
    return -1 * parseUnary(parser);
  }

  if (curr.type == Token::TOKEN_PLUS) {
    parser->readToken(Token::TOKEN_PLUS);
    return parseUnary(parser);
  }

  return parsePrimary(parser);
}

expr parsePrimary(Parser *parser) {
  if (parser->currentToken().type == Token::TOKEN_OPEN_PARENTESIS) {

    parser->readToken(Token::TOKEN_OPEN_PARENTESIS);

    expr expression = parseExpression(parser);

    parser->readToken(Token::TOKEN_CLOSE_PARENTESIS);

    return expression;
  }

  if (parser->currentToken().type == Token::TOKEN_STRING_LITERAL) {
		std::string identifier = parser->currentToken().value;

    parser->readToken(Token::TOKEN_STRING_LITERAL);

    Token curr = parser->currentToken();

    if (curr.type == Token::TOKEN_OPEN_PARENTESIS) {
      parser->readToken(Token::TOKEN_OPEN_PARENTESIS);

      expr call = func_call(identifier.c_str(), {});

      while (true) {

        call.insert(parseExpression(parser));

        if (parser->currentToken().type == Token::TOKEN_CLOSE_PARENTESIS) {
          break;
        }

        parser->readToken(Token::TOKEN_COMMA);
      }

      parser->readToken(Token::TOKEN_CLOSE_PARENTESIS);
    }

    return expr(identifier);
  }

  return parseLiteral(parser);
}

expr parseLiteral(Parser *parser) {

  Token curr = parser->currentToken();

  if (curr.type == Token::TOKEN_INT_LITERAL) {
    expr t = Int::fromString(curr.value.c_str());

    parser->readNumeric();

    return t;
  }

  if (curr.type == Token::TOKEN_FLOAT_LITERAL) {
    double v = strtod(curr.value.c_str(), NULL);

    double integral;

    double fractional = std::modf(v, &integral);

    parser->readNumeric();

    unsigned long long n, d;

    toFraction(fractional, 1000, n, d);

    expr k = Int(integral) + alg::fraction(n, d);

    reduce(&k);

    return k;
  }

  printf("Unexpected expression '%s' on line '%u' at position '%u'",
         curr.value.c_str(), curr.line, curr.position);

  abort();

  return undefined();
}

expr Parser::parse() { return parseExpression(this); }
