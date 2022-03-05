
#include "Parser.hpp"
#include "MathSystem/Algebra/Expression.hpp"
#include "MathSystem/Algebra/Lexer.hpp"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
using namespace alg;

Parser::Parser(string src) : lexer{src} {}

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
  Token token = this->lexer.currentToken();

	if (token.type != tokenType) {
    printf("Unexpected token at line '%u' and position '%u', expecting %u!\n",
           token.line, token.position, tokenType);

		exit(-1);
  }

  this->lexer.getToken();
}

expr parseExpression(Parser *parser);
expr parseTernary(Parser *parser);
expr parseEquality(Parser *parser);
expr parseBooleans(Parser *parser);
expr parseBitwise(Parser *parser);
expr parseComparison(Parser *parser);
expr parseBitwiseShift(Parser *parser);
expr parseTerm(Parser *parser);
expr parseFactor(Parser *parser);
expr parseUnary(Parser *parser);
expr parseArrayAccess(Parser *parser);
expr parseMemberAccess(Parser *parser);
expr parsePostfixSuffixUnary(Parser *parser);
expr parsePrimary(Parser *parser);
expr parseLiteral(Parser *parser);

expr parseExpression(Parser *parser) { return parseTernary(parser); }

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

// FACTOR → UNARY (('/' | '*' | '%' ) FACTOR)*
expr parseFactor(Parser *parser) {

  expr root = parseUnary(parser);

	Token curr = parser->currentToken();

		if(curr.type == Token::TOKEN_DIV) {
			parser->readToken(curr.type);

			return root / parseFactor(parser);
		}

		if(curr.type == Token::TOKEN_TIMES) {
			parser->readToken(curr.type);

			return root * parseFactor(parser);
		}


  return root;
}

// UNARY → ('!' | '-' | '+' | '++' | '--' | '~' )UNARY | UNARY('++' | '--') |
// CALL
expr parseUnary(Parser *parser) {
	Token curr = parser->currentToken();

	if(curr.type == Token::TOKEN_MINUS) {
		parser->readToken(Token::TOKEN_MINUS);
		return -1 * parseUnary(parser);
	}

	if(curr.type == Token::TOKEN_PLUS) {
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
    string identifier = parser->currentToken().value;

		parser->readToken(Token::TOKEN_STRING_LITERAL);

    return expr(identifier);
  }

  return parseLiteral(parser);
}

expr parseLiteral(Parser* parser) {
	Token curr = parser->currentToken();

	if(curr.type == Token::TOKEN_INT_LITERAL) {
		expr t = Int::fromString(curr.value.c_str());

		parser->readNumeric();

		return t;
	}

	if(curr.type == Token::TOKEN_FLOAT_LITERAL) {
		double v = strtod(curr.value.c_str(), NULL);

		parser->readNumeric();

		double integral = 0;

		double fractional = std::modf(v, &integral);

		if (fractional <= 2.7e-17) {
			return alg::expr(Int(integral));
		}

		return alg::fraction(Int(integral), Int(fractional));
	}

	printf("Undefined Expression '%s' on line '%u' at position '%u'", curr.value.c_str(), curr.line, curr.position);

	return undefined();
}

expr Parser::parse() {
	return parseExpression(this);
}
