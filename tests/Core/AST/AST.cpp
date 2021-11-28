#include "Core/AST/AST.hpp"
#include "test.hpp"

using namespace ast;

void should_create_ast_nodes() {
  Expr ast0 = Expr(Kind::Undefined);
  Expr ast1 = Expr(3) + Expr(4) + Expr(5);
  Expr ast2 = Expr(3) + Expr(4) + Expr(5) / Expr(6);
  Expr ast3 = 0;

  assert(ast0.kind() == Kind::Undefined);
  assert(ast1.kind() == Kind::Addition);
  assert(ast1[0].kind() == Kind::Integer);
  assert(ast1[0].value() == 3);
  assert(ast1[1].kind() == Kind::Integer);
  assert(ast1[1].value() == 4);
  assert(ast1[2].kind() == Kind::Integer);
  assert(ast1[2].value() == 5);

  assert(ast2.kind() == Kind::Addition);
  assert(ast2[0].kind() == Kind::Integer);
  assert(ast2[0].value() == 3);
  assert(ast2[1].kind() == Kind::Integer);
  assert(ast2[1].value() == 4);
  assert(ast2[2].kind() == Kind::Fraction);
  assert(ast2[2][0].kind() == Kind::Integer);
  assert(ast2[2][0].value() == 5);
  assert(ast2[2][1].kind() == Kind::Integer);
  assert(ast2[2][1].value() == 6);

  assert(ast3.kind() == Kind::Integer);
  assert(ast3.value() == 0);
  assert(ast3 == 0);
}

void should_match_ast_nodes() {
  Expr ast0 = Expr(3) + Expr(4) + Expr(5);
  Expr ast1 = Expr(3) + Expr(4) + Expr(5);
  Expr ast2 = Expr(3) + Expr(4);

  assert(ast0 == ast1);
  assert(ast0 != ast2);
}

void should_deep_copy_ast_nodes() {
  Expr ast0 = Expr(3) + Expr(4) + Expr(5);
  Expr ast1 = ast0;

  assert(ast0 == ast1);
}

int main() {
  TEST(should_create_ast_nodes)
  TEST(should_match_ast_nodes)
  TEST(should_deep_copy_ast_nodes)
  return 0;
}
