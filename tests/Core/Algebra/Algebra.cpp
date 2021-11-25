#include "Core/Algebra/Algebra.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;
/*
void should_create_algebraic_expressions() {
        Expr e0 = Expr(3);
        Expr e1 = Expr(1) + Expr(2) + Expr(3);
        Expr e2 = Expr(1) - Expr(2) - Expr(3);
        Expr e3 = "x";
        Expr e4 = fraction(1, 2);
        Expr e5 = power("x", 4);
        Expr e7 = Expr(2) * Expr(5);

        assert(e0.kind() == Kind::Integer);
        assert(e0.value() == 3);

        assert(e1.kind() == Kind::Addition);
        assert(e1[0].kind() == Kind::Integer);
        assert(e1[1].kind() == Kind::Integer);
        assert(e1[2].kind() == Kind::Integer);
        assert(e1[0].value() == 1);
        assert(e1[1].value() == 2);
        assert(e1[2].value() == 3);

        assert(e2.kind() == Kind::Subtraction);
        assert(e2[0].kind() == Kind::Integer);
        assert(e2[1].kind() == Kind::Integer);
        assert(e2[2].kind() == Kind::Integer);
        assert(e2[0].value() == 1);
        assert(e2[1].value() == 2);
        assert(e2[2].value() == 3);

        assert(e3.kind() == Kind::Symbol);
        assert(e3.identifier() == "x");

        assert(e4.kind() == Kind::Fraction);
        assert(e4[0].kind() == Kind::Integer);
        assert(e4[1].kind() == Kind::Integer);
        assert(e4[0].value() == 1);
        assert(e4[1].value() == 2);

        assert(e5.kind() == Kind::Power);
        assert(e5[0].kind() == Kind::Symbol);
        assert(e5[1].kind() == Kind::Integer);
        assert(e5[0].identifier() == "x");
        assert(e5[1].value() == 4);

        assert(e7.kind() == Kind::Multiplication);
        assert(e7[0] == 2);
        assert(e7[1] == 5);
}
*/

void should_get_info_of_algebraic_expressions() {
  Expr exp0 = integer(3);
  Expr exp1 = fraction(4, 2);
  Expr exp2 = power(integer(5), integer(2));
  Expr exp3 = base(exp0);
  Expr exp4 = expoent(exp0);
  Expr exp5 = base(exp2);
  Expr exp6 = expoent(exp2);

  assert(exp3 == 3);
  assert(exp4 == 1);
  assert(exp5 == 5);
  assert(exp6 == 2);
}

void should_order_realate_expressions() {
  Expr exp0 = integer(1);
  Expr exp1 = integer(2);
  Expr exp2 = symbol("x");
  Expr exp3 = power(symbol("x"), integer(2));
  Expr exp4 = add({symbol("x"), power(symbol("x"), integer(2))});
  Expr exp5 = add({symbol("x"), power(symbol("x"), integer(3))});
  Expr exp6 = mul({symbol("x"), power(symbol("x"), integer(2))});
  Expr exp7 = mul({symbol("x"), power(symbol("x"), integer(3))});

  assert(orderRelation(exp0, exp1));
  assert(orderRelation(exp1, exp2));
  assert(orderRelation(exp2, exp3));
  assert(orderRelation(exp4, exp5));
  assert(orderRelation(exp6, exp7));
}

int main() {
  //	should_create_algebraic_expressions();
  should_get_info_of_algebraic_expressions();
  should_order_realate_expressions();
  return 0;
}
