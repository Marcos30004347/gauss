#include "Core/Polynomial/Resultant.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "test.hpp"
#include <cstdio>

using namespace ast;
using namespace algebra;
using namespace polynomial;

void should_get_univariate_resultant() {
  Expr ux = add({mul({integer(2), power(symbol("x"), integer(3))}),
                 mul({integer(-3), symbol("x")}), integer(1)});

  Expr vx = add({mul({integer(3), power(symbol("x"), integer(2))}),
                 mul({integer(-4), symbol("x")}), integer(3)});

  Expr x = symbol("x");

  Expr r0 = univariateResultant(ux, vx, x);

  assert(r0.kind() == Kind::Integer);
  assert(r0.value() == 218);
}

void should_get_multivariate_resultants0() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = add({
      mul({
          power(symbol("x"), integer(3)),
          power(symbol("y"), integer(2)),
      }),
      mul({
          integer(6),
          power(symbol("x"), integer(4)),
          symbol("y"),
      }),
      mul({
          integer(9),
          power(symbol("x"), integer(5)),
      }),
      mul({
          integer(4),
          power(symbol("x"), integer(2)),
          power(symbol("y"), integer(2)),
      }),
      mul({
          integer(24),
          power(symbol("x"), integer(3)),
          symbol("y"),
      }),
      mul({
          integer(36),
          power(symbol("x"), integer(3)),
      }),
      mul({
          integer(5),
          symbol("x"),
          power(symbol("y"), integer(2)),
      }),
      mul({
          integer(45),
          power(symbol("x"), integer(3)),
      }),
      mul({
          integer(2),
          power(symbol("y"), integer(2)),
      }),
      mul({
          integer(12),
          symbol("y"),
          symbol("x"),
      }),
      mul({
          integer(18),
          power(symbol("x"), integer(2)),
      }),
  });

  Expr v = add({mul({
                    power(symbol("x"), integer(5)),
                    power(symbol("y"), integer(2)),
                }),
                mul({
                    integer(8),
                    power(symbol("x"), integer(4)),
                    symbol("y"),
                }),
                mul({
                    integer(16),
                    power(symbol("x"), integer(3)),
                }),
                mul({
                    integer(12),
                    power(symbol("x"), integer(4)),
                    power(symbol("y"), integer(2)),
                }),
                mul({
                    integer(96),
                    power(symbol("x"), integer(3)),
                    symbol("y"),
                }),
                mul({
                    integer(192),
                    power(symbol("x"), integer(2)),
                }),
                mul({
                    integer(45),
                    power(symbol("x"), integer(3)),
                    power(symbol("y"), integer(2)),
                }),
                mul({
                    integer(360),
                    symbol("y"),
                    power(symbol("x"), integer(2)),
                }),
                mul({integer(720), symbol("x")}),
                mul({
                    integer(50),
                    power(symbol("x"), integer(2)),
                    power(symbol("y"), integer(2)),
                }),
                mul({
                    integer(400),
                    symbol("y"),
                    symbol("x"),
                }),
                integer(800)});

  Expr K = symbol("Z");
  Expr L = list({symbol("x"), symbol("y")});

  Expr r0 = srPolynomialResultant(u, v, L, K);

  assert(r0.kind() == Kind::Integer);
  assert(r0.value() == 0);
}

void should_get_multivariate_resultants() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = power(x, 3) * power(y, 3) + 6 * power(x, 2) * y +
           5 * x * power(y, 2) + 2 * power(y, 2) + y * x + 3 * power(x, 2);
  Expr v = power(x, 2) * power(y, 2) + 5 * power(x, 3) + 3 * power(x, 3) * y +
           4 * y * x + 8;

  Expr K = Expr("Z");
  Expr L = list({Expr("x"), Expr("y")});

  Expr r0 = multivariateResultant(u, v, L, K);
  Expr r1 = srPolynomialResultant(u, v, L, K);

  assert(r0 == 0);
  assert(r1 == 0);
}

void should_get_remainder_sequence() {
  Expr u = add({power(symbol("x"), integer(8)), power(symbol("x"), integer(6)),
                mul({
                    integer(-3),
                    power(symbol("x"), integer(4)),
                }),
                mul({
                    integer(-3),
                    power(symbol("x"), integer(3)),
                }),
                mul({
                    integer(8),
                    power(symbol("x"), integer(2)),
                }),
                mul({integer(2), symbol("x")}), integer(-5)});

  Expr v = add({mul({
                    integer(3),
                    power(symbol("x"), integer(6)),
                }),
                mul({
                    integer(5),
                    power(symbol("x"), integer(4)),
                }),
                mul({
                    integer(-4),
                    power(symbol("x"), integer(2)),
                }),
                mul({integer(-9), symbol("x")}), integer(21)});

  Expr L = list({symbol("x")});
  Expr Q = symbol("Q");

  // PPRS(u, v, x);
  Expr r = polyRemSeq(u, v, L, Q);

  assert(r[0].kind() == Kind::Integer);
  assert(r[0].value() == 1);

  assert(r[1].kind() == Kind::Integer);
  assert(r[1].value() == 260708);
}

void should_get_remainder_sequence_mv() {
  Expr f = add({
      mul({integer(3), symbol("y"), power(symbol("x"), integer(2))}),
      mul({integer(-1), add({power(symbol("y"), integer(3)), integer(4)})}),
  });

  Expr g =
      add({power(symbol("x"), integer(2)),
           mul({power(symbol("y"), integer(3)), symbol("x")}), integer(-9)});

  Expr L = list({symbol("x"), symbol("y")});
  Expr Z = symbol("Z");

  Expr r = polyRemSeq(f, g, L, Z);

  assert(r[0].kind() == Kind::Integer);
  assert(r[0].value() == 1);

  Expr res = add({mul({integer(-3), power(symbol("y"), integer(10))}),
                  mul({integer(-12), power(symbol("y"), integer(7))}),
                  power(symbol("y"), integer(6)),
                  mul({integer(-54), power(symbol("y"), integer(4))}),
                  mul({integer(8), power(symbol("y"), integer(3))}),
                  mul({integer(729), power(symbol("y"), integer(2))}),
                  mul({integer(-216), symbol("y")}), integer(16)});

  assert(r[1] == (res));
}

void should_get_remainder_sequence_mv1() {
  Expr u = add({power(symbol("x"), integer(8)), power(symbol("x"), integer(6)),
                mul({
                    integer(-3),
                    power(symbol("x"), integer(4)),
                }),
                mul({
                    integer(-3),
                    power(symbol("x"), integer(3)),
                }),
                mul({
                    integer(8),
                    power(symbol("x"), integer(2)),
                }),
                mul({integer(2), symbol("x")}), integer(-5)});

  Expr v = add({mul({
                    integer(3),
                    power(symbol("x"), integer(6)),
                }),
                mul({
                    integer(5),
                    power(symbol("x"), integer(4)),
                }),
                mul({
                    integer(-4),
                    power(symbol("x"), integer(2)),
                }),
                mul({integer(-9), symbol("x")}), integer(21)});

  Expr L = list({symbol("x")});

  Expr Z = symbol("Q");

  Expr r = polyRemSeq(u, v, L, Z);

  assert(r[0].kind() == Kind::Integer);
  assert(r[0].value() == 1);

  assert(r[1].kind() == Kind::Integer);
  assert(r[1].value() == 260708);
}

void should_get_remainder_sequence_mv2() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = power(x, 3) * power(y, 2) + 6 * power(x, 4) * y + 9 * power(x, 5) +
           4 * power(x, 5) + 4 * power(x, 2) * power(y, 2) +
           24 * power(x, 3) * y + 36 * power(x, 4) + 5 * x * power(y, 2) +
           30 * y * power(x, 2) + 45 * power(x, 3) + 2 * power(y, 2) +
           12 * x * y + 18 * power(x, 2);
	/*
  add({
      mul({
          power(symbol("x"), integer(3)),
          power(symbol("y"), integer(2)),
      }),
      mul({
          integer(6),
          power(symbol("x"), integer(4)),
          symbol("y"),
      }),
      mul({
          integer(9),
          power(symbol("x"), integer(5)),
      }),
      mul({
          integer(4),
          power(symbol("x"), integer(2)),
          power(symbol("y"), integer(2)),
      }),
      mul({
          integer(24),
          power(symbol("x"), integer(3)),
          symbol("y"),
      }),
      mul({
          integer(36),
          power(symbol("x"), integer(4)),
      }),
      mul({
          integer(5),
          symbol("x"),
          power(symbol("y"), integer(2)),
      }),
      mul({
          integer(30),
          symbol("y"),
          power(symbol("x"), integer(2)),
      }),
      mul({
          integer(45),
          power(symbol("x"), integer(3)),
      }),
      mul({
          integer(2),
          power(symbol("y"), integer(2)),
      }),
      mul({
          integer(12),
          symbol("y"),
          symbol("x"),
      }),
      mul({
          integer(18),
          power(symbol("x"), integer(2)),
      }),
  });
	*/
  Expr v = power(x, 5)*power(y, 2) + 8*power(x, 4)*y + 16*power(x, 3) + 12*power(x, 4)*power(y, 2) + 96*power(x, 3)*y + 192*power(x, 2) + 45*power(x, 3)*power(y, 2) + 360*y*power(x, 2) + 720*x + 50*power(x,2)*power(y, 2) + 400*x*y + 800;
	/*
		add({mul({
                    power(symbol("x"), integer(5)),
                    power(symbol("y"), integer(2)),
                }),
                mul({
                    integer(8),
                    power(symbol("x"), integer(4)),
                    symbol("y"),
                }),
                mul({
                    integer(16),
                    power(symbol("x"), integer(3)),
                }),
                mul({
                    integer(12),
                    power(symbol("x"), integer(4)),
                    power(symbol("y"), integer(2)),
                }),
                mul({
                    integer(96),
                    power(symbol("x"), integer(3)),
                    symbol("y"),
                }),
                mul({
                    integer(192),
                    power(symbol("x"), integer(2)),
                }),
                mul({
                    integer(45),
                    power(symbol("x"), integer(3)),
                    power(symbol("y"), integer(2)),
                }),
                mul({
                    integer(360),
                    symbol("y"),
                    power(symbol("x"), integer(2)),
                }),
                mul({integer(720), symbol("x")}),
                mul({
                    integer(50),
                    power(symbol("x"), integer(2)),
                    power(symbol("y"), integer(2)),
                }),
                mul({
                    integer(400),
                    symbol("y"),
                    symbol("x"),
                }),
                integer(800)});
	*/
  Expr L = list({ x, y });

  Expr Z = Expr("Z");

  Expr r = polyRemSeq(u, v, L, Z);

	printf("%s\n", r.toString().c_str());

	Expr uv_gcd = add({integer(2), symbol("x")});

  assert(r[0] == (uv_gcd));
  assert(r[1].kind() == Kind::Integer);
  assert(r[1].value() == 0);
}

void should_get_remainder_sequence_mv3() {
  Expr t = add({power(symbol("z"), integer(4)), power(symbol("z"), integer(3)),
                mul({add({integer(2), symbol("x"),
                          mul({
                              integer(-1),
                              power(symbol("x"), integer(2)),
                          })}),
                     power(symbol("z"), integer(2))}),
                mul({add({integer(1), power(symbol("x"), integer(2)),
                          mul({
                              integer(-2),
                              power(symbol("x"), integer(3)),
                          })}),
                     symbol("z")}),
                integer(-2)});

  Expr v = add({power(symbol("x"), integer(4)), integer(-3)});

  Expr u = algebraicExpand(t);

  Expr L = list({symbol("x"), symbol("z")});

  Expr Q = symbol("Q");

  Expr s = polyRemSeq(u, v, L, Q);

  Expr r = add({power(symbol("z"), integer(16)),
                mul({integer(4), power(symbol("z"), integer(15))}),
                mul({integer(14), power(symbol("z"), integer(14))}),
                mul({integer(32), power(symbol("z"), integer(13))}),
                mul({integer(47), power(symbol("z"), integer(12))}),
                mul({integer(92), power(symbol("z"), integer(11))}),
                mul({integer(66), power(symbol("z"), integer(10))}),
                mul({integer(120), power(symbol("z"), integer(9))}),
                mul({integer(-50), power(symbol("z"), integer(8))}),
                mul({integer(-24), power(symbol("z"), integer(7))}),
                mul({integer(-132), power(symbol("z"), integer(6))}),
                mul({integer(-40), power(symbol("z"), integer(5))}),
                mul({integer(-52), power(symbol("z"), integer(4))}),
                mul({integer(-64), power(symbol("z"), integer(3))}),
                mul({integer(-64), power(symbol("z"), integer(2))}),
                mul({integer(-32), symbol("z")}), integer(16)});

  assert(s[0].kind() == Kind::Integer);
  assert(s[0].value() == 1);
  assert(s[1] == r);
}

int main() {
  // TEST(should_get_univariate_resultant)
  // TEST(should_get_multivariate_resultants)
  // TEST(should_get_multivariate_resultants0)
  // TEST(should_get_remainder_sequence)
  // TEST(should_get_remainder_sequence_mv)
  // TEST(should_get_remainder_sequence_mv1)
  TEST(should_get_remainder_sequence_mv2)
  TEST(should_get_remainder_sequence_mv3)
  return 0;
}
