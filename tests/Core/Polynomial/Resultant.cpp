#include "Core/Polynomial/Resultant.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Polynomial/Polynomial.hpp"
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

  Expr r0 = polynomialResultant(ux, vx, list({x}), Expr("Z"));

  assert(r0.kind() == Kind::Integer);
  assert(r0.value() == 218);
}

void should_get_multivariate_resultants0() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = power(x, 3) * power(y, 2) + 6 * power(x, 4) * y + 9 * power(x, 5) +
           4 * power(x, 2) * power(y, 2) + 24 * power(x, 3) * y +
           36 * power(x, 3) + 5 * x * power(y, 2) + 45 * power(x, 3) +
           2 * power(y, 2) + 12 * y * x + 18 * power(x, 2);

  Expr v = power(x, 5) * power(y, 2) + 8 * power(x, 4) * y + 16 * power(x, 3) +
           12 * power(x, 4) * power(y, 2) + 96 * power(x, 3) * y +
           192 * power(x, 2) + 45 * power(x, 3) * power(y, 2) +
           360 * y * power(x, 2) + 720 * x + 50 * power(x, 2) * power(y, 2) +
           400 * x * y + 800;

  Expr K = symbol("Z");
  Expr L = list({symbol("x"), symbol("y")});

  Expr r0 = polynomialResultant(u, v, L, K);

  // TODO: currently big int are not being created from strings, when this gets
  // added, test r == 10734984939700224000*y + 82778463510567321600*(y^2) +
  // 36933286538080419840*(y^3) + 20609600878213595136*(y^4) +
  // 12674699737977323520*(y^5) + 4038186495449235456*(y^6) +
  // 292335413412888576*(y^7) + 133935452101804032*(y^8) +
  // 55974572889440256*(y^9) + 1212185541869568*(y^10) +
  // -1422896046391296*(y^11) + 479307895603200*(y^12) + -38800247132160*(y^13)
  // + -6799984183296*(y^14) + 929350485504*(y^15) + -48693021696*(y^16) +
  // 2346364800*(y^17) + -45950976*(y^18) + 1105920*(y^19) +
  // 104853341271490560000
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
  Expr r1 = polynomialResultant(u, v, L, K);

  assert(r1 == -57024 * y + -133344 * power(y, 2) + -120904 * power(y, 3) +
                   -22656 * power(y, 4) + 23824 * power(y, 5) +
                   49304 * power(y, 6) + 26796 * power(y, 7) +
                   -10328 * power(y, 8) + -4104 * power(y, 9) +
                   1148 * power(y, 10) + 112 * power(y, 11) +
                   -40 * power(y, 12) + 4 * power(y, 13) + -8640);
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

  Expr v = power(x, 5) * power(y, 2) + 8 * power(x, 4) * y + 16 * power(x, 3) +
           12 * power(x, 4) * power(y, 2) + 96 * power(x, 3) * y +
           192 * power(x, 2) + 45 * power(x, 3) * power(y, 2) +
           360 * y * power(x, 2) + 720 * x + 50 * power(x, 2) * power(y, 2) +
           400 * x * y + 800;

  Expr L = list({x, y});

  Expr Z = Expr("Z");

  Expr r = polyRemSeq(u, v, L, Z);

  assert(r[0] == 1);

  // TODO: currently big int are not being created from strings, when this gets
  // added, test r[1] ==  -18572624535748608000*y + 10593940139723980800*(y^2) +
  // -587913046240788480*(y^3) + -2244223197765435392*(y^4) +
  // 1015993222301745152*(y^5) + -20912999483047936*(y^6) +
  // -108919567828385792*(y^7) + 28961535612157952*(y^8) +
  // 1249933094027264*(y^9) + -1831940388552704*(y^10) + 282281363111936*(y^11)
  // + 10696300331008*(y^12) + -8584356855808*(y^13) + 1196684484608*(y^14) +
  // -86404104192*(y^15) + 3759341568*(y^16) + -94371840*(y^17) + 1179648*(y^18)
  // + 12314263137812480000
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
  printf("%s\n", algebraicExpand(t).toString().c_str());

  Expr u = algebraicExpand(t);

  Expr L = list({symbol("x"), symbol("z")});

  Expr Q = symbol("Q");

  Expr s = polyRemSeq(u, v, L, Q);
  printf("-----> %s\n", s.toString().c_str());
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

void should_get_remainder_sequence_mv2_poly_exp() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr L = list({x, y});

  Expr u = polyExpr(
      power(x, 3) * power(y, 2) + 6 * power(x, 4) * y + 13 * power(x, 5) +
          4 * power(x, 2) * power(y, 2) + 24 * power(x, 3) * y +
          36 * power(x, 4) + 5 * x * power(y, 2) + 30 * y * power(x, 2) +
          45 * power(x, 3) + 2 * power(y, 2) + 12 * x * y + 18 * power(x, 2),
      L);

  Expr v =
      polyExpr(power(x, 5) * power(y, 2) + 8 * power(x, 4) * y +
                   16 * power(x, 3) + 12 * power(x, 4) * power(y, 2) +
                   96 * power(x, 3) * y + 192 * power(x, 2) +
                   45 * power(x, 3) * power(y, 2) + 360 * y * power(x, 2) +
                   720 * x + 50 * power(x, 2) * power(y, 2) + 400 * x * y + 800,
               L);

  Expr Z = Expr("Z");

  Expr r = colPolyRemSeq(u, v, L, Z);

  // TODO: currently big int are not being created from strings, when this gets
  // added, test r[1] ==  -18572624535748608000*y + 10593940139723980800*(y^2) +
  // -587913046240788480*(y^3) + -2244223197765435392*(y^4) +
  // 1015993222301745152*(y^5) + -20912999483047936*(y^6) +
  // -108919567828385792*(y^7) + 28961535612157952*(y^8) +
  // 1249933094027264*(y^9) + -1831940388552704*(y^10) + 282281363111936*(y^11)
  // + 10696300331008*(y^12) + -8584356855808*(y^13) + 1196684484608*(y^14) +
  // -86404104192*(y^15) + 3759341568*(y^16) + -94371840*(y^17) + 1179648*(y^18)
  // + 12314263137812480000
}

void should_get_remainder_sequence_mv3_poly_exp() {
  Expr x = Expr("x");
  Expr z = Expr("z");

  Expr L = list({x, z});

  Expr u = polyExpr(power(z, 2) * x + z * power(x, 2) +
                        -1 * power(z, 2) * power(x, 2) + -2 * z * power(x, 3) +
                        z + 2 * power(z, 2) + power(z, 3) + power(z, 4) + -2,
                    L);

  Expr v = polyExpr(power(x, 4) + -3, L);

  Expr Q = Expr("Q");

  Expr s = colPolyRemSeq(u, v, L, Q);

  assert(s == list({add({add({1 * power(z, 0)}) * power(x, 0)}),
                    add({add({mul({1, power(z, 16)}), mul({4, power(z, 15)}),
                              mul({14, power(z, 14)}), mul({32, power(z, 13)}),
                              mul({47, power(z, 12)}), mul({92, power(z, 11)}),
                              mul({66, power(z, 10)}), mul({120, power(z, 9)}),
                              mul({-50, power(z, 8)}), mul({-24, power(z, 7)}),
                              mul({-132, power(z, 6)}), mul({-40, power(z, 5)}),
                              mul({-52, power(z, 4)}), mul({-64, power(z, 3)}),
                              mul({-64, power(z, 2)}), mul({-32, power(z, 1)}),
                              mul({16, power(z, 0)})}) *
                         power(x, 0)})}));
}



void should_get_remainder_sequence_mv_poly_exp() {
	Expr x = Expr("x");
	Expr y = Expr("y");

  Expr L = list({ x, y });

	Expr f = polyExpr(3*y*power(x, 2) + -1*power(y, 3) + -4, L);

  Expr g = polyExpr(power(x, 2) + power(y, 3)*x + -9, L);

  Expr Z = Expr("Z");

  Expr r = colPolyRemSeq(f, g, L, Z);

	assert(r == list({
				add({ add({ 1*power(y, 0) })* power(x, 0) }),
				add({add({
								16*power(y, 0),
								-216*power(y, 1),
								729*power(y, 2),
								8*power(y, 3),
								-54*power(y, 4),
								1*power(y, 6),
								-12*power(y, 7),
								-3*power(y, 10),
				})*power(x, 0)})
	}));
}

void should_get_remainder_sequence_mv1_poly_exp() {
	Expr x = Expr("x");

	Expr L = list({ x });

	Expr u = polyExpr(2*x + 8*power(x, 2) + -3*power(x, 3) + -3*power(x, 4) + power(x, 6) + power(x, 8) + -5, L);


  Expr v = polyExpr(-9*x + -4*power(x, 2) + 5*power(x, 4) + 3*power(x, 6) + 21, L);

  Expr Z = Expr("Q");

	Expr r = colPolyRemSeq(u, v, L, Z);

	assert(r == list({
				add({1*power(x, 0)}),
				add({260708*power(x, 0)}),
			}));
}



int main() {
  // TEST(should_get_univariate_resultant)
  // TEST(should_get_multivariate_resultants)
  // TEST(should_get_multivariate_resultants0)
  // TEST(should_get_remainder_sequence)
  // TEST(should_get_remainder_sequence_mv)
  // TEST(should_get_remainder_sequence_mv1)
  // TEST(should_get_remainder_sequence_mv2)
  // TEST(should_get_remainder_sequence_mv3)

	TEST(should_get_remainder_sequence_mv_poly_exp)
  TEST(should_get_remainder_sequence_mv1_poly_exp)
  TEST(should_get_remainder_sequence_mv2_poly_exp)
  TEST(should_get_remainder_sequence_mv3_poly_exp)
  return 0;
}
