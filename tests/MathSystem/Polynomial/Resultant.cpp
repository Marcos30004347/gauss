#include "MathSystem/Polynomial/Resultant.hpp"
#include "test.hpp"
#include <cstdio>

using namespace alg;
using namespace poly;

// void should_get_multivariate_resultants0() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = pow(x, 3) * pow(y, 2) + 6 * pow(x, 4) * y + 9 * pow(x, 5) +
//            4 * pow(x, 2) * pow(y, 2) + 24 * pow(x, 3) * y +
//            36 * pow(x, 3) + 5 * x * pow(y, 2) + 45 * pow(x, 3) +
//            2 * pow(y, 2) + 12 * y * x + 18 * pow(x, 2);

//   expr v = pow(x, 5) * pow(y, 2) + 8 * pow(x, 4) * y + 16 * pow(x, 3) +
//            12 * pow(x, 4) * pow(y, 2) + 96 * pow(x, 3) * y +
//            192 * pow(x, 2) + 45 * pow(x, 3) * pow(y, 2) +
//            360 * y * pow(x, 2) + 720 * x + 50 * pow(x, 2) * pow(y, 2) +
//            400 * x * y + 800;

//   expr K = symbol("Z");
//   expr L = list({symbol("x"), symbol("y")});

//   expr r0 = polynomialResultant(u, v, L, K);

//   // TODO: currently big int are not being created from strings, when this gets
//   // added, test r == 10734984939700224000*y + 82778463510567321600*(y^2) +
//   // 36933286538080419840*(y^3) + 20609600878213595136*(y^4) +
//   // 12674699737977323520*(y^5) + 4038186495449235456*(y^6) +
//   // 292335413412888576*(y^7) + 133935452101804032*(y^8) +
//   // 55974572889440256*(y^9) + 1212185541869568*(y^10) +
//   // -1422896046391296*(y^11) + 479307895603200*(y^12) + -38800247132160*(y^13)
//   // + -6799984183296*(y^14) + 929350485504*(y^15) + -48693021696*(y^16) +
//   // 2346364800*(y^17) + -45950976*(y^18) + 1105920*(y^19) +
//   // 104853341271490560000
// }

// void should_get_multivariate_resultants() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = pow(x, 3) * pow(y, 3) + 6 * pow(x, 2) * y +
//            5 * x * pow(y, 2) + 2 * pow(y, 2) + y * x + 3 * pow(x, 2);
//   expr v = pow(x, 2) * pow(y, 2) + 5 * pow(x, 3) + 3 * pow(x, 3) * y +
//            4 * y * x + 8;

//   expr K = expr("Z");
//   expr L = list({expr("x"), expr("y")});
//   expr r1 = polynomialResultant(u, v, L, K);

//   assert(r1 == -57024 * y + -133344 * pow(y, 2) + -120904 * pow(y, 3) +
//                    -22656 * pow(y, 4) + 23824 * pow(y, 5) +
//                    49304 * pow(y, 6) + 26796 * pow(y, 7) +
//                    -10328 * pow(y, 8) + -4104 * pow(y, 9) +
//                    1148 * pow(y, 10) + 112 * pow(y, 11) +
//                    -40 * pow(y, 12) + 4 * pow(y, 13) + -8640);
// }

// void should_get_remainder_sequence() {
//   expr u = create(kind::ADD, {pow(symbol("x"), integer(8)), pow(symbol("x"), integer(6)),
//                 create(kind::MUL, {
//                     integer(-3),
//                     pow(symbol("x"), integer(4)),
//                 }),
//                 create(kind::MUL, {
//                     integer(-3),
//                     pow(symbol("x"), integer(3)),
//                 }),
//                 create(kind::MUL, {
//                     integer(8),
//                     pow(symbol("x"), integer(2)),
//                 }),
//                 create(kind::MUL, {integer(2), symbol("x")}), integer(-5)});

//   expr v = create(kind::ADD, {create(kind::MUL, {
//                     integer(3),
//                     pow(symbol("x"), integer(6)),
//                 }),
//                 create(kind::MUL, {
//                     integer(5),
//                     pow(symbol("x"), integer(4)),
//                 }),
//                 create(kind::MUL, {
//                     integer(-4),
//                     pow(symbol("x"), integer(2)),
//                 }),
//                 create(kind::MUL, {integer(-9), symbol("x")}), integer(21)});

//   expr L = list({symbol("x")});
//   expr Q = symbol("Q");

//   // PPRS(u, v, x);
//   expr r = polyRemSeq(u, v, L, Q);

//   assert(r[0].kind() == kind::INT);
//   assert(r[0].value() == 1);

//   assert(r[1].kind() == kind::INT);
//   assert(r[1].value() == 260708);
// }

// void should_get_remainder_sequence_mv() {
//   expr f = create(kind::ADD, {
//       create(kind::MUL, {integer(3), symbol("y"), pow(symbol("x"), integer(2))}),
//       create(kind::MUL, {integer(-1), create(kind::ADD, {pow(symbol("y"), integer(3)), integer(4)})}),
//   });

//   expr g =
//       create(kind::ADD, {pow(symbol("x"), integer(2)),
//            create(kind::MUL, {pow(symbol("y"), integer(3)), symbol("x")}), integer(-9)});

//   expr L = list({symbol("x"), symbol("y")});
//   expr Z = symbol("Z");

//   expr r = polyRemSeq(f, g, L, Z);

//   assert(r[0].kind() == kind::INT);
//   assert(r[0].value() == 1);

//   expr res = create(kind::ADD, {create(kind::MUL, {integer(-3), pow(symbol("y"), integer(10))}),
//                   create(kind::MUL, {integer(-12), pow(symbol("y"), integer(7))}),
//                   pow(symbol("y"), integer(6)),
//                   create(kind::MUL, {integer(-54), pow(symbol("y"), integer(4))}),
//                   create(kind::MUL, {integer(8), pow(symbol("y"), integer(3))}),
//                   create(kind::MUL, {integer(729), pow(symbol("y"), integer(2))}),
//                   create(kind::MUL, {integer(-216), symbol("y")}), integer(16)});

//   assert(r[1] == (res));
// }

// void should_get_remainder_sequence_mv1() {
//   expr u = create(kind::ADD, {pow(symbol("x"), integer(8)), pow(symbol("x"), integer(6)),
//                 create(kind::MUL, {
//                     integer(-3),
//                     pow(symbol("x"), integer(4)),
//                 }),
//                 create(kind::MUL, {
//                     integer(-3),
//                     pow(symbol("x"), integer(3)),
//                 }),
//                 create(kind::MUL, {
//                     integer(8),
//                     pow(symbol("x"), integer(2)),
//                 }),
//                 create(kind::MUL, {integer(2), symbol("x")}), integer(-5)});

//   expr v = create(kind::ADD, {create(kind::MUL, {
//                     integer(3),
//                     pow(symbol("x"), integer(6)),
//                 }),
//                 create(kind::MUL, {
//                     integer(5),
//                     pow(symbol("x"), integer(4)),
//                 }),
//                 create(kind::MUL, {
//                     integer(-4),
//                     pow(symbol("x"), integer(2)),
//                 }),
//                 create(kind::MUL, {integer(-9), symbol("x")}), integer(21)});

//   expr L = list({symbol("x")});

//   expr Z = symbol("Q");

//   expr r = polyRemSeq(u, v, L, Z);

//   assert(r[0].kind() == kind::INT);
//   assert(r[0].value() == 1);

//   assert(r[1].kind() == kind::INT);
//   assert(r[1].value() == 260708);
// }

// void should_get_remainder_sequence_mv2() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = pow(x, 3) * pow(y, 2) + 6 * pow(x, 4) * y + 9 * pow(x, 5) +
//            4 * pow(x, 5) + 4 * pow(x, 2) * pow(y, 2) +
//            24 * pow(x, 3) * y + 36 * pow(x, 4) + 5 * x * pow(y, 2) +
//            30 * y * pow(x, 2) + 45 * pow(x, 3) + 2 * pow(y, 2) +
//            12 * x * y + 18 * pow(x, 2);

//   expr v = pow(x, 5) * pow(y, 2) + 8 * pow(x, 4) * y + 16 * pow(x, 3) +
//            12 * pow(x, 4) * pow(y, 2) + 96 * pow(x, 3) * y +
//            192 * pow(x, 2) + 45 * pow(x, 3) * pow(y, 2) +
//            360 * y * pow(x, 2) + 720 * x + 50 * pow(x, 2) * pow(y, 2) +
//            400 * x * y + 800;

//   expr L = list({x, y});

//   expr Z = expr("Z");

//   expr r = polyRemSeq(u, v, L, Z);

//   assert(r[0] == 1);

//   // TODO: currently big int are not being created from strings, when this gets
//   // added, test r[1] ==  -18572624535748608000*y + 10593940139723980800*(y^2) +
//   // -587913046240788480*(y^3) + -2244223197765435392*(y^4) +
//   // 1015993222301745152*(y^5) + -20912999483047936*(y^6) +
//   // -108919567828385792*(y^7) + 28961535612157952*(y^8) +
//   // 1249933094027264*(y^9) + -1831940388552704*(y^10) + 282281363111936*(y^11)
//   // + 10696300331008*(y^12) + -8584356855808*(y^13) + 1196684484608*(y^14) +
//   // -86404104192*(y^15) + 3759341568*(y^16) + -94371840*(y^17) + 1179648*(y^18)
//   // + 12314263137812480000
// }

// void should_get_remainder_sequence_mv3() {
//   expr t = create(kind::ADD, {pow(symbol("z"), integer(4)), pow(symbol("z"), integer(3)),
//                 create(kind::MUL, {create(kind::ADD, {integer(2), symbol("x"),
//                           create(kind::MUL, {
//                               integer(-1),
//                               pow(symbol("x"), integer(2)),
//                           })}),
//                      pow(symbol("z"), integer(2))}),
//                 create(kind::MUL, {create(kind::ADD, {integer(1), pow(symbol("x"), integer(2)),
//                           create(kind::MUL, {
//                               integer(-2),
//                               pow(symbol("x"), integer(3)),
//                           })}),
//                      symbol("z")}),
//                 integer(-2)});

//   expr v = create(kind::ADD, {pow(symbol("x"), integer(4)), integer(-3)});
//   // printf("%s\n", expand(t).toString().c_str());

//   expr u = expand(t);

//   expr L = list({symbol("x"), symbol("z")});

//   expr Q = symbol("Q");

//   expr s = polyRemSeq(u, v, L, Q);
//   // printf("-----> %s\n", s.toString().c_str());
//   expr r = create(kind::ADD, {pow(symbol("z"), integer(16)),
//                 create(kind::MUL, {integer(4), pow(symbol("z"), integer(15))}),
//                 create(kind::MUL, {integer(14), pow(symbol("z"), integer(14))}),
//                 create(kind::MUL, {integer(32), pow(symbol("z"), integer(13))}),
//                 create(kind::MUL, {integer(47), pow(symbol("z"), integer(12))}),
//                 create(kind::MUL, {integer(92), pow(symbol("z"), integer(11))}),
//                 create(kind::MUL, {integer(66), pow(symbol("z"), integer(10))}),
//                 create(kind::MUL, {integer(120), pow(symbol("z"), integer(9))}),
//                 create(kind::MUL, {integer(-50), pow(symbol("z"), integer(8))}),
//                 create(kind::MUL, {integer(-24), pow(symbol("z"), integer(7))}),
//                 create(kind::MUL, {integer(-132), pow(symbol("z"), integer(6))}),
//                 create(kind::MUL, {integer(-40), pow(symbol("z"), integer(5))}),
//                 create(kind::MUL, {integer(-52), pow(symbol("z"), integer(4))}),
//                 create(kind::MUL, {integer(-64), pow(symbol("z"), integer(3))}),
//                 create(kind::MUL, {integer(-64), pow(symbol("z"), integer(2))}),
//                 create(kind::MUL, {integer(-32), symbol("z")}), integer(16)});

//   assert(s[0].kind() == kind::INT);
//   assert(s[0].value() == 1);
//   assert(s[1] == r);
// }

void should_get_remainder_sequence_mv2_poly_exp() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x, y});

  expr u = polyExpr(
      pow(x, 3) * pow(y, 2) + 6 * pow(x, 4) * y + 13 * pow(x, 5) +
          4 * pow(x, 2) * pow(y, 2) + 24 * pow(x, 3) * y +
          36 * pow(x, 4) + 5 * x * pow(y, 2) + 30 * y * pow(x, 2) +
          45 * pow(x, 3) + 2 * pow(y, 2) + 12 * x * y + 18 * pow(x, 2),
      L);

  expr v =
      polyExpr(pow(x, 5) * pow(y, 2) + 8 * pow(x, 4) * y +
                   16 * pow(x, 3) + 12 * pow(x, 4) * pow(y, 2) +
                   96 * pow(x, 3) * y + 192 * pow(x, 2) +
                   45 * pow(x, 3) * pow(y, 2) + 360 * y * pow(x, 2) +
                   720 * x + 50 * pow(x, 2) * pow(y, 2) + 400 * x * y + 800,
               L);

  expr Z = expr("Z");

  expr r = remSeqPolyExpr(u, v, L, Z);

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
  expr x = expr("x");
  expr z = expr("z");

  expr L = list({x, z});

  expr u = polyExpr(pow(z, 2) * x + z * pow(x, 2) +
                        -1 * pow(z, 2) * pow(x, 2) + -2 * z * pow(x, 3) +
                        z + 2 * pow(z, 2) + pow(z, 3) + pow(z, 4) + -2,
                    L);

  expr v = polyExpr(pow(x, 4) + -3, L);

  expr Q = expr("Q");

  expr s = remSeqPolyExpr(u, v, L, Q);

  assert(s == list({
                  create(kind::ADD, {create(kind::ADD, {create(kind::MUL, {1, pow(z, 16)}), create(kind::MUL, {4, pow(z, 15)}),
                            create(kind::MUL, {14, pow(z, 14)}), create(kind::MUL, {32, pow(z, 13)}),
                            create(kind::MUL, {47, pow(z, 12)}), create(kind::MUL, {92, pow(z, 11)}),
                            create(kind::MUL, {66, pow(z, 10)}), create(kind::MUL, {120, pow(z, 9)}),
                            create(kind::MUL, {-50, pow(z, 8)}), create(kind::MUL, {-24, pow(z, 7)}),
                            create(kind::MUL, {-132, pow(z, 6)}), create(kind::MUL, {-40, pow(z, 5)}),
                            create(kind::MUL, {-52, pow(z, 4)}), create(kind::MUL, {-64, pow(z, 3)}),
                            create(kind::MUL, {-64, pow(z, 2)}), create(kind::MUL, {-32, pow(z, 1)}),
                            create(kind::MUL, {16, pow(z, 0)})}) *
                       pow(x, 0)}),
                  create(kind::ADD, {create(kind::ADD, {1 * pow(z, 0)}) * pow(x, 0)}),

              }));
}

void should_get_remainder_sequence_mv_poly_exp() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x, y});

  expr f = polyExpr(3 * y * pow(x, 2) + -1 * pow(y, 3) + -4, L);

  expr g = polyExpr(pow(x, 2) + pow(y, 3) * x + -9, L);

  expr Z = expr("Z");

  expr r = remSeqPolyExpr(f, g, L, Z);
  assert(r == list({
                  create(kind::ADD, {create(kind::ADD, {
                           16 * pow(y, 0),
                           -216 * pow(y, 1),
                           729 * pow(y, 2),
                           8 * pow(y, 3),
                           -54 * pow(y, 4),
                           1 * pow(y, 6),
                           -12 * pow(y, 7),
                           -3 * pow(y, 10),
                       }) *
                       pow(x, 0)}),
                  create(kind::ADD, {create(kind::ADD, {1 * pow(y, 0)}) * pow(x, 0)}),
              }));
}

void should_get_remainder_sequence_mv1_poly_exp() {
  expr x = expr("x");

  expr L = list({x});

  expr u = polyExpr(2 * x + 8 * pow(x, 2) + -3 * pow(x, 3) +
                        -3 * pow(x, 4) + pow(x, 6) + pow(x, 8) + -5,
                    L);

  expr v = polyExpr(
      -9 * x + -4 * pow(x, 2) + 5 * pow(x, 4) + 3 * pow(x, 6) + 21, L);

  expr Z = expr("Q");

  expr r = remSeqPolyExpr(u, v, L, Z);

  assert(r == list({
                  create(kind::ADD, {260708 * pow(x, 0)}),
                  create(kind::ADD, {1 * pow(x, 0)}),
              }));
}

void should_get_multivariate_resultants0_poly_exp() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x, y});

  expr u = polyExpr(pow(x, 3) * pow(y, 2) + 6 * pow(x, 4) * y +
                        9 * pow(x, 5) + 4 * pow(x, 2) * pow(y, 2) +
                        24 * pow(x, 3) * y + 5 * x * pow(y, 2) +
                        81 * pow(x, 3) + 2 * pow(y, 2) + 12 * y * x +
                        18 * pow(x, 2),
                    L);

  expr v =
      polyExpr(pow(x, 5) * pow(y, 2) + 8 * pow(x, 4) * y +
                   16 * pow(x, 3) + 12 * pow(x, 4) * pow(y, 2) +
                   96 * pow(x, 3) * y + 192 * pow(x, 2) +
                   45 * pow(x, 3) * pow(y, 2) + 360 * y * pow(x, 2) +
                   720 * x + 50 * pow(x, 2) * pow(y, 2) + 400 * x * y + 800,
               L);

  expr K = expr("Z");

  expr r0 = resultantPolyExpr(u, v, L, K);
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

void should_get_multivariate_resultants_poly_exp() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x, y});

  expr u = polyExpr(pow(x, 3) * pow(y, 3) + 6 * pow(x, 2) * y +
                        5 * x * pow(y, 2) + 2 * pow(y, 2) + y * x +
                        3 * pow(x, 2),
                    L);
  expr v = polyExpr(pow(x, 2) * pow(y, 2) + 5 * pow(x, 3) +
                        3 * pow(x, 3) * y + 4 * y * x + 8,
                    L);

  expr K = expr("Z");

	expr r1 = resultantPolyExpr(u, v, L, K);

  assert(r1 == create(kind::ADD, {create(kind::ADD, {-8640 * pow(y, 0), -57024 * pow(y, 1),
                         -133344 * pow(y, 2), -120904 * pow(y, 3),
                         -22656 * pow(y, 4), 23824 * pow(y, 5),
                         49304 * pow(y, 6), 26796 * pow(y, 7),
                         -10328 * pow(y, 8), -4104 * pow(y, 9),
                         1148 * pow(y, 10), 112 * pow(y, 11),
                         -40 * pow(y, 12), 4 * pow(y, 13)}) *
                    pow(x, 0)}));
}

int main() {
  // TEST(should_get_multivariate_resultants0)
  // TEST(should_get_remainder_sequence)
  // TEST(should_get_remainder_sequence_mv)
  // TEST(should_get_remainder_sequence_mv1)
  // TEST(should_get_remainder_sequence_mv2)
  // TEST(should_get_remainder_sequence_mv3)
  // TEST(should_get_multivariate_resultants)

  TEST(should_get_multivariate_resultants_poly_exp)
  TEST(should_get_multivariate_resultants0_poly_exp)
  TEST(should_get_remainder_sequence_mv_poly_exp)
  TEST(should_get_remainder_sequence_mv1_poly_exp)
  TEST(should_get_remainder_sequence_mv2_poly_exp)
  TEST(should_get_remainder_sequence_mv3_poly_exp)
  return 0;
}
