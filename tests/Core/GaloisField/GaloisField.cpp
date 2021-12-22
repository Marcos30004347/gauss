#include "Core/GaloisField/GaloisField.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "test.hpp"
#include <cstdlib>

using namespace ast;
using namespace algebra;
using namespace galoisField;
using namespace polynomial;

void should_project_poly_galois_field() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;

  assert(gf(u, 5, false) == 3 * power(x, 3) + power(x, 2) + 3 * x);
  assert(gf(u, 5, true) == -2 * power(x, 3) + power(x, 2) + -2 * x);

  Expr g = 15 * power(x, 4) * power(y, 3) + 8 * power(x, 3) +
           6 * power(x, 2) * power(y, 2) + 8 * x * y + 17 * y;

  assert(gf(g, 6, false) ==
         3 * power(x, 4) * power(y, 3) + 2 * power(x, 3) + 2 * x * y + 5 * y);
  assert(gf(g, 6, true) ==
         3 * power(x, 4) * power(y, 3) + 2 * power(x, 3) + 2 * x * y + -1 * y);
}

void should_project_poly_expr_galois_field() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = polyExpr(
      15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x, list({x}));

  assert(gfPolyExpr(u, 5, false) ==
         polyExpr(3 * power(x, 3) + power(x, 2) + 3 * x, list({x})));
  assert(gfPolyExpr(u, 5, true) ==
         polyExpr(-2 * power(x, 3) + power(x, 2) + -2 * x, list({x})));

  Expr g = polyExpr(15 * power(x, 4) * power(y, 3) + 8 * power(x, 3) +
                        6 * power(x, 2) * power(y, 2) + 8 * x * y + 17 * y,
                    list({x, y}));

  assert(gfPolyExpr(g, 6, false) ==
         polyExpr(3 * power(x, 4) * power(y, 3) + 2 * power(x, 3) + 2 * x * y +
                      5 * y,
                  list({x, y})));
  assert(gfPolyExpr(g, 6, true) ==
         polyExpr(3 * power(x, 4) * power(y, 3) + 2 * power(x, 3) + 2 * x * y +
                      -1 * y,
                  list({x, y})));
}

void should_ground_project_poly_galois_field() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = (15 * power(x, 4) + 8 * power(x, 3)) * (6 * power(x, 2) + 8 * x);

  assert(groundGf(u, 5, false) == (3 * power(x, 3)) * (power(x, 2) + 3 * x));
  assert(groundGf(u, 5, true) == (-2 * power(x, 3)) * (power(x, 2) + -2 * x));

  Expr g = (15 * power(x, 4) * power(y, 3) + 8 * power(x, 3)) *
           (6 * power(x, 2) * power(y, 2) + 8 * x * y + 17 * y);

  assert(groundGf(g, 6, false) ==
         (3 * power(x, 4) * power(y, 3) + 2 * power(x, 3)) *
             (2 * x * y + 5 * y));
  assert(groundGf(g, 6, true) ==
         (3 * power(x, 4) * power(y, 3) + 2 * power(x, 3)) *
             (2 * x * y + -1 * y));
}

void should_add_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 4 * power(x, 4) + 3 * power(x, 3) + 5 * power(x, 2) + 2 * x;

  assert(addPolyGf(u, v, x, 5) == power(x, 2) + power(x, 3) + -1 * power(x, 4));
}

void should_add_poly_expr_galois_field() {
  Expr x = Expr("x");

  Expr u = polyExpr(
      15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x, list({x}));
  Expr v = polyExpr(4 * power(x, 4) + 3 * power(x, 3) + 5 * power(x, 2) + 2 * x,
                    list({x}));

  assert(addPolyExprGf(u, v, 5) ==
         polyExpr(power(x, 2) + power(x, 3) + -1 * power(x, 4), list({x})));
}

void should_sub_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 4 * power(x, 4) + 3 * power(x, 3) + 5 * power(x, 2) + 2 * x;

  assert(subPolyGf(u, v, x, 5) == x + power(x, 2) + power(x, 4));
}

void should_sub_poly_expr_galois_field() {
  Expr x = Expr("x");

  Expr u = polyExpr(
      15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x, list({x}));
  Expr v = polyExpr(4 * power(x, 4) + 3 * power(x, 3) + 5 * power(x, 2) + 2 * x,
                    list({x}));

  assert(subPolyExprGf(u, v, 5) ==
         polyExpr(x + power(x, 2) + power(x, 4), list({x})));
}

void should_mul_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 4 * power(x, 4) + 3 * power(x, 3) + 5 * power(x, 2) + 2 * x;

  assert(mulPolyGf(u, v, x, 5) ==
         power(x, 2) + 2 * power(x, 3) + -2 * power(x, 6) + 2 * power(x, 7));
}

void should_mul_poly_expr_galois_field() {
  Expr x = Expr("x");

  Expr u = polyExpr(
      15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x, list({x}));
  Expr v = polyExpr(4 * power(x, 4) + 3 * power(x, 3) + 5 * power(x, 2) + 2 * x,
                    list({x}));

  assert(mulPolyExprGf(u, v, 5) ==
         polyExpr(power(x, 2) + 2 * power(x, 3) + -2 * power(x, 6) +
                      2 * power(x, 7),
                  list({x})));
}

void should_div_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 5 * power(x, 2) + 2 * x;

  assert(divPolyGf(u, v, x, 7) == list({-1 * x + 3 * power(x, 2) + 3, 2 * x}));
}

void should_div_poly_expr_galois_field() {
  Expr x = Expr("x");

  Expr u = polyExpr(
      15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x, list({x}));
  Expr v = polyExpr(5 * power(x, 2) + 2 * x, list({x}));

  assert(divPolyExprGf(u, v, x, 7) ==
         list({polyExpr(-1 * x + 3 * power(x, 2) + 3, list({x})),
               polyExpr(2 * x, list({x}))}));

  Expr g = polyExpr(power(x, 2) + -2 * power(x, 3) + 3 * power(x, 4) +
                        -1 * power(x, 6) + 2 * power(x, 7) + power(x, 8),
                    list({x}));
  Expr k = polyExpr(2 * x + 5 * power(x, 2), list({x}));

  assert(
      divPolyExprGf(g, k, x, 7) ==
      list({polyExpr(-2 + 2 * x + power(x, 2) + -1 * power(x, 3) +
                         -1 * power(x, 4) + 2 * power(x, 5) + 3 * power(x, 6),
                     list({x})),
            polyExpr(-3 * power(x, 1), list({x}))}));
}

void should_pow_mod_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 5 * power(x, 2) + 2 * x;

  assert(powModPolyGf(u, v, x, 3, 7) == x);
}

void should_pow_mod_poly_expr_galois_field() {
  Expr x = Expr("x");

  Expr u = polyExpr(
      15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x, list({x}));
  Expr v = polyExpr(5 * power(x, 2) + 2 * x, list({x}));

  assert(powModPolyExprGf(u, v, x, 3, 7) == polyExpr(x, list({x})));
}

void should_get_monic_form_galois_field() {
  Expr x = Expr("x");

  Expr u = 3 * power(x, 4) + 3 * power(x, 3) + -2 * power(x, 2) + x;

  assert(monicPolyGf(u, x, 7) ==
         list({3, -2 * x + -3 * power(x, 2) + power(x, 3) + power(x, 4)}));
}

void should_get_monic_form_poly_expr_galois_field() {
  Expr x = Expr("x");

  Expr u = polyExpr(3 * power(x, 4) + 3 * power(x, 3) + -2 * power(x, 2) + x,
                    list({x}));

  assert(monicPolyExprGf(u, x, 7) ==
         list({polyExpr(3, list({x})),
               polyExpr(-2 * x + -3 * power(x, 2) + power(x, 3) + power(x, 4),
                        list({x}))}));
}

void should_get_gcd_of_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 5 * power(x, 2) + 2 * x;

  Expr gcd = gcdPolyGf(u, v, x, 7);

  assert(gcd == x);
}

void should_perform_extended_euclid_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 5 * power(x, 2) + 2 * x;

  Expr e = extendedEuclidGf(u, v, x, 7);

  assert(e == list({x, -3, -3 * x + 2 * power(x, 2) + 2}));
}

int main() {

  TEST(should_project_poly_galois_field)
  TEST(should_project_poly_expr_galois_field)
  TEST(should_ground_project_poly_galois_field)
  TEST(should_add_poly_galois_field)
  TEST(should_add_poly_expr_galois_field)
  TEST(should_sub_poly_galois_field)
  TEST(should_sub_poly_expr_galois_field)
  TEST(should_mul_poly_galois_field)
  TEST(should_mul_poly_expr_galois_field)
  TEST(should_div_poly_galois_field)
  TEST(should_div_poly_expr_galois_field)
  TEST(should_pow_mod_poly_galois_field)
  TEST(should_pow_mod_poly_expr_galois_field)
  TEST(should_get_monic_form_galois_field)
  TEST(should_get_monic_form_poly_expr_galois_field)
  TEST(should_get_gcd_of_poly_galois_field)
  TEST(should_perform_extended_euclid_galois_field)

  return 0;
}
