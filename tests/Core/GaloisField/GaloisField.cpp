#include "Core/GaloisField/GaloisField.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "test.hpp"
#include <cstdlib>

using namespace ast;
using namespace algebra;
using namespace galoisField;

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

void should_sub_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 4 * power(x, 4) + 3 * power(x, 3) + 5 * power(x, 2) + 2 * x;

  assert(subPolyGf(u, v, x, 5) == x + power(x, 2) + power(x, 4));
}

void should_mul_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 4 * power(x, 4) + 3 * power(x, 3) + 5 * power(x, 2) + 2 * x;

  assert(mulPolyGf(u, v, x, 5) ==
         power(x, 2) + 2 * power(x, 3) + -2 * power(x, 6) + 2 * power(x, 7));
}

void should_div_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 5 * power(x, 2) + 2 * x;

  assert(divPolyGf(u, v, x, 7) == list({-1 * x + 3 * power(x, 2) + 3, 2 * x}));
}

void should_pow_mod_poly_galois_field() {
  Expr x = Expr("x");

  Expr u = 15 * power(x, 4) + 8 * power(x, 3) + 6 * power(x, 2) + 8 * x;
  Expr v = 5 * power(x, 2) + 2 * x;

  assert(powModPolyGf(u, v, x, 3, 7) == x);
}

void should_get_monic_form_galois_field() {
  Expr x = Expr("x");

  Expr u = 3 * power(x, 4) + 3 * power(x, 3) + -2 * power(x, 2) + x;

  assert(monicPolyGf(u, x, 7) ==
         list({3, -2 * x + -3 * power(x, 2) + power(x, 3) + power(x, 4)}));
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
  TEST(should_ground_project_poly_galois_field)
  TEST(should_add_poly_galois_field)
  TEST(should_sub_poly_galois_field)
  TEST(should_mul_poly_galois_field)
  TEST(should_div_poly_galois_field)
  TEST(should_pow_mod_poly_galois_field)
  TEST(should_get_monic_form_galois_field)
  TEST(should_get_gcd_of_poly_galois_field)
  TEST(should_perform_extended_euclid_galois_field)

  return 0;
}
