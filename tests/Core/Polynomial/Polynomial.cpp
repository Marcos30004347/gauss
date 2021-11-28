#include "Core/Polynomial/Polynomial.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Expand/Expand.hpp"

#include "test.hpp"

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace polynomial;

void should_get_polynomial_variable() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr exp0 = 4 * x + power(x, 2) + 5 * power(x, 3);
  Expr exp1 = 4 * x + power(y, 2) + 5 * sin(x);

  Expr vars_exp0 = variables(exp0);
  Expr vars_exp1 = variables(exp1);

  assert(variables(exp0) == list({x}));
  assert(variables(exp1) == list({y, x, sin(x)}));
}

void should_get_degree_of_variables() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr exp0 = 4 * x + power(x, 2) + 5 * power(x, 3);
  Expr exp1 = 4 * x + power(y, 2) + power(sin(x), 5);

  assert(degree(exp0, x) == 3);
  assert(degree(exp1, x) == 1);
  assert(degree(exp1, y) == 2);
  assert(degree(exp1, sin(x)) == 5);
}

void should_get_coefficients() {
  Expr x = Expr("x");
  Expr a = Expr("a");
  Expr b = Expr("b");

  Expr exp0 = 4 * power(x, 2);
  Expr exp1 = a * power(x, 2) + b * power(x, 2);
  Expr exp2 = a * power(x, 2) - b * power(x, 2);

  assert(coeff(exp0, x, 2) == 4);
  assert(coeff(exp1, x, 2) == a + b);
  assert(coeff(exp2, x, 2) == a - b);
}

void should_get_leading_coefficient() {
  Expr x = Expr("x");
  Expr a = Expr("a");
  Expr b = Expr("b");

  Expr exp0 = 4 * power(x, 2);
  Expr exp1 = a * power(x, 2) + b * power(x, 2);
  Expr exp2 = a * power(x, 2) + b * power(x, 3);

  Expr leadcoeff_exp0 = leadCoeff(exp0, x);
  Expr leadcoeff_exp1 = leadCoeff(exp1, x);
  Expr leadcoeff_exp2 = leadCoeff(exp2, x);

  assert(leadCoeff(exp0, x) == 4);
  assert(leadCoeff(exp1, x) == a + b);
  assert(leadCoeff(exp2, x) == b);
}

void should_divided_polynomials() {
  Expr x = Expr("x");
  Expr exp0 = 5 * power(x, 2) + 4 * x + 1;
  Expr exp1 = 2 * x + 3;

  Expr res = divideGPE(exp0, exp1, x);
  assert(res == list({fraction(-7, 4) + fraction(5, 2) * x, fraction(25, 4)}));
}

void should_get_gcd_polynomials() {
  Expr x = Expr("x");
  Expr u = power(x, 7) + -4 * power(x, 5) + -1 * power(x, 2) + 4;
  Expr v = power(x, 5) + -4 * power(x, 3) + -1 * power(x, 2) + 4;

  Expr res = gcdGPE(u, v, x);

  assert(gcdGPE(u, v, x) == 4 + -4 * x + -1 * power(x, 2) + power(x, 3));
}

void should_get_extended_gcd_polynomials() {
  Expr x = Expr("x");
  Expr u = power(x, 7) + -4 * power(x, 5) + -1 * power(x, 2) + 4;
  Expr v = power(x, 5) + -4 * power(x, 3) + -1 * power(x, 2) + 4;

  Expr res = extendedEuclideanAlgGPE(u, v, x);

  assert(res[1] == -1 * x);
  assert(res[2] == 1 + power(x, 3));
  assert(expandAST(res[1] * u + res[2] * v) == res[0]);
}

void should_calculate_monomial_division() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = power(x, 3) + 3 * power(x, 2) * y + 4 * x * power(y, 2);
  Expr v = x * y + 2 * y + 3 * power(y, 2);

  Expr L = list({Expr("x"), Expr("y")});

  assert(monomialPolyDiv(u, v, L) ==
         list({
             -6 + 3 * x + -5 * y,
             power(x, 3) + 12 * y + 28 * power(y, 2) + 15 * power(y, 3),
         }));
}

void should_get_leading_monomial() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 3 * power(x, 2) * y + 4 * x * power(y, 2) + power(y, 3) + x + 3;

  Expr L0 = list({Expr("x"), Expr("y")});
  Expr L1 = list({Expr("y"), Expr("x")});

  assert(leadMonomial(u, L0) == 3 * power(x, 2) * y);
  assert(leadMonomial(u, L1) == power(y, 3));
}

void should_rec_divide_polynomials() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = power(x, 2) * power(y, 2) + x;
  Expr v = x * y + 1;

  Expr L0 = list({x, y});
  Expr L1 = list({y, x});

  Expr Q = Expr("Q");

  Expr R0 = recPolyDiv(u, v, L0, Q);
  Expr R1 = recPolyDiv(u, v, L1, Q);

  assert(R0[0] == x * y);
  assert(R0[1] == x + -1 * x * y);
  assert(R1[0] == -1 + x * y);
  assert(R1[1] == 1 + x);
}

void should_pseudo_divide_polynomials() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 5 * power(x, 4) * power(y, 3) + 3 * x * y + 2;
  Expr v = 2 * power(x, 3) * y + 2 * x + 3;

  Expr R = pseudoDivision(u, v, x);

  assert(R.kind() == Kind::List);
  assert(R.size() == 2);

  assert(R[0] == 10 * x * power(y, 4));
  assert(R[1] == 8 * power(y, 2) + 12 * x * power(y, 3) +
                     -30 * x * power(y, 4) + -20 * power(x, 2) * power(y, 4));
}

void should_normalize_polynomial() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = (2 * y + 1) * x + 6 * y + 3;

  Expr L = list({x, y});
  Expr Q = Expr("Q");

  Expr u_ = normalizePoly(u, L, Q);
  Expr u_res = fraction(3, 2) + fraction(1, 2) * x + 3 * y + x * y;

  assert(u_ == u_res);
}

void should_mv_poly_gcd() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = -1 * y * power(x, 2) + power(y, 3);
  Expr v = y * power(x, 2) + 2 * power(y, 2) * x + power(y, 3);

  Expr L = list({x, y});
  Expr Z = Expr("Z");

  Expr gcd = mvPolyGCD(u, v, L, Z);

  assert(gcd == x * y + power(y, 2));
}

void should_get_coeff_var_parts_of_monomial() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = Expr(4) * Expr(5) * fraction(1, 2) * x * power(x, 2) * power(y, 3);

  Expr S = list({x, y});

  Expr L = coeffVarMonomial(u, S);

  assert(L.kind() == Kind::List);
  assert(L.size() == 2);

  assert(L[0] == Expr(4) * Expr(5) * fraction(1, 2));
  assert(L[1] == x * power(x, 2) * power(y, 3));
}

void should_collect_terms() {
  Expr x = Expr("x");
  Expr y = Expr("y");
  Expr a = Expr("a");
  Expr b = Expr("b");
  Expr c = Expr("c");
  Expr d = Expr("d");

  Expr u = a * x + b * x + c + d;
  Expr S = list({x});
  Expr C = collectTerms(u, S);

  assert(C == (a + b) * x + (c + d));
}

void should_collect_polynomials() {
  Expr x = Expr("x");
  Expr y = Expr("y");
  Expr z = Expr("z");
  Expr u = Expr("u");

  Expr p0 = 4 * power(x, 2) * power(y, 5) * power(z, 3) +
            10 * power(x, 3) * y * z + 2 * power(y, 3) * power(z, 2) +
            power(x, 4) * power(z, 3) + 11 * power(x, 4) * power(y, 3) +
            power(z, 4) * x + y * power(z, 3);

	//printf("%s\n", collect(p0, list({y, x, z})).toString().c_str());
	assert(algebraicExpand(collect(p0, list({x, y, z}))) == algebraicExpand(p0));
}

void should_algebraic_expand_expressions() {
  Expr x = Expr("x");
  Expr y = Expr("y");
  Expr z = Expr("z");

  // Expr u0 = (x*power(y + 1, fraction(3, 2)) + 1)*((x*power(y + 1, fraction(3,
  // 2))) - 1);

  // assert(algebraicExpand(u0) == -1 + power(x, 2) + 3*power(x, 2)*y +
  // 3*power(x, 2)*power(y, 2) + power(x, 2)*power(y, 3));

  Expr u1 = power(x * power(y + 1, fraction(1, 2)) + 1, 4);
  assert(algebraicExpand(u1) ==
         1 + 6 * power(x, 2) + power(x, 4) + 6 * power(x, 2) * y +
             2 * power(x, 4) * y + power(x, 4) * power(y, 2) +
             4 * x * power(1 + y, fraction(1, 2)) +
             4 * power(x, 3) * power(1 + y, fraction(3, 2)));

  Expr u2 = (x + 2) * (x + 3) * (x + 4);

  assert(algebraicExpand(u2) == 24 + 26 * x + 9 * power(x, 2) + power(x, 3));

  Expr u3 = power(x + y + z, 3);
  assert(algebraicExpand(u3) ==
         power(x, 3) + 3 * power(x, 2) * y + 3 * x * power(y, 2) + power(y, 3) +
             3 * power(x, 2) * z + 6 * x * y * z + 3 * power(y, 2) * z +
             3 * x * power(z, 2) + 3 * y * power(z, 2) + power(z, 3));

  Expr u4 = power(x + 1, 2) + power(y + 1, 2);

  assert(algebraicExpand(u4) == 2 + 2 * x + power(x, 2) + 2 * y + power(y, 2));

  Expr u5 = power(power(x + 2, 2) + 3, 2);

  assert(algebraicExpand(u5) ==
         49 + 56 * x + 30 * power(x, 2) + 8 * power(x, 3) + power(x, 4));

  Expr u6 = (-32 * power(z, 3) + 32 * power(z, 4) + 48 * power(z, 5) +
             -24 * power(z, 6) + -48 * power(z, 7) + -36 * power(z, 8) +
             -40 * power(z, 9) + -8 * power(z, 10) + -8 * power(z, 11)) /
            (4 * power(z, 2));

  assert(algebraicExpand(u6) == -8 * z + 8 * power(z, 2) + 12 * power(z, 3) +
                                    -6 * power(z, 4) + -12 * power(z, 5) +
                                    -9 * power(z, 6) + -10 * power(z, 7) +
                                    -2 * power(z, 8) + -2 * power(z, 9));
}

void should_expand_main_operator() {
  Expr x = Expr("x");

  Expr u0 = x * (2 + power(1 + x, 2));

  assert(algebraicExpandRoot(u0) == 2 * x + x * power(1 + x, 2));

  Expr u1 = power(x + power(1 + x, 2), 2);

  assert(algebraicExpandRoot(u1) ==
         power(x, 2) + 2 * x * power(1 + x, 2) + power(1 + x, 4));
}

void should_get_polynomial_content() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 4 * power(x, 2) + -1 * 6 * x;

  Expr L = list({x});

  Expr Z = Expr("Z");
  Expr Q = Expr("Q");

  assert(cont(u, L, Z) == 2);

  Expr t = 2 * x;

  assert(cont(t, L, Z) == 2);

  Expr p = -1 * x;

  assert(cont(p, L, Z) == 1);

  Expr T = list({x, y});

  Expr a = fraction(1, 2) * x * y + 6 * y;

  assert(cont(a, T, Q) == y);

  Expr b = (power(y, 2) + 2 * y + 1) * power(x, 2) + (2 * power(y, 2) - 2) * x +
           (3 * y + 3);

  assert(cont(b, T, Q) == 1 + y);
}

void should_get_polynomial_content_sub_resultant() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 4 * power(x, 2) + -1 * 6 * x;

  Expr L = list({x});

  Expr Z = Expr("Z");
  Expr Q = Expr("Q");

  assert(cont(u, L, Z) == 2);

  Expr t = 2 * x;

  assert(cont(t, L, Z) == 2);

  Expr p = -1 * x;

  assert(cont(p, L, Z) == 1);

  Expr T = list({x, y});

  Expr a = fraction(1, 2) * x * y + 6 * y;

  assert(cont(a, T, Q) == y);

  Expr b = (power(y, 2) + 2 * y + 1) * power(x, 2) + (2 * power(y, 2) - 2) * x +
           (3 * y + 3);

  assert(cont(b, T, Q) == 1 + y);
}

void should_monomial_base_expand_polynomials() {
  Expr x = Expr("x");
  Expr a = Expr("a");
  Expr b = Expr("b");
  Expr t = Expr("t");

  Expr u =
      power(a, 2) * b + 2 * a * power(b, 2) + power(b, 3) + 2 * a + 2 * b + 3;

  Expr v = a + b;

  Expr L = list({a, b});

  Expr r = monomialBasedPolyExpansion(u, v, L, t);

  assert(r == 3 + 2 * t + b * power(t, 2));
}


int main() {
  TEST(should_get_polynomial_variable)
  TEST(should_get_degree_of_variables)
  TEST(should_get_coefficients)
  TEST(should_get_leading_coefficient)
  TEST(should_divided_polynomials)
  TEST(should_get_gcd_polynomials)
  TEST(should_get_extended_gcd_polynomials)
  TEST(should_calculate_monomial_division)
  TEST(should_get_leading_monomial)
  TEST(should_rec_divide_polynomials)
  TEST(should_pseudo_divide_polynomials)
  TEST(should_normalize_polynomial)
  TEST(should_mv_poly_gcd)
  TEST(should_get_coeff_var_parts_of_monomial)
  TEST(should_algebraic_expand_expressions)
  TEST(should_expand_main_operator)
  TEST(should_get_coeff_var_parts_of_monomial)
  TEST(should_collect_terms)
  TEST(should_get_polynomial_content)
  TEST(should_get_polynomial_content_sub_resultant)
  TEST(should_monomial_base_expand_polynomials)
  TEST(should_collect_polynomials)
  return 0;
}
