#include "Core/Polynomial/Polynomial.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Expand/Expand.hpp"

#include "Core/Polynomial/Resultant.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "test.hpp"
#include <climits>
#include <cstdio>
#include <cstdlib>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace polynomial;

void should_get_polynomial_variable() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr exp0 = 4 * x + power(x, 2) + 5 * power(x, 3);
  Expr exp1 = 4 * x + power(y, 2) + 5 * sin(x);

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

  assert(R0[0] == x * y);
  assert(R0[1] == x + -1 * x * y);
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

// void should_mv_poly_gcd() {
//   Expr x = Expr("x");
//   Expr y = Expr("y");

//   Expr u = -1 * y * power(x, 2) + power(y, 3);
//   Expr v = y * power(x, 2) + 2 * power(y, 2) * x + power(y, 3);

//   Expr L = list({x, y});
//   Expr Z = Expr("Z");

//   Expr gcd = mvPolyGCD(u, v, L, Z);

//   assert(gcd == x * y + power(y, 2));
// }

void should_get_poly_gcd_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");
  Expr z = Expr("z");
  Expr Z = Expr("Z");

  Expr L = list({x, y});

  Expr u = polyExpr(-1 * y * power(x, 2) + power(y, 3), L);
  Expr v = polyExpr(y * power(x, 2) + 2 * power(y, 2) * x + power(y, 3), L);

  Expr uv_gcd = gcdPolyExpr(u, v, L, Z);

	assert(uv_gcd == polyExpr(x * y + power(y, 2), L));

  // Expr T = list({x, y, z});

  // Expr a =
  //     polyExpr(3 * z * y * x + 2 * power(z, 3) * power(y, 2) * x +
  //     -2 * z * power(y, 3) * x + -3 * power(z, 3) * x +
  //     -2 * power(z, 2) * y * power(x, 2) + 8 * power(z, 3) * y * power(x, 2) +
  //     2 * power(z, 4) * y * power(x, 2) + -8 * z * power(y, 2) * power(x, 2) +
  //     6 * power(z, 2) * power(y, 2) * power(x, 2) +
  //     -6 * power(y, 3) * power(x, 2) + -12 * power(z, 2) * y * power(x, 3) +
  //     12 * power(z, 3) * y * power(x, 3) + -9 * z * power(y, 2) * power(x, 3) +
  //     2 * power(z, 3) * power(y, 2) * power(x, 3) +
  //     power(z, 5) * power(y, 2) * power(x, 3) +
  //     -1 * power(z, 3) * power(y, 3) * power(x, 3) +
  //     -2 * z * power(y, 4) * power(x, 3) + -3 * power(z, 3) * power(x, 3) +
  //     12 * power(z, 4) * power(x, 3) + 8 * power(z, 3) * y * power(x, 4) +
  //     2 * power(z, 4) * y * power(x, 4) + 4 * power(z, 5) * y * power(x, 4) +
  //     8 * power(z, 2) * power(y, 2) * power(x, 4) +
  //     -4 * power(z, 3) * power(y, 2) * power(x, 4) +
  //     4 * power(z, 4) * power(y, 2) * power(x, 4) +
  //     -8 * z * power(y, 3) * power(x, 4) +
  //     -6 * power(z, 2) * power(y, 3) * power(x, 4) +
  //     -8 * power(y, 4) * power(x, 4) + 12 * power(z, 3) * y * power(x, 5) +
  //     -12 * power(z, 2) * power(y, 2) * power(x, 5) +
  //     power(z, 5) * power(y, 2) * power(x, 5) +
  //     -12 * z * power(y, 3) * power(x, 5) +
  //     -1 * power(z, 3) * power(y, 4) * power(x, 5) +
  //     12 * power(z, 4) * power(x, 5) + 4 * power(z, 5) * y * power(x, 6) +
  //     4 * power(z, 4) * power(y, 2) * power(x, 6) +
  //     -4 * power(z, 3) * power(y, 3) * power(x, 6) +
  //     -4 * power(z, 2) * power(y, 4) * power(x, 6) + -2 * power(z, 2) * y +
	// 						 2 * power(y, 2), T);

  // Expr b =
	// 	polyExpr( -4 * power(z, 2) * y * x + 16 * power(z, 3) * y * x +
  //     4 * power(z, 4) * y * x + -16 * z * power(y, 2) * x +
  //     12 * power(z, 2) * power(y, 2) * x + -12 * power(y, 3) * x +
  //     -36 * power(z, 2) * y * power(x, 2) + 36 * power(z, 3) * y * power(x, 2) +
  //     -27 * z * power(y, 2) * power(x, 2) +
  //     6 * power(z, 3) * power(y, 2) * power(x, 2) +
  //     3 * power(z, 5) * power(y, 2) * power(x, 2) +
  //     -3 * power(z, 3) * power(y, 3) * power(x, 2) +
  //     -6 * z * power(y, 4) * power(x, 2) + -9 * power(z, 3) * power(x, 2) +
  //     36 * power(z, 4) * power(x, 2) + 32 * power(z, 3) * y * power(x, 3) +
  //     8 * power(z, 4) * y * power(x, 3) + 16 * power(z, 5) * y * power(x, 3) +
  //     32 * power(z, 2) * power(y, 2) * power(x, 3) +
  //     -16 * power(z, 3) * power(y, 2) * power(x, 3) +
  //     16 * power(z, 4) * power(y, 2) * power(x, 3) +
  //     -32 * z * power(y, 3) * power(x, 3) +
  //     -24 * power(z, 2) * power(y, 3) * power(x, 3) +
  //     -32 * power(y, 4) * power(x, 3) + 60 * power(z, 3) * y * power(x, 4) +
  //     -60 * power(z, 2) * power(y, 2) * power(x, 4) +
  //     5 * power(z, 5) * power(y, 2) * power(x, 4) +
  //     -60 * z * power(y, 3) * power(x, 4) +
  //     -5 * power(z, 3) * power(y, 4) * power(x, 4) +
  //     60 * power(z, 4) * power(x, 4) + 24 * power(z, 5) * y * power(x, 5) +
  //     24 * power(z, 4) * power(y, 2) * power(x, 5) +
  //     -24 * power(z, 3) * power(y, 3) * power(x, 5) +
  //     -24 * power(z, 2) * power(y, 4) * power(x, 5) + 3 * z * y +
	// 						2 * power(z, 3) * power(y, 2) + -2 * z * power(y, 3) + -3 * power(z, 3), T);


  // Expr ab_gcd = gcdPolyExpr(a, b, T, Z);
  // printf("----> %s\n", ab_gcd.toString().c_str());
}

void should_get_poly_gcd() {
  Expr x = Expr("x");
  Expr y = Expr("y");
  Expr z = Expr("z");
  Expr Z = Expr("Z");

  Expr L = list({x, y});

  Expr u = -1 * y * power(x, 2) + power(y, 3);
  Expr v = y * power(x, 2) + 2 * power(y, 2) * x + power(y, 3);

  Expr uv_gcd = gcdPoly(u, v, L, Z);

  assert(uv_gcd == x * y + power(y, 2));
}

void should_get_heuristic_gcd_of_polys() {
	Expr x = Expr("x");
	Expr y = Expr("y");
	Expr z = Expr("z");

	Expr Z = Expr("Z");
  Expr L = list({x, y});

  Expr u = -1 * y * power(x, 2) + power(y, 3);
  Expr v = y * power(x, 2) + 2 * power(y, 2) * x + power(y, 3);

	Expr U = heuristicGcdPoly(u, v, L, Z);

  assert(U[0] == x * y + power(y, 2));
  assert(U[1] == y + -x);
  assert(U[2] == y + x);

  Expr a =
      3 * z * y * x + 2 * power(z, 3) * power(y, 2) * x +
      -2 * z * power(y, 3) * x + -3 * power(z, 3) * x +
      -2 * power(z, 2) * y * power(x, 2) + 8 * power(z, 3) * y * power(x, 2) +
      2 * power(z, 4) * y * power(x, 2) + -8 * z * power(y, 2) * power(x, 2) +
      6 * power(z, 2) * power(y, 2) * power(x, 2) +
      -6 * power(y, 3) * power(x, 2) + -12 * power(z, 2) * y * power(x, 3) +
      12 * power(z, 3) * y * power(x, 3) + -9 * z * power(y, 2) * power(x, 3) +
      2 * power(z, 3) * power(y, 2) * power(x, 3) +
      power(z, 5) * power(y, 2) * power(x, 3) +
      -1 * power(z, 3) * power(y, 3) * power(x, 3) +
      -2 * z * power(y, 4) * power(x, 3) + -3 * power(z, 3) * power(x, 3) +
      12 * power(z, 4) * power(x, 3) + 8 * power(z, 3) * y * power(x, 4) +
      2 * power(z, 4) * y * power(x, 4) + 4 * power(z, 5) * y * power(x, 4) +
      8 * power(z, 2) * power(y, 2) * power(x, 4) +
      -4 * power(z, 3) * power(y, 2) * power(x, 4) +
      4 * power(z, 4) * power(y, 2) * power(x, 4) +
      -8 * z * power(y, 3) * power(x, 4) +
      -6 * power(z, 2) * power(y, 3) * power(x, 4) +
      -8 * power(y, 4) * power(x, 4) + 12 * power(z, 3) * y * power(x, 5) +
      -12 * power(z, 2) * power(y, 2) * power(x, 5) +
      power(z, 5) * power(y, 2) * power(x, 5) +
      -12 * z * power(y, 3) * power(x, 5) +
      -1 * power(z, 3) * power(y, 4) * power(x, 5) +
      12 * power(z, 4) * power(x, 5) + 4 * power(z, 5) * y * power(x, 6) +
      4 * power(z, 4) * power(y, 2) * power(x, 6) +
      -4 * power(z, 3) * power(y, 3) * power(x, 6) +
      -4 * power(z, 2) * power(y, 4) * power(x, 6) + -2 * power(z, 2) * y +
      2 * power(y, 2);

  Expr b =
      -4 * power(z, 2) * y * x + 16 * power(z, 3) * y * x +
      4 * power(z, 4) * y * x + -16 * z * power(y, 2) * x +
      12 * power(z, 2) * power(y, 2) * x + -12 * power(y, 3) * x +
      -36 * power(z, 2) * y * power(x, 2) + 36 * power(z, 3) * y * power(x, 2) +
      -27 * z * power(y, 2) * power(x, 2) +
      6 * power(z, 3) * power(y, 2) * power(x, 2) +
      3 * power(z, 5) * power(y, 2) * power(x, 2) +
      -3 * power(z, 3) * power(y, 3) * power(x, 2) +
      -6 * z * power(y, 4) * power(x, 2) + -9 * power(z, 3) * power(x, 2) +
      36 * power(z, 4) * power(x, 2) + 32 * power(z, 3) * y * power(x, 3) +
      8 * power(z, 4) * y * power(x, 3) + 16 * power(z, 5) * y * power(x, 3) +
      32 * power(z, 2) * power(y, 2) * power(x, 3) +
      -16 * power(z, 3) * power(y, 2) * power(x, 3) +
      16 * power(z, 4) * power(y, 2) * power(x, 3) +
      -32 * z * power(y, 3) * power(x, 3) +
      -24 * power(z, 2) * power(y, 3) * power(x, 3) +
      -32 * power(y, 4) * power(x, 3) + 60 * power(z, 3) * y * power(x, 4) +
      -60 * power(z, 2) * power(y, 2) * power(x, 4) +
      5 * power(z, 5) * power(y, 2) * power(x, 4) +
      -60 * z * power(y, 3) * power(x, 4) +
      -5 * power(z, 3) * power(y, 4) * power(x, 4) +
      60 * power(z, 4) * power(x, 4) + 24 * power(z, 5) * y * power(x, 5) +
      24 * power(z, 4) * power(y, 2) * power(x, 5) +
      -24 * power(z, 3) * power(y, 3) * power(x, 5) +
      -24 * power(z, 2) * power(y, 4) * power(x, 5) + 3 * z * y +
      2 * power(z, 3) * power(y, 2) + -2 * z * power(y, 3) + -3 * power(z, 3);

  Expr T = list({x, y, z});
	Expr K = heuristicGcdPoly(a, b, T, Z);

  assert(K[0] == 1);
  assert(K[1] == a);
  assert(K[2] == b);


	Expr c = 4*power(z, 5)*y + 4*power(z, 4)*power(y, 2) + -4*power(z, 3)*power(y, 3) + -4*power(z, 2)*power(y, 4);

	Expr d = 12*power(z, 3)*y + -12*power(z, 2)*power(y, 2) + power(z, 5)*power(y, 2) + -12*z*power(y, 3) + -1*power(z, 3)*power(y, 4) + 12*power(z, 4);

	Expr O = list({y, z});

	Expr J = heuristicGcdPoly(c, d, O, Z);

	assert(J[0] == z*power(y, 2) + -1*power(z, 3));
	assert(J[1] == -4*power(z, 2)*y + -4*z*power(y, 2));
	assert(J[2] == -12*y + -1*power(z, 2)*power(y, 2) + -12*z);


	Expr e = power(y, 2);
	Expr f = 2*y;
	Expr P = list({y});
	Expr I = heuristicGcdPoly(e, f, P, Z);

	assert(I[0] == y);
	assert(I[1] == y);
	assert(I[2] == 2);
}


void should_get_heuristic_gcd_of_poly_exprs() {
	Expr x = Expr("x");
	Expr y = Expr("y");
	Expr z = Expr("z");

	Expr Z = Expr("Z");

  Expr T = list({x, y, z});

  Expr L = list({x, y});

  Expr u = polyExpr(-1 * y * power(x, 2) + power(y, 3), L);
  Expr v = polyExpr(y * power(x, 2) + 2 * power(y, 2) * x + power(y, 3), L);

	Expr U = heuristicGcdPolyExpr(u, v, L, Z);

  assert(U[0] == polyExpr(x * y + power(y, 2), L));
  assert(U[1] == polyExpr(y + -x, L));
  assert(U[2] == polyExpr(y + x, L));

	Expr a = polyExpr(
      3 * z * y * x + 2 * power(z, 3) * power(y, 2) * x +
      -2 * z * power(y, 3) * x + -3 * power(z, 3) * x +
      -2 * power(z, 2) * y * power(x, 2) + 8 * power(z, 3) * y * power(x, 2) +
      2 * power(z, 4) * y * power(x, 2) + -8 * z * power(y, 2) * power(x, 2) +
      6 * power(z, 2) * power(y, 2) * power(x, 2) +
      -6 * power(y, 3) * power(x, 2) + -12 * power(z, 2) * y * power(x, 3) +
      12 * power(z, 3) * y * power(x, 3) + -9 * z * power(y, 2) * power(x, 3) +
      2 * power(z, 3) * power(y, 2) * power(x, 3) +
      power(z, 5) * power(y, 2) * power(x, 3) +
      -1 * power(z, 3) * power(y, 3) * power(x, 3) +
      -2 * z * power(y, 4) * power(x, 3) + -3 * power(z, 3) * power(x, 3) +
      12 * power(z, 4) * power(x, 3) + 8 * power(z, 3) * y * power(x, 4) +
      2 * power(z, 4) * y * power(x, 4) + 4 * power(z, 5) * y * power(x, 4) +
      8 * power(z, 2) * power(y, 2) * power(x, 4) +
      -4 * power(z, 3) * power(y, 2) * power(x, 4) +
      4 * power(z, 4) * power(y, 2) * power(x, 4) +
      -8 * z * power(y, 3) * power(x, 4) +
      -6 * power(z, 2) * power(y, 3) * power(x, 4) +
      -8 * power(y, 4) * power(x, 4) + 12 * power(z, 3) * y * power(x, 5) +
      -12 * power(z, 2) * power(y, 2) * power(x, 5) +
      power(z, 5) * power(y, 2) * power(x, 5) +
      -12 * z * power(y, 3) * power(x, 5) +
      -1 * power(z, 3) * power(y, 4) * power(x, 5) +
      12 * power(z, 4) * power(x, 5) + 4 * power(z, 5) * y * power(x, 6) +
      4 * power(z, 4) * power(y, 2) * power(x, 6) +
      -4 * power(z, 3) * power(y, 3) * power(x, 6) +
      -4 * power(z, 2) * power(y, 4) * power(x, 6) + -2 * power(z, 2) * y +
      2 * power(y, 2), T);

  Expr b = polyExpr(
      -4 * power(z, 2) * y * x + 16 * power(z, 3) * y * x +
      4 * power(z, 4) * y * x + -16 * z * power(y, 2) * x +
      12 * power(z, 2) * power(y, 2) * x + -12 * power(y, 3) * x +
      -36 * power(z, 2) * y * power(x, 2) + 36 * power(z, 3) * y * power(x, 2) +
      -27 * z * power(y, 2) * power(x, 2) +
      6 * power(z, 3) * power(y, 2) * power(x, 2) +
      3 * power(z, 5) * power(y, 2) * power(x, 2) +
      -3 * power(z, 3) * power(y, 3) * power(x, 2) +
      -6 * z * power(y, 4) * power(x, 2) + -9 * power(z, 3) * power(x, 2) +
      36 * power(z, 4) * power(x, 2) + 32 * power(z, 3) * y * power(x, 3) +
      8 * power(z, 4) * y * power(x, 3) + 16 * power(z, 5) * y * power(x, 3) +
      32 * power(z, 2) * power(y, 2) * power(x, 3) +
      -16 * power(z, 3) * power(y, 2) * power(x, 3) +
      16 * power(z, 4) * power(y, 2) * power(x, 3) +
      -32 * z * power(y, 3) * power(x, 3) +
      -24 * power(z, 2) * power(y, 3) * power(x, 3) +
      -32 * power(y, 4) * power(x, 3) + 60 * power(z, 3) * y * power(x, 4) +
      -60 * power(z, 2) * power(y, 2) * power(x, 4) +
      5 * power(z, 5) * power(y, 2) * power(x, 4) +
      -60 * z * power(y, 3) * power(x, 4) +
      -5 * power(z, 3) * power(y, 4) * power(x, 4) +
      60 * power(z, 4) * power(x, 4) + 24 * power(z, 5) * y * power(x, 5) +
      24 * power(z, 4) * power(y, 2) * power(x, 5) +
      -24 * power(z, 3) * power(y, 3) * power(x, 5) +
      -24 * power(z, 2) * power(y, 4) * power(x, 5) + 3 * z * y +
      2 * power(z, 3) * power(y, 2) + -2 * z * power(y, 3) + -3 * power(z, 3), T);

	Expr K = heuristicGcdPolyExpr(a, b, T, Z);

  assert(K[0] == polyExpr(1, T));
  assert(K[1] == a);
  assert(K[2] == b);

	Expr O = list({y, z});

	Expr c = polyExpr(4*power(z, 5)*y + 4*power(z, 4)*power(y, 2) + -4*power(z, 3)*power(y, 3) + -4*power(z, 2)*power(y, 4), O);

	Expr d = polyExpr(12*power(z, 3)*y + -12*power(z, 2)*power(y, 2) + power(z, 5)*power(y, 2) + -12*z*power(y, 3) + -1*power(z, 3)*power(y, 4) + 12*power(z, 4), O);

	Expr J = heuristicGcdPolyExpr(c, d, O, Z);

	assert(J[0] == polyExpr(z*power(y, 2) + -1*power(z, 3), O));
	assert(J[1] == polyExpr(-4*power(z, 2)*y + -4*z*power(y, 2), O));
	assert(J[2] == polyExpr(-12*y + -1*power(z, 2)*power(y, 2) + -12*z, O));


	Expr P = list({y});
	Expr e = polyExpr(power(y, 2), P);
	Expr f = polyExpr(2*y, P);
	Expr I = heuristicGcdPolyExpr(e, f, P, Z);

	assert(I[0] == polyExpr(y, P));
	assert(I[1] == polyExpr(y, P));
	assert(I[2] == polyExpr(2, P));
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

  assert(polyExpr(p0, list({y, x, z})) ==
         add({add({add({1 * power(z, 4)}) * power(x, 1),
                   add({1 * power(z, 3)}) * power(x, 4)}) *
                  power(y, 0),
              add({add({1 * power(z, 3)}) * power(x, 0),
                   add({10 * power(z, 1)}) * power(x, 3)}) *
                  power(y, 1),
              add({add({2 * power(z, 2)}) * power(x, 0),
                   add({11 * power(z, 0)}) * power(x, 4)}) *
                  power(y, 3),
              add({add({4 * power(z, 3)}) * power(x, 2)}) * power(y, 5)}));
}

void should_algebraic_expand_expressions() {
  Expr x = Expr("x");
  Expr y = Expr("y");
  Expr z = Expr("z");

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

  Expr u7 =
      power(-1 * (-5926821888 * y + -6440585216L * (power(y, 2)) +
                  -4756602880L * (power(y, 3)) + 2305909760L * (power(y, 4)) +
                  -168882304 * (power(y, 5)) + -268451584 * (power(y, 6)) +
                  31912288 * (power(y, 7)) + 3696960 * (power(y, 8)) +
                  -648480 * (power(y, 9)) + 44472 * (power(y, 10)) +
                  -1456 * (power(y, 11)) + 24 * (power(y, 12)) + -4175495168),
            5 - 4);
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
  Expr z = Expr("z");

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

void should_get_if_poly_expr_is_zero() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 0;

  assert(isZeroPolyExpr(u) == true);

  Expr g = 0 * power(x, 4);

  assert(isZeroPolyExpr(g) == true);

  Expr t = polyExpr(0 * power(x, 3) * y, list({x, y}));

  assert(isZeroPolyExpr(t) == true);

  Expr k = polyExpr(0, list({x, y}));

  assert(isZeroPolyExpr(k) == true);
}

void should_mul_collected_polys() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 2 * x + 4 * power(x, 2) * power(y, 4) + 10 * power(y, 2) +
           8 * power(x, 6) + 7 * power(y, 3) + power(x, 3);

  Expr v = 5 * x + 7 * power(x, 3) * power(y, 2) + 11 * y + 8 * power(x, 4) +
           2 * power(y, 3) + power(x, 5);

  Expr L = list({x, y});

  Expr uv = mulPolyExpr(polyExpr(u, L), polyExpr(v, L));

  assert(uv ==
         add({add({110 * power(y, 3), 77 * power(y, 4), 20 * power(y, 5),
                   14 * power(y, 6)}) *
                  power(x, 0),
              add({22 * power(y, 1), 50 * power(y, 2), 39 * power(y, 3)}) *
                  power(x, 1),
              add({10 * power(y, 0), 44 * power(y, 5), 8 * power(y, 7)}) *
                  power(x, 2),
              add({11 * power(y, 1), 2 * power(y, 3), 90 * power(y, 4),
                   49 * power(y, 5)}) *
                  power(x, 3),
              add({5 * power(y, 0), 94 * power(y, 2), 56 * power(y, 3)}) *
                  power(x, 4),
              add({16 * power(y, 0), 10 * power(y, 2), 7 * power(y, 3),
                   28 * power(y, 6)}) *
                  power(x, 5),
              add({2 * power(y, 0), 88 * power(y, 1), 7 * power(y, 2),
                   16 * power(y, 3), 32 * power(y, 4)}) *
                  power(x, 6),
              add({48 * power(y, 0), 4 * power(y, 4)}) * power(x, 7),
              add({1 * power(y, 0)}) * power(x, 8),
              add({56 * power(y, 2)}) * power(x, 9),
              add({64 * power(y, 0)}) * power(x, 10),
              add({8 * power(y, 0)}) * power(x, 11)}));

  Expr g = 3;
  Expr t = 5;

  assert(mulPolyExpr(polyExpr(g, L), polyExpr(t, L)) == polyExpr(15, L));
  assert(mulPolyExpr(polyExpr(g, list({})), polyExpr(t, list({}))) ==
         polyExpr(15, list({})));

  Expr k = 4 * x + 15 * power(x, 2) + 4 * y * x + 5 * y;
  assert(mulPolyExpr(polyExpr(k, L), 5) ==
         add({
             add({25 * power(y, 1)}) * power(x, 0),
             add({20 * power(y, 0), 20 * power(y, 1)}) * power(x, 1),
             add({75 * power(y, 0)}) * power(x, 2),
         }));

  assert(mulPolyExpr(polyExpr(3, L), polyExpr(4, L)) == polyExpr(12, L));
  assert(mulPolyExpr(polyExpr(4, list({})), polyExpr(3, list({}))) ==
         polyExpr(12, list({})));

  Expr r = 4 * power(x, 3) + 5 * power(x, 2) + 3 * x + 4;
  Expr h = 2 * power(x, 2) + 3 * x + 7;

  assert(mulPolyExpr(polyExpr(r, list({x})), polyExpr(h, list({x}))) ==
         28 * power(x, 0) + 33 * power(x, 1) + 52 * power(x, 2) +
             49 * power(x, 3) + 22 * power(x, 4) + 8 * power(x, 5));
}

void should_add_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 2 * power(x, 2) * y + 3 * power(y, 2) * x + 4 * x;
  Expr v = 3 * power(x, 2) + 4 * power(y, 2) * x + 2 * x;

  Expr L = list({x, y});

  Expr r = addPolyExpr(polyExpr(u, L), polyExpr(v, L));

  assert(r == add({add({6 * power(y, 0), 7 * power(y, 2)}) * power(x, 1),
                   add({3 * power(y, 0), 2 * power(y, 1)}) * power(x, 2)}));

  Expr g = 2 * power(y, 3) * power(x, 5) + 4 * power(x, 4) + 4 * y + 4 * x;
  Expr t = 5 * power(y, 3) * power(x, 5) + 3 * x + y;

  Expr gt = addPolyExpr(polyExpr(g, L), polyExpr(t, L));

  assert(gt == add({add({5 * power(y, 1)}) * power(x, 0),
                    add({7 * power(y, 0)}) * power(x, 1),
                    add({4 * power(y, 0)}) * power(x, 4),
                    add({7 * power(y, 3)}) * power(x, 5)}));

  Expr e = 4 * power(x, 3) + 5 * power(x, 2) + 3 * x + 4;
  Expr h = 2 * power(x, 2) + 3 * x + 7;
  Expr q = 2 * power(x, 2) + -3 * x + 7;

  assert(addPolyExpr(polyExpr(e, list({x})), polyExpr(h, list({x}))) ==
         11 * power(x, 0) + 6 * power(x, 1) + 7 * power(x, 2) +
             4 * power(x, 3));

  assert(addPolyExpr(polyExpr(h, list({x})), polyExpr(q, list({x}))) ==
         4 * power(x, 2) + 14 * power(x, 0));

  assert(addPolyExpr(polyExpr(e, list({x})), polyExpr(0, list({x}))) ==
         polyExpr(e, list({x})));
  assert(addPolyExpr(polyExpr(0, list({x})), polyExpr(e, list({x}))) ==
         polyExpr(e, list({x})));
}

void should_sub_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 2 * y * power(x, 2) + 3 * power(y, 2) * x + 4 * x;
  Expr v = 3 * power(x, 2) + 4 * power(y, 2) * x + 2 * x;

  Expr L = list({x, y});

  assert(subPolyExpr(polyExpr(u, L), polyExpr(v, L)) ==
         add({
             add({2 * power(y, 0), -1 * power(y, 2)}) * power(x, 1),
             add({-3 * power(y, 0), 2 * power(y, 1)}) * power(x, 2),
         }));

  assert(subPolyExpr(polyExpr(3, L), polyExpr(4, L)) == polyExpr(-1, L));
  assert(subPolyExpr(polyExpr(4, list({})), polyExpr(3, list({}))) ==
         polyExpr(1, list({})));

  Expr e = 4 * power(x, 3) + 5 * power(x, 2) + 3 * x + 4;
  Expr h = 2 * power(x, 2) + 3 * x + 7;

  assert(subPolyExpr(polyExpr(e, list({x})), polyExpr(h, list({x}))) ==
         -3 * power(x, 0) + 3 * power(x, 2) + 4 * power(x, 3));

  assert(subPolyExpr(polyExpr(e, list({x})), polyExpr(0, list({x}))) ==
         polyExpr(e, list({x})));
  assert(subPolyExpr(polyExpr(0, list({x})), polyExpr(e, list({x}))) ==
         mulPolyExpr(polyExpr(e, list({x})), -1));

  Expr a = add({add({-1 * power(y, 0), 1 * power(y, 2)}) * power(x, 2)});
  Expr b = add({add({-1 * power(y, 0), 1 * power(y, 2)}) * power(x, 2),
                add({-1 * power(y, 2), 1 * power(y, 3)}) * power(x, 3)});

  assert(subPolyExpr(a, b) ==
         add({add({1 * power(y, 2), -1 * power(y, 3)}) * power(x, 3)}));

  Expr c = -1 * power(y, 0) + -1 * power(y, 1);
  Expr d = add({-1 * power(y, 1)});

  assert(subPolyExpr(c, d) == add({-1 * power(y, 0)}));
}
void should_rec_divide_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = power(x, 2) * power(y, 2) + x;
  Expr v = x * y + 1;

  Expr L0 = list({x, y});
  Expr L1 = list({y, x});

  Expr Q = Expr("Q");

  Expr R0 = divPolyExpr(polyExpr(u, L0), polyExpr(v, L0), L0, Q);

  assert(R0 ==
         list({
             add({add({1 * power(y, 1)}) * power(x, 1)}),
             add({add({1 * power(y, 0), -1 * power(y, 1)}) * power(x, 1)}),
         }));
}

void should_pseudo_divide_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr u = 5 * power(x, 4) * power(y, 3) + 3 * x * y + 2;
  Expr v = 2 * power(x, 3) * y + 2 * x + 3;

  Expr L = list({x, y});

  Expr R = pseudoDivPolyExpr(polyExpr(u, L), polyExpr(v, L), L);

  assert(R ==
         list({add({add({10 * power(y, 4)}) * power(x, 1)}),
               add({
                   add({8 * power(y, 2)}) * power(x, 0),
                   add({12 * power(y, 3), -30 * power(y, 4)}) * power(x, 1),
                   add({-20 * power(y, 4)}) * power(x, 2),
               })}));
}

void should_pow_col_poly() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr L = list({x});

  Expr u = polyExpr(2 * power(x, 2) + 4 * x + 12, L);

  assert(powPolyExpr(u, 2) == add({
                                     144 * power(x, 0),
                                     96 * power(x, 1),
                                     64 * power(x, 2),
                                     16 * power(x, 3),
                                     4 * power(x, 4),
                                 }));

  assert(powPolyExpr(u, 5) == add({
                                     248832 * power(x, 0),
                                     414720 * power(x, 1),
                                     483840 * power(x, 2),
                                     368640 * power(x, 3),
                                     222720 * power(x, 4),
                                     100864 * power(x, 5),
                                     37120 * power(x, 6),
                                     10240 * power(x, 7),
                                     2240 * power(x, 8),
                                     320 * power(x, 9),
                                     32 * power(x, 10),
                                 }));

  Expr T = list({x, y});

  Expr v = polyExpr(2 * power(x, 3) + 4 * power(y, 2) * x + 12 * y + 14, T);
  assert(
      powPolyExpr(v, 3) ==
      add({add({2744 * power(y, 0), 7056 * power(y, 1), 6048 * power(y, 2),
                1728 * power(y, 3)}) *
               power(x, 0),
           add({2352 * power(y, 2), 4032 * power(y, 3), 1728 * power(y, 4)}) *
               power(x, 1),
           add({672 * power(y, 4), 576 * power(y, 5)}) * power(x, 2),
           add({1176 * power(y, 0), 2016 * power(y, 1), 864 * power(y, 2),
                64 * power(y, 6)}) *
               power(x, 3),
           add({672 * power(y, 2), 576 * power(y, 3)}) * power(x, 4),
           add({96 * power(y, 4)}) * power(x, 5),
           add({168 * power(y, 0), 144 * power(y, 1)}) * power(x, 6),
           add({48 * power(y, 2)}) * power(x, 7),
           add({8 * power(y, 0)}) * power(x, 9)}));
}

void should_normalize_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr L = list({x, y});
  Expr Q = Expr("Q");

  Expr u = polyExpr(2 * x * y * x + x + 6 * y + 3, L);

  assert(normalizePolyExpr(u, L, Q) ==
         add({
             add({fraction(3, 2) * power(y, 0), 3 * power(y, 1)}) * power(x, 0),
             add({fraction(1, 2) * power(y, 0), 1 * power(y, 1)}) * power(x, 1),
         }));
}

void should_get_gcd_of_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr L = list({x, y});

  Expr u = polyExpr(-1 * y * power(x, 2) + power(y, 3), L);
  Expr v = polyExpr(y * power(x, 2) + 2 * power(y, 2) * x + power(y, 3), L);

  Expr Z = Expr("Z");

  Expr gcd = gcdPolyExpr(u, v, L, Z);

  assert(gcd == add({
                    add({1 * power(y, 2)}) * power(x, 0),
                    add({1 * power(y, 1)}) * power(x, 1),
                }));
}

void should_get_content_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr L = list({x});

  Expr Z = Expr("Z");
  Expr Q = Expr("Q");

  Expr u = polyExpr(4 * power(x, 2) + -1 * 6 * x, L);

  assert(contPolyExpr(u, L, Z) == 2);

  Expr t = polyExpr(2 * x, L);

  assert(contPolyExpr(t, L, Z) == 2);

  Expr p = polyExpr(-1 * x, L);

  assert(contPolyExpr(p, L, Z) == 1);

  Expr T = list({x, y});

  Expr a = polyExpr(fraction(1, 2) * x * y + 6 * y, T);

  assert(contPolyExpr(a, T, Q) == add({1 * power(y, 1)}));
  Expr b = polyExpr(power(y, 2) * power(x, 2) + 2 * y * power(x, 2) +
                        power(x, 2) + 2 * power(y, 2) * x + -2 * x + 3 * y + 3,
                    T);

  assert(contPolyExpr(b, T, Z) == 1 * power(y, 0) + 1 * power(y, 1));
}

void should_expand_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr L = list({x});
  Expr T = list({x, y});

  Expr u = 10 * power(x, 4) + 4 * power(x, 2) + 10;

  assert(expandPolyExpr(polyExpr(u, L)) == u);

  Expr v = 10*power(x, 4)*power(y, 2) + 6*power(x, 2) + 10*x +
           4*power(x, 4)*y + 15*y + 10*x*y;

	assert(expandPolyExpr(polyExpr(v, T)) == v);
}

void should_diff_poly_expr() {
  Expr x = Expr("x");
  Expr y = Expr("y");

  Expr L = list({x});
  Expr T = list({x, y});

  Expr u = polyExpr(10 * power(x, 4) + 4 * power(x, 2) + 10, L);

  assert(diffPolyExpr(u, x) == 8 * power(x, 1) + 40 * power(x, 3));

  assert(diffPolyExpr(u, y) == Expr(Kind::Addition, {0 * power(x, 0)}));

  Expr v = polyExpr(10 * power(x, 4) * power(y, 2) + 6 * power(x, 2) * y +
                        10 * x + 4 * power(x, 4) * y + 15 * y + 10 * x * y,
                    T);

  assert(diffPolyExpr(v, x) ==
         add({
             add({10 * power(y, 0), 10 * power(y, 1)}) * power(x, 0),
             add({12 * power(y, 1)}) * power(x, 1),
             add({16 * power(y, 1), 40 * power(y, 2)}) * power(x, 3),
         }));
  assert(diffPolyExpr(v, y) ==
         add({
             add({15 * power(y, 0)}) * power(x, 0),
             add({10 * power(y, 0)}) * power(x, 1),
             add({6 * power(y, 0)}) * power(x, 2),
             add({4 * power(y, 0), 20 * power(y, 1)}) * power(x, 4),
         }));
}

void should_remove_denominators_from_polys() {
	Expr x = Expr("x");
	Expr y = Expr("y");
	Expr Z = Expr("Z");

	Expr L = list({x, y});

	Expr u = fraction(1, 2)*x + fraction(1, 3)*y + 1;

	Expr v = removeDenominatorsPoly(u, L, Z);

	assert(v == list({6, 3*x + 2*y + 6}));
}

void should_remove_denominators_from_poly_expr() {
	Expr x = Expr("x");
	Expr y = Expr("y");
	Expr Z = Expr("Z");

	Expr L = list({x, y});

	Expr u = polyExpr(fraction(1, 2)*x + fraction(1, 3)*y + 1, L);

	Expr v = removeDenominatorsPolyExpr(u, L, Z);

	assert(v == list({polyExpr(6, L), polyExpr(3*x + 2*y + 6, L)}));
}


int main() {
  TEST(should_get_polynomial_variable)
  TEST(should_get_degree_of_variables)
  TEST(should_get_coefficients)
  TEST(should_get_leading_coefficient)
  TEST(should_algebraic_expand_expressions)
  TEST(should_divided_polynomials)
  TEST(should_get_gcd_polynomials)
  TEST(should_get_extended_gcd_polynomials)
  TEST(should_calculate_monomial_division)
  TEST(should_get_leading_monomial)
  TEST(should_rec_divide_polynomials)
  TEST(should_pseudo_divide_polynomials)
  TEST(should_normalize_polynomial)
  TEST(should_get_coeff_var_parts_of_monomial)
  TEST(should_expand_main_operator)
  TEST(should_get_coeff_var_parts_of_monomial)
  TEST(should_collect_terms)
	TEST(should_remove_denominators_from_polys);
	TEST(should_remove_denominators_from_poly_expr);
  TEST(should_get_polynomial_content)
  TEST(should_get_polynomial_content_sub_resultant)
  TEST(should_monomial_base_expand_polynomials)
  TEST(should_collect_polynomials)
  TEST(should_get_if_poly_expr_is_zero)
  TEST(should_add_poly_expr)
  TEST(should_mul_collected_polys)
  TEST(should_sub_poly_expr)
  TEST(should_rec_divide_poly_expr)
  TEST(should_pseudo_divide_poly_expr)
  TEST(should_pow_col_poly)
  TEST(should_normalize_poly_expr)
  TEST(should_get_gcd_of_poly_expr)
  TEST(should_expand_poly_expr)
  TEST(should_get_content_poly_expr)
  TEST(should_diff_poly_expr)
  TEST(should_get_poly_gcd_poly_expr)
  TEST(should_get_poly_gcd)
	TEST(should_get_heuristic_gcd_of_polys)
	TEST(should_get_heuristic_gcd_of_poly_exprs)

  return 0;
}
