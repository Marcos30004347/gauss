#include "gauss/Polynomial/Polynomial.hpp"

#include "gauss/Algebra/Expression.hpp"
#include "test.hpp"

#include <climits>
#include <cstdio>
#include <cstdlib>

using namespace alg;
using namespace poly;

// void should_get_polynomial_variable() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr exp0 = 4 * x + pow(x, 2) + 5 * pow(x, 3);
//   expr exp1 = 4 * x + pow(y, 2) + 5 * sin(x);

//   assert(variables(exp0) == set({x}));
//   assert(variables(exp1) == set({y, x, sin(x)}));
// }

// void should_get_degree_of_variables() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr exp0 = 4 * x + pow(x, 2) + 5 * pow(x, 3);
//   expr exp1 = 4 * x + pow(y, 2) + pow(sin(x), 5);

// 	assert(degree(exp0, x) == 3);
//   assert(degree(exp1, x) == 1);
//   assert(degree(exp1, y) == 2);
//   assert(degree(exp1, sin(x)) == 5);
// }

// void should_get_coefficients() {
//   expr x = expr("x");
//   expr a = expr("a");
//   expr b = expr("b");

//   expr exp0 = 4 * pow(x, 2);
//   expr exp1 = a * pow(x, 2) + b * pow(x, 2);
//   expr exp2 = a * pow(x, 2) - b * pow(x, 2);

//   assert(coeff(exp0, x, 2) == 4);
//   assert(coeff(exp1, x, 2) == a + b);
//   assert(coeff(exp2, x, 2) == a - b);
// }

// void should_get_leading_coefficient() {
//   expr x = expr("x");
//   expr a = expr("a");
//   expr b = expr("b");

//   expr exp0 = 4 * pow(x, 2);
//   expr exp1 = a * pow(x, 2) + b * pow(x, 2);
//   expr exp2 = a * pow(x, 2) + b * pow(x, 3);

//   expr leadcoeff_exp0 = leadCoeff(exp0, x);
//   expr leadcoeff_exp1 = leadCoeff(exp1, x);
//   expr leadcoeff_exp2 = leadCoeff(exp2, x);

//   assert(leadCoeff(exp0, x) == 4);
//   assert(leadCoeff(exp1, x) == a + b);
//   assert(leadCoeff(exp2, x) == b);
// }

// void should_divided_polynomials() {
//   expr x = expr("x");
//   expr exp0 = 5 * pow(x, 2) + 4 * x + 1;
//   expr exp1 = 2 * x + 3;

// 	expr res = divideGPE(exp0, exp1, x);

// 	assert(res == list({fraction(-7, 4) + fraction(5, 2) * x, fraction(25, 4)}));
// }

// void should_get_gcd_polynomials() {
//   expr x = expr("x");
//   expr u = pow(x, 7) + -4 * pow(x, 5) + -1 * pow(x, 2) + 4;
//   expr v = pow(x, 5) + -4 * pow(x, 3) + -1 * pow(x, 2) + 4;

//   expr res = gcdGPE(u, v, x);

//   assert(gcdGPE(u, v, x) == 4 + -4 * x + -1 * pow(x, 2) + pow(x, 3));
// }

// void should_calculate_monomial_division() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = pow(x, 3) + 3 * pow(x, 2) * y + 4 * x * pow(y, 2);
//   expr v = x * y + 2 * y + 3 * pow(y, 2);

//   expr L = list({x, y});

//   assert(monomialPolyDiv(u, v, L) ==
//          list({
//              -6 + 3 * x + -5 * y,
//              pow(x, 3) + 12 * y + 28 * pow(y, 2) + 15 * pow(y, 3),
//          }));
// }

// void should_get_leading_monomial() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = 3 * pow(x, 2) * y + 4 * x * pow(y, 2) + pow(y, 3) + x + 3;

//   expr L0 = list({expr("x"), expr("y")});
//   expr L1 = list({expr("y"), expr("x")});

//   assert(leadMonomial(u, L0) == 3 * pow(x, 2) * y);
//   assert(leadMonomial(u, L1) == pow(y, 3));
// }

// void should_rec_divide_polynomials() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = pow(x, 2) * pow(y, 2) + x;
//   expr v = x * y + 1;

//   expr L0 = list({x, y});
//   expr L1 = list({y, x});

//   expr Q = expr("Q");

//   expr R0 = recPolyDiv(u, v, L0, Q);

//   assert(R0[0] == x * y);
//   assert(R0[1] == x + -1 * x * y);
// }

// void should_pseudo_divide_polynomials() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = 5 * pow(x, 4) * pow(y, 3) + 3 * x * y + 2;
//   expr v = 2 * pow(x, 3) * y + 2 * x + 3;

//   expr R = pseudoDivision(u, v, x);

//   assert(R.kind() == kind::LIST);
//   assert(R.size() == 2);

//   assert(R[0] == 10 * x * pow(y, 4));
//   assert(R[1] == 8 * pow(y, 2) + 12 * x * pow(y, 3) +
//                      -30 * x * pow(y, 4) + -20 * pow(x, 2) * pow(y, 4));
// }

// void should_normalize_polynomial() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = (2 * y + 1) * x + 6 * y + 3;

//   expr L = list({x, y});
//   expr Q = expr("Q");

//   expr u_ = normalizePoly(u, L, Q);

// 	expr u_res = fraction(3, 2) + fraction(1, 2) * x + 3 * y + x * y;

//   assert(u_ == u_res);
// }

// void should_mv_poly_gcd() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = -1 * y * pow(x, 2) + pow(y, 3);
//   expr v = y * pow(x, 2) + 2 * pow(y, 2) * x + pow(y, 3);

//   expr L = list({x, y});
//   expr Z = expr("Z");

//   expr gcd = mvPolyGCD(u, v, L, Z);

//   assert(gcd == x * y + pow(y, 2));
// }

void should_get_poly_gcd_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");
  expr z = expr("z");
  expr Z = expr("Z");

  expr L = list({x, y});

  expr u = polyExpr(-1 * y * pow(x, 2) + pow(y, 3), L);
  expr v = polyExpr(y * pow(x, 2) + 2 * pow(y, 2) * x + pow(y, 3), L);

  expr uv_gcd = gcdPolyExpr(u, v, L, Z);

	assert(uv_gcd == polyExpr(x * y + pow(y, 2), L));
}

// void should_get_poly_gcd() {
//   expr x = expr("x");
//   expr y = expr("y");
//   expr z = expr("z");
//   expr Z = expr("Z");

//   expr L = list({x, y});

//   expr u = -1 * y * pow(x, 2) + pow(y, 3);
//   expr v = y * pow(x, 2) + 2 * pow(y, 2) * x + pow(y, 3);

//   expr uv_gcd = gcdPoly(u, v, L, Z);

//   assert(uv_gcd == x * y + pow(y, 2));
// }

// void should_get_heuristic_gcd_of_polys() {
// 	expr x = expr("x");
// 	expr y = expr("y");
// 	expr z = expr("z");

// 	expr Z = expr("Z");
//   expr L = list({x, y});

//   expr u = -1 * y * pow(x, 2) + pow(y, 3);
//   expr v = y * pow(x, 2) + 2 * pow(y, 2) * x + pow(y, 3);

// 	expr U = heuristicGcdPoly(u, v, L, Z);

//   assert(U[0] == x * y + pow(y, 2));
//   assert(U[1] == y + -x);
//   assert(U[2] == y + x);

//   expr a =
//       3 * z * y * x + 2 * pow(z, 3) * pow(y, 2) * x +
//       -2 * z * pow(y, 3) * x + -3 * pow(z, 3) * x +
//       -2 * pow(z, 2) * y * pow(x, 2) + 8 * pow(z, 3) * y * pow(x, 2) +
//       2 * pow(z, 4) * y * pow(x, 2) + -8 * z * pow(y, 2) * pow(x, 2) +
//       6 * pow(z, 2) * pow(y, 2) * pow(x, 2) +
//       -6 * pow(y, 3) * pow(x, 2) + -12 * pow(z, 2) * y * pow(x, 3) +
//       12 * pow(z, 3) * y * pow(x, 3) + -9 * z * pow(y, 2) * pow(x, 3) +
//       2 * pow(z, 3) * pow(y, 2) * pow(x, 3) +
//       pow(z, 5) * pow(y, 2) * pow(x, 3) +
//       -1 * pow(z, 3) * pow(y, 3) * pow(x, 3) +
//       -2 * z * pow(y, 4) * pow(x, 3) + -3 * pow(z, 3) * pow(x, 3) +
//       12 * pow(z, 4) * pow(x, 3) + 8 * pow(z, 3) * y * pow(x, 4) +
//       2 * pow(z, 4) * y * pow(x, 4) + 4 * pow(z, 5) * y * pow(x, 4) +
//       8 * pow(z, 2) * pow(y, 2) * pow(x, 4) +
//       -4 * pow(z, 3) * pow(y, 2) * pow(x, 4) +
//       4 * pow(z, 4) * pow(y, 2) * pow(x, 4) +
//       -8 * z * pow(y, 3) * pow(x, 4) +
//       -6 * pow(z, 2) * pow(y, 3) * pow(x, 4) +
//       -8 * pow(y, 4) * pow(x, 4) + 12 * pow(z, 3) * y * pow(x, 5) +
//       -12 * pow(z, 2) * pow(y, 2) * pow(x, 5) +
//       pow(z, 5) * pow(y, 2) * pow(x, 5) +
//       -12 * z * pow(y, 3) * pow(x, 5) +
//       -1 * pow(z, 3) * pow(y, 4) * pow(x, 5) +
//       12 * pow(z, 4) * pow(x, 5) + 4 * pow(z, 5) * y * pow(x, 6) +
//       4 * pow(z, 4) * pow(y, 2) * pow(x, 6) +
//       -4 * pow(z, 3) * pow(y, 3) * pow(x, 6) +
//       -4 * pow(z, 2) * pow(y, 4) * pow(x, 6) + -2 * pow(z, 2) * y +
//       2 * pow(y, 2);

//   expr b =
//       -4 * pow(z, 2) * y * x + 16 * pow(z, 3) * y * x +
//       4 * pow(z, 4) * y * x + -16 * z * pow(y, 2) * x +
//       12 * pow(z, 2) * pow(y, 2) * x + -12 * pow(y, 3) * x +
//       -36 * pow(z, 2) * y * pow(x, 2) + 36 * pow(z, 3) * y * pow(x, 2) +
//       -27 * z * pow(y, 2) * pow(x, 2) +
//       6 * pow(z, 3) * pow(y, 2) * pow(x, 2) +
//       3 * pow(z, 5) * pow(y, 2) * pow(x, 2) +
//       -3 * pow(z, 3) * pow(y, 3) * pow(x, 2) +
//       -6 * z * pow(y, 4) * pow(x, 2) + -9 * pow(z, 3) * pow(x, 2) +
//       36 * pow(z, 4) * pow(x, 2) + 32 * pow(z, 3) * y * pow(x, 3) +
//       8 * pow(z, 4) * y * pow(x, 3) + 16 * pow(z, 5) * y * pow(x, 3) +
//       32 * pow(z, 2) * pow(y, 2) * pow(x, 3) +
//       -16 * pow(z, 3) * pow(y, 2) * pow(x, 3) +
//       16 * pow(z, 4) * pow(y, 2) * pow(x, 3) +
//       -32 * z * pow(y, 3) * pow(x, 3) +
//       -24 * pow(z, 2) * pow(y, 3) * pow(x, 3) +
//       -32 * pow(y, 4) * pow(x, 3) + 60 * pow(z, 3) * y * pow(x, 4) +
//       -60 * pow(z, 2) * pow(y, 2) * pow(x, 4) +
//       5 * pow(z, 5) * pow(y, 2) * pow(x, 4) +
//       -60 * z * pow(y, 3) * pow(x, 4) +
//       -5 * pow(z, 3) * pow(y, 4) * pow(x, 4) +
//       60 * pow(z, 4) * pow(x, 4) + 24 * pow(z, 5) * y * pow(x, 5) +
//       24 * pow(z, 4) * pow(y, 2) * pow(x, 5) +
//       -24 * pow(z, 3) * pow(y, 3) * pow(x, 5) +
//       -24 * pow(z, 2) * pow(y, 4) * pow(x, 5) + 3 * z * y +
//       2 * pow(z, 3) * pow(y, 2) + -2 * z * pow(y, 3) + -3 * pow(z, 3);

//   expr T = list({x, y, z});
// 	expr K = heuristicGcdPoly(a, b, T, Z);
// 	// printf("K[0] = %s\n", to_string(K[0]).c_str());
// 	// sort(&a);
// 	// printf("a    = %s\n", to_string(a).c_str());
// 	// printf("K[1] = %s\n", to_string(K[1]).c_str());
// 	// printf("K[2] = %s\n", to_string(K[2]).c_str());
//   assert(K[0] == 1);
//   assert(K[1] == a);
//   assert(K[2] == b);


// 	expr c = 4*pow(z, 5)*y + 4*pow(z, 4)*pow(y, 2) + -4*pow(z, 3)*pow(y, 3) + -4*pow(z, 2)*pow(y, 4);

// 	expr d = 12*pow(z, 3)*y + -12*pow(z, 2)*pow(y, 2) + pow(z, 5)*pow(y, 2) + -12*z*pow(y, 3) + -1*pow(z, 3)*pow(y, 4) + 12*pow(z, 4);

// 	expr O = list({y, z});

// 	expr J = heuristicGcdPoly(c, d, O, Z);

// 	assert(J[0] == z*pow(y, 2) + -1*pow(z, 3));
// 	assert(J[1] == -4*pow(z, 2)*y + -4*z*pow(y, 2));
// 	assert(J[2] == -12*y + -1*pow(z, 2)*pow(y, 2) + -12*z);


// 	expr e = pow(y, 2);
// 	expr f = 2*y;
// 	expr P = list({y});
// 	expr I = heuristicGcdPoly(e, f, P, Z);

// 	assert(I[0] == y);
// 	assert(I[1] == y);
// 	assert(I[2] == 2);
// }


void should_get_heuristic_gcd_of_poly_exprs() {
	expr x = expr("x");
	expr y = expr("y");
	expr z = expr("z");

	expr Z = expr("Z");

  expr T = list({x, y, z});

  expr L = list({x, y});

  expr u = polyExpr(-1 * y * pow(x, 2) + pow(y, 3), L);
  expr v = polyExpr(y * pow(x, 2) + 2 * pow(y, 2) * x + pow(y, 3), L);

	expr U = heuristicGcdPolyExpr(u, v, L, Z);

  assert(U[0] == polyExpr(x * y + pow(y, 2), L));
  assert(U[1] == polyExpr(y + -x, L));
  assert(U[2] == polyExpr(y + x, L));

	expr a = polyExpr(
      3 * z * y * x + 2 * pow(z, 3) * pow(y, 2) * x +
      -2 * z * pow(y, 3) * x + -3 * pow(z, 3) * x +
      -2 * pow(z, 2) * y * pow(x, 2) + 8 * pow(z, 3) * y * pow(x, 2) +
      2 * pow(z, 4) * y * pow(x, 2) + -8 * z * pow(y, 2) * pow(x, 2) +
      6 * pow(z, 2) * pow(y, 2) * pow(x, 2) +
      -6 * pow(y, 3) * pow(x, 2) + -12 * pow(z, 2) * y * pow(x, 3) +
      12 * pow(z, 3) * y * pow(x, 3) + -9 * z * pow(y, 2) * pow(x, 3) +
      2 * pow(z, 3) * pow(y, 2) * pow(x, 3) +
      pow(z, 5) * pow(y, 2) * pow(x, 3) +
      -1 * pow(z, 3) * pow(y, 3) * pow(x, 3) +
      -2 * z * pow(y, 4) * pow(x, 3) + -3 * pow(z, 3) * pow(x, 3) +
      12 * pow(z, 4) * pow(x, 3) + 8 * pow(z, 3) * y * pow(x, 4) +
      2 * pow(z, 4) * y * pow(x, 4) + 4 * pow(z, 5) * y * pow(x, 4) +
      8 * pow(z, 2) * pow(y, 2) * pow(x, 4) +
      -4 * pow(z, 3) * pow(y, 2) * pow(x, 4) +
      4 * pow(z, 4) * pow(y, 2) * pow(x, 4) +
      -8 * z * pow(y, 3) * pow(x, 4) +
      -6 * pow(z, 2) * pow(y, 3) * pow(x, 4) +
      -8 * pow(y, 4) * pow(x, 4) + 12 * pow(z, 3) * y * pow(x, 5) +
      -12 * pow(z, 2) * pow(y, 2) * pow(x, 5) +
      pow(z, 5) * pow(y, 2) * pow(x, 5) +
      -12 * z * pow(y, 3) * pow(x, 5) +
      -1 * pow(z, 3) * pow(y, 4) * pow(x, 5) +
      12 * pow(z, 4) * pow(x, 5) + 4 * pow(z, 5) * y * pow(x, 6) +
      4 * pow(z, 4) * pow(y, 2) * pow(x, 6) +
      -4 * pow(z, 3) * pow(y, 3) * pow(x, 6) +
      -4 * pow(z, 2) * pow(y, 4) * pow(x, 6) + -2 * pow(z, 2) * y +
      2 * pow(y, 2), T);

  expr b = polyExpr(
      -4 * pow(z, 2) * y * x + 16 * pow(z, 3) * y * x +
      4 * pow(z, 4) * y * x + -16 * z * pow(y, 2) * x +
      12 * pow(z, 2) * pow(y, 2) * x + -12 * pow(y, 3) * x +
      -36 * pow(z, 2) * y * pow(x, 2) + 36 * pow(z, 3) * y * pow(x, 2) +
      -27 * z * pow(y, 2) * pow(x, 2) +
      6 * pow(z, 3) * pow(y, 2) * pow(x, 2) +
      3 * pow(z, 5) * pow(y, 2) * pow(x, 2) +
      -3 * pow(z, 3) * pow(y, 3) * pow(x, 2) +
      -6 * z * pow(y, 4) * pow(x, 2) + -9 * pow(z, 3) * pow(x, 2) +
      36 * pow(z, 4) * pow(x, 2) + 32 * pow(z, 3) * y * pow(x, 3) +
      8 * pow(z, 4) * y * pow(x, 3) + 16 * pow(z, 5) * y * pow(x, 3) +
      32 * pow(z, 2) * pow(y, 2) * pow(x, 3) +
      -16 * pow(z, 3) * pow(y, 2) * pow(x, 3) +
      16 * pow(z, 4) * pow(y, 2) * pow(x, 3) +
      -32 * z * pow(y, 3) * pow(x, 3) +
      -24 * pow(z, 2) * pow(y, 3) * pow(x, 3) +
      -32 * pow(y, 4) * pow(x, 3) + 60 * pow(z, 3) * y * pow(x, 4) +
      -60 * pow(z, 2) * pow(y, 2) * pow(x, 4) +
      5 * pow(z, 5) * pow(y, 2) * pow(x, 4) +
      -60 * z * pow(y, 3) * pow(x, 4) +
      -5 * pow(z, 3) * pow(y, 4) * pow(x, 4) +
      60 * pow(z, 4) * pow(x, 4) + 24 * pow(z, 5) * y * pow(x, 5) +
      24 * pow(z, 4) * pow(y, 2) * pow(x, 5) +
      -24 * pow(z, 3) * pow(y, 3) * pow(x, 5) +
      -24 * pow(z, 2) * pow(y, 4) * pow(x, 5) + 3 * z * y +
      2 * pow(z, 3) * pow(y, 2) + -2 * z * pow(y, 3) + -3 * pow(z, 3), T);

	expr K = heuristicGcdPolyExpr(a, b, T, Z);

  assert(K[0] == polyExpr(1, T));
  assert(K[1] == a);
  assert(K[2] == b);

	expr O = list({y, z});

	expr c = polyExpr(4*pow(z, 5)*y + 4*pow(z, 4)*pow(y, 2) + -4*pow(z, 3)*pow(y, 3) + -4*pow(z, 2)*pow(y, 4), O);

	expr d = polyExpr(12*pow(z, 3)*y + -12*pow(z, 2)*pow(y, 2) + pow(z, 5)*pow(y, 2) + -12*z*pow(y, 3) + -1*pow(z, 3)*pow(y, 4) + 12*pow(z, 4), O);

	expr J = heuristicGcdPolyExpr(c, d, O, Z);

	assert(J[0] == polyExpr(z*pow(y, 2) + -1*pow(z, 3), O));
	assert(J[1] == polyExpr(-4*pow(z, 2)*y + -4*z*pow(y, 2), O));
	assert(J[2] == polyExpr(-12*y + -1*pow(z, 2)*pow(y, 2) + -12*z, O));


	expr P = list({y});
	expr e = polyExpr(pow(y, 2), P);
	expr f = polyExpr(2*y, P);
	expr I = heuristicGcdPolyExpr(e, f, P, Z);

	assert(I[0] == polyExpr(y, P));
	assert(I[1] == polyExpr(y, P));
	assert(I[2] == polyExpr(2, P));
}



// void should_get_coeff_var_parts_of_monomial() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = expr(4) * expr(5) * fraction(1, 2) * x * pow(x, 2) * pow(y, 3);

//   set S = set({x, y});

//   expr L = coeffVarMonomial(u, S);

//   assert(L.kind() == kind::LIST);
//   assert(L.size() == 2);

//   assert(L[0] == expr(4) * expr(5) * fraction(1, 2));
//   assert(L[1] == x * pow(x, 2) * pow(y, 3));
// }

// void should_collect_terms() {
//   expr x = expr("x");
//   expr y = expr("y");
//   expr a = expr("a");
//   expr b = expr("b");
//   expr c = expr("c");
//   expr d = expr("d");

//   expr u = a * x + b * x + c + d;

// 	set S = set({x});

//   expr C = collectTerms(u, S);

//   assert(C == (a + b) * x + (c + d));
// }

void should_collect_polynomials() {
  expr x = expr("x");
  expr y = expr("y");
  expr z = expr("z");
  expr u = expr("u");

  expr p0 = 4 * pow(x, 2) * pow(y, 5) * pow(z, 3) +
            10 * pow(x, 3) * y * z + 2 * pow(y, 3) * pow(z, 2) +
            pow(x, 4) * pow(z, 3) + 11 * pow(x, 4) * pow(y, 3) +
            pow(z, 4) * x + y * pow(z, 3);

  assert(polyExpr(p0, list({y, x, z})) ==

         create(kind::ADD, {create(kind::ADD, {create(kind::ADD, {1 * pow(z, 4)}) * pow(x, 1),
                   create(kind::ADD, {1 * pow(z, 3)}) * pow(x, 4)}) *
                  pow(y, 0),
              create(kind::ADD, {create(kind::ADD, {1 * pow(z, 3)}) * pow(x, 0),
                   create(kind::ADD, {10 * pow(z, 1)}) * pow(x, 3)}) *
                  pow(y, 1),
              create(kind::ADD, {create(kind::ADD, {2 * pow(z, 2)}) * pow(x, 0),
                   create(kind::ADD, {11 * pow(z, 0)}) * pow(x, 4)}) *
                  pow(y, 3),
              create(kind::ADD, {create(kind::ADD, {4 * pow(z, 3)}) * pow(x, 2)}) * pow(y, 5)}));
}


// void should_expand_main_operator() {
//   expr x = expr("x");

//   expr u0 = x * (2 + pow(1 + x, 2));

//   assert(expandRoot(u0) == 2 * x + x * pow(1 + x, 2));

//   expr u1 = pow(x + pow(1 + x, 2), 2);

//   assert(expandRoot(u1) ==
//          pow(x, 2) + 2 * x * pow(1 + x, 2) + pow(1 + x, 4));
// }

// void should_get_polynomial_content() {
//   expr x = expr("x");
//   expr y = expr("y");
//   expr z = expr("z");

//   expr u = 4 * pow(x, 2) + -1 * 6 * x;

//   expr L = list({x});

//   expr Z = expr("Z");
//   expr Q = expr("Q");

// 	assert(cont(u, L, Z) == 2);

//   expr t = 2 * x;

//   assert(cont(t, L, Z) == 2);

//   expr p = -1 * x;

//   assert(cont(p, L, Z) == 1);

//   expr T = list({x, y});

//   expr a = fraction(1, 2) * x * y + 6 * y;
// 	assert(cont(a, T, Q) == y);

// 	expr b = (pow(y, 2) + 2 * y + 1) * pow(x, 2) + (2 * pow(y, 2) - 2) * x +
//             (3 * y + 3);

//   assert(cont(b, T, Q) == 1 + y);
// }

// void should_get_polynomial_content_sub_resultant() {
//   expr x = expr("x");
//   expr y = expr("y");

//   expr u = 4 * pow(x, 2) + -1 * 6 * x;

//   expr L = list({x});

//   expr Z = expr("Z");
//   expr Q = expr("Q");

//   assert(cont(u, L, Z) == 2);

//   expr t = 2 * x;

//   assert(cont(t, L, Z) == 2);

//   expr p = -1 * x;

//   assert(cont(p, L, Z) == 1);

//   expr T = list({x, y});

//   expr a = fraction(1, 2) * x * y + 6 * y;

//   assert(cont(a, T, Q) == y);

//   expr b = (pow(y, 2) + 2 * y + 1) * pow(x, 2) + (2 * pow(y, 2) - 2) * x +
//            (3 * y + 3);

//   assert(cont(b, T, Q) == 1 + y);
// }

// void should_monomial_base_expand_polynomials() {
//   expr x = expr("x");
//   expr a = expr("a");
//   expr b = expr("b");
//   expr t = expr("t");

//   expr u =
//       pow(a, 2) * b + 2 * a * pow(b, 2) + pow(b, 3) + 2 * a + 2 * b + 3;

//   expr v = a + b;

//   expr L = list({a, b});

//   expr r = monomialBasedPolyExpansion(u, v, L, t);

//   assert(r == 3 + 2 * t + b * pow(t, 2));
// }

void should_get_if_poly_expr_is_zero() {
  expr x = expr("x");
  expr y = expr("y");

  expr u = 0;

  assert(isZeroPolyExpr(u) == true);

  expr g = 0 * pow(x, 4);

  assert(isZeroPolyExpr(g) == true);

  expr t = polyExpr(0 * pow(x, 3) * y, list({x, y}));

  assert(isZeroPolyExpr(t) == true);

  expr k = polyExpr(0, list({x, y}));

  assert(isZeroPolyExpr(k) == true);
}

void should_mul_collected_polys() {
  expr x = expr("x");
  expr y = expr("y");

  expr u = 2 * x + 4 * pow(x, 2) * pow(y, 4) + 10 * pow(y, 2) +
           8 * pow(x, 6) + 7 * pow(y, 3) + pow(x, 3);

  expr v = 5 * x + 7 * pow(x, 3) * pow(y, 2) + 11 * y + 8 * pow(x, 4) +
           2 * pow(y, 3) + pow(x, 5);

  expr L = list({x, y});

  expr uv = mulPolyExpr(polyExpr(u, L), polyExpr(v, L));

  assert(uv ==
         create(kind::ADD, {create(kind::ADD, {110 * pow(y, 3), 77 * pow(y, 4), 20 * pow(y, 5),
                   14 * pow(y, 6)}) *
                  pow(x, 0),
              create(kind::ADD, {22 * pow(y, 1), 50 * pow(y, 2), 39 * pow(y, 3)}) *
                  pow(x, 1),
              create(kind::ADD, {10 * pow(y, 0), 44 * pow(y, 5), 8 * pow(y, 7)}) *
                  pow(x, 2),
              create(kind::ADD, {11 * pow(y, 1), 2 * pow(y, 3), 90 * pow(y, 4),
                   49 * pow(y, 5)}) *
                  pow(x, 3),
              create(kind::ADD, {5 * pow(y, 0), 94 * pow(y, 2), 56 * pow(y, 3)}) *
                  pow(x, 4),
              create(kind::ADD, {16 * pow(y, 0), 10 * pow(y, 2), 7 * pow(y, 3),
                   28 * pow(y, 6)}) *
                  pow(x, 5),
              create(kind::ADD, {2 * pow(y, 0), 88 * pow(y, 1), 7 * pow(y, 2),
                   16 * pow(y, 3), 32 * pow(y, 4)}) *
                  pow(x, 6),
              create(kind::ADD, {48 * pow(y, 0), 4 * pow(y, 4)}) * pow(x, 7),
              create(kind::ADD, {1 * pow(y, 0)}) * pow(x, 8),
              create(kind::ADD, {56 * pow(y, 2)}) * pow(x, 9),
              create(kind::ADD, {64 * pow(y, 0)}) * pow(x, 10),
              create(kind::ADD, {8 * pow(y, 0)}) * pow(x, 11)}));

  expr g = 3;
  expr t = 5;

  assert(mulPolyExpr(polyExpr(g, L), polyExpr(t, L)) == polyExpr(15, L));
  assert(mulPolyExpr(polyExpr(g, list({})), polyExpr(t, list({}))) ==
         polyExpr(15, list({})));

  expr k = 4 * x + 15 * pow(x, 2) + 4 * y * x + 5 * y;
  assert(mulPolyExpr(polyExpr(k, L), 5) ==
         create(kind::ADD, {
             create(kind::ADD, {25 * pow(y, 1)}) * pow(x, 0),
             create(kind::ADD, {20 * pow(y, 0), 20 * pow(y, 1)}) * pow(x, 1),
             create(kind::ADD, {75 * pow(y, 0)}) * pow(x, 2),
         }));

  assert(mulPolyExpr(polyExpr(3, L), polyExpr(4, L)) == polyExpr(12, L));
  assert(mulPolyExpr(polyExpr(4, list({})), polyExpr(3, list({}))) ==
         polyExpr(12, list({})));

  expr r = 4 * pow(x, 3) + 5 * pow(x, 2) + 3 * x + 4;
  expr h = 2 * pow(x, 2) + 3 * x + 7;

  assert(mulPolyExpr(polyExpr(r, list({x})), polyExpr(h, list({x}))) ==
         28 * pow(x, 0) + 33 * pow(x, 1) + 52 * pow(x, 2) +
             49 * pow(x, 3) + 22 * pow(x, 4) + 8 * pow(x, 5));
}

void should_add_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr u = 2 * pow(x, 2) * y + 3 * pow(y, 2) * x + 4 * x;
  expr v = 3 * pow(x, 2) + 4 * pow(y, 2) * x + 2 * x;

  expr L = list({x, y});

  expr r = addPolyExpr(polyExpr(u, L), polyExpr(v, L));

	assert(r == create(kind::ADD, {create(kind::ADD, {6 * pow(y, 0), 7 * pow(y, 2)}) * pow(x, 1),
                   create(kind::ADD, {3 * pow(y, 0), 2 * pow(y, 1)}) * pow(x, 2)}));

  expr g = 2 * pow(y, 3) * pow(x, 5) + 4 * pow(x, 4) + 4 * y + 4 * x;
  expr t = 5 * pow(y, 3) * pow(x, 5) + 3 * x + y;

  expr gt = addPolyExpr(polyExpr(g, L), polyExpr(t, L));

  assert(gt == create(kind::ADD, {create(kind::ADD, {5 * pow(y, 1)}) * pow(x, 0),
                    create(kind::ADD, {7 * pow(y, 0)}) * pow(x, 1),
                    create(kind::ADD, {4 * pow(y, 0)}) * pow(x, 4),
                    create(kind::ADD, {7 * pow(y, 3)}) * pow(x, 5)}));

  expr e = 4 * pow(x, 3) + 5 * pow(x, 2) + 3 * x + 4;
  expr h = 2 * pow(x, 2) + 3 * x + 7;
  expr q = 2 * pow(x, 2) + -3 * x + 7;

  assert(addPolyExpr(polyExpr(e, list({x})), polyExpr(h, list({x}))) ==
         11 * pow(x, 0) + 6 * pow(x, 1) + 7 * pow(x, 2) +
             4 * pow(x, 3));

  assert(addPolyExpr(polyExpr(h, list({x})), polyExpr(q, list({x}))) ==
         4 * pow(x, 2) + 14 * pow(x, 0));

  assert(addPolyExpr(polyExpr(e, list({x})), polyExpr(0, list({x}))) ==
         polyExpr(e, list({x})));
  assert(addPolyExpr(polyExpr(0, list({x})), polyExpr(e, list({x}))) ==
         polyExpr(e, list({x})));
}

void should_sub_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr u = 2 * y * pow(x, 2) + 3 * pow(y, 2) * x + 4 * x;
  expr v = 3 * pow(x, 2) + 4 * pow(y, 2) * x + 2 * x;

  expr L = list({x, y});

  assert(subPolyExpr(polyExpr(u, L), polyExpr(v, L)) ==
         create(kind::ADD, {
             create(kind::ADD, {2 * pow(y, 0), -1 * pow(y, 2)}) * pow(x, 1),
             create(kind::ADD, {-3 * pow(y, 0), 2 * pow(y, 1)}) * pow(x, 2),
         }));

  assert(subPolyExpr(polyExpr(3, L), polyExpr(4, L)) == polyExpr(-1, L));
  assert(subPolyExpr(polyExpr(4, list({})), polyExpr(3, list({}))) ==
         polyExpr(1, list({})));

  expr e = 4 * pow(x, 3) + 5 * pow(x, 2) + 3 * x + 4;
  expr h = 2 * pow(x, 2) + 3 * x + 7;

  assert(subPolyExpr(polyExpr(e, list({x})), polyExpr(h, list({x}))) ==
         -3 * pow(x, 0) + 3 * pow(x, 2) + 4 * pow(x, 3));

  assert(subPolyExpr(polyExpr(e, list({x})), polyExpr(0, list({x}))) ==
         polyExpr(e, list({x})));
  assert(subPolyExpr(polyExpr(0, list({x})), polyExpr(e, list({x}))) ==
         mulPolyExpr(polyExpr(e, list({x})), -1));

  expr a = create(kind::ADD, {create(kind::ADD, {-1 * pow(y, 0), 1 * pow(y, 2)}) * pow(x, 2)});
  expr b = create(kind::ADD, {create(kind::ADD, {-1 * pow(y, 0), 1 * pow(y, 2)}) * pow(x, 2),
                create(kind::ADD, {-1 * pow(y, 2), 1 * pow(y, 3)}) * pow(x, 3)});

  assert(subPolyExpr(a, b) ==
         create(kind::ADD, {create(kind::ADD, {1 * pow(y, 2), -1 * pow(y, 3)}) * pow(x, 3)}));

  expr c = -1 * pow(y, 0) + -1 * pow(y, 1);
  expr d = create(kind::ADD, {-1 * pow(y, 1)});

  assert(subPolyExpr(c, d) == create(kind::ADD, {-1 * pow(y, 0)}));
}

void should_rec_divide_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr u = pow(x, 2) * pow(y, 2) + x;
  expr v = x * y + 1;

  expr L0 = list({x, y});
  // expr L1 = list({y, x});

  expr Q = expr("Q");

  expr R0 = divPolyExpr(polyExpr(u, L0), polyExpr(v, L0), L0, Q);

  assert(R0 ==
         list({
             create(kind::ADD, {create(kind::ADD, {1 * pow(y, 1)}) * pow(x, 1)}),
             create(kind::ADD, {create(kind::ADD, {1 * pow(y, 0), -1 * pow(y, 1)}) * pow(x, 1)}),
         }));
}

void should_pseudo_divide_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr u = 5 * pow(x, 4) * pow(y, 3) + 3 * x * y + 2;
  expr v = 2 * pow(x, 3) * y + 2 * x + 3;

  expr L = list({x, y});

  expr R = pseudoDivPolyExpr(polyExpr(u, L), polyExpr(v, L), L);

  assert(R ==
         list({create(kind::ADD, {create(kind::ADD, {10 * pow(y, 4)}) * pow(x, 1)}),
               create(kind::ADD, {
                   create(kind::ADD, {8 * pow(y, 2)}) * pow(x, 0),
                   create(kind::ADD, {12 * pow(y, 3), -30 * pow(y, 4)}) * pow(x, 1),
                   create(kind::ADD, {-20 * pow(y, 4)}) * pow(x, 2),
               })}));
}

void should_pow_col_poly() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x});

  expr u = polyExpr(2 * pow(x, 2) + 4 * x + 12, L);

  assert(powPolyExpr(u, 2) == create(kind::ADD, {
                                     144 * pow(x, 0),
                                     96 * pow(x, 1),
                                     64 * pow(x, 2),
                                     16 * pow(x, 3),
                                     4 * pow(x, 4),
                                 }));

  assert(powPolyExpr(u, 5) == create(kind::ADD, {
                                     248832 * pow(x, 0),
                                     414720 * pow(x, 1),
                                     483840 * pow(x, 2),
                                     368640 * pow(x, 3),
                                     222720 * pow(x, 4),
                                     100864 * pow(x, 5),
                                     37120 * pow(x, 6),
                                     10240 * pow(x, 7),
                                     2240 * pow(x, 8),
                                     320 * pow(x, 9),
                                     32 * pow(x, 10),
                                 }));

  expr T = list({x, y});

  expr v = polyExpr(2 * pow(x, 3) + 4 * pow(y, 2) * x + 12 * y + 14, T);
  assert(
      powPolyExpr(v, 3) ==
      create(kind::ADD, {create(kind::ADD, {2744 * pow(y, 0), 7056 * pow(y, 1), 6048 * pow(y, 2),
                1728 * pow(y, 3)}) *
               pow(x, 0),
           create(kind::ADD, {2352 * pow(y, 2), 4032 * pow(y, 3), 1728 * pow(y, 4)}) *
               pow(x, 1),
           create(kind::ADD, {672 * pow(y, 4), 576 * pow(y, 5)}) * pow(x, 2),
           create(kind::ADD, {1176 * pow(y, 0), 2016 * pow(y, 1), 864 * pow(y, 2),
                64 * pow(y, 6)}) *
               pow(x, 3),
           create(kind::ADD, {672 * pow(y, 2), 576 * pow(y, 3)}) * pow(x, 4),
           create(kind::ADD, {96 * pow(y, 4)}) * pow(x, 5),
           create(kind::ADD, {168 * pow(y, 0), 144 * pow(y, 1)}) * pow(x, 6),
           create(kind::ADD, {48 * pow(y, 2)}) * pow(x, 7),
           create(kind::ADD, {8 * pow(y, 0)}) * pow(x, 9)}));
}

void should_normalize_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x, y});
  expr Q = expr("Q");

  expr u = polyExpr(2 * x * y * x + x + 6 * y + 3, L);

  assert(normalizePolyExpr(u, L, Q) ==
         create(kind::ADD, {
             create(kind::ADD, {fraction(3, 2) * pow(y, 0), 3 * pow(y, 1)}) * pow(x, 0),
             create(kind::ADD, {fraction(1, 2) * pow(y, 0), 1 * pow(y, 1)}) * pow(x, 1),
         }));
}

void should_get_gcd_of_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x, y});

  expr u = polyExpr(-1 * y * pow(x, 2) + pow(y, 3), L);

	expr v = polyExpr(y * pow(x, 2) + 2 * pow(y, 2) * x + pow(y, 3), L);

  expr Z = expr("Z");

  expr gcd = gcdPolyExpr(u, v, L, Z);

  assert(gcd == create(kind::ADD, {
                    create(kind::ADD, {1 * pow(y, 2)}) * pow(x, 0),
                    create(kind::ADD, {1 * pow(y, 1)}) * pow(x, 1),
                }));
}

void should_get_content_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x});

  expr Z = expr("Z");
  expr Q = expr("Q");

  expr u = polyExpr(4 * pow(x, 2) + -1 * 6 * x, L);

  assert(contPolyExpr(u, L, Z) == 2);

  expr t = polyExpr(2 * x, L);

  assert(contPolyExpr(t, L, Z) == 2);

  expr p = polyExpr(-1 * x, L);

  assert(contPolyExpr(p, L, Z) == 1);

  expr T = list({x, y});

  expr a = polyExpr(fraction(1, 2) * x * y + 6 * y, T);

  assert(contPolyExpr(a, T, Q) == create(kind::ADD, {1 * pow(y, 1)}));
  expr b = polyExpr(pow(y, 2) * pow(x, 2) + 2 * y * pow(x, 2) +
                        pow(x, 2) + 2 * pow(y, 2) * x + -2 * x + 3 * y + 3,
                    T);

  assert(contPolyExpr(b, T, Z) == 1 * pow(y, 0) + 1 * pow(y, 1));
}

void should_expand_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x});
  expr T = list({x, y});

  expr u = 10 * pow(x, 4) + 4 * pow(x, 2) + 10;

  assert(expandPolyExpr(polyExpr(u, L)) == u);

  expr v = 10*pow(x, 4)*pow(y, 2) + 6*pow(x, 2) + 10*x +
           4*pow(x, 4)*y + 15*y + 10*x*y;

	assert(expandPolyExpr(polyExpr(v, T)) == v);
}

void should_diff_poly_expr() {
  expr x = expr("x");
  expr y = expr("y");

  expr L = list({x});
  expr T = list({x, y});

  expr u = polyExpr(10 * pow(x, 4) + 4 * pow(x, 2) + 10, L);

  assert(diffPolyExpr(u, x) == 8 * pow(x, 1) + 40 * pow(x, 3));

  assert(diffPolyExpr(u, y) == expr(kind::ADD, {0 * pow(x, 0)}));

  expr v = polyExpr(10 * pow(x, 4) * pow(y, 2) + 6 * pow(x, 2) * y +
                        10 * x + 4 * pow(x, 4) * y + 15 * y + 10 * x * y,
                    T);

  assert(diffPolyExpr(v, x) ==
         create(kind::ADD, {
             create(kind::ADD, {10 * pow(y, 0), 10 * pow(y, 1)}) * pow(x, 0),
             create(kind::ADD, {12 * pow(y, 1)}) * pow(x, 1),
             create(kind::ADD, {16 * pow(y, 1), 40 * pow(y, 2)}) * pow(x, 3),
         }));
  assert(diffPolyExpr(v, y) ==
         create(kind::ADD, {
             create(kind::ADD, {15 * pow(y, 0)}) * pow(x, 0),
             create(kind::ADD, {10 * pow(y, 0)}) * pow(x, 1),
             create(kind::ADD, {6 * pow(y, 0)}) * pow(x, 2),
             create(kind::ADD, {4 * pow(y, 0), 20 * pow(y, 1)}) * pow(x, 4),
         }));
}

// void should_remove_denominators_from_polys() {
// 	expr x = expr("x");
// 	expr y = expr("y");
// 	expr Z = expr("Z");

// 	expr L = list({x, y});

// 	expr u = fraction(1, 2)*x + fraction(1, 3)*y + 1;

// 	expr v = removeDenominatorsPoly(u, L, Z);

// 	assert(v == list({6, 3*x + 2*y + 6}));
// }

void should_remove_denominators_from_poly_expr() {
	expr x = expr("x");
	expr y = expr("y");
	expr Z = expr("Z");

	expr L = list({x, y});

	expr u = polyExpr(fraction(1, 2)*x + fraction(1, 3)*y + 1, L);

	expr v = removeDenominatorsPolyExpr(u, L, Z);

	assert(v == list({polyExpr(6, L), polyExpr(3*x + 2*y + 6, L)}));
}

void should_get_list_of_variables() {
	expr x = expr("x");
	expr y = expr("y");
	expr z = expr("z");

	expr t = 5*x + 26*pow(x, 4)*z*pow(y, 2);

	list l = getVariableListForPolyExpr(t);

	assert(l == list({x, y, z}));
}

int main() {
  // TEST(should_get_polynomial_variable)
  // TEST(should_get_degree_of_variables)
  // TEST(should_get_coefficients)
  // TEST(should_get_leading_coefficient)
  // TEST(should_divided_polynomials)
  // TEST(should_get_gcd_polynomials)
  // TEST(should_get_leading_monomial)
  // TEST(should_calculate_monomial_division)
  // TEST(should_rec_divide_polynomials)
  // TEST(should_pseudo_divide_polynomials)
  // TEST(should_normalize_polynomial)
  // TEST(should_get_coeff_var_parts_of_monomial)
  // TEST(should_get_coeff_var_parts_of_monomial)
	// TEST(should_collect_terms)
	// TEST(should_remove_denominators_from_polys);
  // TEST(should_get_polynomial_content)
	// TEST(should_get_polynomial_content_sub_resultant)
	// TEST(should_monomial_base_expand_polynomials)
  // TEST(should_get_poly_gcd)
	// TEST(should_get_heuristic_gcd_of_polys)
	TEST(should_remove_denominators_from_poly_expr);
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
	TEST(should_get_heuristic_gcd_of_poly_exprs)
	TEST(should_get_list_of_variables)
  return 0;
}
