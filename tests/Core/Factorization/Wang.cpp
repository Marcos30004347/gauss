#include "Core/Algebra/Expression.hpp"
#include "test.hpp"

#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Factorization/Wang.hpp"
#include <iterator>

using namespace alg;
using namespace polynomial;
using namespace factorization;

void should_get_nondivisors()
{
	expr y = expr("y");
	expr z = expr("z");
	expr K = expr("Z");

	expr F = list({-14, 3, -11, -17 });

	expr L = list({ y, z});

	expr d = nondivisors(4, F, 1, L, K);

	assert(d[0] == 7);
	assert(d[1] == 3);
	assert(d[2] == 11);
	assert(d[3] == 17);
}

void should_get_nondivisors_poly_expr()
{
	expr y = expr("y");
	expr z = expr("z");
	expr K = expr("Z");

	expr F = list({ -14, 3, -11, -17 });

	expr L = list({ y, z});

	expr d = nondivisorsPolyExpr(4, F, 1, L, K);

	assert(d[0] == 7);
	assert(d[1] == 3);
	assert(d[2] == 11);
	assert(d[3] == 17);
}




void should_solve_diophant()
{
	expr x = expr("x");
	expr y = expr("y");

	expr H1 = list({
			44*pow(x, 2) + 42*x + 1,
			126*pow(x, 2) + -9*x + 28,
			187*pow(x, 2) + -23
		});

	expr H2 = list({
			-4*pow(x, 2)*y + -12*pow(x, 2) + -3*x*y + 1,
			-9*pow(x, 2)*y + -9*x + -2*y,
			pow(x, 2)*pow(y, 2) + -9*pow(x, 2) + y + -9
		});

	expr H3 = list({
			-4*pow(x, 2)*y + -12*pow(x, 2) + -3*x*y + 1,
			-9*pow(x, 2)*y + -9*x + -2*y,
			pow(x, 2)*pow(y, 2) + -9*pow(x, 2) + y + -9
		});

	expr c1 = -70686*pow(x, 5) + -5863*pow(x, 4) + -17826*pow(x, 3) + 2009*pow(x, 2) + 5031*x + 74;

	expr c2 =
		9 * pow(x, 5) * pow(y, 4) + 12 * pow(x, 5) * pow(y, 3) +
		-45 * pow(x, 5) * pow(y, 2) + -108 * pow(x, 5) * y +
		-324 * pow(x, 5) + 18 * pow(x, 4) * pow(y, 3) +
		-216 * pow(x, 4) * pow(y, 2) + -810 * pow(x, 4) * y +
		2 * pow(x, 3) * pow(y, 4) + 9 * pow(x, 3) * pow(y, 3) +
		-252 * pow(x, 3) * pow(y, 2) + -288 * pow(x, 3) * y +
		-945 * pow(x, 3) + -30 * pow(x, 2) * pow(y, 2) +
		-414 * pow(x, 2) * y + 2 * x * pow(y, 3) +
		-54 * x * pow(y, 2) + -3 * x * y + 81 * x + 12 * y;

	expr c3 = -36*pow(x, 4)*pow(y, 2) + -108*pow(x, 4)*y + -27*pow(x, 3)*pow(y, 2) + -36*pow(x, 3)*y + -108*pow(x, 3) + -8*pow(x, 2)*pow(y, 2) + -42*pow(x, 2)*y + -6*x*pow(y, 2)+ 9*x + 2*y;

	expr L1 = list({ x });
	expr I1 = list({});

 	expr D1 = multivariateDiophant(H1, c1, L1, I1, 5, 6291469, 1);

	expr R1 = list({ -3*x, -2, 1 });

	assert(D1 == R1);

	expr L2 = list({ x, y });
	expr I2 = list({ -14 });

 	expr D2 = multivariateDiophant(H2, c2, L2, I2, 5, 6291469, 1);
	expr R2 = list({ -1*x*y, -3*x, -6 });

	assert(D2 == R2);

 	expr D3 = multivariateDiophant(H3, c3, L2, I2, 5, 6291469, 1);

	expr R3 = list({ 0, 0, -1 });

	assert(D3 == R3);
}





void should_solve_diophant_poly_expr()
{
	expr x = expr("x");
	expr y = expr("y");
	expr Z = expr("Z");

	expr L1 = list({ x });
	expr L2 = list({ x, y });

	expr H1 = list({
			polyExpr(44*pow(x, 2) + 42*x + 1, L1),
			polyExpr(126*pow(x, 2) + -9*x + 28, L1),
			polyExpr(187*pow(x, 2) + -23, L1)
		});

	expr H2 = list({
			polyExpr(-4*pow(x, 2)*y + -12*pow(x, 2) + -3*x*y + 1, L2),
			polyExpr(-9*pow(x, 2)*y + -9*x + -2*y, L2),
			polyExpr(pow(x, 2)*pow(y, 2) + -9*pow(x, 2) + y + -9, L2)
		});

	expr H3 = list({
			polyExpr(-4*pow(x, 2)*y + -12*pow(x, 2) + -3*x*y + 1, L2),
			polyExpr(-9*pow(x, 2)*y + -9*x + -2*y, L2),
			polyExpr(pow(x, 2)*pow(y, 2) + -9*pow(x, 2) + y + -9, L2)
		});

	expr c1 = polyExpr(-70686*pow(x, 5) + -5863*pow(x, 4) + -17826*pow(x, 3) + 2009*pow(x, 2) + 5031*x + 74, L1);

	expr c2 =
		polyExpr(
		9 * pow(x, 5) * pow(y, 4) + 12 * pow(x, 5) * pow(y, 3) +
		-45 * pow(x, 5) * pow(y, 2) + -108 * pow(x, 5) * y +
		-324 * pow(x, 5) + 18 * pow(x, 4) * pow(y, 3) +
		-216 * pow(x, 4) * pow(y, 2) + -810 * pow(x, 4) * y +
		2 * pow(x, 3) * pow(y, 4) + 9 * pow(x, 3) * pow(y, 3) +
		-252 * pow(x, 3) * pow(y, 2) + -288 * pow(x, 3) * y +
		-945 * pow(x, 3) + -30 * pow(x, 2) * pow(y, 2) +
		-414 * pow(x, 2) * y + 2 * x * pow(y, 3) +
		-54 * x * pow(y, 2) + -3 * x * y + 81 * x + 12 * y, L2);

	expr c3 = polyExpr(-36*pow(x, 4)*pow(y, 2) + -108*pow(x, 4)*y + -27*pow(x, 3)*pow(y, 2) + -36*pow(x, 3)*y + -108*pow(x, 3) + -8*pow(x, 2)*pow(y, 2) + -42*pow(x, 2)*y + -6*x*pow(y, 2)+ 9*x + 2*y, L2);

	expr I1 = list({});

	expr D1 = multivariateDiophantPolyExpr(H1, c1, L1, I1, 5, 6291469, 1, Z);

	expr R1 = list({ polyExpr(-3*x, L1), polyExpr(-2, L1), polyExpr(1, L1) });

	assert(D1 == R1);

	expr I2 = list({ -14 });

	expr D2 = multivariateDiophantPolyExpr(H2, c2, L2, I2, 5, 6291469, 1, Z);

	expr R2 = list({ polyExpr(-1*x*y, L2), polyExpr(-3*x, L2), polyExpr(-6, L2) });

	assert(D2 == R2);

 	expr D3 = multivariateDiophantPolyExpr(H3, c3, L2, I2, 5, 6291469, 1, Z);

	expr R3 = list({ polyExpr(0, L2), polyExpr(0, L2), polyExpr(-1, L2) });

	assert(D3 == R3);
}



void should_get_lead_coeffs()
{
  expr x = expr("x");
  expr y = expr("y");
  expr z = expr("z");

  expr f =
      4 * pow(x, 6) * pow(y, 4) * pow(z, 2) +
      4 * pow(x, 6) * pow(y, 3) * pow(z, 3) +
      -4 * pow(x, 6) * pow(y, 2) * pow(z, 4) +
      -4 * pow(x, 6) * y * pow(z, 5) +
      pow(x, 5) * pow(y, 4) * pow(z, 3) +
      12 * pow(x, 5) * pow(y, 3) * z +
      -1 * pow(x, 5) * pow(y, 2) * pow(z, 5) +
      12 * pow(x, 5) * pow(y, 2) * pow(z, 2) +
      -12 * pow(x, 5) * y * pow(z, 3) + -12 * pow(x, 5) * pow(z, 4) +
      8 * pow(x, 4) * pow(y, 4) +
      6 * pow(x, 4) * pow(y, 3) * pow(z, 2) +
      8 * pow(x, 4) * pow(y, 3) * z +
      -4 * pow(x, 4) * pow(y, 2) * pow(z, 4) +
      4 * pow(x, 4) * pow(y, 2) * pow(z, 3) +
      -8 * pow(x, 4) * pow(y, 2) * pow(z, 2) +
      -4 * pow(x, 4) * y * pow(z, 5) + -2 * pow(x, 4) * y * pow(z, 4) +
      -8 * pow(x, 4) * y * pow(z, 3) + 2 * pow(x, 3) * pow(y, 4) * z +
      pow(x, 3) * pow(y, 3) * pow(z, 3) +
      -1 * pow(x, 3) * pow(y, 2) * pow(z, 5) +
      -2 * pow(x, 3) * pow(y, 2) * pow(z, 3) +
      9 * pow(x, 3) * pow(y, 2) * z + -12 * pow(x, 3) * y * pow(z, 3) +
      12 * pow(x, 3) * y * pow(z, 2) + -12 * pow(x, 3) * pow(z, 4) +
      3 * pow(x, 3) * pow(z, 3) + 6 * pow(x, 2) * pow(y, 3) +
      -6 * pow(x, 2) * pow(y, 2) * pow(z, 2) +
      8 * pow(x, 2) * pow(y, 2) * z + -2 * pow(x, 2) * y * pow(z, 4) +
      -8 * pow(x, 2) * y * pow(z, 3) + 2 * pow(x, 2) * y * pow(z, 2) +
      2 * x * pow(y, 3) * z + -2 * x * pow(y, 2) * pow(z, 3) +
      -3 * x * y * z + 3 * x * pow(z, 3) + -2 * pow(y, 2) +
      2 * y * pow(z, 2);

  expr K = expr("Z");

  expr L = list({ x, y, z });
  expr a = list({ -14, 3 });

  expr p = replace(f, L[1], a[0]);
  expr t = replace(p, L[2], a[1]);

	expr k = reduce(t);

  expr d = cont(k, L, K);
  expr s = pp(k, d, L, K);

  expr F = list({
			y, z, y + z, y + -z
  });

  expr sF = list({
      -14,
      3,
      -11,
      -17,
  });

  expr sqf = sqfFactors(s, L[0], K);

  expr wlc = wangLeadingCoeff(f, d, sqf[1], F, sF, a, L, K);
  assert(wlc[0] == f);

  expr S = set({ 187*pow(x, 2) + -23, 44*pow(x, 2) + 42*x + 1, 126*pow(x, 2) + -9*x + 28});
  expr Q = set({
      wlc[1][0],
      wlc[1][1],
      wlc[1][2],
  });

  assert(S == Q);

	expr q = set({
			-4*y + -4*z,
			(y + z)*(y + -1*z),
			-1*y*pow(z, 2)
		});

  expr E = set({
      wlc[2][0],
      wlc[2][1],
      wlc[2][2],
  });

  assert(E == q);
}




void should_get_lead_coeffs_poly_expr()
{
  expr x = expr("x");
  expr y = expr("y");
  expr z = expr("z");

  expr K = expr("Z");

  expr L = list({ x, y, z });
  expr a = list({ -14, 3 });

  expr f = polyExpr(
      4 * pow(x, 6) * pow(y, 4) * pow(z, 2) +
      4 * pow(x, 6) * pow(y, 3) * pow(z, 3) +
      -4 * pow(x, 6) * pow(y, 2) * pow(z, 4) +
      -4 * pow(x, 6) * y * pow(z, 5) +
      pow(x, 5) * pow(y, 4) * pow(z, 3) +
      12 * pow(x, 5) * pow(y, 3) * z +
      -1 * pow(x, 5) * pow(y, 2) * pow(z, 5) +
      12 * pow(x, 5) * pow(y, 2) * pow(z, 2) +
      -12 * pow(x, 5) * y * pow(z, 3) + -12 * pow(x, 5) * pow(z, 4) +
      8 * pow(x, 4) * pow(y, 4) +
      6 * pow(x, 4) * pow(y, 3) * pow(z, 2) +
      8 * pow(x, 4) * pow(y, 3) * z +
      -4 * pow(x, 4) * pow(y, 2) * pow(z, 4) +
      4 * pow(x, 4) * pow(y, 2) * pow(z, 3) +
      -8 * pow(x, 4) * pow(y, 2) * pow(z, 2) +
      -4 * pow(x, 4) * y * pow(z, 5) + -2 * pow(x, 4) * y * pow(z, 4) +
      -8 * pow(x, 4) * y * pow(z, 3) + 2 * pow(x, 3) * pow(y, 4) * z +
      pow(x, 3) * pow(y, 3) * pow(z, 3) +
      -1 * pow(x, 3) * pow(y, 2) * pow(z, 5) +
      -2 * pow(x, 3) * pow(y, 2) * pow(z, 3) +
      9 * pow(x, 3) * pow(y, 2) * z + -12 * pow(x, 3) * y * pow(z, 3) +
      12 * pow(x, 3) * y * pow(z, 2) + -12 * pow(x, 3) * pow(z, 4) +
      3 * pow(x, 3) * pow(z, 3) + 6 * pow(x, 2) * pow(y, 3) +
      -6 * pow(x, 2) * pow(y, 2) * pow(z, 2) +
      8 * pow(x, 2) * pow(y, 2) * z + -2 * pow(x, 2) * y * pow(z, 4) +
      -8 * pow(x, 2) * y * pow(z, 3) + 2 * pow(x, 2) * y * pow(z, 2) +
      2 * x * pow(y, 3) * z + -2 * x * pow(y, 2) * pow(z, 3) +
      -3 * x * y * z + 3 * x * pow(z, 3) + -2 * pow(y, 2) +
      2 * y * pow(z, 2), L);

  expr p = evalPolyExpr(f, L[1], a[0]);
  expr t = evalPolyExpr(p, L[2], a[1]);

	expr T = list({ x });

	expr d = contPolyExpr(t, T, K);
  expr s = ppPolyExpr(t, T, K);

	expr R = rest(L);

  expr F = list({
			polyExpr(y, R), polyExpr(z, R), polyExpr(y + z, R), polyExpr(y + -z, R)
  });

  expr sF = list({ -14, 3, -11, -17 });

	expr sqf = sqfFactorsPolyExpr(s, T, K);
	assert(sqf[1].size() == 3);
	expr wlc = wangLeadingCoeffPolyExpr(f, d, sqf[1], F, sF, a, L, K);

	assert(wlc[0] == f);

	assert(wlc[1].size() == 3);
  expr S = set({ polyExpr(187*pow(x, 2) + -23, T), polyExpr(44*pow(x, 2) + 42*x + 1, T), polyExpr(126*pow(x, 2) + -9*x + 28, T)});
	expr Q = set({
      wlc[1][0],
      wlc[1][1],
      wlc[1][2],
  });

  assert(S == Q);

	expr q = set({
			polyExpr(-4*y + -4*z, R),
			polyExpr(pow(y,2) + -1*pow(z, 2), R),
			polyExpr(-1*y*pow(z, 2), R)
		});

  expr E = set({
      wlc[2][0],
      wlc[2][1],
      wlc[2][2],
  });


	assert(E == q);
}




void should_factorize_multivariate_polynomials()
{
	expr x = expr("x");
	expr y = expr("y");
	expr z = expr("z");
	expr u = expr("u");

	expr L = list({ x });

	expr K = expr("Z");

	//expr U0 = factors(0, L, K);

	//printf("\n\n-----> U0 = %s\n\n", U0.toString().c_str());

	//expr U1 = factors(3, L, K);

	//printf("\n\n-----> U1 = %s\n\n", U1.toString().c_str());

	//expr U2 = factors(-8, L, K);

	//printf("\n\n-----> U2 = %s\n\n", U2.toString().c_str());

	//expr U3 = factors(pow(x, 2) + -9, L, K);

	//printf("\n\n-----> U3 = %s\n\n", U3.toString().c_str());

	expr R = list({ x , y });

	//expr U4 = factors(pow(x, 2)*pow(y, 2) + 6*pow(x, 2)*y + 9*pow(x, 2) + -1, R, K);

	//printf("\n\n-----> U4 = %s\n\n", U4.toString().c_str());
	expr T = list({ x, y, z });

  expr U5 = factors(pow(x, 2)*pow(y, 2)*pow(z, 2) + -9, T, K);

	printf("\n\n-----> U5 = %s\n\n", to_string(U5).c_str());
		return;
	expr E = list({ x, y, z, u });

	expr U6 = factors(pow(x, 2)*pow(y, 2)*pow(z, 2)*pow(u, 2) + -9, E, K);

	printf("\n\n-----> U6 = %s\n\n", to_string(U6).c_str());

	// printf("U6 = %s\n", U6->toString().c_str());


  expr u7 =
      4 * pow(x, 6) * pow(y, 4) * pow(z, 2) +
      4 * pow(x, 6) * pow(y, 3) * pow(z, 3) +
      -4 * pow(x, 6) * pow(y, 2) * pow(z, 4) +
      -4 * pow(x, 6) * y * pow(z, 5) +
      pow(x, 5) * pow(y, 4) * pow(z, 3) +
      12 * pow(x, 5) * pow(y, 3) * z +
      -1 * pow(x, 5) * pow(y, 2) * pow(z, 5) +
      12 * pow(x, 5) * pow(y, 2) * pow(z, 2) +
      -12 * pow(x, 5) * y * pow(z, 3) + -12 * pow(x, 5) * pow(z, 4) +
      8 * pow(x, 4) * pow(y, 4) +
      6 * pow(x, 4) * pow(y, 3) * pow(z, 2) +
      8 * pow(x, 4) * pow(y, 3) * z +
      -4 * pow(x, 4) * pow(y, 2) * pow(z, 4) +
      4 * pow(x, 4) * pow(y, 2) * pow(z, 3) +
      -8 * pow(x, 4) * pow(y, 2) * pow(z, 2) +
      -4 * pow(x, 4) * y * pow(z, 5) + -2 * pow(x, 4) * y * pow(z, 4) +
      -8 * pow(x, 4) * y * pow(z, 3) + 2 * pow(x, 3) * pow(y, 4) * z +
      pow(x, 3) * pow(y, 3) * pow(z, 3) +
      -1 * pow(x, 3) * pow(y, 2) * pow(z, 5) +
      -2 * pow(x, 3) * pow(y, 2) * pow(z, 3) +
      9 * pow(x, 3) * pow(y, 2) * z + -12 * pow(x, 3) * y * pow(z, 3) +
      12 * pow(x, 3) * y * pow(z, 2) + -12 * pow(x, 3) * pow(z, 4) +
      3 * pow(x, 3) * pow(z, 3) + 6 * pow(x, 2) * pow(y, 3) +
      -6 * pow(x, 2) * pow(y, 2) * pow(z, 2) +
      8 * pow(x, 2) * pow(y, 2) * z + -2 * pow(x, 2) * y * pow(z, 4) +
      -8 * pow(x, 2) * y * pow(z, 3) + 2 * pow(x, 2) * y * pow(z, 2) +
      2 * x * pow(y, 3) * z + -2 * x * pow(y, 2) * pow(z, 3) +
      -3 * x * y * z + 3 * x * pow(z, 3) + -2 * pow(y, 2) +
      2 * y * pow(z, 2);

	expr U7 = factors(u7, T, K);

	printf("\n\n-----> U7 = %s\n\n", to_string(U7).c_str());
}




void should_factorize_multivariate_polynomials_poly_expr()
{
	expr x = expr("x");
	expr y = expr("y");
	expr z = expr("z");
	expr u = expr("u");

	expr L = list({ x });

	expr K = expr("Z");

	//expr U0 = factorsPolyExpr(polyExpr(0, L), L, K);

	//printf("\n\n-----> U0 = %s\n\n", U0.toString().c_str());

	//expr U1 = factorsPolyExpr(polyExpr(3, L), L, K);

	//printf("\n\n-----> U1 = %s\n\n", U1.toString().c_str());

	//expr U2 = factorsPolyExpr(polyExpr(-8, L), L, K);

	//printf("\n\n-----> U2 = %s\n\n", U2.toString().c_str());

	//expr U3 = factorsPolyExpr(polyExpr(pow(x, 2) + -9, L), L, K);

	//printf("\n\n-----> U3 = %s\n\n", U3.toString().c_str());

	expr R = list({ x , y });

	expr U4 = factorsPolyExpr(polyExpr(pow(x, 2)*pow(y, 2) + 6*pow(x, 2)*y + 9*pow(x, 2) + -1, R), R, K);

	printf("\n\n-----> U4 = %s\n\n", to_string(U4).c_str());

	expr T = list({ x, y, z });

	expr U5 = factorsPolyExpr(polyExpr(pow(x, 2)*pow(y, 2)*pow(z, 2) + -9, T), T, K);

	printf("\n\n-----> U5 = %s\n\n", to_string(U5).c_str());
	expr E = list({ x, y, z, u });

	expr U6 = factorsPolyExpr(polyExpr(pow(x, 2)*pow(y, 2)*pow(z, 2)*pow(u, 2) + -9, E), E, K);

	printf("\n\n-----> U6 = %s\n\n", to_string(U6).c_str());

	// printf("U6 = %s\n", U6->toString().c_str());


  expr u7 =
      4 * pow(x, 6) * pow(y, 4) * pow(z, 2) +
      4 * pow(x, 6) * pow(y, 3) * pow(z, 3) +
      -4 * pow(x, 6) * pow(y, 2) * pow(z, 4) +
      -4 * pow(x, 6) * y * pow(z, 5) +
      pow(x, 5) * pow(y, 4) * pow(z, 3) +
      12 * pow(x, 5) * pow(y, 3) * z +
      -1 * pow(x, 5) * pow(y, 2) * pow(z, 5) +
      12 * pow(x, 5) * pow(y, 2) * pow(z, 2) +
      -12 * pow(x, 5) * y * pow(z, 3) + -12 * pow(x, 5) * pow(z, 4) +
      8 * pow(x, 4) * pow(y, 4) +
      6 * pow(x, 4) * pow(y, 3) * pow(z, 2) +
      8 * pow(x, 4) * pow(y, 3) * z +
      -4 * pow(x, 4) * pow(y, 2) * pow(z, 4) +
      4 * pow(x, 4) * pow(y, 2) * pow(z, 3) +
      -8 * pow(x, 4) * pow(y, 2) * pow(z, 2) +
      -4 * pow(x, 4) * y * pow(z, 5) + -2 * pow(x, 4) * y * pow(z, 4) +
      -8 * pow(x, 4) * y * pow(z, 3) + 2 * pow(x, 3) * pow(y, 4) * z +
      pow(x, 3) * pow(y, 3) * pow(z, 3) +
      -1 * pow(x, 3) * pow(y, 2) * pow(z, 5) +
      -2 * pow(x, 3) * pow(y, 2) * pow(z, 3) +
      9 * pow(x, 3) * pow(y, 2) * z + -12 * pow(x, 3) * y * pow(z, 3) +
      12 * pow(x, 3) * y * pow(z, 2) + -12 * pow(x, 3) * pow(z, 4) +
      3 * pow(x, 3) * pow(z, 3) + 6 * pow(x, 2) * pow(y, 3) +
      -6 * pow(x, 2) * pow(y, 2) * pow(z, 2) +
      8 * pow(x, 2) * pow(y, 2) * z + -2 * pow(x, 2) * y * pow(z, 4) +
      -8 * pow(x, 2) * y * pow(z, 3) + 2 * pow(x, 2) * y * pow(z, 2) +
      2 * x * pow(y, 3) * z + -2 * x * pow(y, 2) * pow(z, 3) +
      -3 * x * y * z + 3 * x * pow(z, 3) + -2 * pow(y, 2) +
      2 * y * pow(z, 2);

	expr U7 = factorsPolyExpr(polyExpr(u7, T), T, K);

	printf("\n\n-----> U7 = %s\n\n", to_string(U7).c_str());

	return;
}

void should_lift_factors()
{
	expr x = expr("x");
	expr y = expr("y");
	expr z = expr("z");

	expr f =
      4 * pow(x, 6) * pow(y, 4) * pow(z, 2) +
      4 * pow(x, 6) * pow(y, 3) * pow(z, 3) +
      -4 * pow(x, 6) * pow(y, 2) * pow(z, 4) +
      -4 * pow(x, 6) * y * pow(z, 5) +
      pow(x, 5) * pow(y, 4) * pow(z, 3) +
      12 * pow(x, 5) * pow(y, 3) * z +
      -1 * pow(x, 5) * pow(y, 2) * pow(z, 5) +
      12 * pow(x, 5) * pow(y, 2) * pow(z, 2) +
      -12 * pow(x, 5) * y * pow(z, 3) + -12 * pow(x, 5) * pow(z, 4) +
      8 * pow(x, 4) * pow(y, 4) +
      6 * pow(x, 4) * pow(y, 3) * pow(z, 2) +
      8 * pow(x, 4) * pow(y, 3) * z +
      -4 * pow(x, 4) * pow(y, 2) * pow(z, 4) +
      4 * pow(x, 4) * pow(y, 2) * pow(z, 3) +
      -8 * pow(x, 4) * pow(y, 2) * pow(z, 2) +
      -4 * pow(x, 4) * y * pow(z, 5) + -2 * pow(x, 4) * y * pow(z, 4) +
      -8 * pow(x, 4) * y * pow(z, 3) + 2 * pow(x, 3) * pow(y, 4) * z +
      pow(x, 3) * pow(y, 3) * pow(z, 3) +
      -1 * pow(x, 3) * pow(y, 2) * pow(z, 5) +
      -2 * pow(x, 3) * pow(y, 2) * pow(z, 3) +
      9 * pow(x, 3) * pow(y, 2) * z + -12 * pow(x, 3) * y * pow(z, 3) +
      12 * pow(x, 3) * y * pow(z, 2) + -12 * pow(x, 3) * pow(z, 4) +
      3 * pow(x, 3) * pow(z, 3) + 6 * pow(x, 2) * pow(y, 3) +
      -6 * pow(x, 2) * pow(y, 2) * pow(z, 2) +
      8 * pow(x, 2) * pow(y, 2) * z + -2 * pow(x, 2) * y * pow(z, 4) +
      -8 * pow(x, 2) * y * pow(z, 3) + 2 * pow(x, 2) * y * pow(z, 2) +
      2 * x * pow(y, 3) * z + -2 * x * pow(y, 2) * pow(z, 3) +
      -3 * x * y * z + 3 * x * pow(z, 3) + -2 * pow(y, 2) +
      2 * y * pow(z, 2);



	expr U = list({
			44*pow(x, 2) + 42*x + 1,
			126*pow(x, 2) + -9*x + 28,
			187*pow(x, 2) + -23
		});


	expr LC = list({
			-4*y + -4*z,
			-1*y*pow(z, 2),
			pow(y, 2) + -1*pow(z, 2)
		});

	expr a = list({ -14, 3 });

	expr L = list({ x, y, z });

	expr K = expr("Z");

	expr r = wangEEZ(f, U, LC, a, 6291469, L, K);

	assert(r.size() == 3);

	assert(expand(r[0] * r[1] * r[2]) == f);
}



void should_lift_factors_poly_expr()
{
	expr x = expr("x");
	expr y = expr("y");
	expr z = expr("z");

	expr L = list({ x, y, z });

	expr f = polyExpr(
      4 * pow(x, 6) * pow(y, 4) * pow(z, 2) +
      4 * pow(x, 6) * pow(y, 3) * pow(z, 3) +
      -4 * pow(x, 6) * pow(y, 2) * pow(z, 4) +
      -4 * pow(x, 6) * y * pow(z, 5) +
      pow(x, 5) * pow(y, 4) * pow(z, 3) +
      12 * pow(x, 5) * pow(y, 3) * z +
      -1 * pow(x, 5) * pow(y, 2) * pow(z, 5) +
      12 * pow(x, 5) * pow(y, 2) * pow(z, 2) +
      -12 * pow(x, 5) * y * pow(z, 3) + -12 * pow(x, 5) * pow(z, 4) +
      8 * pow(x, 4) * pow(y, 4) +
      6 * pow(x, 4) * pow(y, 3) * pow(z, 2) +
      8 * pow(x, 4) * pow(y, 3) * z +
      -4 * pow(x, 4) * pow(y, 2) * pow(z, 4) +
      4 * pow(x, 4) * pow(y, 2) * pow(z, 3) +
      -8 * pow(x, 4) * pow(y, 2) * pow(z, 2) +
      -4 * pow(x, 4) * y * pow(z, 5) + -2 * pow(x, 4) * y * pow(z, 4) +
      -8 * pow(x, 4) * y * pow(z, 3) + 2 * pow(x, 3) * pow(y, 4) * z +
      pow(x, 3) * pow(y, 3) * pow(z, 3) +
      -1 * pow(x, 3) * pow(y, 2) * pow(z, 5) +
      -2 * pow(x, 3) * pow(y, 2) * pow(z, 3) +
      9 * pow(x, 3) * pow(y, 2) * z + -12 * pow(x, 3) * y * pow(z, 3) +
      12 * pow(x, 3) * y * pow(z, 2) + -12 * pow(x, 3) * pow(z, 4) +
      3 * pow(x, 3) * pow(z, 3) + 6 * pow(x, 2) * pow(y, 3) +
      -6 * pow(x, 2) * pow(y, 2) * pow(z, 2) +
      8 * pow(x, 2) * pow(y, 2) * z + -2 * pow(x, 2) * y * pow(z, 4) +
      -8 * pow(x, 2) * y * pow(z, 3) + 2 * pow(x, 2) * y * pow(z, 2) +
      2 * x * pow(y, 3) * z + -2 * x * pow(y, 2) * pow(z, 3) +
      -3 * x * y * z + 3 * x * pow(z, 3) + -2 * pow(y, 2) +
      2 * y * pow(z, 2), L);

	expr U = list({
			polyExpr(44*pow(x, 2) + 42*x + 1, list({x})),
			polyExpr(126*pow(x, 2) + -9*x + 28, list({x})),
			polyExpr(187*pow(x, 2) + -23, list({x}))
		});

	expr LC = list({
			polyExpr(-4*y + -4*z, rest(L)),
			polyExpr(-1*y*pow(z, 2), rest(L)),
			polyExpr(pow(y, 2) + -1*pow(z, 2), rest(L))
		});

	expr a = list({ -14, 3 });

	//	expr L = list({ x, y, z });

	expr K = expr("Z");

	expr r = wangEEZPolyExpr(f, U, LC, a, 6291469, L, K);

	assert(r.size() == 3);

	expr g = mulPolyExpr(r[0], r[1]);
	g = mulPolyExpr(g, r[2]);

	assert(g == f);
}



int main()
{

	  // TEST(should_get_nondivisors)
	TEST(should_get_nondivisors_poly_expr)

		//TEST(should_solve_diophant)
	TEST(should_solve_diophant_poly_expr)

		//	TEST(should_get_lead_coeffs)
	TEST(should_get_lead_coeffs_poly_expr)

		//TEST(should_lift_factors)
	TEST(should_lift_factors_poly_expr)

		//TEST(should_factorize_multivariate_polynomials)
	TEST(should_factorize_multivariate_polynomials_poly_expr)

	return 0;
}
