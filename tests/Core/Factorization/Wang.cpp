#include "Core/AST/AST.hpp"
#include "Core/AST/Int.hpp"
#include "test.hpp"

#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Factorization/Wang.hpp"
#include "Core/Simplification/Simplification.hpp"
#include <iterator>

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace simplification;
using namespace factorization;

void should_get_nondivisors()
{
	Expr y = Expr("y");
	Expr z = Expr("z");
	Expr K = Expr("Z");

	Expr F = list({-14, 3, -11, -17 });

	Expr L = list({ y, z});

	Expr d = nondivisors(4, F, 1, L, K);

	assert(d[0] == 7);
	assert(d[1] == 3);
	assert(d[2] == 11);
	assert(d[3] == 17);
}

void should_get_nondivisors_poly_expr()
{
	Expr y = Expr("y");
	Expr z = Expr("z");
	Expr K = Expr("Z");

	Expr F = list({ -14, 3, -11, -17 });

	Expr L = list({ y, z});

	Expr d = nondivisorsPolyExpr(4, F, 1, L, K);

	assert(d[0] == 7);
	assert(d[1] == 3);
	assert(d[2] == 11);
	assert(d[3] == 17);
}




void should_solve_diophant()
{
	Expr x = Expr("x");
	Expr y = Expr("y");

	Expr H1 = list({
			44*power(x, 2) + 42*x + 1,
			126*power(x, 2) + -9*x + 28,
			187*power(x, 2) + -23
		});

	Expr H2 = list({
			-4*power(x, 2)*y + -12*power(x, 2) + -3*x*y + 1,
			-9*power(x, 2)*y + -9*x + -2*y,
			power(x, 2)*power(y, 2) + -9*power(x, 2) + y + -9
		});

	Expr H3 = list({
			-4*power(x, 2)*y + -12*power(x, 2) + -3*x*y + 1,
			-9*power(x, 2)*y + -9*x + -2*y,
			power(x, 2)*power(y, 2) + -9*power(x, 2) + y + -9
		});

	Expr c1 = -70686*power(x, 5) + -5863*power(x, 4) + -17826*power(x, 3) + 2009*power(x, 2) + 5031*x + 74;

	Expr c2 =
		9 * power(x, 5) * power(y, 4) + 12 * power(x, 5) * power(y, 3) +
		-45 * power(x, 5) * power(y, 2) + -108 * power(x, 5) * y +
		-324 * power(x, 5) + 18 * power(x, 4) * power(y, 3) +
		-216 * power(x, 4) * power(y, 2) + -810 * power(x, 4) * y +
		2 * power(x, 3) * power(y, 4) + 9 * power(x, 3) * power(y, 3) +
		-252 * power(x, 3) * power(y, 2) + -288 * power(x, 3) * y +
		-945 * power(x, 3) + -30 * power(x, 2) * power(y, 2) +
		-414 * power(x, 2) * y + 2 * x * power(y, 3) +
		-54 * x * power(y, 2) + -3 * x * y + 81 * x + 12 * y;

	Expr c3 = -36*power(x, 4)*power(y, 2) + -108*power(x, 4)*y + -27*power(x, 3)*power(y, 2) + -36*power(x, 3)*y + -108*power(x, 3) + -8*power(x, 2)*power(y, 2) + -42*power(x, 2)*y + -6*x*power(y, 2)+ 9*x + 2*y;

	Expr L1 = list({ x });
	Expr I1 = list({});

 	Expr D1 = multivariateDiophant(H1, c1, L1, I1, 5, 6291469, 1);

	Expr R1 = list({ -3*x, -2, 1 });

	assert(D1 == R1);

	Expr L2 = list({ x, y });
	Expr I2 = list({ -14 });

 	Expr D2 = multivariateDiophant(H2, c2, L2, I2, 5, 6291469, 1);

	Expr R2 = list({ -1*x*y, -3*x, -6 });

	assert(D2 == R2);

 	Expr D3 = multivariateDiophant(H3, c3, L2, I2, 5, 6291469, 1);

	Expr R3 = list({ 0, 0, -1 });

	assert(D3 == R3);
}





void should_solve_diophant_poly_expr()
{
	Expr x = Expr("x");
	Expr y = Expr("y");
	Expr Z = Expr("Z");

	Expr L1 = list({ x });
	Expr L2 = list({ x, y });

	Expr H1 = list({
			polyExpr(44*power(x, 2) + 42*x + 1, L1),
			polyExpr(126*power(x, 2) + -9*x + 28, L1),
			polyExpr(187*power(x, 2) + -23, L1)
		});

	Expr H2 = list({
			polyExpr(-4*power(x, 2)*y + -12*power(x, 2) + -3*x*y + 1, L2),
			polyExpr(-9*power(x, 2)*y + -9*x + -2*y, L2),
			polyExpr(power(x, 2)*power(y, 2) + -9*power(x, 2) + y + -9, L2)
		});

	Expr H3 = list({
			polyExpr(-4*power(x, 2)*y + -12*power(x, 2) + -3*x*y + 1, L2),
			polyExpr(-9*power(x, 2)*y + -9*x + -2*y, L2),
			polyExpr(power(x, 2)*power(y, 2) + -9*power(x, 2) + y + -9, L2)
		});

	Expr c1 = polyExpr(-70686*power(x, 5) + -5863*power(x, 4) + -17826*power(x, 3) + 2009*power(x, 2) + 5031*x + 74, L1);

	Expr c2 =
		polyExpr(
		9 * power(x, 5) * power(y, 4) + 12 * power(x, 5) * power(y, 3) +
		-45 * power(x, 5) * power(y, 2) + -108 * power(x, 5) * y +
		-324 * power(x, 5) + 18 * power(x, 4) * power(y, 3) +
		-216 * power(x, 4) * power(y, 2) + -810 * power(x, 4) * y +
		2 * power(x, 3) * power(y, 4) + 9 * power(x, 3) * power(y, 3) +
		-252 * power(x, 3) * power(y, 2) + -288 * power(x, 3) * y +
		-945 * power(x, 3) + -30 * power(x, 2) * power(y, 2) +
		-414 * power(x, 2) * y + 2 * x * power(y, 3) +
		-54 * x * power(y, 2) + -3 * x * y + 81 * x + 12 * y, L2);

	Expr c3 = polyExpr(-36*power(x, 4)*power(y, 2) + -108*power(x, 4)*y + -27*power(x, 3)*power(y, 2) + -36*power(x, 3)*y + -108*power(x, 3) + -8*power(x, 2)*power(y, 2) + -42*power(x, 2)*y + -6*x*power(y, 2)+ 9*x + 2*y, L2);

	Expr I1 = list({});

 	Expr D1 = multivariateDiophantPolyExpr(H1, c1, L1, I1, 5, 6291469, 1, Z);

	Expr R1 = list({ polyExpr(-3*x, L1), polyExpr(-2, L1), polyExpr(1, L1) });

	assert(D1 == R1);

	Expr I2 = list({ -14 });

 	Expr D2 = multivariateDiophantPolyExpr(H2, c2, L2, I2, 5, 6291469, 1, Z);

	Expr R2 = list({ polyExpr(-1*x*y, L2), polyExpr(-3*x, L2), polyExpr(-6, L2) });
	assert(D2 == R2);

 	Expr D3 = multivariateDiophant(H3, c3, L2, I2, 5, 6291469, 1);

	Expr R3 = list({ polyExpr(0, L2), polyExpr(0, L2), polyExpr(-1, L2) });

	assert(D3 == R3);
}



void should_get_lead_coeffs()
{
  Expr x = Expr("x");
  Expr y = Expr("y");
  Expr z = Expr("z");

  Expr f =
      4 * power(x, 6) * power(y, 4) * power(z, 2) +
      4 * power(x, 6) * power(y, 3) * power(z, 3) +
      -4 * power(x, 6) * power(y, 2) * power(z, 4) +
      -4 * power(x, 6) * y * power(z, 5) +
      power(x, 5) * power(y, 4) * power(z, 3) +
      12 * power(x, 5) * power(y, 3) * z +
      -1 * power(x, 5) * power(y, 2) * power(z, 5) +
      12 * power(x, 5) * power(y, 2) * power(z, 2) +
      -12 * power(x, 5) * y * power(z, 3) + -12 * power(x, 5) * power(z, 4) +
      8 * power(x, 4) * power(y, 4) +
      6 * power(x, 4) * power(y, 3) * power(z, 2) +
      8 * power(x, 4) * power(y, 3) * z +
      -4 * power(x, 4) * power(y, 2) * power(z, 4) +
      4 * power(x, 4) * power(y, 2) * power(z, 3) +
      -8 * power(x, 4) * power(y, 2) * power(z, 2) +
      -4 * power(x, 4) * y * power(z, 5) + -2 * power(x, 4) * y * power(z, 4) +
      -8 * power(x, 4) * y * power(z, 3) + 2 * power(x, 3) * power(y, 4) * z +
      power(x, 3) * power(y, 3) * power(z, 3) +
      -1 * power(x, 3) * power(y, 2) * power(z, 5) +
      -2 * power(x, 3) * power(y, 2) * power(z, 3) +
      9 * power(x, 3) * power(y, 2) * z + -12 * power(x, 3) * y * power(z, 3) +
      12 * power(x, 3) * y * power(z, 2) + -12 * power(x, 3) * power(z, 4) +
      3 * power(x, 3) * power(z, 3) + 6 * power(x, 2) * power(y, 3) +
      -6 * power(x, 2) * power(y, 2) * power(z, 2) +
      8 * power(x, 2) * power(y, 2) * z + -2 * power(x, 2) * y * power(z, 4) +
      -8 * power(x, 2) * y * power(z, 3) + 2 * power(x, 2) * y * power(z, 2) +
      2 * x * power(y, 3) * z + -2 * x * power(y, 2) * power(z, 3) +
      -3 * x * y * z + 3 * x * power(z, 3) + -2 * power(y, 2) +
      2 * y * power(z, 2);

  Expr K = Expr("Z");

  Expr L = list({ x, y, z });
  Expr a = list({ -14, 3 });

  Expr p = deepReplace(f, L[1], a[0]);
  Expr t = deepReplace(p, L[2], a[1]);

	Expr k = reduceAST(t);

  Expr d = cont(k, L, K);
  Expr s = pp(k, d, L, K);

  Expr F = list({
			y, z, y + z, y + -z
  });

  Expr sF = list({
      -14,
      3,
      -11,
      -17,
  });

  Expr sqf = sqfFactors(s, L[0], K);

  Expr wlc = wangLeadingCoeff(f, d, sqf[1], F, sF, a, L, K);
  assert(wlc[0] == f);

  Expr S = set({ 187*power(x, 2) + -23, 44*power(x, 2) + 42*x + 1, 126*power(x, 2) + -9*x + 28});
  Expr Q = set({
      wlc[1][0],
      wlc[1][1],
      wlc[1][2],
  });

  assert(S == Q);

	Expr q = set({
			-4*y + -4*z,
			(y + z)*(y + -1*z),
			-1*y*power(z, 2)
		});

  Expr E = set({
      wlc[2][0],
      wlc[2][1],
      wlc[2][2],
  });

  assert(E == q);
}

void should_factorize_multivariate_polynomials()
{
	Expr x = Expr("x");
	Expr y = Expr("y");
	Expr z = Expr("z");
	Expr u = Expr("u");

	Expr L = list({ x });

	Expr K = Expr("Z");

	Expr U0 = factors(0, L, K);

	printf("\n\n-----> U0 = %s\n\n", U0.toString().c_str());

	Expr U1 = factors(3, L, K);

	printf("\n\n-----> U1 = %s\n\n", U1.toString().c_str());

	Expr U2 = factors(-8, L, K);

	printf("\n\n-----> U2 = %s\n\n", U2.toString().c_str());

	Expr U3 = factors(power(x, 2) + -9, L, K);

	printf("\n\n-----> U3 = %s\n\n", U3.toString().c_str());

	Expr R = list({ x , y });

	Expr U4 = factors(power(x, 2)*power(y, 2) + 6*power(x, 2)*y + 9*power(x, 2) + -1, R, K);

	printf("\n\n-----> U4 = %s\n\n", U4.toString().c_str());

	Expr T = list({ x, y, z });

	Expr U5 = factors(power(x, 2)*power(y, 2)*power(z, 2) + -9, T, K);

	printf("\n\n-----> U5 = %s\n\n", U5.toString().c_str());
	Expr E = list({ x, y, z, u });

	Expr U6 = factors(power(x, 2)*power(y, 2)*power(z, 2)*power(u, 2) + -9, E, K);

	printf("\n\n-----> U6 = %s\n\n", U6.toString().c_str());

	// printf("U6 = %s\n", U6->toString().c_str());


  Expr u7 =
      4 * power(x, 6) * power(y, 4) * power(z, 2) +
      4 * power(x, 6) * power(y, 3) * power(z, 3) +
      -4 * power(x, 6) * power(y, 2) * power(z, 4) +
      -4 * power(x, 6) * y * power(z, 5) +
      power(x, 5) * power(y, 4) * power(z, 3) +
      12 * power(x, 5) * power(y, 3) * z +
      -1 * power(x, 5) * power(y, 2) * power(z, 5) +
      12 * power(x, 5) * power(y, 2) * power(z, 2) +
      -12 * power(x, 5) * y * power(z, 3) + -12 * power(x, 5) * power(z, 4) +
      8 * power(x, 4) * power(y, 4) +
      6 * power(x, 4) * power(y, 3) * power(z, 2) +
      8 * power(x, 4) * power(y, 3) * z +
      -4 * power(x, 4) * power(y, 2) * power(z, 4) +
      4 * power(x, 4) * power(y, 2) * power(z, 3) +
      -8 * power(x, 4) * power(y, 2) * power(z, 2) +
      -4 * power(x, 4) * y * power(z, 5) + -2 * power(x, 4) * y * power(z, 4) +
      -8 * power(x, 4) * y * power(z, 3) + 2 * power(x, 3) * power(y, 4) * z +
      power(x, 3) * power(y, 3) * power(z, 3) +
      -1 * power(x, 3) * power(y, 2) * power(z, 5) +
      -2 * power(x, 3) * power(y, 2) * power(z, 3) +
      9 * power(x, 3) * power(y, 2) * z + -12 * power(x, 3) * y * power(z, 3) +
      12 * power(x, 3) * y * power(z, 2) + -12 * power(x, 3) * power(z, 4) +
      3 * power(x, 3) * power(z, 3) + 6 * power(x, 2) * power(y, 3) +
      -6 * power(x, 2) * power(y, 2) * power(z, 2) +
      8 * power(x, 2) * power(y, 2) * z + -2 * power(x, 2) * y * power(z, 4) +
      -8 * power(x, 2) * y * power(z, 3) + 2 * power(x, 2) * y * power(z, 2) +
      2 * x * power(y, 3) * z + -2 * x * power(y, 2) * power(z, 3) +
      -3 * x * y * z + 3 * x * power(z, 3) + -2 * power(y, 2) +
      2 * y * power(z, 2);

	Expr U7 = factors(u7, T, K);

	printf("\n\n-----> U7 = %s\n\n", U7.toString().c_str());
}

void should_get_evaluation_points()
{
	Expr U = add({
		mul({
			integer(4),
			power(symbol("x"), integer(6)),
			power(symbol("y"), integer(4)),
			power(symbol("z"), integer(2))
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(6)),
			power(symbol("y"), integer(3)),
			power(symbol("z"), integer(3))
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(6)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(4))
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(6)),
			symbol("y"),
			power(symbol("z"), integer(5))
		}),
		mul({
			power(symbol("x"), integer(5)),
			power(symbol("y"), integer(4)),
			power(symbol("z"), integer(3))
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(5)),
			power(symbol("y"), integer(3)),
			symbol("z")
		}),
		mul({
			integer(-1),
			power(symbol("x"), integer(5)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(5)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(5)),
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(5)),
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)),
			power(symbol("y"), integer(4)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(4)),
			power(symbol("y"), integer(3)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(4)),
			power(symbol("y"), integer(3)),
			symbol("z")
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(4)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(4),
			power(symbol("x"), integer(4)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(4)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-4),
			power(symbol("x"), integer(4)),
			symbol("y"),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(4)),
			symbol("y"),
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(4)),
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(3)),
			power(symbol("y"), integer(4)),
			symbol("z")
		}),
		mul({
			power(symbol("x"), integer(3)),
			power(symbol("y"), integer(3)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-1),
			power(symbol("x"), integer(3)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(5)),
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(3)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(9),
			power(symbol("x"), integer(3)),
			power(symbol("y"), integer(2)),
			symbol("z")
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(3)),
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(12),
			power(symbol("x"), integer(3)),
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(-12),
			power(symbol("x"), integer(3)),
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(3),
			power(symbol("x"), integer(3)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(6),
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(3)),
		}),
		mul({
			integer(-6),
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(8),
			power(symbol("x"), integer(2)),
			power(symbol("y"), integer(2)),
			symbol("z")
		}),
		mul({
			integer(-2),
			power(symbol("x"), integer(2)),
			symbol("y"),
			power(symbol("z"), integer(4)),
		}),
		mul({
			integer(-8),
			power(symbol("x"), integer(2)),
			symbol("y"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(2),
			power(symbol("x"), integer(2)),
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("x"),
			power(symbol("y"), integer(3)),
			symbol("z"),
		}),
		mul({
			integer(-2),
			symbol("x"),
			power(symbol("y"), integer(2)),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-3),
			symbol("x"),
			symbol("y"),
			symbol("z"),
		}),
		mul({
			integer(3),
			symbol("x"),
			power(symbol("z"), integer(3)),
		}),
		mul({
			integer(-2),
			power(symbol("y"), integer(2)),
		}),
		mul({
			integer(2),
			symbol("y"),
			power(symbol("z"), integer(2)),
		}),
	});

	Expr V = list({
		symbol("y"),
		symbol("z"),
		add({symbol("y"), symbol("z")}),
		add({symbol("y"), mul({ integer(-1), symbol("z")}) }),
	});

	Expr G = integer(4);
	Expr L = list({ symbol("x"), symbol("y"), symbol("z") });
	Expr K = symbol("Z");

	Expr E = set({});

	E = getEvaluationPoints(U, G, V, L, K, 3, E);

	printf("%s\n", E.toString().c_str());

}

void should_lift_factors()
{
	Expr x = Expr("x");
	Expr y = Expr("y");
	Expr z = Expr("z");

	Expr f =
      4 * power(x, 6) * power(y, 4) * power(z, 2) +
      4 * power(x, 6) * power(y, 3) * power(z, 3) +
      -4 * power(x, 6) * power(y, 2) * power(z, 4) +
      -4 * power(x, 6) * y * power(z, 5) +
      power(x, 5) * power(y, 4) * power(z, 3) +
      12 * power(x, 5) * power(y, 3) * z +
      -1 * power(x, 5) * power(y, 2) * power(z, 5) +
      12 * power(x, 5) * power(y, 2) * power(z, 2) +
      -12 * power(x, 5) * y * power(z, 3) + -12 * power(x, 5) * power(z, 4) +
      8 * power(x, 4) * power(y, 4) +
      6 * power(x, 4) * power(y, 3) * power(z, 2) +
      8 * power(x, 4) * power(y, 3) * z +
      -4 * power(x, 4) * power(y, 2) * power(z, 4) +
      4 * power(x, 4) * power(y, 2) * power(z, 3) +
      -8 * power(x, 4) * power(y, 2) * power(z, 2) +
      -4 * power(x, 4) * y * power(z, 5) + -2 * power(x, 4) * y * power(z, 4) +
      -8 * power(x, 4) * y * power(z, 3) + 2 * power(x, 3) * power(y, 4) * z +
      power(x, 3) * power(y, 3) * power(z, 3) +
      -1 * power(x, 3) * power(y, 2) * power(z, 5) +
      -2 * power(x, 3) * power(y, 2) * power(z, 3) +
      9 * power(x, 3) * power(y, 2) * z + -12 * power(x, 3) * y * power(z, 3) +
      12 * power(x, 3) * y * power(z, 2) + -12 * power(x, 3) * power(z, 4) +
      3 * power(x, 3) * power(z, 3) + 6 * power(x, 2) * power(y, 3) +
      -6 * power(x, 2) * power(y, 2) * power(z, 2) +
      8 * power(x, 2) * power(y, 2) * z + -2 * power(x, 2) * y * power(z, 4) +
      -8 * power(x, 2) * y * power(z, 3) + 2 * power(x, 2) * y * power(z, 2) +
      2 * x * power(y, 3) * z + -2 * x * power(y, 2) * power(z, 3) +
      -3 * x * y * z + 3 * x * power(z, 3) + -2 * power(y, 2) +
      2 * y * power(z, 2);



	Expr U = list({
			44*power(x, 2) + 42*x + 1,
			126*power(x, 2) + -9*x + 28,
			187*power(x, 2) + -23
		});


	Expr LC = list({
			-4*y + -4*z,
			-1*y*power(z, 2),
			power(y, 2) + -1*power(z, 2)
		});

	Expr a = list({ -14, 3 });

	Expr L = list({ x, y, z });

	Expr K = Expr("Z");

	Expr r = wangEEZ(f, U, LC, a, 6291469, L, K);

	assert(r.size() == 3);
	assert(algebraicExpand(r[0] * r[1] * r[2]) == f);

	printf("r ========= %s\n", r.toString().c_str());
	// assert(r == list({
	// 			z*y*x + 4*(y + z)*power(x, 2) + -1,
	// 			3*z*x + power(z, 2)*y*power(x, 2) + 2*y,
	// 			(power(y, 2) + -1*power(z, 2))*power(x, 2) + y + -1*power(z, 2)
	// 		}));

}

int main()
{

	TEST(should_get_nondivisors)
	TEST(should_get_nondivisors_poly_expr)

	TEST(should_solve_diophant)
	TEST(should_solve_diophant_poly_expr)

	TEST(should_get_lead_coeffs)

	TEST(should_get_evaluation_points)

	TEST(should_lift_factors)

	TEST(should_factorize_multivariate_polynomials)

	return 0;
}
