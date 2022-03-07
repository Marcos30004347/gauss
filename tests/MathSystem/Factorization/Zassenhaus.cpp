#include "MathSystem/Algebra/Expression.hpp"
#include "test.hpp"

#include "MathSystem/Factorization/Utils.hpp"
#include "MathSystem/Factorization/Zassenhaus.hpp"
#include "MathSystem/Polynomial/Polynomial.hpp"

using namespace alg;
using namespace poly;
using namespace factorization;


void should_factorize_zassenhaus_poly_expr() {
  expr x = expr("x");
  expr Z = expr("Z");

	expr L = list({x});

  expr f = polyExpr(pow(x, 4) + -1, L);

  expr z = zassenhausPolyExpr(f, L, Z);

  assert(z == list({
				polyExpr(x + -1, L),
				polyExpr(x + 1, L),
				polyExpr(pow(x, 2) + 1, L),
			}));

  expr g = polyExpr(6*pow(x, 4) + 5*pow(x, 3) + 15*pow(x, 2) + 5*x + 4, L);

	expr r = zassenhausPolyExpr(g, L, Z);

	assert(r == list({
				polyExpr(x + 3*pow(x, 2) + 1, L),
				polyExpr(x + 2*pow(x, 2) + 4, L),
			}));
}


void should_distinct_degree_factorize_poly_expr() {
  expr x = expr("x");
	expr L = list({x});

  expr f = polyExpr(pow(x, 8) + pow(x, 7) + -1 * pow(x, 6) + pow(x, 5) +
										-1 * pow(x, 3) + -1 * pow(x, 2) + -1 * x, L);

  expr d = cantorZassenhausDDFPolyExpr(f, L, 3);

	assert(d == list({
				list({polyExpr(x, L), 1}),
				list({polyExpr(x + pow(x, 3) + pow(x, 4) + -1, L), 2}),
				list({polyExpr(-1 * x + pow(x, 3) + 1, L), 3}),
			}));
}



void should_equal_degree_factorize_poly_expr() {
  expr x = expr("x");
	expr L = list({x});

  expr a1 = polyExpr(pow(x, 5) + -1, L);
  expr a2 = polyExpr(pow(x, 10) + pow(x, 5) + 1, L);

  expr t = cantorZassenhausEDFPolyExpr(a1, L, 1, 11);
	expr h = cantorZassenhausEDFPolyExpr(a2, L, 2, 11);

  // because equal degree factoring is a probabilistic method,
  // the terms will come in  a random order, so we sort them
  // to make shure that the comparison will works
  t = sortTerms(t);
  h = sortTerms(h);

	assert(t == list({
				polyExpr(x + -5, L),
				polyExpr(x + -4, L),
				polyExpr(x + -3, L),
				polyExpr(x + -1, L),
				polyExpr(x + 2, L),
			}));

  assert(h == list({
				polyExpr(pow(x, 2) + -2 * x + 4, L),
				polyExpr(pow(x, 2) + x + 1, L),
				polyExpr(pow(x, 2) + 3 * x + -2, L),
				polyExpr(pow(x, 2) + 4 * x + 5, L),
				polyExpr(pow(x, 2) + 5 * x + 3, L),
			}));
}

void should_factorize_cantor_zassenhaus_poly_expr() {
  expr x = expr("x");
	expr L = list({ x });

  expr a = polyExpr(pow(x, 15) + -1, L);

  expr f = cantorZassenhausPolyExpr(a, L, 11);

	assert(f == list({polyExpr(1, L),
				list({
						polyExpr(pow(x, 2) + -2*x + 4, L),
						polyExpr(pow(x, 2) + x + 1, L),
						polyExpr(pow(x, 2) + 3*x + -2, L),
						polyExpr(pow(x, 2) + 4*x + 5, L),
						polyExpr(pow(x, 2) + 5*x + 3, L),
						polyExpr(x + -5, L),
						polyExpr(x + -4, L),
						polyExpr(x + -3, L),
						polyExpr(x + -1, L),
						polyExpr(x + 2, L),
					})}));
}


int main() {
	TEST(should_distinct_degree_factorize_poly_expr)
  TEST(should_equal_degree_factorize_poly_expr)
	TEST(should_factorize_cantor_zassenhaus_poly_expr)
	TEST(should_factorize_zassenhaus_poly_expr)

	return 0;
}
