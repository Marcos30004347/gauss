#include "Core/Factorization/Zassenhaus.hpp"
#include "Core/AST/AST.hpp"
#include "Core/Algebra/Algebra.hpp"
#include "Core/Algebra/List.hpp"
#include "Core/Factorization/Utils.hpp"
#include "test.hpp"

using namespace ast;
using namespace algebra;
using namespace polynomial;
using namespace factorization;

void should_factorize_zassenhaus() {
  Expr x = symbol("x");
  Expr Z = symbol("Z");

  Expr f = power(x, 4) + -1;


  Expr z = zassenhaus(f, x, Z);

  assert(z == list({
				x + -1,
				x + 1,
				power(x, 2) + 1
			}));

  Expr g = 6*power(x, 4) + 5*power(x, 3) + 15*power(x, 2) + 5*x + 4;

	Expr r = zassenhaus(g, x, Z);

	assert(r == list({
				x + 3*power(x, 2) + 1,
				x + 2*power(x, 2) + 4
			}));
}

void should_distinct_degree_factorize() {
  Expr x = Expr("x");

  Expr f = power(x, 8) + power(x, 7) + -1 * power(x, 6) + power(x, 5) +
           -1 * power(x, 3) + -1 * power(x, 2) + -1 * x;

  Expr d = cantorZassenhausDDF(f, x, 3);

  assert(d == list({
                  list({x, 1}),
                  list({x + power(x, 3) + power(x, 4) + -1, 2}),
                  list({-1 * x + power(x, 3) + 1, 3}),
              }));
}

void should_equal_degree_factorize() {
  Expr x = Expr("x");

  Expr a1 = power(x, 5) + -1;
  Expr a2 = power(x, 10) + power(x, 5) + 1;

  Expr t = cantorZassenhausEDF(a1, x, 1, 11);
  Expr h = cantorZassenhausEDF(a2, x, 2, 11);

  // because equal degree factoring is a probabilistic method,
  // the terms will come in  a random order, so we sort them
  // to make shure that the comparison will works
  t = sortTerms(t);
  h = sortTerms(h);

  assert(t == list({
                  x + -1,
                  x + -3,
                  x + -4,
                  x + -5,
                  x + 2,
              }));

  assert(h == list({
                  power(x, 2) + 3 * x + -2,
                  power(x, 2) + x + 1,
                  power(x, 2) + 5 * x + 3,
                  power(x, 2) + -2 * x + 4,
                  power(x, 2) + 4 * x + 5,
              }));
}

void should_factorize_cantor_zassenhaus() {
  Expr x = Expr("x");

  Expr a = power(x, 15) - 1;

  Expr f = cantorZassenhaus(a, x, 11);

  assert(f == list({1, list({x + -2, power(x, 2) + -3 * x + -2, x + 1,
                             power(x, 2) + -1 * x + 1, x + 3,
                             power(x, 2) + -5 * x + 3, x + 4,
                             power(x, 2) + 2 * x + 4, x + 5,
                             power(x, 2) + -4 * x + 5})}));
}

int main() {
  TEST(should_distinct_degree_factorize)
  TEST(should_equal_degree_factorize)
	TEST(should_factorize_cantor_zassenhaus)
	TEST(should_factorize_zassenhaus)

	return 0;
}
