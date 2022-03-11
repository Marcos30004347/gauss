#include <cstdlib>

#define TEST_TIME_REPORT_UNIT TEST_TIME_REPORT_NS

#include "test.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>

#include "gauss/Algebra/Expression.hpp"

using namespace alg;

void should_construct_expr() {
  expr expr0 = create(kind::ADD);

  assert(is(&expr0, kind::ADD));

  expr expr1 = create(kind::ADD, {integer(3), integer(4), integer(5)});

  assert(is(&expr1, kind::ADD));

  assert(is(operand(&expr1, 0), kind::INT));
  assert(is(operand(&expr1, 1), kind::INT));
  assert(is(operand(&expr1, 2), kind::INT));

  assert(get_val(operand(&expr1, 0)) == 3);
  assert(get_val(operand(&expr1, 1)) == 4);
  assert(get_val(operand(&expr1, 2)) == 5);

  expr expr2 = symbol("x");

  assert(is(&expr2, kind::SYM));

  assert(strcmp(get_id(&expr2), "x") == 0);
}

void should_insert_and_remove_from_expr() {
  expr a = create(kind::ADD);

  a.insert(integer(0), 0);
  a.insert(integer(1), 1);
  a.insert(integer(2), 2);
  a.insert(integer(3), 3);
  a.insert(integer(4), 4);
  a.insert(integer(5), 5);
  a.insert(integer(6), 6);
  a.insert(integer(7), 7);
  a.insert(integer(8), 8);
  a.insert(integer(9), 9);

  assert(size_of(&a) == 10);

  assert(kind_of(operand(&a, 0)) == kind::INT);
  assert(kind_of(operand(&a, 1)) == kind::INT);
  assert(kind_of(operand(&a, 2)) == kind::INT);
  assert(kind_of(operand(&a, 3)) == kind::INT);
  assert(kind_of(operand(&a, 4)) == kind::INT);
  assert(kind_of(operand(&a, 5)) == kind::INT);
  assert(kind_of(operand(&a, 6)) == kind::INT);
  assert(kind_of(operand(&a, 7)) == kind::INT);
  assert(kind_of(operand(&a, 8)) == kind::INT);
  assert(kind_of(operand(&a, 9)) == kind::INT);

  assert(get_val(operand(&a, 0)) == 0);
  assert(get_val(operand(&a, 1)) == 1);
  assert(get_val(operand(&a, 2)) == 2);
  assert(get_val(operand(&a, 3)) == 3);
  assert(get_val(operand(&a, 4)) == 4);
  assert(get_val(operand(&a, 5)) == 5);
  assert(get_val(operand(&a, 6)) == 6);
  assert(get_val(operand(&a, 7)) == 7);
  assert(get_val(operand(&a, 8)) == 8);
  assert(get_val(operand(&a, 9)) == 9);

  a.insert(integer(10), 2);

  assert(kind_of(operand(&a, 1)) == kind::INT);
  assert(kind_of(operand(&a, 2)) == kind::INT);
  assert(kind_of(operand(&a, 3)) == kind::INT);

  assert(get_val(operand(&a, 1)) == 1);
  assert(get_val(operand(&a, 2)) == 10);
  assert(get_val(operand(&a, 3)) == 2);

  assert(size_of(&a) == 11);

  a.remove(2);

  assert(size_of(&a) == 10);

  assert(get_val(operand(&a, 1)) == 1);
  assert(get_val(operand(&a, 2)) == 2);
  assert(get_val(operand(&a, 3)) == 3);
}

void should_eval_exprs() {
  expr x = expr("x");
  expr y = expr("y");
  expr z = expr("z");

  // expr a =
  //     (3 * z * y * x + 2 * pow(z, 3) * pow(y, 2) * x + -2 * z * pow(y, 3) * x +
  //      -3 * pow(z, 3) * x + -2 * pow(z, 2) * y * pow(x, 2) +
  //      8 * pow(z, 3) * y * pow(x, 2) + 2 * pow(z, 4) * y * pow(x, 2) +
  //      -8 * z * pow(y, 2) * pow(x, 2) + 6 * pow(z, 2) * pow(y, 2) * pow(x, 2) +
  //      -6 * pow(y, 3) * pow(x, 2) + -12 * pow(z, 2) * y * pow(x, 3) +
  //      12 * pow(z, 3) * y * pow(x, 3) + -9 * z * pow(y, 2) * pow(x, 3) +
  //      2 * pow(z, 3) * pow(y, 2) * pow(x, 3) +
  //      pow(z, 5) * pow(y, 2) * pow(x, 3) +
  //      -1 * pow(z, 3) * pow(y, 3) * pow(x, 3) + -2 * z * pow(y, 4) * pow(x, 3) +
  //      -3 * pow(z, 3) * pow(x, 3) + 12 * pow(z, 4) * pow(x, 3) +
  //      8 * pow(z, 3) * y * pow(x, 4) + 2 * pow(z, 4) * y * pow(x, 4) +
  //      4 * pow(z, 5) * y * pow(x, 4) + 8 * pow(z, 2) * pow(y, 2) * pow(x, 4) +
  //      -4 * pow(z, 3) * pow(y, 2) * pow(x, 4) +
  //      4 * pow(z, 4) * pow(y, 2) * pow(x, 4) + -8 * z * pow(y, 3) * pow(x, 4) +
  //      -6 * pow(z, 2) * pow(y, 3) * pow(x, 4) + -8 * pow(y, 4) * pow(x, 4) +
  //      12 * pow(z, 3) * y * pow(x, 5) +
  //      -12 * pow(z, 2) * pow(y, 2) * pow(x, 5) +
  //      pow(z, 5) * pow(y, 2) * pow(x, 5) + -12 * z * pow(y, 3) * pow(x, 5) +
  //      -1 * pow(z, 3) * pow(y, 4) * pow(x, 5) + 12 * pow(z, 4) * pow(x, 5) +
  //      4 * pow(z, 5) * y * pow(x, 6) + 4 * pow(z, 4) * pow(y, 2) * pow(x, 6) +
  //      -4 * pow(z, 3) * pow(y, 3) * pow(x, 6) +
  //      -4 * pow(z, 2) * pow(y, 4) * pow(x, 6) + -2 * pow(z, 2) * y +
  //      2 * pow(y, 2));

  // expr b =
  //     (-4 * pow(x, 6) * pow(y, 4) * pow(z, 2) +
  //      -4 * pow(x, 6) * pow(y, 3) * pow(z, 3) +
  //      4 * pow(x, 6) * pow(y, 2) * pow(z, 4) + 4 * pow(x, 6) * y * pow(z, 5));

	// expr c = reduce(a - b);

  expr _1 = expr(1);
  expr _4 = expr(4);
  expr _7 = expr(7);

  expr a2 = y * z * x * pow(y, 2) * 4 * pow(x, 2) * z + (_1 + _1 + _1) +
            x * pow(y, 2) * 4 * pow(x, 2) * y + (x + 4 + pow(z, 3)) +
            (2 + z + _1) + 4 + (_4 * _1 * _7);

  reduce(&a2);

  assert(a2 == 4 * pow(x, 3) * pow(y, 3) * pow(z, 2) +
                   4 * pow(x, 3) * pow(y, 3) + pow(z, 3) + x + z + 42);
}

void should_expand_expr() {
  expr x = expr("x");
  expr y = expr("y");
  expr z = expr("z");

  expr a = (x + 2) * (x + 3) * (x + 4);

  expand(&a);

  assert(a == pow(x, 3) + 9 * pow(x, 2) + 26 * x + 24);

  expr b = pow(x * sqrt(y + 1) + 1, 4);

  expand(&b);

	assert(b == pow(x, 4) * pow(y, 2) +
                  4 * pow(x, 3) * pow(y + 1, fraction(3, 2)) + 6 * pow(x, 2) +
                  4 * x * pow(y + 1, fraction(1, 2)) + 2 * pow(x, 4) * y +
                  6 * pow(x, 2) * y + pow(x, 4) + 1);

  expr c = pow(x + y + z, 3);

  expand(&c);

  assert(c == 3 * x * pow(y, 2) + 3 * x * pow(z, 2) + 3 * y * pow(z, 2) +
                  3 * pow(x, 2) * y + 3 * pow(x, 2) * z + 3 * pow(y, 2) * z +
                  6 * x * y * z + pow(x, 3) + pow(y, 3) + pow(z, 3));

  expr d = pow(x + 1, 2);

  expand(&d);
  assert(d == pow(x, 2) + 2 * x + 1);

  expr e = pow(pow(x + 2, 2) + 3, 2);

  expand(&e);

  assert(e == pow(x, 4) + 8 * pow(x, 3) + 30 * pow(x, 2) + 56 * x + 49);

  expr f = (-32 * pow(z, 3) + 32 * pow(z, 4) + 48 * pow(z, 5) +
            -24 * pow(z, 6) + -48 * pow(z, 7) + -36 * pow(z, 8) +
            -40 * pow(z, 9) + -8 * pow(z, 10) + -8 * pow(z, 11)) /
           (4 * pow(z, 2));

  expand(&f);

  assert(f == -2 * pow(z, 9) + -2 * pow(z, 8) + -10 * pow(z, 7) +
                  -9 * pow(z, 6) + -12 * pow(z, 5) + -6 * pow(z, 4) +
                  12 * pow(z, 3) + 8 * pow(z, 2) + -8 * z);
}

void should_eval_equality() {
  expr x = expr("x");
  expr a = 3 + 4 * x;
  expr b = 4 * x + 3;
  expr c = 4 * x + 4;

  expr d = 4 * pow(x, 3) + 4 * pow(x, 2) + 6 * x + 7;
  expr e = 4 * pow(x, 2) + 4 * pow(x, 3) + 7 + 6 * x;
  assert(a == b);

	assert(a != c);
  assert(d == e);
  assert(a != d);
}

void should_perform_list_operations() {
  list t = {1, 2, 3};

  assert(t.size() == 3);

  assert(t[0] == 1);
  assert(t[1] == 2);
  assert(t[2] == 3);

  t.insert(4);

  assert(t.size() == 4);

  assert(t[3] == 4);


  list g = insert(t, 5);

  assert(g.size() == 5);

  assert(g[0] == 1);
  assert(g[1] == 2);
  assert(g[2] == 3);
  assert(g[3] == 4);
  assert(g[4] == 5);

  list k = remove(g, 3);

  assert(k.size() == 4);

  assert(k[0] == 1);
  assert(k[1] == 2);
  assert(k[2] == 3);
  assert(k[3] == 5);

	list u = remove(k, { expr(2) });

  assert(u.size() == 3);

  assert(u[0] == 1);
  assert(u[1] == 3);
  assert(u[2] == 5);

  list f = join(t, k);

  assert(f.size() == 8);

  assert(f[0] == 1);
  assert(f[1] == 2);
  assert(f[2] == 3);
  assert(f[3] == 4);
  assert(f[4] == 1);
  assert(f[5] == 2);
  assert(f[6] == 3);
  assert(f[7] == 5);
}

void should_perform_set_operations() {
  set t = {1, 2, 3, 2};

  assert(t.size() == 3);
  assert(t[0] == 1);
  assert(t[1] == 2);
  assert(t[2] == 3);

  set k = {4, 4, 5};

  set u = unification(t, k);

  assert(u.size() == 5);

  assert(u[0] == 1);
  assert(u[1] == 2);
  assert(u[2] == 3);
  assert(u[3] == 4);
  assert(u[4] == 5);

  set h = difference(u, k);
  assert(h.size() == 3);

  assert(h[0] == 1);
  assert(h[1] == 2);
  assert(h[2] == 3);

  set l = {2, 1};

  set v = intersection(h, l);

  assert(v.size() == 2);

  assert(v[0] == 1);
  assert(v[1] == 2);
}

void should_simplify_additions() {
  expr x = expr("x");
  expr y = expr("y");
  expr z = expr("z");

  expr exp0 = expr(2) + expr(2);
  expr exp1 = expr(2) + expr(4) + expr(5) + expr(6);
  expr exp2 = expr(3) + expr(5) + expr(6);
  expr exp3 = expr(5) + expr(6) + expr(3);
  expr exp4 = expr(2) + expr(3) + expr(4) + expr(5);
  expr exp5 = x + x;
  expr exp6 = 2 + x + 2 + x;
  expr exp7 = x + 3 + y + 5;
  expr exp8 = x + 1 + 2 + x + 3 + y + 4 + z + 5;
  expr exp9 = expr(1) + expr(2) + expr(3) + expr(4) + expr(5) + expr(6) +
              expr(7) + expr(8) + expr(9);
  expr exp10 = create(kind::ADD, {2 * (x + 1) + 3 * x + 4 * (x + 1) + 4,
                                  create(kind::ADD, {3 + x + y, 3 * x + 4})});

  expr exp11 = x + y + z + -x;
  expr exp12 = -x + 2 * x + x + -x + y + -y;
  expr exp13 = inf() + x + y + 10;
  expr exp14 = -inf() + x + y + z + 14;
  expr exp15 = x + y + inf() + -14 + -inf();

  assert(reduce(exp0) == 4);
  assert(reduce(exp1) == 17);
  assert(reduce(exp2) == 14);
  assert(reduce(exp3) == 14);
  assert(reduce(exp4) == 14);
  assert(reduce(exp5) == 2 * x);
  assert(reduce(exp6) == 2 * x + 4);
  assert(reduce(exp7) == x + y + 8);
  assert(reduce(exp8) == 2 * x + y + z + 15);
  assert(reduce(exp9) == 45);
  assert(reduce(exp10) == 7 * x + y + 6 * (x + 1) + 11);
  assert(reduce(exp11) == y + z);
  assert(reduce(exp12) == x);
  assert(reduce(exp13) == inf());
  assert(reduce(exp14) == -inf());
  assert(reduce(exp15) == undefined());
}

void should_simplify_products() {
  expr x = expr("x");
  expr y = expr("y");

  expr exp0 = expr(2) * expr(2);
  expr exp1 = expr(kind::MUL, {expr(2) * expr(3), expr(4) * expr(5)});
  expr exp2 = x * x;
  expr exp3 = expr(kind::MUL, {2 * x, 2 * x});
  expr exp4 = 2 * x * 4 * y * 6 * x * 14 * y * 10 * x;

  assert(reduce(exp0) == 4);
  assert(reduce(exp1) == 120);
  assert(reduce(exp2) == pow(expr("x"), 2));
  assert(reduce(exp3) == 4 * pow(expr("x"), 2));
  assert(reduce(exp4) == 6720 * pow(x, 3) * pow(y, 2));
}

void should_simplify_subtractions() {
  expr x = expr("x");
  expr a = expr("a");
  expr b = expr("b");
  expr c = expr("c");
  expr d = expr("d");
  expr e = expr("e");
  expr f = expr("f");
  expr g = expr("g");

  expr exp0 = expr(3) - expr(2);
  expr exp1 = expr(2) - expr(3);
  expr exp2 = x - x;
  expr exp3 = (a + b + c) - (d - e);
  expr exp4 = expr(kind::SUB, {
                                          a + b + c,
                                          d - e,
                                          f - g,
                                      });

  expr exp5 = expr(kind::SUB, {
                                          a - b - c,
                                          d - e,
                                          f - g,
		});
  expr exp6 = expr(1) - expr(2) - expr(3) - expr(5) - expr(7) - x - expr(4) -
              expr(6) - a;

  assert(reduce(exp0) == 1);
  assert(reduce(exp1) == -1);
  assert(reduce(exp2) == 0);
  assert(reduce(exp3) == a + b + c + -d + e);
  assert(reduce(exp4) == a + b + c + -d + e + -f + g);
  assert(reduce(exp5) == a + -b + -c + -d + e + -f + g);
  assert(reduce(exp6) == -a + -x + -26);
}

void should_simplify_pows() {
  expr exp0 = pow(integer(2), integer(2));
  expr exp1 = pow(fraction(1, 2), integer(2));
  expr exp2 = pow(integer(2), integer(0));
  expr exp3 = pow(integer(2), integer(1));
  expr exp4 = pow(symbol("x"), integer(0));
  expr exp5 = pow(integer(0), integer(0));
  expr exp6 = pow(integer(0), integer(3));

  expr res_exp0 = reduce(exp0);
  expr res_exp1 = reduce(exp1);
  expr res_exp2 = reduce(exp2);
  expr res_exp3 = reduce(exp3);
  expr res_exp4 = reduce(exp4);
  expr res_exp5 = reduce(exp5);
  expr res_exp6 = reduce(exp6);

  assert(res_exp0 == 4);
  assert(res_exp1 == fraction(1, 4));
  assert(res_exp2 == 1);
  assert(res_exp3 == 2);
  assert(res_exp4 == 1);
  assert(res_exp5 == undefined());
  assert(res_exp6 == 0);
}

void should_simplify_divisions() {
  expr exp0 = expr(4) / expr(2);
  expr exp1 = expr("x") / expr("x");
  expr exp2 = (expr(2) * expr("x")) / expr("x");
  expr exp3 = pow(expr("x"), 2) / expr("x");

  expr res_exp0 = reduce(exp0);
  expr res_exp1 = reduce(exp1);
  expr res_exp2 = reduce(exp2);
  expr res_exp3 = reduce(exp3);

  assert(res_exp0.kind() == kind::INT);
  assert(res_exp0.value() == 2);

  assert(res_exp1.kind() == kind::INT);
  assert(res_exp1.value() == 1);

  assert(res_exp2.kind() == kind::INT);
  assert(res_exp2.value() == 2);

  assert(res_exp3.kind() == kind::SYM);
  assert(res_exp3.identifier() == "x");
}

void should_simplify_expressions_matrix() {
	expr A = mat(3, 3, {1, 1, 1, 2, 2, 2, 3, 3, 3});
	expr B = mat(3, 3,
							 {1, 1, 1, 2, 2, 2, 3, 3, 3} );

	expr C = mat(5, 4, {1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5});

	expr D = mat(4, 5, {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4});

	expr E = mat(5, 5, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });
	expr F = mat(4, 4);

	expr G = mat(5, 5, {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });

	expr t = (A + C + D + B + E + F + G);

	expr k = reduce(t);

	printf("%s\n", to_string(k).c_str());

	// expr H = mat(4, 2);
	// expr I = mat(2, 4);
	// expr J = mat(3, 3);
	// expr K = mat(7, 3);
	// expr L = mat(2, 2);
	// expr M = mat(3, 3);
	// expr N = mat(1, 1);

	// expr o = H * J * I * L * K * N * M;

	// expr j = reduce(o);


	// printf("%s\n", to_string(j).c_str());

	expr l = svd_matrix(A);
	printf("%s\n", to_string(l).c_str());
}

int main() {
  TEST(should_construct_expr)
  TEST(should_eval_equality)
  TEST(should_insert_and_remove_from_expr)
  TEST(should_eval_exprs)
  TEST(should_expand_expr)
  TEST(should_perform_list_operations);
  TEST(should_perform_set_operations)
  TEST(should_simplify_additions)
  TEST(should_simplify_products)
  TEST(should_simplify_subtractions)
  TEST(should_simplify_pows)
  TEST(should_simplify_divisions)
	TEST(should_simplify_expressions_matrix)
  return 0;
}
