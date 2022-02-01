#include <cstdlib>

#define TEST_TIME_REPORT_UNIT TEST_TIME_REPORT_NS

#include "test.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>

#include "Core/AST/AST3.hpp"

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

  expr f = (-32 * pow(z, 3) + 32 * pow(z, 4) + 48 * pow(z, 5) + -24 * pow(z, 6) +
           -48 * pow(z, 7) + -36 * pow(z, 8) + -40 * pow(z, 9) +
           -8 * pow(z, 10) + -8 * pow(z, 11)) /
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

	t.append(4);

	assert(t.size() == 4);

	assert(t[3] == 4);

	list g = append(t,  5);

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

	assert(t[0] == 3);
	assert(t[1] == 2);
	assert(t[2] == 1);

	set k = { 4, 4, 5 };

	set u = unnification(t, k);

	assert(u.size() == 5);

	assert(u[0] == 5);
	assert(u[1] == 4);
	assert(u[2] == 3);
	assert(u[3] == 2);
	assert(u[4] == 1);

	set h = difference(u, k);

	assert(h.size() == 3);

	assert(h[0] == 3);
	assert(h[1] == 2);
	assert(h[2] == 1);

	set l = { 2, 1 };

	set v = intersection(h, l);

	assert(v.size() == 2);

	assert(v[0] == 2);
	assert(v[1] == 1);
}


int main() {
  TEST(should_construct_expr)
  TEST(should_eval_equality)
  TEST(should_insert_and_remove_from_expr)
  TEST(should_eval_exprs)
  TEST(should_expand_expr)
	TEST(should_perform_list_operations);
	TEST(should_perform_set_operations)

	return 0;
}
