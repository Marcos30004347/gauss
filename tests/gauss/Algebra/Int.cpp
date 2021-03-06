#include "gauss/Algebra/Int.hpp"
#include "test.hpp"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

void should_get_quotient_of_div_by_powers_of_two() {
  assert(quoPow2(2, 1) == 1);
  assert(quoPow2(4, 1) == 2);
  assert(quoPow2(4, 2) == 1);
  assert(quoPow2(4, 3) == 0);
  assert(quoPow2(8, 2) == 2);
  assert(quoPow2(8, 3) == 1);
  assert(quoPow2(9, 1) == 4);
  assert(quoPow2(7, 2) == 1);
  assert(quoPow2(7, 3) == 0);
}

void should_get_remainder_of_div_by_powers_of_two() {
  assert(modPow2(2, 1) == 0);
  assert(modPow2(4, 1) == 0);
  assert(modPow2(4, 2) == 0);
  assert(modPow2(4, 3) == 4);
  assert(modPow2(8, 2) == 0);
  assert(modPow2(8, 3) == 0);
  assert(modPow2(9, 1) == 1);
  assert(modPow2(7, 2) == 3);
  assert(modPow2(7, 3) == 7);
}

void should_create_bints_from_types() {
  bint<1, short>* v0 = bint<1, short>::from<int>(123);

  assert(v0->digit[0] == 1);
	assert(v0->digit[1] == 1);
	assert(v0->digit[2] == 0);
	assert(v0->digit[3] == 1);
	assert(v0->digit[4] == 1);
	assert(v0->digit[5] == 1);
	assert(v0->digit[6] == 1);

	delete v0;

	bint<30>* v1 = bint<30>::from<int>(123);
	assert(v1->digit[0] == 123);

	delete v1;
}

void should_abs_add_digits_bints() {
	bint<1>* v0 = bint<1>::from<int>(2);
  bint<1>* v1 = bint<1>::from<int>(2);
	bint<1>* v2 = new bint<1>();

	bint<1>::abs_add_digits(v0, v1, v2);

	assert(v2->digit[0] == 0);
	assert(v2->digit[1] == 0);
	assert(v2->digit[2] == 1);

	delete v0;
	delete v1;
	delete v2;

	bint<3>* v3 = bint<3>::from<int>(923456);
  bint<3>* v4 = bint<3>::from<int>(187654);
	bint<3>* v5 = new bint<3>();

	bint<3>::abs_add_digits(v3, v4, v5);

	assert(v5->digit[0] == 6);
	assert(v5->digit[1] == 0);
	assert(v5->digit[2] == 1);
	assert(v5->digit[3] == 2);
	assert(v5->digit[4] == 7);
	assert(v5->digit[5] == 1);
	assert(v5->digit[6] == 4);

	delete v3;
	delete v4;
	delete v5;
}

void should_abs_sub_digits_bints() {
	bint<1>* v0 = bint<1>::from<int>(3);
  bint<1>* v1 = bint<1>::from<int>(2);
	bint<1>* v2 = new bint<1>();

	bint<1>::abs_sub_digits(v0, v1, v2);

	assert(v2->digit[0] == 1);

	delete v0;
	delete v1;
	delete v2;

	bint<3>* v3 = bint<3>::from<int>(923456);
  bint<3>* v4 = bint<3>::from<int>(187654);
	bint<3>* v5 = new bint<3>();

	bint<3>::abs_sub_digits(v3, v4, v5);

	assert(v5->digit[0] == 2);
	assert(v5->digit[1] == 7);
	assert(v5->digit[2] == 0);
	assert(v5->digit[3] == 5);
	assert(v5->digit[4] == 3);
	assert(v5->digit[5] == 6);
	assert(v5->digit[6] == 2);

	delete v3;
	delete v4;
	delete v5;
}

void should_abs_mul_digits_bints() {
	bint<1>* v0 = bint<1>::from<int>(3);
  bint<1>* v1 = bint<1>::from<int>(2);
	bint<1>* v2 = new bint<1>();

	bint<1>::abs_mul_digits(v0, v1, v2);

	assert(v2->digit[0] == 0);
	assert(v2->digit[1] == 1);
	assert(v2->digit[2] == 1);

	delete v0;
	delete v1;
	delete v2;

	bint<3>* v3 = bint<3>::from<int>(923456);
  bint<3>* v4 = bint<3>::from<int>(187654);
	bint<3>* v5 = new bint<3>();

	bint<3>::abs_mul_digits(v3, v4, v5);

	assert(v5->digit[0] == 0);
	assert(v5->digit[1] == 0);
	assert(v5->digit[2] == 6);
	assert(v5->digit[3] == 5);
	assert(v5->digit[4] == 4);
	assert(v5->digit[5] == 5);
	assert(v5->digit[6] == 1);
	assert(v5->digit[7] == 7);
	assert(v5->digit[8] == 0);
	assert(v5->digit[9] == 3);
	assert(v5->digit[10] == 1);
	assert(v5->digit[11] == 4);
	assert(v5->digit[12] == 2);

	delete v3;
	delete v4;
	delete v5;
}

void should_abs_square_digits_bints() {

	bint<1>* v0 = bint<1>::from<int>(3);
	bint<1>* v1 = new bint<1>();

	bint<1>::abs_square_digits(v0, v1);

	assert(v1->digit[0] == 1);
	assert(v1->digit[1] == 0);
	assert(v1->digit[2] == 0);
	assert(v1->digit[3] == 1);

	delete v0;
	delete v1;

	bint<3>* v2 = bint<3>::from<int>(923456);
	bint<3>* v3 = new bint<3>();

	bint<3>::abs_square_digits(v2, v3);

	assert(v3->digit[0] == 0);
	assert(v3->digit[1] == 0);
	assert(v3->digit[2] == 0);
	assert(v3->digit[3] == 0);
	assert(v3->digit[4] == 1);
	assert(v3->digit[5] == 1);
	assert(v3->digit[6] == 7);
	assert(v3->digit[7] == 0);
	assert(v3->digit[8] == 5);
	assert(v3->digit[9] == 1);
	assert(v3->digit[10] == 2);
	assert(v3->digit[11] == 3);
	assert(v3->digit[12] == 4);
	assert(v3->digit[13] == 1);

	delete v2;
	delete v3;
}

void should_add_bints() {
	bint<3>* v0 = bint<3>::from(10);
	bint<3>* v1 = bint<3>::from(10);
	bint<3>* v2 = new bint<3>();

	bint<3>::add(v0, v1, v2);

	bint<3>* r0 = bint<3>::from(20);
	assert(bint<3>::compare(v2, r0) == 0);

	delete v0;
	delete v1;
	delete v2;
	delete r0;

	bint<3>* v3 = bint<3>::from(-2);
	bint<3>* v4 = bint<3>::from(-2);
	bint<3>* v5 = new bint<3>();

	bint<3>::add(v3, v4, v5);
	bint<3>* r1 = bint<3>::from(-4);

	assert(bint<3>::compare(v5, r1) == 0);

	delete v3;
	delete v4;
	delete v5;
	delete r1;

	bint<3>* v6 = bint<3>::from(5);
	bint<3>* v7 = bint<3>::from(-2);
	bint<3>* v8 = new bint<3>();

	bint<3>::add(v6, v7, v8);
	bint<3>* r2 = bint<3>::from(3);
	assert(bint<3>::compare(v8, r2) == 0);

	delete v6;
	delete v7;
	delete v8;
	delete r2;

	bint<3>* v9 = bint<3>::from(-2);
	bint<3>* v10 = bint<3>::from(5);
	bint<3>* v11 = new bint<3>();

	bint<3>::add(v9, v10, v11);
	bint<3>* r3 = bint<3>::from(3);

	assert(bint<3>::compare(v11, r3) == 0);

	delete v9;
	delete v10;
	delete v11;
	delete r3;

	bint<3>* v12 = bint<3>::from(10);
	bint<3>* v13 = bint<3>::from(3);
	bint<3>* v14 = new bint<3>();

	bint<3>::add(v12, v13, v14);
	bint<3>* r4 = bint<3>::from(13);
	assert(bint<3>::compare(v14, r4) == 0);

	delete v12;
	delete v13;
	delete v14;
	delete r4;

	bint<3>* v15 = bint<3>::from(3);
	bint<3>* v16 = bint<3>::from(10);
	bint<3>* v17 = new bint<3>();

	bint<3>::add(v15, v16, v17);
	bint<3>* r5 = bint<3>::from(13);
	assert(bint<3>::compare(v17, r5) == 0);

	delete v15;
	delete v16;
	delete v17;
	delete r5;
}

void should_digits_shift_big_ints() {
	bint<1>::digit_t carry;
	bint<1>* c = new bint<1>();

	bint<1>* b = bint<1>::from(5);

	c->resize(b->size);

	carry = bint<1>::digits_lshift(b->digit, b->size, 1, c->digit);

	assert(carry == 1);
	assert(c->digit[0] == 0);
	assert(c->digit[1] == 1);
	assert(c->digit[2] == 0);

	carry = bint<1>::digits_rshift(b->digit, b->size, 2, c->digit);

	assert(carry == 1);
	assert(c->digit[0] == 1);
	assert(c->digit[1] == 0);
	assert(c->digit[2] == 0);

	delete b;
	delete c;
}

void should_divide_big_ints() {
	bint<1>* a = bint<1>::from(10);
	bint<1>* b = bint<1>::from(4);
	bint<1>* q0 = new bint<1>();
	bint<1> *r0 = new bint<1>();

	bint<1>::div(a, b, q0, r0);

	assert(q0->digit[0] == 0);
	assert(q0->digit[1] == 1);
	assert(r0->digit[0] == 0);
	assert(r0->digit[1] == 1);

	delete a;
	delete b;
	delete q0;
	delete r0;

	bint<2>* c = bint<2>::from(100);
	bint<2>* d = bint<2>::from(20);
	bint<2>* q1 = new bint<2>();
	bint<2>* r1 = new bint<2>();

	bint<2>::div(c, d, q1, r1);

	assert(q1->digit[0] == 1);
	assert(q1->digit[1] == 1);
	assert(r1->size == 0);

	delete c;
	delete d;
	delete q1;
	delete r1;

	bint<30>* e = bint<30>::from(1021);
	bint<30>* f = bint<30>::from(100);
	bint<30>* q2 = new bint<30>();
	bint<30>* r2 = new bint<30>();

	bint<30>::div(e, f, q2, r2);

	assert(q2->digit[0] == 10);
	assert(r2->digit[0] == 21);

	delete e;
	delete f;
	delete q2;
	delete r2;

	bint<2>* g = bint<2>::from(100);
	bint<2>* h = bint<2>::from(3);
	bint<2>* q3 = new bint<2>();
	bint<2>* r3 = new bint<2>();

	bint<2>::div(g, h, q3, r3);

	assert(q3->digit[0] == 1);
	assert(q3->digit[1] == 0);
	assert(q3->digit[2] == 2);
	assert(r3->digit[0] == 1);

	delete g;
	delete h;
	delete q3;
	delete r3;

	bint<30>* i = bint<30>::from(3331244669);
	bint<30>* j = bint<30>::from(1801814625);
  bint<30>* q4 = new bint<30>();
	bint<30>* r4 = new bint<30>();

	bint<30>::div(i, j, q4, r4);

	assert(q4->digit[0] == 1);
	assert(r4->digit[0] == 455688220);
	assert(r4->digit[1] == 1);

	delete i;
	delete j;
	delete q4;
	delete r4;
}

#include <float.h>
#include <cfloat>
#include <cmath>

void should_convert_big_int_to_double() {
	double v;
	bint<1>* a = bint<1>::from(3);
	bint<1>::to_double(a, &v);
	assert(std::abs(v - 3.0) <= std::numeric_limits<double>::epsilon());

	delete a;

	bint<1>* b = bint<1>::from(10002312);
	bint<1>::to_double(b, &v);
	assert(std::abs(v - 10002312.0) <= std::numeric_limits<double>::epsilon());
	delete b;

	bint<5>* c = bint<5>::from(123218312);
	bint<5>::to_double(c, &v);
	assert(std::abs(v - 123218312.0) <= std::numeric_limits<double>::epsilon());
	delete c;

	bint<30>* d = bint<30>::from(std::numeric_limits<long long>::max());
	bint<30>::to_double(d, &v);
	assert(std::abs(v - (double)std::numeric_limits<long long>::max()) <= std::numeric_limits<double>::epsilon());
	delete d;
}

void should_get_frexp_of_bints() {
	int e;
	int e0;

	bool overflow;

	bint<1>* b0 = bint<1>::from(1230);
	double fr0 = bint<1>::frexp(b0, &e, &overflow);

	assert(e == 11);
	assert(overflow == false);
	assert(std::abs(fr0 - std::frexp(1230, &e0)) <= std::numeric_limits<double>::epsilon());

	delete b0;
}

void should_convert_big_int_to_long_long() {
	long long v;
	bint<1>* a = bint<1>::from(3);

	bint<1>::to_long(a, &v);
	assert(v == 3);
	delete a;

	bint<1>* b = bint<1>::from(10002312);
	bint<1>::to_long(b, &v);
	assert(v == 10002312);
	delete b;

	bint<5>* c = bint<5>::from(123218312);
	bint<5>::to_long(c, &v);
	assert(v == 123218312);
	delete c;
}

void should_elevate_big_int_to_expoent() {
	bint<1>* a = bint<1>::from(4);
	bint<1>* b = bint<1>::from(2);
	bint<1>* c = bint<1>::pow(a, b);

	assert(c->digit[0] == 0);
	assert(c->digit[1] == 0);
	assert(c->digit[2] == 0);
	assert(c->digit[3] == 0);
	assert(c->digit[4] == 1);

	delete a;
	delete b;
	delete c;

	bint<4>* d = bint<4>::from(101);
	bint<4>* e = bint<4>::from(3);
	bint<4>* f = bint<4>::pow(d, e);

	assert(f->digit[0] == 13);
	assert(f->digit[1] == 9);
	assert(f->digit[2] == 8);
	assert(f->digit[3] == 11);
	assert(f->digit[4] == 15);

	delete d;
	delete e;
	delete f;
}

void should_convert_numbers_to_string() {
	bint<1>* a = bint<1>::from(10);
	assert(a->to_string() ==  "10");
	delete a;

	bint<1>* b = bint<1>::from(101231232);
	assert(b->to_string() ==  "101231232");
	delete b;

	bint<1>* c = bint<1>::from(101212312312331232);
	assert(c->to_string() ==  "101212312312331232");
	delete c;

	bint<4>* d = bint<4>::from(91322312354552);
	assert(d->to_string() ==  "91322312354552");
	delete d;
}

void should_get_ceil_log2() {
	bint<1>* a = bint<1>::from((1 << 4));
	bint<1>* b = bint<1>::from((1 << 4) + 1);
	bint<30>* c = bint<30>::from(1234215);
	bint<30>* d = bint<30>::from(std::numeric_limits<long long>::max());

	bint<1>* al2 = bint<1>::ceil_log2(a);
	bint<1>* bl2 = bint<1>::ceil_log2(b);
	bint<30>* cl2 = bint<30>::ceil_log2(c);
	bint<30>* dl2 = bint<30>::ceil_log2(d);

	assert(al2->digit[0] == 0);
	assert(al2->digit[1] == 0);
	assert(al2->digit[2] == 1);

	assert(bl2->digit[0] == 0);
	assert(bl2->digit[1] == 0);
	assert(bl2->digit[2] == 1);

	assert(cl2->digit[0] == 21);
	assert(dl2->digit[0] == 63);

	delete a;
	delete b;
	delete c;
	delete d;
	delete al2;
	delete bl2;
	delete cl2;
	delete dl2;
}

void should_shift_bints() {
	bint<1>* a = bint<1>::from((1 << 12) + (1 << 7));

	bint<1>* b = bint<1>::lshift(a, 5);

	assert(b->digit[0] == 0);
	assert(b->digit[1] == 0);
	assert(b->digit[2] == 0);
	assert(b->digit[3] == 0);
	assert(b->digit[4] == 0);
	assert(b->digit[5] == 0);
	assert(b->digit[6] == 0);
	assert(b->digit[7] == 0);
	assert(b->digit[8] == 0);
	assert(b->digit[9] == 0);
	assert(b->digit[10] == 0);
	assert(b->digit[11] == 0);
	assert(b->digit[12] == 1);
	assert(b->digit[13] == 0);
	assert(b->digit[14] == 0);
	assert(b->digit[15] == 0);
	assert(b->digit[16] == 0);
	assert(b->digit[17] == 1);


	bint<1>* c = bint<1>::rshift(b, 5);

	assert(c->digit[0] == 0);
	assert(c->digit[1] == 0);
	assert(c->digit[2] == 0);
	assert(c->digit[3] == 0);
	assert(c->digit[4] == 0);
	assert(c->digit[5] == 0);
	assert(c->digit[6] == 0);
	assert(c->digit[7] == 1);
	assert(c->digit[8] == 0);
	assert(c->digit[9] == 0);
	assert(c->digit[10] == 0);
	assert(c->digit[11] == 0);
	assert(c->digit[12] == 1);

	bint<30>* d = bint<30>::from(4);

	bint<30>* e = bint<30>::lshift(d, 4);

	assert(e->digit[0] == (1 << 6));

	bint<30>* f = bint<30>::rshift(e, 4);

	assert(f->digit[0] == (1 << 2));

	delete a;
	delete b;
	delete c;
	delete d;
	delete e;
	delete f;
}

void should_get_sqrt_of_bints() {
	bint<30>* v = bint<30>::from(100042);
	bint<30>* s = new bint<30>();
	bint<30>* r = new bint<30>();

	bint<30>::isqrt(v, s, r);

	assert(s->digit[0] == 316);
	assert(r->digit[0] == 186);

	delete v;
	delete s;
	delete r;
}


int main() {
	TEST(should_get_quotient_of_div_by_powers_of_two)
  TEST(should_get_remainder_of_div_by_powers_of_two)
  TEST(should_create_bints_from_types)
	TEST(should_abs_add_digits_bints)
	TEST(should_abs_sub_digits_bints)
	TEST(should_abs_mul_digits_bints)
	TEST(should_abs_square_digits_bints)
	TEST(should_add_bints)
	TEST(should_digits_shift_big_ints)
	TEST(should_divide_big_ints)
	TEST(should_convert_big_int_to_double)
	TEST(should_convert_big_int_to_long_long)
	TEST(should_get_frexp_of_bints)
	TEST(should_elevate_big_int_to_expoent)
	TEST(should_convert_numbers_to_string)
	TEST(should_get_ceil_log2)
	TEST(should_shift_bints)
	TEST(should_get_sqrt_of_bints)
}
