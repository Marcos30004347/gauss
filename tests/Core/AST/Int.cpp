
#include "Core/AST/Int.hpp"
#include "test.hpp"
#include <cstdint>

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
  bint<1, short> v0 = bint<1, short>::from<int>(123);

  assert(v0.digit[0] == 1);
	assert(v0.digit[1] == 1);
	assert(v0.digit[2] == 0);
	assert(v0.digit[3] == 1);
	assert(v0.digit[4] == 1);
	assert(v0.digit[5] == 1);
	assert(v0.digit[6] == 1);

	bint<30> v1 = bint<30>::from<int>(123);
	assert(v1.digit[0] == 123);
}


void should_add_bints() {
	bint<1> v0 = bint<1>::from<int>(2);
  bint<1> v1 = bint<1>::from<int>(2);

	bint<1>::digit_t* a = new bint<1>::digit_t[v0.size + 1];
	bint<1>::add(v0.digit, v0.size, v1.digit, v1.size, a);

	assert(a[0] == 0);
	assert(a[1] == 0);
	assert(a[2] == 1);

	bint<3> v2 = bint<3>::from<int>(923456);
  bint<3> v3 = bint<3>::from<int>(187654);

	bint<3>::digit_t *b = new bint<3>::digit_t[v2.size + 1];
	bint<3>::add(v2.digit, v2.size, v3.digit, v3.size, b);

	assert(b[0] == 6);
	assert(b[1] == 0);
	assert(b[2] == 1);
	assert(b[3] == 2);
	assert(b[4] == 7);
	assert(b[5] == 1);
	assert(b[6] == 4);

	delete[] a;
	delete[] b;
}

void should_sub_bints() {
	bint<1> v0 = bint<1>::from<int>(3);
  bint<1> v1 = bint<1>::from<int>(2);

	bint<1>::digit_t* a = new bint<1>::digit_t[v0.size + 1];
	bint<1>::sub(v0.digit, v0.size, v1.digit, v1.size, a);

	assert(a[0] == 1);
	assert(a[1] == 0);
	assert(a[2] == 0);

	bint<3> v2 = bint<3>::from<int>(923456);
  bint<3> v3 = bint<3>::from<int>(187654);

	bint<3>::digit_t *b = new bint<3>::digit_t[v2.size];
	bint<3>::sub(v2.digit, v2.size, v3.digit, v3.size, b);

	assert(b[0] == 2);
	assert(b[1] == 7);
	assert(b[2] == 0);
	assert(b[3] == 5);
	assert(b[4] == 3);
	assert(b[5] == 6);
	assert(b[6] == 2);

	delete[] a;
	delete[] b;
}

void should_mul_bints() {
	bint<1> v0 = bint<1>::from<int>(3);
  bint<1> v1 = bint<1>::from<int>(2);

	bint<1>::digit_t* a = new bint<1>::digit_t[v0.size + v1.size];
	bint<1>::mul(v0.digit, v0.size, v1.digit, v1.size, a);

	assert(a[0] == 0);
	assert(a[1] == 1);
	assert(a[2] == 1);
	assert(a[3] == 0);
	assert(a[4] == 0);
	assert(a[5] == 0);

	delete[] a;

	bint<3> v2 = bint<3>::from<int>(923456);
  bint<3> v3 = bint<3>::from<int>(187654);

	bint<3>::digit_t *b = new bint<3>::digit_t[v2.size + v3.size];
	bint<3>::mul(v2.digit, v2.size, v3.digit, v3.size, b);

	assert(b[0] == 0);
	assert(b[1] == 0);
	assert(b[2] == 6);
	assert(b[3] == 5);
	assert(b[4] == 4);
	assert(b[5] == 5);
	assert(b[6] == 1);
	assert(b[7] == 7);
	assert(b[8] == 0);
	assert(b[9] == 3);
	assert(b[10] == 1);
	assert(b[11] == 4);
	assert(b[12] == 2);

	delete[] b;
}

void should_sqr_bints() {

	bint<1> v0 = bint<1>::from<int>(3);
	bint<1>::digit_t* a = new bint<1>::digit_t[2*v0.size];
	bint<1>::square(v0.digit, v0.size, a);

	assert(a[0] == 1);
	assert(a[1] == 0);
	assert(a[2] == 0);
	assert(a[3] == 1);

	delete[] a;

	bint<3> v2 = bint<3>::from<int>(923456);
	bint<3>::digit_t* b = new bint<3>::digit_t[2*v2.size];
	bint<3>::square(v2.digit, v2.size, b);

	assert(b[0] == 0);
	assert(b[1] == 0);
	assert(b[2] == 0);
	assert(b[3] == 0);
	assert(b[4] == 1);
	assert(b[5] == 1);
	assert(b[6] == 7);
	assert(b[7] == 0);
	assert(b[8] == 5);
	assert(b[9] == 1);
	assert(b[10] == 2);
	assert(b[11] == 3);
	assert(b[12] == 4);
	assert(b[13] == 1);

	delete[] b;
}

int main() {
  TEST(should_get_quotient_of_div_by_powers_of_two)
  TEST(should_get_remainder_of_div_by_powers_of_two)
  TEST(should_create_bints_from_types)
	TEST(should_add_bints)
	TEST(should_sub_bints)
	TEST(should_mul_bints)
	TEST(should_sqr_bints)
}
