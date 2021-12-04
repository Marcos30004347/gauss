
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

void should_abs_add_digits_bints() {
	bint<1> v0 = bint<1>::from<int>(2);
  bint<1> v1 = bint<1>::from<int>(2);
	bint<1> v2;

	bint<1>::abs_add_digits(&v0, &v1, &v2);

	assert(v2.digit[0] == 0);
	assert(v2.digit[1] == 0);
	assert(v2.digit[2] == 1);

	bint<3> v3 = bint<3>::from<int>(923456);
  bint<3> v4 = bint<3>::from<int>(187654);
	bint<3> v5;
	bint<3>::abs_add_digits(&v3, &v4, &v5);

	assert(v5.digit[0] == 6);
	assert(v5.digit[1] == 0);
	assert(v5.digit[2] == 1);
	assert(v5.digit[3] == 2);
	assert(v5.digit[4] == 7);
	assert(v5.digit[5] == 1);
	assert(v5.digit[6] == 4);

}

void should_abs_sub_digits_bints() {
	bint<1> v0 = bint<1>::from<int>(3);
  bint<1> v1 = bint<1>::from<int>(2);
	bint<1> v2;

	bint<1>::abs_sub_digits(&v0, &v1, &v2);

	assert(v2.digit[0] == 1);
	assert(v2.digit[1] == 0);
	assert(v2.digit[2] == 0);

	bint<3> v3 = bint<3>::from<int>(923456);
  bint<3> v4 = bint<3>::from<int>(187654);
	bint<3> v5;

	bint<3>::abs_sub_digits(&v3, &v4, &v5);

	assert(v5.digit[0] == 2);
	assert(v5.digit[1] == 7);
	assert(v5.digit[2] == 0);
	assert(v5.digit[3] == 5);
	assert(v5.digit[4] == 3);
	assert(v5.digit[5] == 6);
	assert(v5.digit[6] == 2);
}

void should_abs_mul_digits_bints() {
	bint<1> v0 = bint<1>::from<int>(3);
  bint<1> v1 = bint<1>::from<int>(2);
	bint<1> v2;

	bint<1>::abs_mul_digits(&v0, &v1, &v2);

	assert(v2.digit[0] == 0);
	assert(v2.digit[1] == 1);
	assert(v2.digit[2] == 1);
	assert(v2.digit[3] == 0);
	assert(v2.digit[4] == 0);
	assert(v2.digit[5] == 0);


	bint<3> v3 = bint<3>::from<int>(923456);
  bint<3> v4 = bint<3>::from<int>(187654);
	bint<3> v5;

	bint<3>::abs_mul_digits(&v3, &v4, &v5);

	assert(v5.digit[0] == 0);
	assert(v5.digit[1] == 0);
	assert(v5.digit[2] == 6);
	assert(v5.digit[3] == 5);
	assert(v5.digit[4] == 4);
	assert(v5.digit[5] == 5);
	assert(v5.digit[6] == 1);
	assert(v5.digit[7] == 7);
	assert(v5.digit[8] == 0);
	assert(v5.digit[9] == 3);
	assert(v5.digit[10] == 1);
	assert(v5.digit[11] == 4);
	assert(v5.digit[12] == 2);
}

void should_abs_square_digits_bints() {

	bint<1> v0 = bint<1>::from<int>(3);
	bint<1> v1;

	bint<1>::abs_square_digits(&v0, &v1);

	assert(v1.digit[0] == 1);
	assert(v1.digit[1] == 0);
	assert(v1.digit[2] == 0);
	assert(v1.digit[3] == 1);

	bint<3> v2 = bint<3>::from<int>(923456);
	bint<3> v3;

	bint<3>::abs_square_digits(&v2, &v3);

	assert(v3.digit[0] == 0);
	assert(v3.digit[1] == 0);
	assert(v3.digit[2] == 0);
	assert(v3.digit[3] == 0);
	assert(v3.digit[4] == 1);
	assert(v3.digit[5] == 1);
	assert(v3.digit[6] == 7);
	assert(v3.digit[7] == 0);
	assert(v3.digit[8] == 5);
	assert(v3.digit[9] == 1);
	assert(v3.digit[10] == 2);
	assert(v3.digit[11] == 3);
	assert(v3.digit[12] == 4);
	assert(v3.digit[13] == 1);
}

void should_add_bints() {
	bint<3> v0 = bint<3>::from(10);
	bint<3> v1 = bint<3>::from(10);
	bint<3> v2;

	bint<3>::add(&v0, &v1, &v2);
	bint<3> r0 = bint<3>::from(20);
	assert(bint<3>::compare(&v2, &r0) == 0);

	bint<3> v3 = bint<3>::from(-2);
	bint<3> v4 = bint<3>::from(-2);
	bint<3> v5;

	bint<3>::add(&v3, &v4, &v5);
	bint<3> r1 = bint<3>::from(-4);
	assert(bint<3>::compare(&v5, &r1) == 0);

	bint<3> v6 = bint<3>::from(5);
	bint<3> v7 = bint<3>::from(-2);
	bint<3> v8;

	bint<3>::add(&v6, &v7, &v8);
	bint<3> r2 = bint<3>::from(3);
	assert(bint<3>::compare(&v8, &r2) == 0);

	bint<3> v9 = bint<3>::from(-2);
	bint<3> v10 = bint<3>::from(5);
	bint<3> v11;

	bint<3>::add(&v9, &v10, &v11);
	v11.printRep();
	bint<3> r3 = bint<3>::from(3);
	assert(bint<3>::compare(&v11, &r3) == 0);

	bint<3> v12 = bint<3>::from(10);
	bint<3> v13 = bint<3>::from(3);
	bint<3> v14;

	bint<3>::add(&v12, &v13, &v14);
	bint<3> r4 = bint<3>::from(13);
	assert(bint<3>::compare(&v14, &r4) == 0);

	bint<3> v15 = bint<3>::from(3);
	bint<3> v16 = bint<3>::from(10);
	bint<3> v17;

	bint<3>::add(&v15, &v16, &v17);
	bint<3> r5 = bint<3>::from(13);
	assert(bint<3>::compare(&v17, &r5) == 0);
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
}
