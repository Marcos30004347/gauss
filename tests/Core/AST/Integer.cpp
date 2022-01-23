#include "Core/AST/Integer.hpp"
#include "test.hpp"
#include <cassert>

void should_add_ints() {
	assert(Int(1) + Int(0) == Int(1));
	assert(Int(-1) + Int(1) == Int(0));
	assert(Int(1) + Int(-1) == Int(0));
	assert(Int(2000) + Int(2000) == Int(4000));
	assert(Int(-2000) + Int(2000) == Int(0));
	assert(Int(2000) + Int(-2000) == Int(0));
	assert(Int(4000) + Int(-200) == Int(3800));
	assert(Int(-2000) + Int(-2000) == Int(-4000));
	assert(Int(-2000) + Int(-1000) == Int(-3000));
	assert(Int(2000) + Int(-2000) == Int(0));
}

void should_sub_ints() {
	assert(Int(1) - Int(0) == Int(1));
	assert(Int(-1) - Int(1) == Int(-2));
	assert(Int(1) - Int(-1) == Int(2));
	assert(Int(2000) - Int(2000) == Int(0));
	assert(Int(2000) - Int(-2000) == Int(4000));
	assert(Int(-2000) - Int(2000) == Int(-4000));
	assert(Int(4000) - Int(-200) == Int(4200));
	assert(Int(-2000) - Int(-2000) == Int(0));
	assert(Int(-2000) - Int(-1000) == Int(-1000));
}

void should_mul_ints() {
	assert(Int(1) * Int(1) == Int(1));
	assert(Int(-1) * Int(1) == Int(-1));
	assert(Int(1) * Int(-1) == Int(-1));
	assert(Int(-1) * Int(-1) == Int(1));
	assert(Int(2) * Int(2) == Int(4));
	assert(Int(4) * Int(-3) == Int(-12));
	assert(Int(10023) * Int(1212) == Int(12147876));
}

void should_div_ints() {
	assert(Int(4) / Int(2) == Int(2));
	assert(Int(-12) / Int(3) == Int(-4));
	assert(Int(4) / Int(-2) == Int(-2));
	assert(Int(10) / Int(5) == Int(2));
	assert(Int(412343212341232342) / Int(2) == Int(412343212341232342/2));
	assert(Int(123434515134212) / Int(12343124) == Int(123434515134212/12343124));
}

void should_rem_ints() {
	assert(Int(5) % Int(2) == Int(1));
	assert(Int(-5) % Int(2) == Int(-1));
	assert(Int(5) % Int(-2) == Int(1));
	assert(Int(12342341234) % Int(2341) == Int(12342341234 % 2341));
	assert(Int(51234) % Int(2332) == Int(51234 % 2332));
	assert(Int(5123) % Int(3) == Int(5123 % 3));
}

void should_pow_ints() {
	assert(pow(Int(4), Int(2)) == Int(16));
	assert(pow(Int(2), Int(3)) == Int(8));
	assert(pow(Int(1), Int(5)) == Int(1));
	assert(pow(Int(10), Int(3)) == Int(1000));
	assert(pow(Int(4), 0.5) == 2.0);
	assert(pow(Int(-3), Int(2)) == Int(9));
	assert(pow(Int(-3), Int(3)) == Int(-27));
}

void should_gcd_ints() {
	assert(gcd(Int(4), Int(2)) == Int(2));
	assert(gcd(Int(9), Int(6)) == Int(3));
	assert(gcd(Int(-4), Int(2)) == Int(2));
	assert(gcd(Int(-4), Int(-2)) == Int(-2));
}

void should_lcm_ints() {
	assert(lcm(Int(4), Int(2)) == Int(4));
	assert(lcm(Int(4), Int(3)) == Int(12));
}


void should_copy_ints() {
	Int u;
  Int i = 2;
	Int j = 4;
	Int t = Int(i);
	Int k = Int(j);
	Int r = Int(Int(3));

	assert(u == 0);
	assert(i == 2);
	assert(j == 4);
	assert(t == 2);
	assert(k == 4);
	assert(r == 3);
}

void should_increment_and_decrement_ints() {
	Int i = 0;

	assert(i++ == 0);
	assert(++i == 2);
	assert(i == 2);
	assert(i-- == 2);
	assert(--i == 0);
	assert(i == 0);
}

void should_invert_ints() {
	Int i = 1;

	assert(-i = -1);
	assert(+i = 1);
}

int main() {
	TEST(should_add_ints)
	TEST(should_sub_ints)
	TEST(should_mul_ints)
	TEST(should_div_ints)
	TEST(should_rem_ints)
	TEST(should_pow_ints)
	TEST(should_gcd_ints)
	TEST(should_lcm_ints)
	TEST(should_copy_ints)
	TEST(should_increment_and_decrement_ints)
	TEST(should_invert_ints)
}
