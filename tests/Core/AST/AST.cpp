//#include "Core/AST/AST.hpp"
#include <cstdlib>
#define TEST_TIME_REPORT_UNIT TEST_TIME_REPORT_NS

#include "test.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>

// using namespace ast;

// void should_create_ast_nodes() {
//   Expr ast0 = Expr(Kind::Undefined);
//   Expr ast1 = Expr(3) + Expr(4) + Expr(5);
//   Expr ast2 = Expr(3) + Expr(4) + Expr(5) / Expr(6);
//   Expr ast3 = 0;

//   assert(ast0.kind() == Kind::Undefined);
//   assert(ast1.kind() == Kind::Addition);
//   assert(ast1[0].kind() == kind::INT);
//   assert(ast1[0].value() == 3);
//   assert(ast1[1].kind() == kind::INT);
//   assert(ast1[1].value() == 4);
//   assert(ast1[2].kind() == kind::INT);
//   assert(ast1[2].value() == 5);

//   assert(ast2.kind() == Kind::Addition);
//   assert(ast2[0].kind() == kind::INT);
//   assert(ast2[0].value() == 3);
//   assert(ast2[1].kind() == kind::INT);
//   assert(ast2[1].value() == 4);
//   assert(ast2[2].kind() == Kind::Fraction);
//   assert(ast2[2][0].kind() == kind::INT);
//   assert(ast2[2][0].value() == 5);
//   assert(ast2[2][1].kind() == kind::INT);
//   assert(ast2[2][1].value() == 6);

//   assert(ast3.kind() == kind::INT);
//   assert(ast3.value() == 0);
//   assert(ast3 == 0);
// }

// void should_match_ast_nodes() {
//   Expr ast0 = Expr(3) + Expr(4) + Expr(5);
//   Expr ast1 = Expr(3) + Expr(4) + Expr(5);
//   Expr ast2 = Expr(3) + Expr(4);

//   assert(ast0 == ast1);
//   assert(ast0 != ast2);
// }

// void should_deep_copy_ast_nodes() {
//   Expr ast0 = Expr(3) + Expr(4) + Expr(5);
//   Expr ast1 = ast0;

//   assert(ast0 == ast1);
// }

// #include <stdlib.h>
// #include <string.h>

// class buffer {
//   buffer(uint64_t *b, int64_t *r, uint64_t *p) {
//     buff = b;
//     refs = r;
//     prnt = p;
//     incref();
//   }

//   uint64_t *prnt;

// public:
//   uint64_t *buff;
//   int64_t *refs;

//   buffer(const uint64_t s) {
//     refs = new int64_t(1);
//     prnt = buff = new uint64_t[s];
//   }

//   buffer() {
//     refs = nullptr;
//     buff = nullptr;
//     prnt = nullptr;
//   }

//   buffer(buffer &b) {
//     refs = b.refs;
//     buff = b.buff;
//     prnt = b.prnt;

//     incref();
//   }

//   buffer(buffer &&b) {
//     buff = b.buff;
//     refs = b.refs;
//     prnt = b.prnt;
//   }

//   ~buffer() {
//     if (refs == nullptr || *refs == 0)
//       return;

//     decref();

//     if (*refs == 0) {
//       delete refs;
//       delete prnt;
//     }
//   }

//   void view(uint64_t idx, buffer &b) {
//     b.prnt = prnt;
//     b.refs = refs;
//     b.buff = buff + idx;

//     incref();
//   }

//   inline int64_t refcnt() { return *refs; }
//   inline int64_t incref() { return *refs += 1; }
//   inline int64_t decref() {
//     if (*refs == 1) {
//       delete prnt;
//       delete refs;
//       return 0;
//     }
//     return *refs -= 1;
//   }

//   inline uint64_t &operator[](uint64_t i) { return buff[i]; }

//   inline buffer operator+(uint64_t i) { return buffer(buff + i, refs, prnt);
//   }

//   inline buffer operator-(uint64_t i) { return buffer(buff - i, refs, prnt);
//   }

//   inline buffer &operator=(buffer &a) {
//     if (*refs == 1) {
//       delete prnt;
//       delete refs;
//     }

//     buff = a.buff;
//     refs = a.refs;
//     prnt = a.prnt;

//     incref();

//     return *this;
//   }

//   inline buffer &operator=(buffer &&a) {

//     if (*refs == 1) {
//       delete prnt;
//       delete refs;
//     }

//     decref();

//     buff = a.buff;
//     refs = a.refs;
//     prnt = a.prnt;

//     incref();

//     return *this;
//   }

//   inline uint64_t *operator&() const { return buff; }
// };

// class expr {
//   // TODO: move desc to start after the heading, this will allow us to use
//   one
//   // less uint64 in the buffer.
//   // TODO: adding new elements will append its content to the end of the
//   buffer,
//   // but desc will contain right index.
//   // TODO: fix inserting, the size of the element being inserted needs to be
//   // added to all successive desc element.

//   expr() {}

// public:
//   buffer buff;

//   const static uint64_t kind_idx = 0;
//   const static uint64_t size_idx = 1;
//   const static uint64_t desc_idx = 2;

//   const static size_t head_size = 3;
//   const static size_t desc_marg = 1;

//   inline Kind kind() { return (Kind)buff[kind_idx]; }

//   inline uint64_t size() { return buff[size_idx]; }
//   inline uint64_t desc() { return buff[desc_idx]; }
//   inline uint64_t memory_size() { return buff[buff[desc_idx] +
//   buff[size_idx]]; }

//   inline void emplace_header(Kind k, uint64_t s) {
//     buff[kind_idx] = k;
//     buff[size_idx] = s;
//     buff[desc_idx] = head_size + s;

// 		// first elements starts after the head
// 		buff[buff[desc_idx]] = head_size;
//     // last element if desc holds the size of the buffer
//     buff[buff[desc_idx] + buff[size_idx]] = buff[desc_idx] + buff[size_idx] +
//     1;
// 	}

//   expr(Kind k) : buff(head_size + desc_marg) { emplace_header(k, 0); }

//   expr(int v) : buff(head_size + desc_marg + 2) {
//    emplace_header(kind::INT, /*sizeof int storage*/ 1);
// 		buff[head_size] = v;  }

//   expr(buffer &&b) : buff(0) { buff = std::move(b); }

//   inline uint64_t value() { return buff[buff[desc()]]; }

//   inline void insert(expr a) {
//     uint64_t t_msize = memory_size();
//     uint64_t a_msize = a.memory_size();

//     uint64_t res_size = t_msize + a_msize + 1;

// 		buffer bold = buff;

//     buff = buffer(res_size);

//     uint64_t *bptr = &buff;

//     memcpy(bptr, &bold, sizeof(uint64_t) * bold[desc_idx]);

//     bptr = bptr + bold[desc_idx];

//     memcpy(bptr, &a.buff, sizeof(uint64_t) * a_msize);

//     bptr = bptr + a_msize;

//     memcpy(bptr, &bold + bold[desc_idx], sizeof(uint64_t) * size());

//     buff[size_idx] = size() + 1;
//     buff[desc_idx] = bold[desc_idx] + a_msize;
//     buff[buff[desc_idx] + bold[size_idx]] = bold[desc_idx];
//     buff[buff[desc_idx] + buff[size_idx]] = buff[desc_idx] + buff[size_idx] +
//     1;
//   }

// 	// append a's content to the end of the operand array, but access will
// 	// be indexed by idx
//   inline void insert(expr a, uint64_t idx) {
//     assert(idx < size());

//     const uint64_t t_msize = memory_size();
//     const uint64_t a_msize = a.memory_size();

//     buffer bold = buff;

//     buff = buffer(t_msize + a_msize + 1);

//     uint64_t *bptr = &buff;

//     memcpy(bptr, &bold, sizeof(uint64_t) * bold[desc_idx]);

//     bptr = bptr + bold[desc_idx];

//     memcpy(bptr, &a.buff, sizeof(uint64_t) * a_msize);

//     bptr = bptr + a_msize;

//     memcpy(bptr, &bold + bold[desc_idx], sizeof(uint64_t) * idx);

//     bptr = bptr + idx + 1;

//     memcpy(bptr, &bold + bold[desc_idx] + idx, sizeof(uint64_t) * (size() -
//     idx));

// 	  buff[size_idx] = size() + 1;
//     buff[desc_idx] = bold[desc_idx] + a_msize;
// 		buff[buff[desc_idx] + idx] = bold[desc_idx];
//     buff[buff[desc_idx] + buff[size_idx]] = buff[desc_idx] + buff[size_idx] +
//     1;
//   }

//   inline uint64_t size_of_operand(uint64_t idx) {
//     return buff[desc() + idx + 1] - buff[desc() + idx];
//   }

//   inline void remove(uint64_t idx) {
// 		assert(idx < size());

// 		buffer bold = buff;

// 		uint64_t* dlt = &buff + buff[buff[desc_idx] + idx];

// 		uint64_t t = memory_size();
// 		uint64_t s = dlt[dlt[desc_idx] + dlt[size_idx]];

// 		buff = buffer(t - s);

// 		uint64_t* bptr = &buff;

// 		memcpy(bptr, &bold, sizeof(uint64_t) * bold[bold[desc_idx] +
// idx]);

// 		bptr += bold[bold[desc_idx] + idx];

// 		memcpy(bptr, &bold + bold[bold[desc_idx] + idx] + s,
// sizeof(uint64_t) * (bold[desc_idx] - bold[bold[desc_idx] + idx]));

// 		bptr += bold[desc_idx] - bold[bold[desc_idx] + idx];

// 		memcpy(bptr, &bold + bold[desc_idx], sizeof(uint64_t) * idx);

// 		bptr += idx;

// 		memcpy(bptr, &bold + bold[desc_idx] + idx + 1, sizeof(uint64_t)
// * (bold[size_idx] - idx + 1));

// 		buff[size_idx] -= 1;
// 		buff[desc_idx] -= s;
// 		buff[buff[desc_idx] + buff[size_idx]] -= s;

// 		for(size_t i = 0; i < buff[size_idx]; i++) {
// 			if(buff[buff[desc_idx] + i] > bold[bold[desc_idx] +
// idx]) { 				buff[buff[desc_idx] + i] -= s;
// 			}
// 		}
// 	}

//   expr operator[](uint64_t idx) {
//     expr p;
//     buff.view(buff[desc() + idx], p.buff);
//     return p;
//   }

//   void printRec() {
//     if (kind() == kind::INT) {
//       std::cout << value();
//       return;
//     }

//     for (uint64_t i = 0; i < size(); i++) {
//       operator[](i).printRec();
//       if (i < size() - 1)
//         std::cout << " + ";
//     }
//   }

//   void print() {
//     printRec();
//     std::cout << "\n";
//   }

//   void printBuffer() {
//     uint64_t s = memory_size();
//     for (size_t i = 0; i < 30; i++) {
//       std::cout << buff[i] << " ";
//     }
//     std::cout << "\n";
//   }
// };

// void test_buffer() {
//   buffer c;
//   {
//     buffer a(10);

//     assert(a.refcnt() == 1);

//     {
//       buffer b = a;
//       assert(a.refcnt() == 2);
//     }

//     assert(a.refcnt() == 1);

//     a[0] = 3;
//     a[1] = 4;

//     assert(a[0] == 3);
//     assert(a[1] == 4);

//     a.view(1, c);

//     assert(a.refcnt() == 2);
//     assert(a.refcnt() == c.refcnt());
//   }

//   assert(c.refcnt() == 1);
//   assert(c[0] == 4);

//   buffer d(10);
//   assert(d.refcnt() == 1);
//   buffer e = d;

//   assert(d.refcnt() == 2);
//   assert(e.refcnt() == 2);

//   d = buffer(4);

//   assert(d.refcnt() == 1);
//   assert(e.refcnt() == 1);
// }

#include "Core/AST/AST3.hpp"
//#include "Core/AST/AST4.hpp"

using namespace ast_teste;

void should_construct_ast() {
  ast ast0 = create(kind::ADD);

  assert(is(&ast0, kind::ADD));

  ast ast1 = create(kind::ADD, {integer(3), integer(4), integer(5)});

  assert(is(&ast1, kind::ADD));

  assert(is(operand(&ast1, 0), kind::INT));
  assert(is(operand(&ast1, 1), kind::INT));
  assert(is(operand(&ast1, 2), kind::INT));

  assert(get_val(operand(&ast1, 0)) == 3);
  assert(get_val(operand(&ast1, 1)) == 4);
  assert(get_val(operand(&ast1, 2)) == 5);

  ast ast2 = symbol("x");

  assert(is(&ast2, kind::SYM));

  assert(strcmp(get_id(&ast2), "x") == 0);
}

void should_insert_and_remove_from_ast() {
  ast a = create(kind::ADD);

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

void should_eval_asts() {
  ast x = ast("x");
  ast y = ast("y");
  ast z = ast("z");

  ast _1 = ast(1);
  ast _4 = ast(4);
  ast _7 = ast(7);

  ast a2 = y * z * x * pow(y, 2) * 4 * pow(x, 2) * z + (_1 + _1 + _1) +
           x * pow(y, 2) * 4 * pow(x, 2) * y + (x + 4 + pow(z, 3)) +
           (2 + z + _1) + 4 + (_4 * _1 * _7);

  reduce(&a2);

  assert(a2 == 4 * pow(x, 3) * pow(y, 3) * pow(z, 2) +
                   4 * pow(x, 3) * pow(y, 3) + pow(z, 3) + x + z + 42);
}

void should_expand_ast() {
  ast x = ast("x");
  ast y = ast("y");
  ast z = ast("z");

  ast a = (x + 2) * (x + 3) * (x + 4);

  expand(&a);

  assert(a == pow(x, 3) + 9 * pow(x, 2) + 26 * x + 24);

  ast b = pow(x * sqrt(y + 1) + 1, 4);

  expand(&b);
  assert(b == pow(x, 4) * pow(y, 2) +
                  4 * pow(x, 3) * pow(y + 1, fraction(3, 2)) + 6 * pow(x, 2) +
                  4 * x * pow(y + 1, fraction(1, 2)) + 2 * pow(x, 4) * y +
                  6 * pow(x, 2) * y + pow(x, 4) + 1);

  ast c = pow(x + y + z, 3);

  expand(&c);

  assert(c == 3 * x * pow(y, 2) + 3 * x * pow(z, 2) + 3 * y * pow(z, 2) +
                  3 * pow(x, 2) * y + 3 * pow(x, 2) * z + 3 * pow(y, 2) * z +
                  6 * x * y * z + pow(x, 3) + pow(y, 3) + pow(z, 3));

  ast d = pow(x + 1, 2);

  expand(&d);
  assert(d == pow(x, 2) + 2 * x + 1);

  ast e = pow(pow(x + 2, 2) + 3, 2);

  expand(&e);

  assert(e == pow(x, 4) + 8 * pow(x, 3) + 30 * pow(x, 2) + 56 * x + 49);

  ast f = (-32 * pow(z, 3) + 32 * pow(z, 4) + 48 * pow(z, 5) + -24 * pow(z, 6) +
           -48 * pow(z, 7) + -36 * pow(z, 8) + -40 * pow(z, 9) +
           -8 * pow(z, 10) + -8 * pow(z, 11)) /
          (4 * pow(z, 2));

  expand(&f);

  assert(f == -2 * pow(z, 9) + -2 * pow(z, 8) + -10 * pow(z, 7) +
                  -9 * pow(z, 6) + -12 * pow(z, 5) + -6 * pow(z, 4) +
                  12 * pow(z, 3) + 8 * pow(z, 2) + -8 * z);
}

void should_eval_equality() {
  ast x = ast("x");
  ast a = 3 + 4 * x;
  ast b = 4 * x + 3;
  ast c = 4 * x + 4;

  ast d = 4 * pow(x, 3) + 4 * pow(x, 2) + 6 * x + 7;
  ast e = 4 * pow(x, 2) + 4 * pow(x, 3) + 7 + 6 * x;

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

	list u = remove(k, { ast(2) });

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
  TEST(should_construct_ast)
  TEST(should_eval_equality)
  TEST(should_insert_and_remove_from_ast)
  TEST(should_eval_asts)
  TEST(should_expand_ast)
	TEST(should_perform_list_operations);
	TEST(should_perform_set_operations)

	return 0;
}
