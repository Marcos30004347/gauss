//#include "Core/AST/AST.hpp"
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
//   assert(ast1[0].kind() == Kind::Integer);
//   assert(ast1[0].value() == 3);
//   assert(ast1[1].kind() == Kind::Integer);
//   assert(ast1[1].value() == 4);
//   assert(ast1[2].kind() == Kind::Integer);
//   assert(ast1[2].value() == 5);

//   assert(ast2.kind() == Kind::Addition);
//   assert(ast2[0].kind() == Kind::Integer);
//   assert(ast2[0].value() == 3);
//   assert(ast2[1].kind() == Kind::Integer);
//   assert(ast2[1].value() == 4);
//   assert(ast2[2].kind() == Kind::Fraction);
//   assert(ast2[2][0].kind() == Kind::Integer);
//   assert(ast2[2][0].value() == 5);
//   assert(ast2[2][1].kind() == Kind::Integer);
//   assert(ast2[2][1].value() == 6);

//   assert(ast3.kind() == Kind::Integer);
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
//    emplace_header(Kind::Integer, /*sizeof int storage*/ 1);
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
//     if (kind() == Kind::Integer) {
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

using namespace ast_teste;

void should_construct_ast() {
  ast *ast0 = ast_create(ast::add);

  assert(ast_is_kind(ast0, ast::add));

  ast_delete(ast0);

  ast *ast1 =
      ast_create(ast::add, {ast_integer(3), ast_integer(4), ast_integer(5)});

  assert(ast_is_kind(ast1, ast::add));

  assert(ast_is_kind(ast_operand(ast1, 0), ast::integer));
  assert(ast_is_kind(ast_operand(ast1, 1), ast::integer));
  assert(ast_is_kind(ast_operand(ast1, 2), ast::integer));

  assert(ast_value(ast_operand(ast1, 0)) == 3);
  assert(ast_value(ast_operand(ast1, 1)) == 4);
  assert(ast_value(ast_operand(ast1, 2)) == 5);

  ast_delete(ast1);

  ast *ast2 = ast_symbol("x");

  assert(ast_is_kind(ast2, ast::symbol));

	assert(strcmp(ast_id(ast2), "x") == 0);

  ast_delete(ast2);
}

void should_ast_insert_and_ast_remove_from_ast() {
  ast *a = ast_create(ast::add);

  ast_insert(a, ast_integer(0), 0);
  ast_insert(a, ast_integer(1), 1);
  ast_insert(a, ast_integer(2), 2);
  ast_insert(a, ast_integer(3), 3);
  ast_insert(a, ast_integer(4), 4);
  ast_insert(a, ast_integer(5), 5);
  ast_insert(a, ast_integer(6), 6);
  ast_insert(a, ast_integer(7), 7);
  ast_insert(a, ast_integer(8), 8);
  ast_insert(a, ast_integer(9), 9);

  assert(ast_size(a) == 10);

  assert(ast_kind(ast_operand(a, 0)) == ast::integer);
  assert(ast_kind(ast_operand(a, 1)) == ast::integer);
  assert(ast_kind(ast_operand(a, 2)) == ast::integer);
  assert(ast_kind(ast_operand(a, 3)) == ast::integer);
  assert(ast_kind(ast_operand(a, 4)) == ast::integer);
  assert(ast_kind(ast_operand(a, 5)) == ast::integer);
  assert(ast_kind(ast_operand(a, 6)) == ast::integer);
  assert(ast_kind(ast_operand(a, 7)) == ast::integer);
  assert(ast_kind(ast_operand(a, 8)) == ast::integer);
  assert(ast_kind(ast_operand(a, 9)) == ast::integer);

  assert(ast_value(ast_operand(a, 0)) == 0);
  assert(ast_value(ast_operand(a, 1)) == 1);
  assert(ast_value(ast_operand(a, 2)) == 2);
  assert(ast_value(ast_operand(a, 3)) == 3);
  assert(ast_value(ast_operand(a, 4)) == 4);
  assert(ast_value(ast_operand(a, 5)) == 5);
  assert(ast_value(ast_operand(a, 6)) == 6);
  assert(ast_value(ast_operand(a, 7)) == 7);
  assert(ast_value(ast_operand(a, 8)) == 8);
  assert(ast_value(ast_operand(a, 9)) == 9);

  ast_insert(a, ast_integer(10), 2);

  assert(ast_kind(ast_operand(a, 1)) == ast::integer);
  assert(ast_kind(ast_operand(a, 2)) == ast::integer);
  assert(ast_kind(ast_operand(a, 3)) == ast::integer);

  assert(ast_value(ast_operand(a, 1)) == 1);
  assert(ast_value(ast_operand(a, 2)) == 10);
  assert(ast_value(ast_operand(a, 3)) == 2);

  assert(ast_size(a) == 11);

  ast_remove(a, 2);

  assert(ast_size(a) == 10);

  assert(ast_value(ast_operand(a, 1)) == 1);
  assert(ast_value(ast_operand(a, 2)) == 2);
  assert(ast_value(ast_operand(a, 3)) == 3);

  ast_delete(a);
}

void should_sort_asts() {
  ast *a0 =
      ast_create(ast::add, {ast_integer(1), ast_integer(3), ast_integer(2)});

	//printf("%s\n", ast_to_string(a0).c_str());
  ast_sort(a0);
	//printf("%s\n", ast_to_string(a0).c_str());


  // assert(ast_is_kind(ast_operand(a0, 0), ast::integer));
  // assert(ast_is_kind(ast_operand(a0, 1), ast::integer));
  // assert(ast_is_kind(ast_operand(a0, 2), ast::integer));

  // assert(ast_value(ast_operand(a0, 0)) == 1);
  // assert(ast_value(ast_operand(a0, 1)) == 2);
  // assert(ast_value(ast_operand(a0, 2)) == 3);

  ast_delete(a0);

  ast *a1 = ast_create(
      ast::mul, {
                    ast_symbol("x"),
                    ast_create(ast::pow, {ast_symbol("y"), ast_integer(2)}),
                    ast_integer(4),

                    ast_create(ast::pow, {ast_symbol("x"), ast_integer(2)}),
                    ast_symbol("y"),
                });
	//printf("%s\n", ast_to_string(a1).c_str());
  ast_sort(a1);
	//printf("%s\n", ast_to_string(a1).c_str());

	// assert(ast_is_kind(ast_operand(a1, 0), ast::integer));
  // assert(ast_value(ast_operand(a1, 0)) == 4);

	// assert(ast_is_kind(ast_operand(a1, 1), ast::symbol));
  // assert(ast_id(ast_operand(a1, 1)) == "x");

	// assert(ast_is_kind(ast_operand(a1, 2), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(a1, 2), 0), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a1, 2), 0)) == "x");
  // assert(ast_is_kind(ast_operand(ast_operand(a1, 2), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a1, 2), 1)) == 2);

	// assert(ast_is_kind(ast_operand(a1, 3), ast::symbol));
  // assert(ast_id(ast_operand(a1, 3)) == "y");
  // assert(ast_is_kind(ast_operand(a1, 4), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(a1, 4), 0), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a1, 4), 0)) == "y");
  // assert(ast_is_kind(ast_operand(ast_operand(a1, 4), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a1, 4), 1)) == 2);

  ast_delete(a1);

	TIMED_SECTION_START("AST2 alloc");
  ast *a2 = ast_create(
      ast::add,
      {
          ast_create(
              ast::mul,
              {
                  ast_symbol("y"),
                  ast_symbol("z"),
                  ast_symbol("x"),
                  ast_create(ast::pow, {ast_symbol("y"), ast_integer(2)}),
                  ast_integer(4),
                  ast_create(ast::pow, {ast_symbol("x"), ast_integer(2)}),
                  ast_symbol("z"),
              }),

          ast_create(ast::add,
                      {ast_integer(1), ast_integer(1), ast_integer(1)}),
          ast_create(
              ast::mul,
              {
                  ast_symbol("x"),
                  ast_create(ast::pow, {ast_symbol("y"), ast_integer(2)}),
                  ast_integer(4),

                  ast_create(ast::pow, {ast_symbol("x"), ast_integer(2)}),
                  ast_symbol("y"),
              }),

          ast_create(
              ast::add,
              {
                  ast_symbol("x"),
                  ast_integer(4),
                  ast_create(ast::pow, {ast_symbol("z"), ast_integer(3)}),
              }),
					ast_symbol("x"),
					ast_symbol("z"),
					ast_symbol("y"),
					ast_create(ast::mul , {ast_symbol("x"), ast_integer(2)}),
					ast_create(ast::mul , {ast_symbol("z"), ast_integer(2)}),
          ast_create(ast::add,
                     {ast_integer(2), ast_integer(4), ast_integer(1)}),
          ast_create(ast::add,
                     {ast_integer(2), ast_symbol("z"), ast_integer(1)}),

          ast_create(
              ast::mul,
              {
                  ast_symbol("w"),
                  ast_create(ast::pow, {ast_symbol("x"), ast_integer(2)}),
                  ast_integer(4),

                  ast_create(ast::pow, {ast_symbol("z"), ast_integer(2)}),
                  ast_symbol("y"),
              }),
					ast_integer(4),
          ast_create(ast::add,
                     {ast_integer(2), ast_symbol("z"),
                      ast_create(ast::pow, {ast_symbol("x"), ast_integer(3)})}),

          ast_create(ast::mul,
                     {
                         ast_integer(4),
                         ast_integer(1),
                         ast_integer(7),

                     }),

          ast_create(ast::add,
                     {ast_integer(2), ast_symbol("x"),
                      ast_create(ast::pow, {ast_symbol("x"), ast_integer(2)})}),
      });

  TIMED_SECTION_START("AST2");
  ast_sort(a2);
  TIMED_SECTION_STOP("AST2");
  printf("%s\n", ast_to_string(a2).c_str());

  TIMED_SECTION_STOP("AST2 alloc");
  // printf("%s\n", ast_to_string(a2).c_str());

  return;

  // assert(ast_is_kind(ast_operand(a2, 0), ast::mul));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 0), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 0), 0)) == 1);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 0), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 0), 1)) == 4);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 0), 2), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 0), 2)) == 7);

  // assert(ast_is_kind(ast_operand(a2, 1), ast::mul));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 1), 0)) == 4);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 1), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 1), 1)) == "x");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 2), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 1), 2), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 1),
  // 2), 0)) == "x"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 1), 2), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 1), 2), 1)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 3), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 1), 3)) == "y");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 4), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 1), 4), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 1),
  // 4), 0)) == "y"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 1), 4), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 1), 4), 1)) == 2);

  // assert(ast_is_kind(ast_operand(a2, 1), ast::mul));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 1), 0)) == 4);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 1), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 1), 1)) == "x");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 2), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 1), 2), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 1),
  // 2), 0)) == "x"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 1), 2), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 1), 2), 1)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 3), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 1), 3)) == "y");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 1), 4), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 1), 4), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 1),
  // 4), 0)) == "y"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 1), 4), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 1), 4), 1)) == 2);

  // assert(ast_is_kind(ast_operand(a2, 2), ast::mul));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 2), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 2), 0)) == 4);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 2), 1), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 2), 1), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 2),
  // 1), 0)) == "x"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 2), 1), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 2), 1), 1)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 2), 2), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 2), 2), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 2),
  // 2), 0)) == "y"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 2), 2), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 2), 2), 1)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 2), 3), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 2), 3)) == "z");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 2), 4), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 2), 4)) == "z");

  // assert(ast_is_kind(ast_operand(a2, 3), ast::mul));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 3), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 3), 0)) == 4);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 3), 1), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 3), 1)) == "w");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 3), 2), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 3), 2), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 3),
  // 2), 0)) == "x"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 3), 2), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 3), 2), 1)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 3), 3), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 3), 3)) == "y");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 3), 4), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 3), 4), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 3),
  // 4), 0)) == "z"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 3), 4), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 3), 4), 1)) == 2);

  // assert(ast_is_kind(ast_operand(a2, 4), ast::add));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 4), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 4), 0)) == 1);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 4), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 4), 1)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 4), 2), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 4), 2)) == 4);

  // assert(ast_is_kind(ast_operand(a2, 5), ast::add));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 5), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 5), 0)) == 1);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 5), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 5), 1)) == 3);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 5), 2), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 5), 2)) == 4);

  // assert(ast_is_kind(ast_operand(a2, 6), ast::add));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 6), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 6), 0)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 6), 1), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 6), 1)) == "x");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 6), 2), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 6), 2), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 6),
  // 2), 0)) == "x"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 6), 2), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 6), 2), 1)) == 2);

  // assert(ast_is_kind(ast_operand(a2, 7), ast::add));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 7), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 7), 0)) == 1);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 7), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 7), 1)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 7), 2), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 7), 2)) == "z");

  // assert(ast_is_kind(ast_operand(a2, 8), ast::add));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 8), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 8), 0)) == 2);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 8), 1), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 8), 1), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 8),
  // 1), 0)) == "x"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 8), 1), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 8), 1), 1)) == 3);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 8), 2), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 8), 2)) == "z");

  // assert(ast_is_kind(ast_operand(a2, 9), ast::add));
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 9), 0), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(a2, 9), 0)) == 4);
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 9), 1), ast::symbol));
  // assert(ast_id(ast_operand(ast_operand(a2, 9), 1)) == "x");
  // assert(ast_is_kind(ast_operand(ast_operand(a2, 9), 2), ast::pow));
  // assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2, 9), 2), 0),
  // ast::symbol)); assert(ast_id(ast_operand(ast_operand(ast_operand(a2, 9),
  // 2), 0)) == "z"); assert(ast_is_kind(ast_operand(ast_operand(ast_operand(a2,
  // 9), 2), 1), ast::integer));
  // assert(ast_value(ast_operand(ast_operand(ast_operand(a2, 9), 2), 1)) == 3);




  ast *a3 = ast_create(
      ast::add,
      {ast_create(ast::pow, {ast_symbol("x"), ast_integer(2)}),
       ast_create(ast::pow, {ast_integer(2), ast_integer(3)}),
       ast_create(ast::pow, {ast_symbol("z"), ast_integer(3)}), ast_symbol("x"),
       ast_symbol("z"), ast_create(ast::fact, {ast_symbol("x")}),
       ast_create(ast::mul,
                  {
                      ast_create(ast::pow, {ast_symbol("y"), ast_integer(3)}),
                      ast_create(ast::pow, {ast_symbol("x"), ast_integer(3)}),
                  })});

  printf("%s\n", ast_to_string(a3).c_str());
  ast_sort(a3);
  printf("%s\n", ast_to_string(a3).c_str());

  ast_delete(a3);
}



void should_eval_consts_ast() {
	ast* u0 = ast_create(ast::add, { ast_symbol("x") });
	ast* v0 = ast_create(ast::add, { ast_symbol("x") });

	ast* r0 = eval_add_mulsym(u0, 0, v0, 0, true);

	assert(ast_operand(v0, 0) == nullptr);

	printf("%s\n", ast_to_string(r0).c_str());

	ast* u1 = ast_create(ast::add, { ast_create(ast::mul, { ast_integer(2), ast_symbol("x") }) });
	ast* v1 = ast_create(ast::add, { ast_create(ast::mul, { ast_integer(3), ast_symbol("x") }) });

	ast* r1 = eval_add_mulsym(u1, 0, v1, 0, true);

	printf("%s\n", ast_to_string(r1).c_str());

	ast* u2 = ast_create(ast::add, { ast_create(ast::mul, { ast_symbol("x"), ast_symbol("y")}) });
	ast* v2 = ast_create(ast::add, { ast_create(ast::mul, { ast_integer(3), ast_symbol("x"), ast_symbol("y") }) });

	ast* r2 = eval_add_mulsym(u2, 0, v2, 0, true);

	printf("%s\n", ast_to_string(r2).c_str());

	ast* u3 = ast_create(ast::add, { ast_create(ast::mul, { ast_integer(3), ast_symbol("x") }) });
	ast* v3 = ast_create(ast::add, { ast_symbol("x") });

	ast* r3 = eval_add_mulsym(u3, 0, v3, 0, true);

	printf("%s\n", ast_to_string(r3).c_str());
}

int main() {
  TEST(should_construct_ast)
  TEST(should_ast_insert_and_ast_remove_from_ast)
  TEST(should_sort_asts)

	TEST(should_eval_consts_ast)

  return 0;
}
