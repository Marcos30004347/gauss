#include "Core/AST/AST.hpp"
#include "test.hpp"
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>

using namespace ast;

void should_create_ast_nodes() {
  Expr ast0 = Expr(Kind::Undefined);
  Expr ast1 = Expr(3) + Expr(4) + Expr(5);
  Expr ast2 = Expr(3) + Expr(4) + Expr(5) / Expr(6);
  Expr ast3 = 0;

  assert(ast0.kind() == Kind::Undefined);
  assert(ast1.kind() == Kind::Addition);
  assert(ast1[0].kind() == Kind::Integer);
  assert(ast1[0].value() == 3);
  assert(ast1[1].kind() == Kind::Integer);
  assert(ast1[1].value() == 4);
  assert(ast1[2].kind() == Kind::Integer);
  assert(ast1[2].value() == 5);

  assert(ast2.kind() == Kind::Addition);
  assert(ast2[0].kind() == Kind::Integer);
  assert(ast2[0].value() == 3);
  assert(ast2[1].kind() == Kind::Integer);
  assert(ast2[1].value() == 4);
  assert(ast2[2].kind() == Kind::Fraction);
  assert(ast2[2][0].kind() == Kind::Integer);
  assert(ast2[2][0].value() == 5);
  assert(ast2[2][1].kind() == Kind::Integer);
  assert(ast2[2][1].value() == 6);

  assert(ast3.kind() == Kind::Integer);
  assert(ast3.value() == 0);
  assert(ast3 == 0);
}

void should_match_ast_nodes() {
  Expr ast0 = Expr(3) + Expr(4) + Expr(5);
  Expr ast1 = Expr(3) + Expr(4) + Expr(5);
  Expr ast2 = Expr(3) + Expr(4);

  assert(ast0 == ast1);
  assert(ast0 != ast2);
}

void should_deep_copy_ast_nodes() {
  Expr ast0 = Expr(3) + Expr(4) + Expr(5);
  Expr ast1 = ast0;

  assert(ast0 == ast1);
}

#include <stdlib.h>
#include <string.h>

class buffer {
  buffer(uint64_t *b, int64_t *r, uint64_t *p) {
    buff = b;
    refs = r;
    prnt = p;
    incref();
  }

  uint64_t *prnt;



public:
  uint64_t *buff;
  int64_t *refs;

  buffer(const uint64_t s) {
    refs = new int64_t(1);
    prnt = buff = new uint64_t[s];
  }

  buffer() {
    refs = nullptr;
    buff = nullptr;
    prnt = nullptr;
  }

	buffer(buffer &b) {
		refs = b.refs;
    buff = b.buff;
    prnt = b.prnt;

		incref();
  }

	buffer(buffer &&b) {
    buff = b.buff;
    refs = b.refs;
    prnt = b.prnt;
  }

  ~buffer() {
    if (refs == nullptr || *refs == 0)
      return;

    decref();

    if (*refs == 0) {
      delete refs;
      delete prnt;
    }
  }

  void view(uint64_t idx, buffer &b) {
    b.prnt = prnt;
    b.refs = refs;
    b.buff = buff + idx;

    incref();
  }

  inline int64_t refcnt() { return *refs; }
  inline int64_t incref() { return *refs += 1; }
  inline int64_t decref() {
    if (*refs == 1) {
      delete prnt;
      delete refs;
      return 0;
    }
    return *refs -= 1;
  }

  inline uint64_t &operator[](uint64_t i) { return buff[i]; }

  inline buffer operator+(uint64_t i) { return buffer(buff + i, refs, prnt); }

  inline buffer operator-(uint64_t i) { return buffer(buff - i, refs, prnt); }

  inline buffer &operator=(buffer &a) {
    if (*refs == 1) {
      delete prnt;
      delete refs;
    }

    buff = a.buff;
    refs = a.refs;
    prnt = a.prnt;

    incref();

    return *this;
  }

  inline buffer &operator=(buffer &&a) {

		if (*refs == 1) {
      delete prnt;
      delete refs;
    }

    decref();

    buff = a.buff;
    refs = a.refs;
    prnt = a.prnt;

		incref();

    return *this;
  }

  inline uint64_t *operator&() const { return buff; }
};

class expr {
	// TODO: move desc to start after the heading, this will allow us to use one less uint64 in the buffer.
	// TODO: adding new elements will append its content to the end of the buffer, but desc will contain right index.
	// TODO: fix inserting, the size of the element being inserted needs to be added to all successive desc element.

  expr() {}

public:
  buffer buff;

  const static uint64_t kind_idx = 0;
  const static uint64_t size_idx = 1;
  const static uint64_t desc_idx = 2;

  const static size_t head_size = 2;
  const static size_t desc_marg = 1;

  inline Kind kind() { return (Kind)buff[kind_idx]; }

  inline uint64_t size() { return buff[size_idx]; }
  inline uint64_t desc() { return 2; }
  inline uint64_t memory_size() { return buff[desc_idx + buff[size_idx]]; }

  inline void emplace_header(Kind k, uint64_t s) {
    buff[kind_idx] = k;
    buff[size_idx] = s;
		buff[desc_idx + s] = buff[size_idx] + buff[desc_idx] + 1;
	}

  expr(Kind k) : buff(head_size + desc_marg) { emplace_header(k, 0); }

  expr(int v) : buff(head_size + desc_marg + 2) {
    emplace_header(Kind::Integer, 1);

		buff[desc_idx + 0] = head_size + buff[size_idx] + desc_marg;
		buff[desc_idx + 1] = buff[desc_idx + buff[size_idx]];

		buff[buff[desc_idx]] = v;
	}

  expr(buffer &&b) : buff(0) { buff = std::move(b); }

	inline uint64_t value() {
		//TODO: assert kind is Integer is worth ?
		return buff[buff[desc_idx]];
	}

  inline void insert(expr a) {
    const uint64_t t_msize = memory_size();
    const uint64_t a_msize = a.memory_size();

    buffer bold = buff;

    buff = buffer(t_msize + a_msize + 1);

    uint64_t *bptr = &buff;

		uint64_t data_marg = head_size + bold[size_idx] + 1;

		// copy old header and descriptor
		memcpy(bptr, &bold, sizeof(uint64_t) * data_marg);

		// set new header and desc
		buff[desc_idx + buff[size_idx]] = t_msize + 1;
		buff[size_idx] += 1;
		buff[desc_idx + buff[size_idx]] += a_msize;

		// copy rest of data
		uint64_t data_size = t_msize - data_marg;

		bptr += head_size + buff[size_idx] + 1;
		memcpy(bptr, &bold + data_marg, sizeof(uint64_t) * data_size);

		// copy a's content
		bptr += data_size;
		memcpy(bptr, &a.buff, sizeof(uint64_t) * a_msize);

		return;

		// const uint64_t t_msize = memory_size();
    // const uint64_t a_msize = a.memory_size();

    // buffer bold = buff;

    // buff = buffer(t_msize + a_msize + 1);

    // uint64_t *bptr = &buff;

		// memcpy(bptr, &bold, sizeof(uint64_t) * bold[desc_idx]);

    // bptr = &buff + bold[desc_idx];

    // memcpy(bptr, &a.buff, sizeof(uint64_t) * a_msize);

    // bptr = bptr + a_msize;

    // memcpy(bptr, &bold + bold[desc_idx], sizeof(uint64_t) * size());

    // buff[size_idx] = size() + 1;
    // buff[desc_idx] = bold[desc_idx] + a_msize;

    // buff[buff[desc_idx] + bold[size_idx]] = bold[desc_idx];
    // buff[buff[desc_idx] + buff[size_idx]] = buff[desc_idx] + buff[size_idx] + 1;
  }

  inline void insert(expr a, uint64_t idx) {
		//TODO: fix
		//TODO: desc needs to be increased by the size of a
		//TODO: suprisingly this is faster that the first method, discover why

		assert(idx < size());

    const uint64_t t_msize = memory_size();
    const uint64_t a_msize = a.memory_size();

    buffer bold = buff;

    buff = buffer(t_msize + a_msize + 1);

    uint64_t *bptr = &buff;

		// copy old header and descriptor
		memcpy(bptr, &bold, sizeof(uint64_t) * (head_size + bold[size_idx] + 1));

		// set new header and desc
		buff[desc_idx + buff[size_idx]] = t_msize + 1;
		buff[size_idx] += 1;
		buff[desc_idx + buff[size_idx]] += a_msize;

		// copy rest of data
		uint64_t data_marg = head_size + bold[size_idx] + 1;
		uint64_t data_size = t_msize - data_marg;

		bptr += head_size + buff[size_idx] + 1;
		memcpy(bptr, &bold + data_marg, sizeof(uint64_t) * data_size);

		// copy a's content
		bptr += data_size;
		memcpy(bptr, &a.buff, sizeof(uint64_t) * a_msize);

		return;
                /*
                // copy content up to index
                memcpy(bptr, &bold, sizeof(uint64_t) * bold[desc_idx + idx]);

    bptr = &buff + bold[desc_idx + idx];

                // copy a into its idx place
    memcpy(bptr, &a.buff, sizeof(uint64_t) * a_msize);

    bptr = bptr + a_msize;

                // copy rest of data
                memcpy(bptr, &bold + bold[desc_idx + idx], sizeof(uint64_t) *
(desc_idx - bold[desc_idx + idx])); bptr = bptr + desc_idx - bold[desc_idx +
idx];

                // copy desc
    memcpy(bptr, &bold + desc_idx, sizeof(uint64_t) * idx);

                bptr = bptr + idx + 1;

                memcpy(bptr, &bold + desc_idx + idx, sizeof(uint64_t) * (size()
- idx));

    buff[size_idx] = size() + 1;
    buff[desc_idx] = desc_idx + a_msize;

    buff[buff[desc_idx] + bold[size_idx]] = desc_idx;
    buff[buff[desc_idx] + buff[size_idx]] = buff[desc_idx] + buff[size_idx] + 1;

    buff[buff[desc_idx] + idx] = bold[desc_idx + idx];
		 */
		}

	inline uint64_t size_of_operand(uint64_t idx) {
		return buff[desc() + idx + 1] - buff[desc() + idx];
	}

	inline void remove(uint64_t idx) {
		assert(idx < size());

		uint64_t tail = buff[desc() + size() - 1] - buff[desc() + idx];
		uint64_t rsize = size_of_operand(idx);

		buffer bold = buff;
		printf("A\n");
		buff = buffer(memory_size() - rsize - 1);

		uint64_t* bptr = &buff;

		printf("B\n");
		memcpy(bptr, &bold, sizeof(uint64_t) * bold[desc() + idx]);

		bptr = bptr + bold[desc() + idx];

		printf("C\n");
		memcpy(bptr, &bold + bold[desc() + idx + 1], sizeof(uint64_t) * tail);

		bptr = bptr + tail;

		printf("D\n");
		memcpy(bptr, &bold + desc(), sizeof(uint64_t)*idx);

		bptr = bptr + idx;

		printf("E\n");
		memcpy(bptr, &bold + desc() + idx + 1, sizeof(uint64_t) * size() - idx - 1);
		printf("F\n");

		buff[size_idx] -= 1;
		buff[desc_idx] -= rsize;



		for(int i =0; i < 20; i++) {
			std::cout << buff[i] << " ";
		}
		std::cout << std::endl;
	}

  expr operator[](uint64_t idx) {
    expr p;
    buff.view(buff[desc() + idx], p.buff);
    return p;
  }

  void printRec() {
    if (kind() == Kind::Integer) {
      std::cout << buff[head_size];
      return;
    }

    for (uint64_t i = 0; i < size(); i++) {
      operator[](i).printRec();
      if (i < size() - 1)
        std::cout << " + ";
    }
  }

  void print() {
    printRec();
    std::cout << "\n";
  }

  void printBuffer() {
		uint64_t s = memory_size();
    for (size_t i = 0; i < s; i++) {
      std::cout << buff[i] << " ";
    }
    std::cout << "\n";
  }
};

void test_buffer() {
	buffer c;
	{
		buffer a(10);

		assert(a.refcnt() == 1);

		{
			buffer b = a;
			assert(a.refcnt() == 2);
		}

		assert(a.refcnt() == 1);

		a[0] = 3;
		a[1] = 4;

		assert(a[0] == 3);
		assert(a[1] == 4);

		a.view(1, c);

		assert(a.refcnt() == 2);
		assert(a.refcnt() == c.refcnt());
	}

	assert(c.refcnt() == 1);
	assert(c[0] == 4);

	buffer d(10);
	assert(d.refcnt() == 1);
	buffer e = d;

	assert(d.refcnt() == 2);
	assert(e.refcnt() == 2);

	d = buffer(4);

	assert(d.refcnt() == 1);
	assert(e.refcnt() == 1);
}

int main() {
  TEST(should_create_ast_nodes)
  TEST(should_match_ast_nodes)
  TEST(should_deep_copy_ast_nodes)
	TEST(test_buffer)

  TIMED_SECTION_START("NEW AST")
  expr a(Kind::Addition);
  a.insert(expr(7));
  a.insert(expr(2));
  a.insert(expr(3));
  TIMED_SECTION_STOP("NEW AST")

  TIMED_SECTION_START("OLD AST")
  Expr a(Kind::Addition);
  a.insert(Expr(7));
  a.insert(Expr(2));
  a.insert(Expr(3));
  TIMED_SECTION_STOP("OLD AST")

	TIMED_SECTION_START("NEW AST")
  expr a(Kind::Addition);
  a.insert(expr(7));
  a.insert(expr(2), 0);
  a.insert(expr(3), 1);
  TIMED_SECTION_STOP("NEW AST")

	TIMED_SECTION_START("OLD AST")
  Expr a(Kind::Addition);
  a.insert(Expr(7));
  a.insert(Expr(2), 0);
  a.insert(Expr(3), 1);
  TIMED_SECTION_STOP("OLD AST")

	TIMED_SECTION_START("NEW AST")
  expr a(Kind::Addition);
  a.insert(expr(7));
  a.insert(expr(2), 0);
  a.insert(expr(3), 1);
	a.print();
	a.printBuffer();
	a.remove(0);
	a.printBuffer();
  TIMED_SECTION_STOP("NEW AST")


  return 0;
}
