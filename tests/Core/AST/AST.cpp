#include "Core/AST/AST.hpp"
#include "test.hpp"
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

#include <string.h>
#include <stdlib.h>

class buffer {
	buffer(uint64_t* b, uint64_t* r) {
		buff = b;
		ref_count = r;
	}

public:
	uint64_t* buff;
	uint64_t* ref_count;

	buffer(const uint64_t size) {
		ref_count = new uint64_t(1);
		buff = new uint64_t[size];
	}

	buffer(buffer& b) {
		ref_count = b.ref_count;
		buff = b.buff;

		*ref_count = *ref_count + 1;
	}

	buffer(buffer&& b) {
		buff = b.buff;
		ref_count = b.ref_count;
	}

	~buffer() {
		*ref_count = *ref_count - 1;
		//std::cout <<"ref count: " << *ref_count << "\n";

		if(*ref_count == 0) {
			delete ref_count;
			delete buff;
		}
	}

	buffer view(uint64_t idx) {
		*ref_count = *ref_count + 2;
		return buffer(buff + idx, ref_count);
	}

	inline uint64_t& operator[](uint64_t i) {
		return buff[i];
	}

	inline buffer operator+(uint64_t i) {
		*ref_count = *ref_count + 1;
		return buffer(buff + i, ref_count);
	}

	inline buffer operator-(uint64_t i) {
		*ref_count = *ref_count + 1;
		return buffer(buff - i, ref_count);
	}

	inline buffer& operator=(buffer& a) {
		if(*ref_count == 1) {
			delete buff;
			delete ref_count;
		}

		buff = a.buff;
		ref_count = a.ref_count;

		*ref_count = *ref_count + 1;

		return *this;
	}

	inline buffer& operator=(buffer&& a) {
		if(*ref_count == 1) {
			delete buff;
			delete ref_count;
		}

		buff = a.buff;
		ref_count = a.ref_count;
		return *this;
	}

};

class expression {
public:
  buffer buff;

  const static uint64_t kind_idx = 0;
  const static uint64_t size_idx = 1;
  const static uint64_t desc_idx = 2;

  const static size_t head_size = 3;
  const static size_t desc_marg = 1;

  inline Kind kind() { return (Kind)buff[kind_idx]; }

  inline uint64_t size() { return buff[size_idx]; }

  inline uint64_t desc() { return buff[desc_idx]; }

  inline uint64_t memory_size() {
		return buff[desc() + size()];
	}

  inline void emplace_header(Kind k, uint64_t s) {
    buff[kind_idx] = k;
    buff[size_idx] = s;
    buff[desc_idx] = head_size + s;

		// first elements starts after the head
		buff[buff[desc_idx]] = head_size;
    // last element if desc holds the size of the buffer
    buff[buff[desc_idx] + buff[size_idx]] = buff[desc_idx] + buff[size_idx] + 1;
  }

  expression(Kind k): buff(head_size + desc_marg) {
    emplace_header(k, 0);
  }

  expression(int v): buff(head_size + desc_marg + 2) {
    emplace_header(Kind::Integer, /*sizeof int storage*/ 1);
		buff[head_size] = v;
  }

	expression(buffer&& b): buff(0) {
		buff = std::move(b);
	}

	inline void insert(expression a) {
    uint64_t t_msize =   memory_size();
		uint64_t a_msize = a.memory_size();

    uint64_t res_size = t_msize + a_msize + 1;

    // allocate a new buffer that can hold current contents plus a
    buffer b = buffer(res_size);

		uint64_t* ptr;

		ptr = b.buff;

		// copy previous content
		memcpy(ptr, buff.buff, sizeof(uint64_t) * buff[desc_idx]);

		ptr = ptr + buff[desc_idx];

		// copy a content of a
    memcpy(ptr, a.buff.buff, sizeof(uint64_t) * a_msize);

		ptr = ptr + a_msize;

		// copy old desc content
		memcpy(ptr, buff.buff + buff[desc_idx], sizeof(uint64_t) * size());

		b[kind_idx] = kind();
    b[size_idx] = size() + 1;
 		b[desc_idx] = buff[desc_idx] + a_msize;

		b[b[desc_idx] + buff[size_idx]] = buff[desc_idx];
		b[b[desc_idx] + b[size_idx]] = b[desc_idx] + b[size_idx] + 1;

		buff = b;
  }

	expression operator[](uint64_t idx) {
		return expression(buff.view(buff[desc() + idx]));
	}

	void printRec() {
		if(kind() == Kind::Integer) {
			std::cout << buff[head_size];
			return;
		}

		printf("+");
		for(uint64_t i = 0; i < size(); i++) {
			operator[](i).print() ;
			std::cout << " ";
		}
	}
	void print() {
		printRec();
		std::cout << "\n";
	}

	void printBuffer() {
		for(size_t i = 0; i < memory_size(); i++) {
			std::cout << buff[i] << " ";
		}
		std::cout << "\n";
	}
};


int main() {
  TEST(should_create_ast_nodes)
  TEST(should_match_ast_nodes)
  TEST(should_deep_copy_ast_nodes)

	expression a(Kind::Addition);

	a.printBuffer();
	a.insert(expression(7));
	a.printBuffer();
	a.insert(expression(2));
	a.printBuffer();
	//a.print();
	// a.print();
	// a.insert(expression(3));
	// a.print();

	return 0;
}
