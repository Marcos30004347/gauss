
#include <cassert>
#include <cstdint>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <utility>
#include <stdio.h>

class buffer {
  buffer(uint32_t *b, int64_t *r, uint32_t *p) {
    buff = b;
    refs = r;
    prnt = p;

    incref();
  }

  uint32_t *prnt;

public:
  uint32_t *buff;
  int64_t *refs;

  buffer(const uint32_t s) {
    refs = new int64_t(1);
    prnt = buff = new uint32_t[s];
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

  inline uint32_t &operator[](uint64_t i) { return buff[i]; }

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

  inline uint32_t *operator&() const { return buff; }
};

enum Type {
	UNDEFINED = 0,
	INTEGER  = (1 << 0),
	ADDITION = (1 << 1),
};

class expr {
	static const uint32_t null_idx = -1;
	// idx where kind is stored in the data buffer
	static const uint32_t kind_idx = 0;
	// idx where count is stored in the data buffer
	static const uint32_t size_idx = 1;
	// idx where indx_buffer start is stored in the data buffer
  static const uint32_t indx_idx = 2;
	// idx where addr_buffer start is stored in the data buffer
  static const uint32_t addr_idx = 3;
	// idx where total size of the buffer is stored in the data buffer
  static const uint64_t memo_idx = 4;

	static const uint32_t head_size = 5;

  buffer data;

	inline expr(buffer &&b) { data = std::move(b); }
	inline expr(): data(1) { data[kind_idx] = Type::UNDEFINED; };

  public:

	inline expr(Type k): data(head_size) {
		data[kind_idx] = k;
		data[size_idx] = 0;
		data[indx_idx] = null_idx;
		data[addr_idx] = null_idx;
		data[memo_idx] = head_size;
	}

	inline expr(uint32_t i): data(head_size + 3) {
		// size of a uint32_t
		uint32_t s = 1;

		// kind and size
		data[kind_idx] = Type::INTEGER;
		data[size_idx] = 1;

		// indx_buff starts after the data buff
		data[indx_idx] = head_size + s;
		// addr_buff starts after the indx buffer by count units
		data[addr_idx] = data[indx_idx] + data[size_idx];

		// total memory usage is the size of the
		// header(kind, size, indx_idx, addr_idx, memory, data_idx) +
		// the size of the data + the 2*size
		data[memo_idx] = head_size + s + 2*data[size_idx];

		// first element of addr_buff points to where the data starts
		data[data[addr_idx]] = head_size;
		// first element of indx_buff points to the first element of addr_buff
		data[data[indx_idx]] = data[addr_idx];

		// set the integer content
		data[data[data[data[indx_idx]]]] = i;
	}

	inline uint32_t size() { return data[size_idx]; }
	inline uint32_t memsize() { return data[memo_idx]; }
	inline uint32_t addr(uint64_t idx) { return data[data[data[indx_idx] + idx]]; }
	inline uint32_t* indx_buffer() { return &data + data[indx_idx]; };
	inline uint32_t* addr_buffer() { return &data + data[addr_idx]; };
	inline uint32_t* data_buffer() { return &data + data[data[indx_idx]]; };

	inline void fill_with_op(uint64_t idx, expr& e) {
		data.view(data[addr(idx)] , e.data);
	}

	inline void insert(expr e, int64_t idx) {
		assert(idx <= size());

		buffer odata = data;

		data = buffer(memsize() + e.memsize() + 2);

		uint32_t* nbuff = &data;
		uint32_t* obuff = &odata;

		// copy header
		memcpy(nbuff, &odata, head_size);

		// update where idx_buff and addr_buff starts by summing e.memsize()
		data[size_idx] += 1;
		data[memo_idx] += e.memsize();
		data[indx_idx] += e.memsize();
		data[addr_idx] += e.memsize();

		// copy current data to new buffer
		nbuff += head_size;
		obuff += head_size;

		size_t data_size = e.memsize() - (e.size() << 1) - head_size;

		memcpy(nbuff, obuff, data_size);

		// append e.data to new buffer data
		nbuff += data_size;
		obuff += data_size;

		memcpy(nbuff, &e.data, e.memsize());

		// copy idx_buff up to idx
		nbuff += e.memsize();
		memcpy(nbuff, obuff, idx);

		// copy idx to its place
		nbuff += idx;
		nbuff[0] = odata[memo_idx] + 1;
		nbuff += 1;

		// copy rest of idx_buff and addr_buffer
		memcpy(nbuff, obuff, (odata[size_idx] << 1) - idx);

		// append new addr to addr buff
		nbuff += 2*odata[size_idx] - idx;
		nbuff[0] = head_size + data_size;
	}

	expr operator[](uint64_t idx) {
		expr e;
		fill_with_op(idx, e);
		return e;
	}

	void printBuffer() {
		for(uint64_t i = 0; i < memsize(); i++) {
			printf("%u ", data[i]);
		}
		printf("\n");
	}

};
