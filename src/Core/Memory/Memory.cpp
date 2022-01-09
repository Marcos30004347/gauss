#include "Memory.hpp"


struct mem_header {
	size_t size;
	size_t refs;
	short bucket;
};

struct mem_chunck {
	void* data;
	char size;

};


void* bucket[8];


void *malloc(size_t size) {
}
void *realloc(void *ptr, size_t size) {
}

void free(void *ptr) {
}

void incref(void *ptr) {
	((mem_header*)((char*)ptr - sizeof(mem_header)))->refs += 1;
}

void decref(void *ptr) {
	mem_header* h = (mem_header*)((char*)ptr - sizeof(mem_header));

	if(h->refs == 1) {
		free(ptr);
	} else {
		h->refs -= 1;
	}
}
