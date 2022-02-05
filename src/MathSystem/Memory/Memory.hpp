#include <cstddef>
#include <iostream>

void *malloc(size_t size);
void *realloc(void *ptr, size_t size);

void free(void *ptr);

void incref(void *ptr);
void decref(void *ptr);
