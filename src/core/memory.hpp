#ifndef CORE_MEMORY_H
#define CORE_MEMORY_H

namespace core {
namespace memory {

void* alloc(unsigned long size);
void* copy(void *dst, void *src, long count);
void free(void* ptr);

} // namespace memory
} // namespace core

#endif