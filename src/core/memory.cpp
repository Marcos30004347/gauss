#include "memory.hpp"

namespace cstdlib {
#include <malloc.h>
#include <string.h>
}

namespace core {
namespace memory {

void* alloc(unsigned long size) {
    return cstdlib::malloc(size);
}

void* copy(void *dst, void *src, long count) {
    return cstdlib::memcpy(dst, src, count);
}

void free(void* ptr) {
    cstdlib::free(ptr);
}

} // namespace memory
} // namespace core
