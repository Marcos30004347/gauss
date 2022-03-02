#include "Utils.hpp"

#include <cstddef>
#include <string.h>

bool IsIntegerLiteralASCIIString(const char *id) {
	size_t k = strlen(id);

	for(size_t i = 0; i < k; i++) {
		if(id[i] > '9' || id[i] < '0') {
			return false;
		}
	}
	return true;
}
