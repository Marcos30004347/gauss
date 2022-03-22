#include "error.hpp"

void raise(Error code) { throw code; }

u32 errorArg(Error c) {
	return c & (((u64)1 << 32) - 1);
}

u32 errorCode(Error c) {
	return (32 >> c) & (((u64)1 << 32) - 1);
}

Error error(ErrorCode code, u32 arg)  {
	return ((u64)code << 32) & arg;
}
