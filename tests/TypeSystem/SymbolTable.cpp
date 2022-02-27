#include "TypeSystem/SymbolTable.hpp"

#include <cassert>
#include <string.h>

int main() {
	symbol_registry* t = symbol_registry_create();

	ull i = symbol_registry_set_entry(t, "test");
	ull j = symbol_registry_set_entry(t, "test");
	ull k = symbol_registry_set_entry(t, "testtest");

	assert(strcmp(symbol_registry_get_symbol(t, i), "test") == 0);
	assert(strcmp(symbol_registry_get_symbol(t, j), "test") == 0);
	assert(strcmp(symbol_registry_get_symbol(t, k), "testtest") == 0);

	symbol_registry_destroy(t);

  return 0;
}
