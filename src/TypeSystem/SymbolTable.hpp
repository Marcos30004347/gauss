#ifndef SYMBOL_TABLE
#define SYMBOL_TABLE

#include <cstddef>

#define ull unsigned long long

#define UNDEFINED_SYMBOL ((ull)-1)

struct symbol_registry;

symbol_registry *symbol_registry_create();

void symbol_registry_destroy(symbol_registry *table);

ull symbol_registry_set_entry(symbol_registry *table, const char *id);

// ull symbol_registry_get_entry(symbol_registry *table, const char *id);

const char* symbol_registry_get_symbol(symbol_registry* registry,  ull key);

#endif
