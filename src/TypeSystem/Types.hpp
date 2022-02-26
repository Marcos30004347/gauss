#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>
#include <cstddef>

#include "SymbolTable.hpp"

typedef ull key;

#define INVALID_KEY ((key)-1)

struct name {
	ull symbol_id;
};

struct term {
  #define type_idx(id) (1 << 31 | id)
  #define term_idx(id) (~(1 << 31) & id)

  enum kind {
    // types
    TYPE_UNS = type_idx(0),
    TYPE_ANY = type_idx(1),
    TYPE_ARR = type_idx(2),
    TYPE_SYM = type_idx(3),

    // terms
		TERM_SYM = term_idx(1),
		TERM_LAM = term_idx(2),
		TERM_FUN = term_idx(3),
		TERM_APP = term_idx(4),
		TERM_VAR = term_idx(5),
		TERM_ARG = term_idx(6),
	};

	// the kind of the term
	kind term_kind;

	// key to the context of this term on the storage.
	key ctx_key;

	// key to the members of this term on the storage.
	key mem_key;
};

struct storage;
struct context;

storage* storage_create();

void storage_destroy(storage *s);

key context_create(storage *strg, key prev);

void context_destroy(context *ctx);

key context_add_binding(storage* strg, key ctx_key, key term_key, key type_key);

key context_get_binding_type_key(storage* strg, key ctx_key, key decl);

context* storage_get_context(storage* strg, key ctx);

term &storage_get_term_from_key(storage *ctx, key idx);

term &storage_get_child_term(storage *ctx, term &t, size_t idx);

key  storage_insert(storage*, key ctx, term::kind, key*, size_t);

key type_kind(storage *ctx);

key type_unspecified(storage *ctx);

key term_symbol(storage* ctx,const char* id);

void print_term(storage* strg, key term);

#endif
