#ifndef TYPES_HPP
#define TYPES_HPP

#include <vector>
#include <cstddef>

#include "SymbolTable.hpp"

typedef unsigned long long key;

#define INVALID_KEY ((key)-1)

// holds a reference to a symbol that can be
// accessed using the storage object.
union name_ref { key value; };

// holds a reference to a context that can be
// accessed using the storage object.
union ctx_ref { key value; };

// holds a reference to a term that can be
// accessed using the storage object.
union term_ref { key value; };

// holds a reference to a members list that can be
// accessed using the storage object.
union members_ref { key value; };

// holds a reference to a binding on some context
union binding_ref { key value; };

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

	// The kind of the term
	kind term_kind;

	// Reference to the context of the term.
	// The context object itself can only
	// be accessed using the storage.
	ctx_ref ctx_key;

	// Reference to the members of the term.
	// The term members/childs can only be
	// accessed using the storage.
	members_ref mem_key;
};

struct storage;
struct context;

ctx_ref root_context(storage*);

storage* storage_create();

void storage_destroy(storage *s);

ctx_ref context_create(storage *strg, ctx_ref prev);

void context_destroy(context *ctx);

binding_ref context_add_binding(storage* strg, ctx_ref ctx_key, term_ref term_key, term_ref type_key);

context* storage_get_context(storage* strg, ctx_ref ctx);

term_ref context_get_binding_type_key(storage* strg, ctx_ref ctx_key, term_ref decl);

term_ref storage_get_child_key(storage* strg, key t, size_t idx);

term_ref storage_get_child_key(storage* strg, term& t, size_t idx);

term &storage_get_term_from_key(storage *strg, term_ref idx);

term &storage_get_child_term(storage *strg, term &t, size_t idx);

term_ref storage_insert(storage*, ctx_ref ctx, term::kind, term_ref*, size_t);

term_ref type_kind(storage *storage);

term_ref type_unspecified(storage *strg);

term_ref term_symbol(storage* strg, ctx_ref ctx, const char* id);

void print_term(storage* strg, term_ref term);

#endif
