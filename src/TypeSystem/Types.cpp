#include "Types.hpp"
#include "TypeSystem/SymbolTable.hpp"

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#define BUCKET_SIZE_LOG2 7
#define CONTEXT_BUCKET_SIZE_LOG2 5

#define BUCKET_SIZE (1 << BUCKET_SIZE_LOG2)
#define CONTEXT_BUCKET_SIZE (1 << BUCKET_SIZE_LOG2)

#define TERM_INDEX_MASK (((key)1 << BUCKET_SIZE_LOG2) - 1)
#define TERM_BUCKET_MASK ~TERM_INDEX_MASK

#define CONTEXT_INDEX_MASK (((key)1 << CONTEXT_BUCKET_SIZE_LOG2) - 1)
#define CONTEXT_TERM_BUCKET_MASK ~CONTEXT_INDEX_MASK

#define build_key(bucket, index) (index & (bucket << BUCKET_SIZE_LOG2))
#define build_ctx_key(bucket, index) (index & (bucket << CONTEXT_BUCKET_SIZE_LOG2))

#define key_bucket(k) ((k & TERM_BUCKET_MASK) >> BUCKET_SIZE_LOG2)
#define key_index(k) (k & TERM_INDEX_MASK)

#define ctx_key_bucket(k) ((k & CONTEXT_BUCKET_MASK) >> CONTEXT_BUCKET_SIZE_LOG2)
#define ctx_key_index(k) (k & CONTEXT_INDEX_MASK)


inline static size_t ceil(size_t x, size_t y) {
  return x == 0 ? 0 : 1 + ((x - 1) / y);
}

/*
 * @brief Context structure, it links declarations terms to types terms.
 * Internally all data is stored lenarly on byckets of BUCKET_SIZE size.
 */
struct context {
	// List of buckets, each bucket holds BUCKET_SIZE keys to storage terms
	key** types;
	// List of buckets, each bucket holds BUCKET_SIZE keys to storage terms.
	key** terms;

	// Total number of bindings on this context
	size_t size;

	// key to the previous context on the storage
	key prev_ctx;

};

struct storage {
  // Bucket list of terms, each bucket have BUCKET_SIZE terms
  // contiguous in memory.
  term **terms;

	// Bucket list of keys, each bucket have BUCKET_SIZE keys
	// contiguous in memory, each term 't' have a 'mem_key' that
	// is a key to a element on this bucket list, the first
	// element 'n' is the number of members on the term 't', and
	// the following 'n' keys are keys of the members of 't', what
	// the keys points to are kind dependent, for instance, it could
	// be a term on 'terms' bucket list, or a key in the 'symbols' for
	// a identifier.
  key **childs;

	context** contexts;

	size_t ctx_count;

  size_t t_count;
  size_t c_count;

  symbol_registry *symbols;

	key root_ctx = -1;

  key term_kind_key = -1;
  key term_unsp_key = -1;
};

key storage_new_empty_context(storage* strg, key prev) {
	bool inc = strg->ctx_count > 0 && ((strg->ctx_count % CONTEXT_BUCKET_SIZE) == 0);

	if(inc) {
		context** old_contexts = strg->contexts;

		strg->contexts = new context*[1 + ceil(strg->ctx_count, CONTEXT_BUCKET_SIZE)];

		for(size_t i = 0; i < 1 + ceil(strg->ctx_count, CONTEXT_BUCKET_SIZE); i++) {
			strg->contexts[i] = old_contexts[i];
		}

		delete[] old_contexts;
	}

	key b = strg->ctx_count / CONTEXT_BUCKET_SIZE;
	key i = strg->ctx_count % CONTEXT_BUCKET_SIZE;

	context* ctx = &strg->contexts[b][i];

	ctx->size = 0;

	ctx->terms = new key*[1];
  ctx->terms[0] = new key[CONTEXT_BUCKET_SIZE];

	ctx->types = new key*[1];
  ctx->types[0] = new key[CONTEXT_BUCKET_SIZE];

	ctx->prev_ctx = prev;

	return build_ctx_key(b, i);
}

key context_add_binding(storage* strg, key ctx_key, key term_key, key type_key) {
	// TODO: check if term_key is already inserted on the current ctx

	context* ctx = storage_get_context(strg, ctx_key);

	bool inc = ctx->size > 0 && ((ctx->size % CONTEXT_BUCKET_SIZE) == 0);

	if(inc) {
		key** old_terms = ctx->terms;
		key** old_types = ctx->types;

		ctx->terms = new key*[1 + ceil(ctx->size, CONTEXT_BUCKET_SIZE)];
		ctx->types = new key*[1 + ceil(ctx->size, CONTEXT_BUCKET_SIZE)];

		for(size_t i = 0; i < 1 + ceil(ctx->size, CONTEXT_BUCKET_SIZE); i++) {
			ctx->terms[i] = old_terms[i];
			ctx->types[i] = old_types[i];
		}

		delete[] old_terms;
		delete[] old_types;
	}

	key b = ctx->size / CONTEXT_BUCKET_SIZE;
	key i = ctx->size % CONTEXT_BUCKET_SIZE;

	ctx->terms[b][i] = term_key;
	ctx->types[b][i] = type_key;

	ctx->size = ctx->size + 1;

	return build_ctx_key(b, i);
}

void context_destroy(context* ctx) {
	for(size_t i = 0; i <= ctx->size / CONTEXT_BUCKET_SIZE; i++) {
		delete[] ctx->types[i];
		delete[] ctx->terms[i];
	}

	delete[] ctx->types;
	delete[] ctx->terms;
}

key context_get_binding_type_key(storage* strg, key ctx_key, key decl) {
	if(ctx_key == INVALID_KEY) {
		return INVALID_KEY;
	}

	context* ctx = storage_get_context(strg, ctx_key);

	long long i = ctx->size / CONTEXT_BUCKET_SIZE;

	for(; i >= 0; i--) {
		long long j = CONTEXT_BUCKET_SIZE - 1;

		for(; j >= 0; j--) {
			if(ctx->terms[i][j] == decl) {
				return ctx->types[i][j];
			}
		}
	}

	return context_get_binding_type_key(strg, ctx->prev_ctx, decl);
}

storage *storage_create() {
  storage *strg = new storage();

  strg->terms = new term *[1];
  strg->terms[0] = new term[BUCKET_SIZE];

  strg->childs = new key *[1];
  strg->childs[0] = new key[BUCKET_SIZE];

  strg->contexts = new context*[1];
  strg->contexts[0] = new context[CONTEXT_BUCKET_SIZE];

  strg->t_count = 0;
  strg->c_count = 0;

	strg->ctx_count = 0;

  strg->symbols = symbol_registry_create();

	strg->root_ctx = storage_new_empty_context(strg, -1);

	// Here we add all the base types and terms to the root context.

	// Add Unspecified Type, this is used on declaration that types are not
	// specified directly, it will be the responsability of the type checker
	// to infer this type.
	strg->term_unsp_key = storage_insert(strg, strg->root_ctx, term::TYPE_UNS, 0, 0);

	// Add Any Type, that is the kind of Types.
  strg->term_kind_key = storage_insert(strg, strg->root_ctx, term::TYPE_ANY, 0, 0);

  return strg;
}

void storage_destroy(storage *strg) {
  symbol_registry_destroy(strg->symbols);

  for (size_t i = 0; i <= strg->t_count / BUCKET_SIZE; i++) {
    delete[] strg->terms[i];
  }

  for (size_t i = 0; i <= strg->c_count / BUCKET_SIZE; i++) {
    delete[] strg->childs[i];
  }

	for (size_t i = 0; i <= strg->ctx_count / BUCKET_SIZE; i++) {
		for(size_t j = 0; j < BUCKET_SIZE; j++) {
			context_destroy(&strg->contexts[i][j]);
		}

		delete[] strg->contexts[i];
  }


  delete[] strg->terms;
  delete[] strg->childs;
  delete[] strg->contexts;

  delete strg;
}

name storage_name_create(storage *strg, const char *id) {
  name n;

  ull k = symbol_registry_get_entry(strg->symbols, id);

  if (k != UNDEFINED_SYMBOL) {

    n.symbol_id = k;

    return n;
  }

  n.symbol_id = symbol_registry_set_entry(strg->symbols, id);

  return n;
}

key storage_insert(storage *strg, key ctx, term::kind kind, key *childs, size_t count) {
  bool inc = (strg->t_count > 0) && ((strg->t_count % BUCKET_SIZE) == 0);

  if (inc) {
    // increase strg->terms
    size_t n = 1 + ceil(strg->t_count, BUCKET_SIZE);

    term **old_terms = strg->terms;
    // name** old_names = strg->names;

    strg->terms = new term *[n];
    // strg->names = new name*[n];

    size_t i = 0;

    for (; i <= strg->t_count / BUCKET_SIZE; i++) {
      strg->terms[i] = old_terms[i];
      // strg->names[i] = old_names[i];
    }

    for (; i <= (strg->c_count + count) / BUCKET_SIZE; i++) {
      strg->terms[i] = new term[BUCKET_SIZE];
      // strg->names[i] = new name[BUCKET_SIZE];
    }

    delete[] old_terms;
    // delete[] old_names;
  }

  key t_b = strg->t_count / BUCKET_SIZE;
  key t_i = strg->t_count % BUCKET_SIZE;

  // strg->names[t_b][t_i] = nm;
  strg->terms[t_b][t_i].term_kind = kind;

  inc = (strg->c_count % BUCKET_SIZE) + count < BUCKET_SIZE;

  if (inc) {
    // increase strg->childs
    size_t n = ceil(strg->c_count + count, BUCKET_SIZE);

    key **old_childs = strg->childs;

    strg->childs = new key *[n];

    size_t i = 0;

    for (; i <= strg->c_count / BUCKET_SIZE; i++) {
      strg->childs[i] = old_childs[i];
    }

    for (; i <= (strg->c_count + count) / BUCKET_SIZE; i++) {
      strg->childs[i] = new key[BUCKET_SIZE];
    }

    delete[] old_childs;
  }

  key c_b = strg->c_count / BUCKET_SIZE;
  key c_i = strg->c_count % BUCKET_SIZE;

  strg->childs[c_b][c_i] = count;

  for (size_t b = 0; b <= (c_i + count) / BUCKET_SIZE; b++) {
    for (size_t i = 0; i < count; i++) {
      strg->childs[c_b + b][c_i + i + 1] = childs[i];
    }
  }

  strg->c_count = strg->c_count + count + 1;
  strg->t_count = strg->t_count + 1;

  strg->terms[t_b][t_i].mem_key = build_key(c_b, c_i);
	strg->terms[t_b][t_i].ctx_key = ctx;

  return build_key(t_b, t_i);
}

key storage_get_child_key(storage *strg, term &t, size_t i) {
  key c_key = t.mem_key;

  size_t bck = key_bucket(c_key);
  size_t idx = key_index(c_key);

  assert(strg->childs[bck][idx] > i /* out of bounds */);

  // first element is the number of childs
  idx += 1;

  bck += (idx + i) / BUCKET_SIZE;

  return strg->childs[bck][(idx + i) % BUCKET_SIZE];
}

term &storage_get_term_from_key(storage *strg, key idx) {
  return strg->terms[key_bucket(idx)][key_index(idx)];
}

term &storage_get_child_term(storage *strg, term &t, size_t i) {
  return storage_get_term_from_key(strg, storage_get_child_key(strg, t, i));
}

void print_term_rec(storage *strg, key term_key) {
	term &t = storage_get_term_from_key(strg, term_key);

	key ctx = t.ctx_key;

	key app_func, app_argu;

	key arg_name, arg_type;

	key fun_name, fun_argu, fun_body, fun_type;

	key lam_argu, lam_body;

	key var_name, var_type, var_valu;

	key arr_arg0, arr_arg1;

	switch (t.term_kind) {

	case term::TERM_APP:
		app_func = storage_get_child_key(strg, t, 0);
	  app_argu = storage_get_child_key(strg, t, 1);

		print_term(strg, app_func);

    printf(" ");

    print_term(strg, app_argu);

    break;

	case term::TERM_ARG:
		arg_name = storage_get_child_key(strg, t, 0);
    arg_type = context_get_binding_type_key(strg, ctx, arg_name);

		print_term(strg, arg_name);
		printf(" : ");
		print_term(strg, arg_type);

		break;

	case term::TERM_FUN:
    fun_name = storage_get_child_key(strg, t, 0);
    fun_argu = storage_get_child_key(strg, t, 1);
    fun_body = storage_get_child_key(strg, t, 2);

    fun_type = context_get_binding_type_key(strg, ctx, fun_name);

    printf("func ");

    print_term(strg, fun_name);

    printf(" ( ");

    print_term(strg, fun_argu);

    printf(" ) ");

    printf(" = ( ");
    print_term(strg, fun_body);
    printf(" ) ");

    printf(" : ");

    print_term(strg, fun_type);

    printf(";");

		break;

	case term::TERM_LAM:
		lam_argu = storage_get_child_key(strg, t, 0);
		lam_body = storage_get_child_key(strg, t, 1);

		printf("λ");

    print_term(strg, lam_argu);

		printf(".");

    print_term(strg, lam_body);

    printf(";");

		break;

  case term::TERM_VAR:
		var_name = storage_get_child_key(strg, t, 0);
		var_valu = storage_get_child_key(strg, t, 1);
		var_type = context_get_binding_type_key(strg, ctx, term_key);

		printf("let ");

    print_term(strg, var_name);

		printf(" : ");

		print_term(strg, var_type);

		printf(" = ");

    print_term(strg, var_valu);

    printf(";");

    break;

	case term::TERM_SYM:
    printf("%s", symbol_registry_get_symbol(strg->symbols, t.mem_key));
    break;

	case term::TYPE_SYM:
    print_term(strg, storage_get_child_key(strg, t, 0));
    break;

	case term::TYPE_ANY:
    printf("*");
    break;

	case term::TYPE_UNS:
    printf("?");
    break;

	case term::TYPE_ARR:
		arr_arg0 = storage_get_child_key(strg, t, 0);
		arr_arg1 = storage_get_child_key(strg, t, 1);

		printf("(");
    print_term(strg, arr_arg0);
		printf(" -> ");
    print_term(strg, arr_arg1);
    printf(")");

    break;
  }
}

key root_context(storage *strg) { return strg->root_ctx; }

key type_kind(storage *strg) { return strg->term_kind_key; }

key type_unspecified(storage *strg) { return strg->term_unsp_key; }

key term_symbol(storage *strg, key ctx, const char *id) {
  name n = storage_name_create(strg, id);

  key c[1] = {n.symbol_id};

  return storage_insert(strg, ctx, term::TERM_SYM, c, 1);
}

key term_var_decl(storage *strg, key ctx, const char *id, key value_key, key type_key) {
  key sym_key = term_symbol(strg, id);

  key c[] = { sym_key, value_key };

	key var = storage_insert(strg, ctx, term::TERM_VAR, c, 2);

	context_add_binding(strg, ctx, sym_key, type_key);

	return var;
}


key term_lambda_decl(storage *strg, key ctx, key arg_symbol, key arg_type, key body, key body_type) {}

key term_function_decl(storage *strg, key ctx, key arg_symbol, key arg_type, key body, key body_type) {
	key found = context_get_binding_type_key(strg, ctx, arg_type);

	if(found == INVALID_KEY) {
		// TODO: handle error because type is not defined
		exit(1);
	}

	key c[2] = { arg_symbol, body };

	key func_ctx = storage_new_empty_context(strg, ctx);

	context_add_binding(strg, func_ctx, arg_symbol, arg_type);

	return storage_insert(strg, ctx, term::TERM_FUN, c, 2);
}

key term_call(storage *strg, key ctx, key fun, key arg) {
  key c[2] = {fun, arg};

  return storage_insert(strg, ctx, term::TERM_APP, c, 2);
}

// void create_call() {}
