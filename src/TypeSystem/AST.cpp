#include "TypeSystem/AST.hpp"
#include "TypeSystem/Utils.hpp"
#include "TypeSystem/SymbolTable.hpp"

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <cstring>
// #include <cstdio>
#include <initializer_list>


#define BUCKET_SIZE_LOG2 7
#define CONTEXT_BUCKET_SIZE_LOG2 5

#define BUCKET_SIZE ((key)1 << BUCKET_SIZE_LOG2)
#define CONTEXT_BUCKET_SIZE ((key)1 << BUCKET_SIZE_LOG2)

#define TERM_INDEX_MASK (((key)1 << BUCKET_SIZE_LOG2) - 1)
#define TERM_BUCKET_MASK ~TERM_INDEX_MASK

#define CONTEXT_INDEX_MASK (((key)1 << CONTEXT_BUCKET_SIZE_LOG2) - 1)
#define CONTEXT_BUCKET_MASK ~CONTEXT_INDEX_MASK

#define build_key(bucket, index) ((key)index | ((key)bucket << BUCKET_SIZE_LOG2))
#define build_ctx_key(bucket, index) ((key)index | ((key)bucket << CONTEXT_BUCKET_SIZE_LOG2))

#define key_bucket(k) ((k & TERM_BUCKET_MASK) >> BUCKET_SIZE_LOG2)
#define key_index(k) (k & TERM_INDEX_MASK)

#define ctx_key_bucket(k) ((k & CONTEXT_BUCKET_MASK) >> CONTEXT_BUCKET_SIZE_LOG2)
#define ctx_key_index(k) (k & CONTEXT_INDEX_MASK)

#define min(a, b) a > b ? b : a

inline static size_t ceil(size_t x, size_t y) {
  return x == 0 ? 0 : 1 + ((x - 1) / y);
}


struct ASTManager {
	// Nodes bucket list.
  ASTNode **terms;

	// ASTNodeKeys bucket list. This holds the keys to the
	// members of any ASTNode.
  ASTNodeKey **childs;

	// Total terms count
  size_t t_count;

	// Total child count
  size_t c_count;

	// The Symbol registry, it stores all the symbols
	// referenced by any ASTNode of kind AST_..._SYMBOL.
  symbol_registry *symbols;

  ASTNodeKey ASTAnyTypeKey = {};
  ASTNodeKey ASTUnspecifiedTypeKey = {};
};

ASTNodeKey ASTNodeKeyFromUnsignedLongLong(key idx) {
	ASTNodeKey ref = {};
	ref.value = idx;
	return ref;
}

ASTMembersKey ASTMembersKeyFRomUnsignedLongLong(key idx) {
	ASTMembersKey ref = {};

	ref.value = idx;

	return ref;
}

ASTNodeKey ASTCreateNode(ASTManager *strg, ASTNode::ASTKind kind,
                         ASTNodeKey *childs, size_t count);

ASTManager *ASTManagerCreate() {
  ASTManager *strg = new ASTManager();

  strg->terms = new ASTNode *[1];
  strg->terms[0] = new ASTNode[BUCKET_SIZE];

  strg->childs = new ASTNodeKey *[1];
  strg->childs[0] = new ASTNodeKey[BUCKET_SIZE];

  strg->t_count = 0;
  strg->c_count = 0;

  strg->symbols = symbol_registry_create();

	strg->ASTAnyTypeKey = ASTCreateNode(strg, ASTNode::AST_TYPE_UNSPECIFIED, 0, 0);
  strg->ASTUnspecifiedTypeKey = ASTCreateNode(strg, ASTNode::AST_TYPE_ANY, 0, 0);

  return strg;
}

void ASTManagerDestroy(ASTManager *strg) {
  for (size_t i = 0; i <= strg->t_count / BUCKET_SIZE; i++) {
    delete[] strg->terms[i];
  }

  for (size_t i = 0; i <= strg->c_count / BUCKET_SIZE; i++) {
    delete[] strg->childs[i];
  }

  delete[] strg->terms;
  delete[] strg->childs;

  symbol_registry_destroy(strg->symbols);

	delete strg;
}

SymbolKey ASTManagerRegisterSymbol(ASTManager *strg, const char *id) {
  SymbolKey n;

  n.value = symbol_registry_set_entry(strg->symbols, id);

  return n;
}

ASTNodeKey ASTCreateNode(ASTManager *strg,  ASTNode::ASTKind kind, ASTNodeKey *childs, size_t count) {
  bool inc = (strg->t_count > 0) && ((strg->t_count % BUCKET_SIZE) == 0);

	if (inc) {
		// increase strg->terms
    size_t n = 1 + ceil(strg->t_count, BUCKET_SIZE);

		ASTNode **old_terms = strg->terms;

    strg->terms = new ASTNode *[n];

    size_t i = 0;

    for (; i < n - 1; i++) {
      strg->terms[i] = old_terms[i];
    }

		strg->terms[n - 1] = new ASTNode[BUCKET_SIZE];

    delete[] old_terms;
  }

	key t_b = strg->t_count / BUCKET_SIZE;
  key t_i = strg->t_count % BUCKET_SIZE;

  strg->terms[t_b][t_i].node_kind = kind;

	inc = (strg->c_count / BUCKET_SIZE) < ((strg->c_count + count + 1) / BUCKET_SIZE);

	if (inc) {
    // increase strg->childs
		size_t q = 1 + (strg->c_count / BUCKET_SIZE);
		size_t k = 1 + ceil((strg->c_count + count + 1), BUCKET_SIZE);

		// size_t n = ceil(strg->c_count + count + 1, BUCKET_SIZE);
    ASTNodeKey **old_childs = strg->childs;

    strg->childs = new ASTNodeKey *[k];

    size_t i = 0;

    for (; i < q; i++) {
      strg->childs[i] = old_childs[i];
    }

    for (; i < k; i++) {
      strg->childs[i] = new ASTNodeKey[BUCKET_SIZE];
    }

    delete[] old_childs;
  }

	key c_b = strg->c_count / BUCKET_SIZE;
	key c_i = strg->c_count % BUCKET_SIZE;

  strg->childs[c_b][c_i] = ASTNodeKeyFromUnsignedLongLong(count);

	for (size_t i = 0; i < count; i++) {
	 	size_t b = (c_i + 1 + i) / BUCKET_SIZE;
	 	strg->childs[c_b + b][(c_i + 1 + i) % BUCKET_SIZE] = childs[i];
	}

  strg->c_count = strg->c_count + count + 1;
  strg->t_count = strg->t_count + 1;

  strg->terms[t_b][t_i].node_members_key = ASTMembersKeyFRomUnsignedLongLong(build_key(c_b, c_i));

  return ASTNodeKeyFromUnsignedLongLong(build_key(t_b, t_i));
}

ASTNodeKey ASTGetChildNodeKey(ASTManager *strg, ASTNodeKey trm, size_t i) {
	ASTNode& t = ASTGetNode(strg, trm);

  ASTMembersKey c_key = t.node_members_key;

  size_t bck = key_bucket(c_key.value);
  size_t idx = key_index(c_key.value);

  assert(strg->childs[bck][idx].value > i /* out of bounds */);

  idx += 1;

  bck += (idx + i) / BUCKET_SIZE;

  return strg->childs[bck][(idx + i) % BUCKET_SIZE];
}

ASTNodeKey ASTGetChildNodeKey(ASTManager *strg, ASTNode &t, size_t i) {
  ASTMembersKey c_key = t.node_members_key;

  size_t bck = key_bucket(c_key.value);

  size_t idx = key_index(c_key.value);

	assert(strg->childs[bck][idx].value > i);

  idx += 1;

  bck += (idx + i) / BUCKET_SIZE;

  return strg->childs[bck][(idx + i) % BUCKET_SIZE];
}

ASTNode &ASTGetNode(ASTManager *strg, ASTNodeKey idx) {
  return strg->terms[key_bucket(idx.value)][key_index(idx.value)];
}

ASTNode &ASTGetChildNode(ASTManager *strg, ASTNode &t, size_t i) {
  return ASTGetNode(strg, ASTGetChildNodeKey(strg, t, i));
}

void ASTPrintRec(ASTManager *strg, ASTNodeKey term_key) {
	// ASTNode &t = ASTGetNode(strg, term_key);

	// // ctx_ref ctx = t.ctx_key;

	// ASTNodeKey app_func, app_argu;

	// ASTNodeKey arg_name, arg_type;

	// ASTNodeKey fun_name, fun_argu, fun_body, fun_type;

	// ASTNodeKey lam_argu, lam_body;

	// ASTNodeKey var_name, var_type, var_valu;

	// ASTNodeKey arr_arg0, arr_arg1;

	// ASTNodeKey sym_symb;

	// ASTNodeKey int_symb;

	// ASTNodeKey bin_symb, bin_left, bin_righ;

	// ASTNodeKey arg_list, nxt_list;

	// ASTNodeKey if_cond, if_body, if_else;

	// switch (t.node_kind) {
	// case ASTNode::AST_TERM_FLOW_IF_ELSE:
	//   if_cond = ASTGetChildNodeKey(strg, t, 0);
	// 	if_body = ASTGetChildNodeKey(strg, t, 1);
	//   if_else = ASTGetChildNodeKey(strg, t, 2);

	// 	printf("if ");

	// 	ASTPrintRec(strg, if_cond);

	// 	printf(" then ");

	// 	ASTPrintRec(strg, if_body);

	// 	if(if_else.value != INVALID_KEY) {
	// 		ASTPrintRec(strg, if_else);
	// 	}

	// 	printf(" endif");

	// 	break;
  // case ASTNode::AST_TERM_INTEGER_LITERAL:
	// 	int_symb = ASTGetChildNodeKey(strg, t, 0);

	// 	ASTPrintRec(strg, int_symb);

	// 	break;
  // case ASTNode::AST_TERM_BINARY_OPERATION:

	// 	bin_symb = ASTGetChildNodeKey(strg, t, 0);
	// 	bin_left = ASTGetChildNodeKey(strg, t, 1);
	//   bin_righ = ASTGetChildNodeKey(strg, t, 2);

	// 	printf("(");
	// 	ASTPrintRec(strg, bin_left);
	// 	printf(" ");
	// 	ASTPrintRec(strg, bin_symb);
	// 	printf(" ");
	// 	ASTPrintRec(strg, bin_righ);
	// 	printf(")");

	// 	break;
  // case ASTNode::AST_TERMS_SEQUENCE:
	// 	arg_list = ASTGetChildNodeKey(strg, t, 0);
	// 	nxt_list = ASTGetChildNodeKey(strg, t, 1);

	// 	ASTPrintRec(strg, arg_list);

	// 	printf(";\n");

	// 	if(nxt_list.value != INVALID_KEY) {
	// 		ASTPrintRec(strg, nxt_list);
	// 	}

	// 	break;

	// case ASTNode::AST_TERM_FUNCTION_CALL:
	// 	app_func = ASTGetChildNodeKey(strg, t, 0);
	//   app_argu = ASTGetChildNodeKey(strg, t, 1);

	// 	ASTPrintRec(strg, app_func);

  //   printf(" ");

  //   ASTPrintRec(strg, app_argu);

  //   break;

	// case ASTNode::AST_TERM_ARGUMENT:
	// 	arg_name = ASTGetChildNodeKey(strg, t, 0);
  //   // arg_type = context_get_binding_type_key(strg, ctx, arg_name);

	// 	ASTPrintRec(strg, arg_name);
	// 	printf(" : ");
	// 	ASTPrintRec(strg, arg_type);

	// 	break;

	// case ASTNode::AST_TERM_FUNCTION_DECLARATION:
  //   fun_name = ASTGetChildNodeKey(strg, t, 0);
  //   fun_argu = ASTGetChildNodeKey(strg, t, 1);
  //   fun_body = ASTGetChildNodeKey(strg, t, 2);

  //   // fun_type = context_get_binding_type_key(strg, ctx, fun_name);

  //   printf("func ");

  //   ASTPrintRec(strg, fun_name);

  //   printf(" ( ");

  //   ASTPrintRec(strg, fun_argu);

  //   printf(" ) ");

  //   printf(" = ( ");
  //   ASTPrintRec(strg, fun_body);
  //   printf(" ) ");

  //   printf(" : ");

  //   ASTPrintRec(strg, fun_type);

  //   printf(";");

	// 	break;

	// case ASTNode::AST_TERM_LAMBDA:
	// 	lam_argu = ASTGetChildNodeKey(strg, t, 0);
	// 	lam_body = ASTGetChildNodeKey(strg, t, 1);

	// 	printf("λ");

  //   ASTPrintRec(strg, lam_argu);

	// 	printf(".");

  //   ASTPrintRec(strg, lam_body);

  //   printf(";");

	// 	break;

  // case ASTNode::AST_TERM_VARIABLE_DECLARATION:
	// 	var_name = ASTGetChildNodeKey(strg, t, 0);
	// 	var_valu = ASTGetChildNodeKey(strg, t, 1);
	// 	// var_type = context_get_binding_type_key(strg, ctx, term_key);

	// 	printf("let ");

  //   ASTPrintRec(strg, var_name);

	// 	printf(" : ");

	// 	ASTPrintRec(strg, var_type);

	// 	printf(" = ");

  //   ASTPrintRec(strg, var_valu);

  //   printf(";");

  //   break;

	// case ASTNode::AST_TERM_SYMBOL:
	// 	sym_symb = ASTGetChildNodeKey(strg, t, 0);

  //   printf("%s", symbol_registry_get_symbol(strg->symbols, sym_symb.value));

	// 	break;

	// case ASTNode::AST_TYPE_SYMBOL:
	// 	sym_symb = ASTGetChildNodeKey(strg, t, 0);

	// 	printf("%s", symbol_registry_get_symbol(strg->symbols, sym_symb.value));

	// 	break;

	// case ASTNode::AST_TYPE_ANY:
  //   printf("Any");
  //   break;

	// case ASTNode::AST_TYPE_UNSPECIFIED:
  //   printf("Unspecified");
  //   break;

	// case ASTNode::AST_TYPE_ARROW:
	// 	arr_arg0 = ASTGetChildNodeKey(strg, t, 0);
	// 	arr_arg1 = ASTGetChildNodeKey(strg, t, 1);

	// 	printf("(");
  //   ASTPrintRec(strg, arr_arg0);
	// 	printf(" -> ");
  //   ASTPrintRec(strg, arr_arg1);
  //   printf(")");

  //   break;
  // }
}

void ASTPrint(ASTManager *strg, ASTNodeKey term) {
	ASTPrintRec(strg, term);
}

ASTNodeKey ASTSymbolTermNode(ASTManager *manager, const char *id) {
	SymbolKey sym = ASTManagerRegisterSymbol(manager, id);

	ASTNodeKey child[] = { ASTNodeKeyFromUnsignedLongLong(sym.value) };

	return ASTCreateNode(manager, ASTNode::AST_TERM_SYMBOL, child, 1);
}

ASTNodeKey ASTSymbolTypeNode(ASTManager *manager, const char *id) {
	SymbolKey sym = ASTManagerRegisterSymbol(manager, id);

	ASTNodeKey child[] = { ASTNodeKeyFromUnsignedLongLong(sym.value) };

	return ASTCreateNode(manager, ASTNode::AST_TYPE_SYMBOL, child, 1);
}

ASTNodeKey ASTArrowTypeNode(ASTManager *manager, ASTNodeKey from, ASTNodeKey to) {
	ASTNodeKey child[] = { from, to };

	return ASTCreateNode(manager, ASTNode::AST_TYPE_ARROW, child, 2);
}

ASTNodeKey ASTUnspecifiedTypeNode(ASTManager* manager) {
	return manager->ASTUnspecifiedTypeKey;
}

ASTNodeKey ASTAnyTypeNode(ASTManager* manager) {
	return manager->ASTAnyTypeKey;
}

ASTNodeKey ASTIntegerLiteralNode(ASTManager * manager, const char* int_sym) {
	assert(IsIntegerLiteralASCIIString(int_sym));

	SymbolKey sym = ASTManagerRegisterSymbol(manager, int_sym);

	ASTNodeKey child[] = { ASTNodeKeyFromUnsignedLongLong(sym.value) };

	return ASTCreateNode(manager, ASTNode::AST_TERM_INTEGER_LITERAL, child, 1);
}

ASTNodeKey ASTBinaryOperationNode(ASTManager *manager, ASTNodeKey op, ASTNodeKey left, ASTNodeKey right) {

	ASTNodeKey child[] = { op, left, right };

	return ASTCreateNode(manager, ASTNode::AST_TERM_BINARY_OPERATION, child, 3);
}

ASTNodeKey ASTLambdaDeclarationNode(ASTManager* manager, ASTNodeKey arg, ASTNodeKey body) {
	ASTNodeKey child[] = { arg, body };

	return ASTCreateNode(manager, ASTNode::AST_TERM_LAMBDA, child, 2);
}

ASTNodeKey ASTFunctionDeclarationNode(ASTManager* manager, ASTNodeKey funcname, ASTNodeKey arguments, ASTNodeKey return_type, ASTNodeKey body) {

	ASTNodeKey child[] = { funcname, arguments, return_type, body };

	return ASTCreateNode(manager, ASTNode::AST_TERM_FUNCTION_DECLARATION, child, 4);
}

ASTNodeKey ASTVariableDeclarationNode(ASTManager* manager, ASTNodeKey varname, ASTNodeKey value, ASTNodeKey type) {
	ASTNodeKey child[] = { varname, value, type };

	return ASTCreateNode(manager, ASTNode::AST_TERM_VARIABLE_DECLARATION, child, 3);
}

ASTNodeKey ASTSequencialTermsNode(ASTManager* manager, ASTNodeKey val, ASTNodeKey seq) {
	assert(ASTGetNode(manager, seq).node_kind == ASTNode::AST_TERMS_SEQUENCE);

	ASTNodeKey child[] = { val, seq };

	return ASTCreateNode(manager, ASTNode::AST_TERMS_SEQUENCE, child, 2);
}


ASTNodeKey ASTIfNode(ASTManager* manager, ASTNodeKey cond, ASTNodeKey body, ASTNodeKey elif) {

	ASTNodeKey child[] = { cond, body, elif };

	return ASTCreateNode(manager, ASTNode::AST_TERM_FLOW_IF_ELSE, child, 3);

}
