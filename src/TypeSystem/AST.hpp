#ifndef AST_HPP
#define AST_HPP

#include <vector>
#include <cstddef>

#include "SymbolTable.hpp"

typedef unsigned long long key;

#define INVALID_KEY ((key)-1)

// holds a reference to a symbol that can be
// accessed using the ASTManager object.
union SymbolKey { key value; };

// holds a reference to a term that can be
// accessed using the ASTManager object.
union ASTNodeKey { key value; };

// holds a reference to a members list that can be
// accessed using the ASTManager object.
union ASTMembersKey { key value; };

struct ASTNode {
  #define type_idx(id) (1 << 31 | id)
  #define term_idx(id) (~(1 << 31) & id)

  enum ASTKind {
    // types
    AST_TYPE_UNSPECIFIED          = type_idx(0),
    AST_TYPE_ANY                  = type_idx(1),
    AST_TYPE_ARROW                = type_idx(2),
    AST_TYPE_SYMBOL               = type_idx(3),

    // terms
		AST_TERM_SYMBOL               = term_idx(1),
		AST_TERM_LAMBDA               = term_idx(2),
		AST_TERM_FUNCTION_CALL        = term_idx(3),
		AST_TERM_FUNCTION_DECLARATION = term_idx(4),
		AST_TERM_VARIABLE_DECLARATION = term_idx(5),
		AST_TERM_ARGUMENT             = term_idx(6),
	};

	// The kind of the term
	ASTKind node_kind;

	// Key to the child members of the node
	ASTMembersKey node_members_key;
};

struct ASTManager;

// Creates a ASTManager
ASTManager *ASTManagerCreate();

// Destroy a ASTManager
void ASTManagerDestroy(ASTManager *s);

// Get the key of the child of the ASTNode with key 'root'
ASTNodeKey ASTGetChildNodeKey(ASTManager *, ASTNodeKey root, size_t child_idx);

// Get the key of the child of the ASTNode 'root'
ASTNodeKey ASTGetChildNodeKey(ASTManager *, ASTNode &root, size_t child_idx);

// Get the ASTNode from a key
ASTNode &ASTGetNode(ASTManager *, ASTNodeKey key);

// Get the child_idx'ith child ASTNode from the ASTNode 'root'
ASTNode &ASTGetChildNode(ASTManager *, ASTNode &root, size_t child_idx);

// Print the AST from the root of key 'k'
void ASTPrint(ASTManager *, ASTNodeKey k);

// Create a ASTNode of kind Symbol Term : id
// It have only one child, that is a key to the SymbolRegistry.
ASTNodeKey ASTSymbolTermNode(ASTManager*, const char* id);

// Create a ASTNode of kind Symbol Type.
// It have only one child, that is a key to the SymbolRegistry.
ASTNodeKey ASTSymbolTypeNode(ASTManager*, const char* id);

// Create a ASTNode of kind Unspecified Type.
// It have no children.
ASTNodeKey ASTUnspecifiedTypeNode(ASTManager*);

// Create a ASTNode of kind Any Type.
// It have no children.
ASTNodeKey ASTAnyTypeNode(ASTManager*);

// Create a ASTNode of kind Arrow Type : from -> to.
// It has to children, the ASTNode of key 'from', and
// the ASTNode of key 'to'.
ASTNodeKey ASTArrowTypeNode(ASTManager*, ASTNodeKey from, ASTNodeKey to);


#endif
