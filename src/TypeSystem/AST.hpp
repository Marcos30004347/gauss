#ifndef AST_HPP
#define AST_HPP

#include <vector>
#include <cstddef>

#include "SymbolTable.hpp"

/*
 * This is the file containing the Abstract Syntax Tree(AST) definitions,
 * it defines a Node class, and a Manager class for AST entities. All nodes
 * are allocated and managed by the Manager, this means that a used should
 * request the manager for a node with a given ASTKind and children and it is
 * the responability of the Manager to allocate and deallocate resources for
 * that Nodes, the Manager will return a ASTNodeKey object, that is the object
 * that can be passed to the Manager to access the Node and its children.
 *
 * Internally all the Nodes are stored in Bulked Lists of a constant length,
 * this makes possible for a fast and cache friendly store of nodes. Assuming
 * that all Nodes get parsed and allocated linearly, children nodes will be
 * close in memory.
 *
 * References:
 * [1] Herb Sutter's talk on machine architecture - https://www.youtube.com/watch?v=L7zSU9HI-6I
 * [2] Software optimization resources - https://www.agner.org/optimize/
 * [3] What every programmer should know about memory - https://lwn.net/Articles/250967/
 */



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

struct Kind {
	char data;
};

struct ASTNode {

  #define type_idx(id) ((char)1 << 6 | id)
  #define term_idx(id) ((char)3 << 6 | id)
  #define util_idx(id) ((char)2 << 6 | id)

	enum ASTKind : unsigned char {
    // types
    AST_TYPE_UNSPECIFIED          = type_idx(1),
    AST_TYPE_ANY                  = type_idx(2),
    AST_TYPE_ARROW                = type_idx(3),
    AST_TYPE_SYMBOL               = type_idx(4),

		// TODO: AST_TYPE_TUPLE
		// TODO: AST_TYPE_UNION
		// TODO: AST_TYPE_DEPENDENT_TUPLE
		// TODO: AST_TYPE_DEPENDENT_UNION

		// terms
    AST_TERM_SYMBOL               = term_idx(1),
    AST_TERM_LAMBDA               = term_idx(2),
    AST_TERM_FUNCTION_CALL        = term_idx(3),
    AST_TERM_FUNCTION_DECLARATION = term_idx(4),
    AST_TERM_VARIABLE_DECLARATION = term_idx(5),
    AST_TERM_ARGUMENT             = term_idx(6),
    AST_TERM_BINARY_OPERATION     = term_idx(7),
    AST_TERM_INTEGER_LITERAL      = term_idx(9),
    AST_TERM_FLOW_IF_ELSE         = term_idx(10),
    // AST_TERM_LOOP_WHILE           = term_idx(12),
    // AST_TERM_LOOP_DO_WHILE        = term_idx(13),
		// AST_TERM_LOOP_FOR             = term_idx(14),

		// TODO: AST_TERM_TUPLE
		// TODO: AST_TERM_UNION
		// TODO: AST_TERM_DEPENDENT_TUPLE
		// TODO: AST_TERM_DEPENDENT_UNION

    // utils
		AST_TERMS_SEQUENCE             = util_idx(1),
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

ASTNodeKey ASTIntegerLiteralNode(ASTManager *, const char* int_sym);

ASTNodeKey ASTBinaryOperationNode(ASTManager *, ASTNodeKey op, ASTNodeKey left, ASTNodeKey right);

ASTNodeKey ASTLambdaDeclarationNode(ASTManager*, ASTNodeKey arg, ASTNodeKey body);

ASTNodeKey ASTFunctionDeclarationNode(ASTManager*, ASTNodeKey funcname, ASTNodeKey arguments, ASTNodeKey return_type, ASTNodeKey body);

ASTNodeKey ASTVariableDeclarationNode(ASTManager*, ASTNodeKey varname, ASTNodeKey value, ASTNodeKey type_anotation);

ASTNodeKey ASTSequencialTermsNode(ASTManager*, ASTNodeKey val, ASTNodeKey seq);

ASTNodeKey ASTIfNode(ASTManager*, ASTNodeKey cond, ASTNodeKey body, ASTNodeKey elif);
// ASTNodeKey ASTDoWhileLoopNode(ASTManager*, ASTNodeKey cond, ASTNodeKey body, ASTNodeKey elif);
// ASTNodeKey ASTWhileLoopNode(ASTManager*, ASTNodeKey cond, ASTNodeKey body, ASTNodeKey elif);
// ASTNodeKey ASTForLoopNode(ASTManager*, ASTNodeKey cond, ASTNodeKey body, ASTNodeKey elif);

#endif
