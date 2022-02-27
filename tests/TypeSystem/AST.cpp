#include "TypeSystem/AST.hpp"

#include <cassert>
#include <cstddef>
#include <cstdio>

#include <string>



int main() {
	ASTManager* manager = ASTManagerCreate();

	for(size_t j = 0; j < 256; j++) {
		ASTNodeKey x = ASTSymbolTypeNode(manager, ("x" + std::to_string(j)).c_str());

		ASTNodeKey y = ASTSymbolTypeNode(manager, ("y" + std::to_string(j)).c_str());

		ASTNodeKey c = ASTArrowTypeNode(manager, x, y);

		ASTNode &arr_type = ASTGetNode(manager, c);

		assert(arr_type.node_kind == ASTNode::AST_TYPE_ARROW);

		printf("%lu: ", 3*j);

		ASTPrint(manager, c);

		printf("\n");
	}

	ASTManagerDestroy(manager);

	return 0;
}
