#include "Tensor.hpp"
#include "Core/Algebra/List.hpp"

using namespace ast;
using namespace algebra;

namespace tensor {

int count_list_rank(AST* r) {
	if(r->kind() != Kind::List) return 0;
	return 1 + count_list_rank(r->operand(0));
}

int assert_tensor_data(AST* r) {
	if(r->kind() != Kind::List) return 0;
	for(unsigned int i=0; i<r->numberOfOperands(); i++) {

	}
	return 1 + count_list_rank(r->operand(0));
}

int kroneckerDelta(unsigned int i, unsigned int j) {
	return i == j;
}

AST* tensor(AST* v) {
	int rank = count_list_rank(v);
	return new AST(Kind::Tensor, { integer(rank), v });
}

}
