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

void t() {
	AST* tensor2d = tensor(list({
		list({ integer(0), integer(1), integer(2) }),
		list({ integer(3), integer(4), integer(5) }),
		list({ integer(6), integer(7), integer(8) }),
	}));

	// tensor2d->operand(0)->operand(0) = 0
	// tensor2d->operand(0)->operand(1) = 1
	// tensor2d->operand(0)->operand(2) = 2
	// tensor2d->operand(1)->operand(0) = 3
	// tensor2d->operand(1)->operand(1) = 4
	// tensor2d->operand(1)->operand(2) = 5
	// tensor2d->operand(2)->operand(0) = 6
	// tensor2d->operand(2)->operand(1) = 7
	// tensor2d->operand(2)->operand(2) = 8

	AST* tensor3d = tensor(list({
		list({
			list({ integer(0), integer(1), integer(2) }),
			list({ integer(3), integer(4), integer(5) }),
			list({ integer(6), integer(7), integer(8) }),
		}),
		list({
			list({ integer(9), integer(10), integer(11) }),
			list({ integer(12), integer(13), integer(14) }),
			list({ integer(15), integer(16), integer(17) }),
		}),
		list({
			list({ integer(18), integer(19), integer(20) }),
			list({ integer(21), integer(22), integer(23) }),
			list({ integer(24), integer(25), integer(26) }),
		}),
	}));

	// tensor2d->operand(0)->operand(0)->operand(0) = 0
	// tensor2d->operand(0)->operand(0)->operand(1) = 1
	// tensor2d->operand(0)->operand(0)->operand(2) = 2
	// tensor2d->operand(0)->operand(1)->operand(0) = 3
	// tensor2d->operand(0)->operand(1)->operand(1) = 4
	// tensor2d->operand(0)->operand(1)->operand(2) = 5
	// tensor2d->operand(0)->operand(2)->operand(0) = 6
	// tensor2d->operand(0)->operand(2)->operand(1) = 7
	// tensor2d->operand(0)->operand(2)->operand(2) = 8
	// tensor2d->operand(1)->operand(0)->operand(0) = 9
	// tensor2d->operand(1)->operand(0)->operand(1) = 10
	// tensor2d->operand(1)->operand(0)->operand(2) = 11
	// tensor2d->operand(1)->operand(1)->operand(0) = 12
	// tensor2d->operand(1)->operand(1)->operand(1) = 13
	// tensor2d->operand(1)->operand(1)->operand(2) = 14
	// tensor2d->operand(1)->operand(2)->operand(0) = 15
	// tensor2d->operand(1)->operand(2)->operand(1) = 16
	// tensor2d->operand(1)->operand(2)->operand(2) = 17
	// tensor2d->operand(2)->operand(0)->operand(0) = 18
	// tensor2d->operand(2)->operand(0)->operand(1) = 19
	// tensor2d->operand(2)->operand(0)->operand(2) = 20
	// tensor2d->operand(2)->operand(1)->operand(0) = 21
	// tensor2d->operand(2)->operand(1)->operand(1) = 22
	// tensor2d->operand(2)->operand(1)->operand(2) = 23
	// tensor2d->operand(2)->operand(2)->operand(0) = 24
	// tensor2d->operand(2)->operand(2)->operand(1) = 25
	// tensor2d->operand(2)->operand(2)->operand(2) = 26

}

}
