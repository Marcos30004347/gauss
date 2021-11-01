#include <assert.h>
#include "Polynomial.hpp"

#include "Core/Simplification/Simplification.hpp"
#include "Core/Simplification/Subtraction.hpp"

using namespace ast;
using namespace algebra;
using namespace simplification;

namespace expand {


void getBinomialsParameters(Int n, Int sum, Int* out, Int index, std::vector<std::vector<Int>>& res) {
	if (index > n || sum < 0)
		return;

	if (index == n) {
		if(sum == 0) {

			res.push_back(std::vector<Int>());

			for(long long i=0; i<n; i++) {
				res[res.size() - 1].push_back(out[i]);
			}
		}
		return;
	}

	for (Int i = 0; i <= 9; i++) {
		out[index.longValue()] = i;
		getBinomialsParameters(n, sum - i, out, index + 1, res);
	}
}
 
std::vector<std::vector<Int>> findNDigitNums(Int n, Int sum) {
	std::vector<std::vector<Int>> res;
    Int out[n.longValue() + 1];

    for (Int i = 0; i <= 9; i++) {
        out[0] = i;
        getBinomialsParameters(n, sum - i, out, 1, res);
    }

	return res;
}

AST* expandMultinomial(AST* b, AST* e) {
	assert(b->kind() == Kind::Addition || b->kind() == Kind::Subtraction);
	assert(e->kind() == Kind::Integer && e->value() > 0);

	AST* k;

	if(b->kind() == Kind::Subtraction) {
		k = reduceSubtractionAST(b);
	} else {
		k = b->copy();
	}

	Int m = k->numberOfOperands();
	Int n = e->value();

	AST* r = new AST(Kind::Addition);

	Int s[m.longValue() + 1] = {0};
	s[m.longValue()] = n;

	std::vector<std::vector<Int>> binomials = findNDigitNums(m, n);
	for(std::vector<Int> bi : binomials) {
		AST* p = new AST(Kind::Multiplication);

		p->includeOperand(binomial(n, bi));

		for(int i=0; i<m; i++) {
			p->includeOperand(power(k->operand(i)->copy(), integer(bi[i])));
		}
		r->includeOperand(p);
	}

	AST* res = reduceAST(r);

	delete k;
	delete r;

	return res;
}

AST* expandMultinomialAST(AST* u) {
	AST* b = base(u);
	AST* e = expoent(u);

	if(
		(b->kind() == Kind::Addition || b->kind() == Kind::Subtraction) &&
		b->numberOfOperands() > 1 && e->kind() == Kind::Integer && e->value() > 0
	) {
		AST* res = expandMultinomial(b, e);
		delete b; delete e;
		return res;
	}

	delete b; delete e;
	return u->copy();
}


}
