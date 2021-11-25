#include "Division.hpp"

using namespace ast;
using namespace algebra;

namespace expand {

Expr expandDivision(Expr n, Expr d) {
	if(n.size() > 1) {
		Expr e = Expr(n.kind());
		for(unsigned int i=0; i<n.size(); i++)
			e.insert(
				div(
					n[i],
					d
				)
			);			
		return e;
	}
	return div(n, d);
}

Expr expandDivisionAST(Expr n) {
	return expandDivision(n[0], n[1]);
}

}
