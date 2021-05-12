#include "Trigonometry.hpp"
#include "Core/Algebra/List.hpp"

using namespace ast;
using namespace algebra;

namespace simplification {

AST* substituteTrig(AST* u) {
	if(
		u->kind() == Kind::Integer ||
		u->kind() == Kind::Fraction ||
		u->kind() == Kind::Symbol
	) return u->deepCopy();

	AST* g = mapUnaryAST(u, substituteTrig);

	if(g->kind() == Kind::FunctionCall) {
		
		if(g->funName() == "tan") {
			AST* k = div(
				funCall("sin", { g->operand(0)->deepCopy() }),
				funCall("cos", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "cot") {
			AST* k = div(
				funCall("cos", { g->operand(0)->deepCopy() }),
				funCall("sin", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "sec") {
			AST* k = div(
				integer(1),
				funCall("cos", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}

		if(g->funName() == "csc") {
			AST* k = div(
				integer(1),
				funCall("sin", { g->operand(0)->deepCopy() })
			);
			delete g;
			return k;
		}
	}

	return g;
}



}
