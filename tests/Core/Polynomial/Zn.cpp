#include "Core/Polynomial/Zp.hpp"
#include "Core/Polynomial/Polynomial.hpp"
#include "Core/Expand/Expand.hpp"
#include "Core/Simplification/Simplification.hpp"
#include "Core/Simplification/Addition.hpp"

#include <assert.h>

using namespace ast;
using namespace expand;
using namespace algebra;
using namespace simplification;
using namespace polynomial;

void should_project_u_in_Z() {
	AST* x = symbol("x");

	AST* u = add({
		power(symbol("x"), integer(2)),
		mul({integer(5), symbol("x")}),
		integer(8)
	});

	AST* u_Tnn5 = Tnn(u, x, 5);
	AST* u_Ts5 = Ts(u, x, 5);

	assert(u_Tnn5->kind() == Kind::Addition);
	assert(u_Tnn5->operand(0)->kind() == Kind::Integer);
	assert(u_Tnn5->operand(0)->value() == 3);
	assert(u_Tnn5->operand(1)->kind() == Kind::Power);
	assert(u_Tnn5->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(u_Tnn5->operand(1)->operand(0)->identifier() == "x");
	assert(u_Tnn5->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(u_Tnn5->operand(1)->operand(1)->value() == 2);

	assert(u_Ts5->kind() == Kind::Addition);
	assert(u_Ts5->operand(0)->kind() == Kind::Integer);
	assert(u_Ts5->operand(0)->value() == -2);
	assert(u_Ts5->operand(1)->kind() == Kind::Power);
	assert(u_Ts5->operand(1)->operand(0)->kind() == Kind::Symbol);
	assert(u_Ts5->operand(1)->operand(0)->identifier() == "x");
	assert(u_Ts5->operand(1)->operand(1)->kind() == Kind::Integer);
	assert(u_Ts5->operand(1)->operand(1)->value() == 2);

	delete x;
	delete u;
	delete u_Tnn5;
	delete u_Ts5;
}

void should_get_extended_euclidean_alg_gpe_in_Zp() {
	AST* x = symbol("x");

	AST* u = add({
		power(symbol("x"), integer(5)),
		mul({integer(2), power(symbol("x"), integer(3))}),
		mul({integer(2), power(symbol("x"), integer(2))}),
		symbol("x"),
		integer(2)
	});

	AST* v = add({
		power(symbol("x"), integer(4)),
		mul({integer(2), power(symbol("x"), integer(3))}),
		mul({integer(2), power(symbol("x"), integer(2))})
	});

	AST* euc = extendedEuclideanAlgGPE_Zp(u, v, x, 3);

	assert(euc->operand(0)->kind() == Kind::Addition);
	assert(euc->operand(0)->operand(0)->kind() == Kind::Integer);
	assert(euc->operand(0)->operand(0)->value() == 2);
	assert(euc->operand(0)->operand(1)->kind() == Kind::Multiplication);
	assert(euc->operand(0)->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(euc->operand(0)->operand(1)->operand(0)->value() == 2);
	assert(euc->operand(0)->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(euc->operand(0)->operand(1)->operand(1)->identifier() == "x");
	assert(euc->operand(0)->operand(2)->kind() == Kind::Power);
	assert(euc->operand(0)->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(euc->operand(0)->operand(2)->operand(0)->identifier() == "x");
	assert(euc->operand(0)->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(euc->operand(0)->operand(2)->operand(1)->value() == 2);

	assert(euc->operand(1)->kind() == Kind::Addition);
	assert(euc->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(euc->operand(1)->operand(0)->value() == 1);
	assert(euc->operand(1)->operand(1)->kind() == Kind::Multiplication);
	assert(euc->operand(1)->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(euc->operand(1)->operand(1)->operand(0)->value() == 2);
	assert(euc->operand(1)->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(euc->operand(1)->operand(1)->operand(1)->identifier() == "x");

	assert(euc->operand(2)->kind() == Kind::Power);
	assert(euc->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(euc->operand(2)->operand(0)->identifier() == "x");
	assert(euc->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(euc->operand(2)->operand(1)->value() == 2);

	AST* k_ = add({
		mul({euc->operand(1)->deepCopy(), u->deepCopy()}),
		mul({euc->operand(2)->deepCopy(), v->deepCopy()})
	});

	AST* k = Tnn(k_, x, 3);

	assert(k->match(euc->operand(0)));

	delete x;
	delete u;
	delete v;
	delete euc;
	delete k_;
	delete k;
}

int main() {
	should_project_u_in_Z();
	// should_get_extended_euclidean_alg_gpe_in_Zp();
	return 0;
}

