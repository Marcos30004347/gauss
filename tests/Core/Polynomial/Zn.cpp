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
	AST* x = symb("x");

	AST* u = add({
		pow(symb("x"), inte(2)),
		mul({inte(5), symb("x")}),
		inte(8)
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
	AST* x = symb("x");

	AST* u = add({
		pow(symb("x"), inte(5)),
		mul({inte(2), pow(symb("x"), inte(3))}),
		mul({inte(2), pow(symb("x"), inte(2))}),
		symb("x"),
		inte(2)
	});

	AST* v = add({
		pow(symb("x"), inte(4)),
		mul({inte(2), pow(symb("x"), inte(3))}),
		mul({inte(2), pow(symb("x"), inte(2))})
	});

	std::vector<AST*> euc = extendedEuclideanAlgGPE_Zp(u, v, x, 3);

	assert(euc[0]->kind() == Kind::Addition);
	assert(euc[0]->operand(0)->kind() == Kind::Integer);
	assert(euc[0]->operand(0)->value() == 2);
	assert(euc[0]->operand(1)->kind() == Kind::Multiplication);
	assert(euc[0]->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(euc[0]->operand(1)->operand(0)->value() == 2);
	assert(euc[0]->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(euc[0]->operand(1)->operand(1)->identifier() == "x");
	assert(euc[0]->operand(2)->kind() == Kind::Power);
	assert(euc[0]->operand(2)->operand(0)->kind() == Kind::Symbol);
	assert(euc[0]->operand(2)->operand(0)->identifier() == "x");
	assert(euc[0]->operand(2)->operand(1)->kind() == Kind::Integer);
	assert(euc[0]->operand(2)->operand(1)->value() == 2);

	assert(euc[1]->kind() == Kind::Addition);
	assert(euc[1]->operand(0)->kind() == Kind::Integer);
	assert(euc[1]->operand(0)->value() == 1);
	assert(euc[1]->operand(1)->kind() == Kind::Multiplication);
	assert(euc[1]->operand(1)->operand(0)->kind() == Kind::Integer);
	assert(euc[1]->operand(1)->operand(0)->value() == 2);
	assert(euc[1]->operand(1)->operand(1)->kind() == Kind::Symbol);
	assert(euc[1]->operand(1)->operand(1)->identifier() == "x");

	assert(euc[2]->kind() == Kind::Power);
	assert(euc[2]->operand(0)->kind() == Kind::Symbol);
	assert(euc[2]->operand(0)->identifier() == "x");
	assert(euc[2]->operand(1)->kind() == Kind::Integer);
	assert(euc[2]->operand(1)->value() == 2);

	AST* k_ = add({
		mul({euc[1]->deepCopy(), u->deepCopy()}),
		mul({euc[2]->deepCopy(), v->deepCopy()})
	});

	AST* k = Tnn(k_, x, 3);

	assert(k->match(euc[0]));

	delete x;
	delete u;
	delete v;
	delete euc[0];
	delete euc[1];
	delete euc[2];
	delete k_;
	delete k;
}

int main() {
	should_project_u_in_Z();
	should_get_extended_euclidean_alg_gpe_in_Zp();
	return 0;
}

