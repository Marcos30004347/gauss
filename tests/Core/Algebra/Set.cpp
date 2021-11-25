#include "Core/Algebra/Set.hpp"
#include "Core/Algebra/Algebra.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;

void should_create_sets() {
  Expr S = set({0, 1, 2});

  assert(S.kind() == Kind::Set);
  assert(S[0].kind() == Kind::Integer);
  assert(S[0].value() == 0);
  assert(S[1].kind() == Kind::Integer);
  assert(S[1].value() == 1);
  assert(S[2].kind() == Kind::Integer);
  assert(S[2].value() == 2);
}

void should_not_include_duplicates_on_set() {
  Expr S = set({0, 1, 2});
  S.insert(0);
  assert(S.size() == 3);
}

void should_get_set_difference() {
  Expr S = set({0, 1, 2});
  Expr A = set({0, 2});

  assert(difference(S, A) == set({1}));
}

void should_get_set_unnification() {
  Expr S = set({0, 1, 2});
  Expr A = set({0, 3});

  Expr M = unification(S, A);

  assert(M == set({0, 1, 2, 3}));
}

void should_get_set_intersection() {
  Expr S = set({0, 1, 2});
  Expr A = set({0, 3});
  Expr M = intersection(S, A);

  assert(M == set({0}));
}

void should_get_if_element_exists_inside_set() {
  Expr S = set({ 0, 1, 2});

  assert(exists(S, 1) == true);
}

void should_get_combinations_of_elements() {
  Expr S = set({0, 1, 3});

  Expr M = combination(S, 2);

  assert(M.kind() == Kind::Set);
  assert(M.size() == 3);

  assert(M[0] == set({0, 1}));
  assert(M[1] == set({0, 3}));
  assert(M[2] == set({1, 3}));
}

int main() {
  should_create_sets();
	should_not_include_duplicates_on_set();
	should_get_set_difference();
	should_get_set_unnification();
	should_get_set_intersection();
	should_get_if_element_exists_inside_set();
	should_get_combinations_of_elements();
  return 0;
}
