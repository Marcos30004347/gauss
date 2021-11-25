#include "Core/Algebra/List.hpp"
#include "Core/Algebra/Algebra.hpp"

#include <assert.h>

using namespace ast;
using namespace algebra;

void should_create_lists() {
  Expr L = list({0, 1, 2});
  assert(L == list({0, 1, 2}));
}

void should_remove_elements_from_list() {
  Expr L = list({0, 1, 2});
  Expr M = list({0, 2});

  assert(remove(L, M) == list({1}));

  Expr A = list({0, 0, 1, 2});
  Expr B = list({0, 2});

  assert(remove(A, B) == list({0, 1}));
}

void should_join_elements_from_list() {
  Expr L = list({0, 1, 2});
  Expr M = list({3, 4});

  assert(join(L, M) == list({0, 1, 2, 3, 4}));
}

void should_adjoin_elements_from_list() {
  Expr L = list({1, 2});
  Expr K = adjoin(0, L);

  assert(K == list({0, 1, 2}));
}

void should_get_rest_elements_from_list() {
  Expr L = list({0, 1, 2});
  assert(rest(L) == list({1, 2}));
}

void should_get_first_elements_from_list() {
  Expr L = list({0, 1, 2});
  Expr f = first(L);

  assert(f.kind() == Kind::Integer);
  assert(f.value() == 0);
}

void should_compare_lists() {
  Expr L = list({0, 1, 2});
  Expr A = list({0, 1, 2});
  Expr B = list({2, 1, 0});

  assert(L == A);
  assert(L != B);
}

int main() {
  should_create_lists();
  should_remove_elements_from_list();
  should_join_elements_from_list();
  should_adjoin_elements_from_list();
  should_get_rest_elements_from_list();
  should_get_first_elements_from_list();
  should_compare_lists();
  return 0;
}
