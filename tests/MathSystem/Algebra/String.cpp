#include <cassert>
#include <stdio.h>

#include "MathSystem/Algebra/String.hpp"
int main() {
  string t = "abacaxi";

  assert(t[0] == 'a');
  assert(t[1] == 'b');
  assert(t[2] == 'a');
  assert(t[3] == 'c');
  assert(t[4] == 'a');
  assert(t[5] == 'x');
  assert(t[6] == 'i');
  assert(t[7] == '\0');

  string h = t + "aba";

  assert(t == "abacaxi");

  assert(h[0] == 'a');
  assert(h[1] == 'b');
  assert(h[2] == 'a');
  assert(h[3] == 'c');
  assert(h[4] == 'a');
  assert(h[5] == 'x');
  assert(h[6] == 'i');
  assert(h[7] == 'a');
  assert(h[8] == 'b');
  assert(h[9] == 'a');
  assert(h[10] == '\0');

  assert(string::to_string(123) == "123");
}
