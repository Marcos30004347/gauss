#include "TypeSystem/Types.hpp"

#include <cassert>
#include <cstddef>
#include <cstdio>

#include <string>


int main() {
	storage* s = storage_create();

	for(size_t j = 0; j < 256; j++) {
		term_ref x = term_symbol(s, root_context(s), ("x" + std::to_string(j)).c_str());
		term_ref y = term_symbol(s, root_context(s), ("y" + std::to_string(j)).c_str());

		term_ref childs[2] = {x, y};

		term_ref c = storage_insert(s, root_context(s), term::TYPE_ARR, childs, 2);

		term &C = storage_get_term_from_key(s, c);

		assert(C.term_kind == term::TYPE_ARR);
		printf("%lu: ", 3*j);
		print_term(s, c);
		printf("\n");
	}

	storage_destroy(s);

	return 0;
}
