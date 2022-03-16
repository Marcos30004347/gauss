#include "gauss/Polynomial/Roots.hpp"
#include "gauss/Algebra/Reduction.hpp"
#include "gauss/Algebra/Expression.hpp"
#include "test.hpp"

#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>

using namespace alg;
using namespace poly;

void asssert_is_small(expr t) {
	assert(is(&t, kind::MUL | kind::INT | kind::FRAC));

	if(is(&t, kind::MUL)) {
		assert(is(&t[0],  kind::INT | kind::FRAC));

		if(is(&t[0], kind::FRAC)) {
			assert(numerator(t[0]).value() / denominator(t[0]).value() == 0);
		}

		if(is(&t[0], kind::INT)) {
			assert(t[0].value() == 0);
		}

	} else {

		if(is(&t, kind::FRAC)) {
			assert(numerator(t).value() / denominator(t).value() == 0);
		}

		if(is(&t, kind::INT)) {
			assert(t.value() == 0);
		}
	}
}

int main() {
	expr x = symbol("x");

	expr p = 5*pow(x, 4) + 3*pow(x, 3) + x + -3;

	expr R = realPolyRoots(p);

	for(size_t i = 0; i < R.size(); i++) {
		expr k = replace(p, x, R[i]);

		expr t = expand(k);

		if(t.kind() == kind::ADD) {
			for(size_t i = 0; i < t.size(); i++) {
				asssert_is_small(t[i]);
			}
		} else {
			asssert_is_small(t);
		}
	}


	return 0;
}
