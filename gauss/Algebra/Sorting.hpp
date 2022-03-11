#include "Expression.hpp"

namespace alg {

int compare(expr *const a, expr *const b, kind ctx);

void sort(expr *a, kind k);

void sort_childs(expr *a, kind k, long int l, long int r);

void sort_childs(expr *a, long int l, long int r);

void sort(expr *a);


}
