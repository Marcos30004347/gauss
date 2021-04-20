#include "core/collection.hpp"

#include <cstdlib>
#include <cstdio>

using namespace core;
using namespace containers;

int* new_int(int i) {
    return new int(i);
}

void should_store_elements() {
    collection<int*> coll0;

    coll0.insert(new int(3));
    coll0.insert(new int(4));
    coll0.insert(new int(5));
    coll0.insert(new int(6));

    coll0.have(new int(3));
    coll0.have(new_int(5));

    int* i = new int(3);
    coll0.have(i);

    collection<int*> coll1 = { new_int(1), new_int(2), new_int(3) };
}

int main() {
    should_store_elements();
    return 0;
}