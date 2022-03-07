#include "TypeSystem/AST.hpp"
#include <cstddef>
#include <stack>

typedef unsigned u32;
typedef unsigned long long u64;

#define MANAGER_BUCKET_SIZE_LOG2 7

#define MANAGER_BUCKET_SIZE ((u64)1 << MANAGER_BUCKET_SIZE_LOG2)

#define MANAGER_INDEX_MASK (((u64)1 << MANAGER_BUCKET_SIZE_LOG2) - 1)

#define MANAGER_BUCKET_MASK ~MANAGER_INDEX_MASK

#define build_key(bucket, index)                                               \
  ((u64)index | ((u64)bucket << MANAGER_BUCKET_SIZE_LOG2))

#define key_bucket(k) (((u64)k & MANAGER_BUCKET_MASK) >> MANAGER_BUCKET_SIZE_LOG2)
#define key_index(k) ((u64)k & MANAGER_INDEX_MASK)

template<typename T, unsigned size>
struct Bucket {
	// data array
	T data[size];
};

template <typename T>
struct BucketSlot {
	// Actual data of this slot
	T data;

	// the total count of 'data' u64 used.
	u64 length;

	// The number of references to this node.
	u64 ref_count;

	// reference to the 'data' reserver array.
	u64 data_ref;
};

template<typename T>
struct Manager {
	unsigned buckets_count;

	// stack with free slots
	std::stack<u64> slots;

	// head of the bucket list
	Bucket<BucketSlot<T>, MANAGER_BUCKET_SIZE>* buckets;

	// array of members;
	Bucket<u64, MANAGER_BUCKET_SIZE>* members;
};


template <typename T>
u64 ManagerReserve(Manager<T>* manager) {
	u64 ref = manager->slots.pop();

	if(manager->slots.size() == 0) {
		Bucket<BucketSlot<T>, MANAGER_BUCKET_SIZE>* old = manager->buckets;

		manager->buckets = new Bucket<BucketSlot<T>, MANAGER_BUCKET_SIZE>[manager->buckets_count + 1];

		for(int i = 0; i < manager->buckets_count; i++) {
			manager->buckets[i] = old[i];
		}

		delete[] old;

		manager->buckets_count += 1;
	}

	Bucket<BucketSlot<T>, MANAGER_BUCKET_SIZE>* b = manager->head;

	u64 buc = key_bucket(ref);
	u64 idx = key_index(ref);

	while (buc) {
		if(b->next == nullptr) {
			// TODO: error
		}

		buc -= 1;

		b = b->next;
	}

	b->data[idx].~T();

	return ref;
}

template <typename T>
void ManagerFree(Manager<T>* manager, u64 key) {
	manager->slots.push(key);
}
