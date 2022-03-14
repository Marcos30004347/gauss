#ifndef PRIMES_H
#define PRIMES_H

#include <vector>

class Primes {
private:
	// array of primes
	std::vector<int> primes;

	// array of least prime factor
	int* lp;

	// caches all prime numbers up to N
	void cachePrimesUpTo(int N);
public:

	Primes(unsigned int N);

	~Primes();

	unsigned int count();

	std::vector<int> factorsOf(int i);

	int operator[](unsigned int idx);
};

extern Primes primes;

#endif
