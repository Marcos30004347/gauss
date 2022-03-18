#ifndef PRIMES_H
#define PRIMES_H

#include <vector>

class Primes {
private:
	// array of primes
	std::vector<unsigned long long> primes;

	// array of least prime factor
	std::vector<unsigned long long> lp;

	void cacheMorePrimes();

public:

	Primes();

	~Primes();

	unsigned int count();

	std::vector<long long> factorsOf(long long i);

	int operator[](unsigned int idx);
};

extern Primes primes;

#endif
