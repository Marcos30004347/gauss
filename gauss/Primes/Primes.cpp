#include "Primes.hpp"
#include <cstddef>
#include <limits>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

Primes primes = Primes();

Primes::Primes()
{
	lp = std::vector<unsigned long long>();

	primes = std::vector<unsigned long long>();

	cacheMorePrimes();
}

Primes::~Primes() {}

int Primes::operator[](unsigned int idx) {

	while(idx >= this->primes.size()) {
		cacheMorePrimes();
	}

	return this->primes[idx];
}

unsigned int Primes::count() { return this->primes.size(); }

void Primes::cacheMorePrimes() {
	unsigned long long n = lp.size();

	if(lp.size() > 2147483647) {
		printf("trying to store to many primes!\n");
		abort();
	}

	for(size_t i = 0; i < 5000; i++) {
		this->lp.push_back(0);
	}

	if(primes.size()) {
		for (unsigned long long i = primes[primes.size() - 1]; i < n; ++i) {
			for (
					 size_t j=0;
					 (j < primes.size()) &&
						 (primes[j] <= lp[i]) &&
						 (i * primes[j] <= lp.size());
					 ++j) {

				lp[i * primes[j]] = primes[j];
			}
		}
	}

	for (unsigned long long i= (n == 0 ? 2 : n); i < lp.size(); ++i) {
		if (lp[i] == 0) {
			lp[i] = i;

			primes.push_back(i);
		}

		for (size_t j=0; j < primes.size() && primes[j]<=lp[i] && i*primes[j] < lp.size(); ++j)
			lp[i * primes[j]] = primes[j];
	}
}

std::vector<unsigned long long> Primes::factorsOf(unsigned long long i) {
	unsigned long long k = i;

	std::vector<unsigned long long> f;

	while (k != 1) {
		f.push_back(lp[k]);

		while(k % f.back() == 0) {
			k = k / f.back();
		}
	}

	return f;
}
