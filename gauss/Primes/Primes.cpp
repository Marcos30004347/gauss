#include "Primes.hpp"

Primes primes = Primes(50000000);

Primes::Primes(unsigned int N)
{
	this->lp = new int[N+1]{0};
	this->cachePrimesUpTo(N);
}

Primes::~Primes() {
	delete []this->lp;
}

int Primes::operator[](unsigned int idx) {
	return this->primes[idx];
}

unsigned int Primes::count() {
	return this->primes.size();
}

void Primes::cachePrimesUpTo(int N) {
	this->primes.clear();

	for (int i=2; i<=N; ++i) {
		if (lp[i] == 0) {
			lp[i] = i;
			primes.push_back(i);
		}

		for (int j=0; j<(int)primes.size() && primes[j]<=lp[i] && i*primes[j]<=N; ++j)
			lp[i * primes[j]] = primes[j];
	}

}

std::vector<int> Primes::factorsOf(int i) {
	int k = i;

	std::vector<int> f;

	while (k != 1) {
		f.push_back(lp[k]);

		while(k % f.back() == 0) {
			k = k/f.back();
		}
	}

	return f;
}
