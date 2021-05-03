#ifndef PRIMES_H
#define PRIMES_H

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <filesystem>


unsigned long nth_prime(unsigned long idx) {
	if(idx <= 10000000) {
		std::ifstream in;

		in.open("../Assets/primes_0.txt");

		std::string s;

		s.reserve(120);    

		for(int i = 0; i < idx; ++i) {
				std::getline(in, s);
		}

		std::getline(in,s);
		
		unsigned long p = strtoul(s.c_str(), nullptr, 0);
	
		return p; 
	}

	return 0;
}

#endif
