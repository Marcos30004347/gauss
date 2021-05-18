#!/bin/bash

if [ -f ./__primes.hpp ]; then rm -rf ./__primes.hpp; fi

echo "#ifndef __PRIMES_FILE_H" >> __primes.hpp
echo "#define __PRIMES_FILE_H" >> __primes.hpp
echo "#include <stdio.h>" >> __primes.hpp
echo "#include <stdlib.h>" >> __primes.hpp

a=0

code="unsigned long get_nth_prime(unsigned long i) {"

for filename in `ls *.hpp | sort -V`; do
    up_to_prime=$(echo "$filename" | sed -nr 's/.*primes_(.*).hpp.*/\1/p')
    
    # Using 500 includes, what let us to include 500*10000 = 5.000.000 primes
    # on our list
    if [ $a -gt 500 ]; then
        break
    fi

    a=$((a+1))

    echo "#include \"$filename\"" >> __primes.hpp
    
    cond="i >= $up_to_prime && i < $up_to_prime + 10000"
    idx="i - $up_to_prime"

    if [ "$up_to_prime" = "0" ]; then
        cond="i < 10000"
        idx="i"
    fi

    code="${code} if($cond) { return primes_$up_to_prime[$idx]; }"
done

code="${code}; printf(\"Error\"); abort(); return 0; }"

echo $code >> __primes.hpp

echo "#endif" >> __primes.hpp
