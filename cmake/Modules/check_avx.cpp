//
// Created by germelcar on 10/3/17.
//

#include <iostream>
#include <immintrin.h>

int main() {
	__m256 a, b;
	float vals[8] = {1, 2, 3, 4, 5, 6, 7, 8};
	const int mask = 123;

	a = _mm256_loadu_ps(vals);
	b = a;
	b = _mm256_dp_ps (a, a, 3);
	_mm256_storeu_ps(vals,b);

	std::cout << "AVX";
	return 0;
}

