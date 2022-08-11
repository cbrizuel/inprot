//
// Created by germelcar on 10/3/17.
//

#include <iostream>
#include <immintrin.h>

int main() {

	volatile __m256i a, b;
	a = _mm256_set1_epi8 (1);
	b = _mm256_add_epi8 (a,a);

	std::cout << "AVX2";
	return 0;
}

