//
// Created by germelcar on 10/3/17.
//

#include <iostream>
#include <immintrin.h>

int main() {

	long long a[2] = {  1, 2 };
	long long b[2] = { -1, 3 };
	long long c[2];
	__m128i va = _mm_loadu_si128 ((__m128i*)a);
	__m128i vb = _mm_loadu_si128 ((__m128i*)b);
	__m128i vc = _mm_cmpgt_epi64 (va, vb);
	_mm_storeu_si128 ((__m128i*)c, vc);

	std::cout << "SSE42";
	return 0;
}

