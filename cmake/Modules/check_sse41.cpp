//
// Created by germelcar on 10/3/17.
//

#include <iostream>
#include <smmintrin.h>

int main() {

	volatile __m128 a, b;
	float vals[4] = {1, 2, 3, 4};
	const int mask = 123;
	a = _mm_loadu_ps (vals);
	b = a;
	b = _mm_dp_ps (a, a, 4);
	_mm_storeu_ps (vals,b);

	std::cout << "SSE41";
	return 0;
}

