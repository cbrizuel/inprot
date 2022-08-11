//
// Created by germelcar on 10/3/17.
//

#include <iostream>
#include <emmintrin.h>

int main () {
	volatile __m128d a, b;
	double vals[2] = {0};

	a = _mm_loadu_pd (vals);
	b = _mm_add_pd (a,a);
	_mm_storeu_pd (vals,b);

	std::cout << "SSE2";
	return 0;
}

