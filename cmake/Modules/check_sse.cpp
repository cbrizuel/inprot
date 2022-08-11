//
// Created by germelcar on 10/3/17.
//

#include <iostream>
#include <xmmintrin.h>

int main () {
	volatile __m128 a, b;
	float vals[4] = {0};

	a = _mm_loadu_ps (vals);
	b = a;
	b = _mm_add_ps (a,b);
	_mm_storeu_ps (vals,b);

	std::cout << "SSE";
	return 0;
}

