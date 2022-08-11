//
// Created by germelcar on 8/31/17.
//

#ifndef INPROT_HYDROPHOBICITY_SCALE_H
#define INPROT_HYDROPHOBICITY_SCALE_H

#include <unordered_map>

template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;


namespace md {

	namespace scales {

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> CID2 {
				{'A', -0.08},
				{'L', 1.24},
				{'R', -0.09},
				{'K', -0.09},
				{'N', -0.70},
				{'M', 1.27},
				{'D', -0.71},
				{'F', 1.53},
				{'C', 0.76},
				{'P', -0.01},
				{'Q', -0.40},
				{'S', -0.93},
				{'E', -1.31},
				{'T', -0.59},
				{'G', -0.84},
				{'W', 2.25},
				{'H', 0.43},
				{'Y', 1.53},
				{'I', 1.39},
				{'V', 1.09}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> CID4 {
				{'A', 0.17},
				{'L', 0.96},
				{'R', -0.70},
				{'K', -0.62},
				{'N', -0.90},
				{'M', 0.60},
				{'D', -1.05},
				{'F', 1.29},
				{'C', 1.24},
				{'P', -0.21},
				{'Q', -1.20},
				{'S', -0.83},
				{'E', -1.19},
				{'T', -0.62},
				{'G', -0.57},
				{'W', 1.51},
				{'H', -0.25},
				{'Y', 0.66},
				{'I', 2.06},
				{'V', 1.21}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> CID5 {
				{'A', 0.02},
				{'L', 1.14},
				{'R', -0.42},
				{'K', -0.41},
				{'N', -0.77},
				{'M', 1.00},
				{'D', -1.04},
				{'F', 1.35},
				{'C', 0.77},
				{'P', -0.09},
				{'Q', -1.10},
				{'S', -0.97},
				{'E', -1.14},
				{'T', -0.77},
				{'G', -0.80},
				{'W', 1.71},
				{'H', 0.26},
				{'Y', 1.11},
				{'I', 1.81},
				{'V', 1.13}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> NormalizedEisenberg {
				{'A',  0.62},
				{'L',  1.10},
				{'R', -2.50},
				{'K', -1.50},
				{'N', -0.78},
				{'M',  0.64},
				{'D', -0.90},
				{'F',  1.20},
				{'C',  0.29},
				{'P',  0.12},
				{'Q', -0.85},
				{'S', -0.18},
				{'E', -0.74},
				{'T', -0.05},
				{'G',  0.48},
				{'W',  0.81},
				{'H', -0.40},
				{'Y',  0.26},
				{'I',  1.40},
				{'V',  1.10}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> ManavalanPonnuswamy {
				{'A', 12.97},
				{'L', 14.90},
				{'R', 11.72},
				{'K', 11.36},
				{'N', 11.42},
				{'M', 14.39},
				{'D', 10.85},
				{'F', 14.00},
				{'C', 14.63},
				{'P', 11.37},
				{'Q', 11.76},
				{'S', 11.23},
				{'E', 11.89},
				{'T', 11.69},
				{'G', 12.43},
				{'W', 13.93},
				{'H', 12.16},
				{'Y', 13.42},
				{'I', 15.67},
				{'V', 15.71}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> Ponnuswamy5 {
				{'A', 14.60},
				{'L', 16.49},
				{'R', 13.24},
				{'K', 13.28},
				{'N', 11.79},
				{'M', 16.23},
				{'D', 13.78},
				{'F', 14.18},
				{'C', 15.90},
				{'P', 14.10},
				{'Q', 12.02},
				{'S', 13.36},
				{'E', 13.59},
				{'T', 14.50},
				{'G', 14.18},
				{'W', 13.90},
				{'H', 15.35},
				{'Y', 14.76},
				{'I', 14.10},
				{'V', 16.30}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> Prabhakaran {
				{'A', -06.70},
				{'L', -11.70},
				{'R', 51.50},
				{'K', 36.80},
				{'N', 20.10},
				{'M', -14.20},
				{'D', 38.50},
				{'F', -15.50},
				{'C', -08.40},
				{'P', 00.80},
				{'Q', 17.20},
				{'S', -02.50},
				{'E', 34.30},
				{'T', -05.00},
				{'G', -04.20},
				{'W', -07.90},
				{'H', 12.60},
				{'Y', 2.90},
				{'I', -13.00},
				{'V', -10.90}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> SweetEisenberg {
				{'A', -0.40},
				{'L', 1.22},
				{'R', -0.59},
				{'K', -0.67},
				{'N', -0.92},
				{'M', 1.02},
				{'D', -1.31},
				{'F', 1.92},
				{'C', 0.17},
				{'P', -0.49},
				{'Q', -0.91},
				{'S', -0.55},
				{'E', -1.22},
				{'T', -0.28},
				{'G', -0.67},
				{'W', 0.50},
				{'H', -0.64},
				{'Y', 1.67},
				{'I', 1.25},
				{'V', 0.91}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> Zimmerman {
				{'A', 0.83},
				{'L', 2.52},
				{'R', 0.83},
				{'K', 1.60},
				{'N', 0.09},
				{'M', 1.40},
				{'D', 0.64},
				{'F', 2.75},
				{'C', 1.48},
				{'P', 2.70},
				{'Q', 0.00},
				{'S', 0.14},
				{'E', 0.65},
				{'T', 0.54},
				{'G', 0.10},
				{'W', 0.31},
				{'H', 1.10},
				{'Y', 2.97},
				{'I', 3.07},
				{'V', 1.79}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> Wolfenden {
				{'A', 1.12},
				{'L', 1.18},
				{'R', -2.55},
				{'K', -0.80},
				{'N', -0.83},
				{'M', 0.55},
				{'D', -0.83},
				{'F', 0.67},
				{'C', 0.59},
				{'P', 0.54},
				{'Q', -0.78},
				{'S', -0.05},
				{'E', -0.92},
				{'T', -0.02},
				{'G', 1.20},
				{'W', -0.19},
				{'H', -0.93},
				{'Y', -0.23},
				{'I', 1.16},
				{'V', 1.13}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> CasariSippl {
				{'A', 0.20},
				{'L', 0.50},
				{'R', -0.70},
				{'K', -1.60},
				{'N', -0.50},
				{'M', 0.50},
				{'D', -1.40},
				{'F', 1.00},
				{'C', 1.90},
				{'P', -1.00},
				{'Q', -1.10},
				{'S', -0.70},
				{'E', -1.30},
				{'T', -0.40},
				{'G', -0.10},
				{'W', 1.60},
				{'H', 0.40},
				{'Y', 0.50},
				{'I', 1.40},
				{'V', 0.70}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> Tossi {
				{'A', -1.1},
				{'L', 9.7},
				{'R', -10.0},
				{'K', -9.9},
				{'N', -7.1},
				{'M', 4.6},
				{'D', -8.3},
				{'F', 10.00},
				{'C', -2.3},
				{'P', -0.20},
				{'Q', -6.0},
				{'S', -4.3},
				{'E', -8.3},
				{'T', -3.8},
				{'G', -2.4},
				{'W', 9.7},
				{'H', -3.8},
				{'Y', 2.5},
				{'I', 8.7},
				{'V', 4.1}
		};

	}

}

#endif //INPROT_HYDROPHOBICITY_SCALE_H
