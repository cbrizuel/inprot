//
// Created by germelcar on 8/31/17.
//

#ifndef INPROT_HYDROPHILICITY_SCALE_H
#define INPROT_HYDROPHILICITY_SCALE_H

#include <unordered_map>

template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;


namespace md {

	namespace scales {

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> KuhnHydrov {
				{'A', 0.78},
				{'L', 0.56},
				{'R', 1.58},
				{'K', 1.10},
				{'N', 1.20},
				{'M', 0.66},
				{'D', 1.35},
				{'F', 0.47},
				{'C', 0.55},
				{'P', 0.69},
				{'Q', 1.19},
				{'S', 1},
				{'E', 1.45},
				{'T', 1.05},
				{'G', 0.68},
				{'W', 0.70},
				{'H', 0.99},
				{'Y', 1},
				{'I', 0.47},
				{'V', 0.51}

		};

	}

}

#endif //INPROT_HYDROPHILICITY_SCALE_H
