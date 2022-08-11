//
// Created by germelcar on 8/31/17.
//

#ifndef INPROT_PKA_SCALE_H
#define INPROT_PKA_SCALE_H


#include <unordered_map>

template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;


namespace md {

	namespace scales {

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> IPC = {
				{'C', 7.555},
				{'D', 3.872},
				{'E', 4.412},
				{'H', 5.637},
				{'K', 9.052},
				{'R', 11.84},
				{'Y', 10.85},
				{'#', 9.094}, // This is equivalent to: NH2
				{'@', 2.869} // This is equivalent to: COOH
		};

	}
}


#endif //INPROT_PKA_SCALE_H
