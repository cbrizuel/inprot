//
// Created by germelcar on 8/31/17.
//

#ifndef INPROT_CHARGE_SCALE_H
#define INPROT_CHARGE_SCALE_H

#include <unordered_map>

template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;


namespace md {

	namespace scales {

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> Klein {
				{'A', 0},
				{'L', 0},
				{'R', 1},
				{'K', 1},
				{'N', 0},
				{'M', 0},
				{'D', -1},
				{'F', 0},
				{'C', 0},
				{'P', 0},
				{'Q', 0},
				{'S', 0},
				{'E', -1},
				{'T', 0},
				{'G', 0},
				{'W', 0},
				{'H', 0},
				{'Y', 0},
				{'I', 0},
				{'V', 0}
		};

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::unordered_map<char, T> ChartonCTDC {
				{'D', 0},
				{'E', 0},
				{'K', 1},
				{'R', 1},
				{'H', 1},
				{'Y', 1},
				{'W', 1},
				{'F', 1},
				{'C', 1},
				{'M', 1},
				{'S', 0},
				{'T', 0},
				{'N', 1},
				{'Q', 1},
				{'G', 0},
				{'A', 0},
				{'V', 0},
				{'L', 0},
				{'I', 0},
				{'P', 0}
		};

	}

}


#endif //INPROT_CHARGE_SCALE_H
