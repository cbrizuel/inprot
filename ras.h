//
// Created by germelcar on 9/1/17.
//

#ifndef INPROT_REDUCED_ALPHABETS_H
#define INPROT_REDUCED_ALPHABETS_H

#include <cstdint>
#include <unordered_map>


namespace md {

	namespace ra {


		//
		// Enumerations - Alphabets
		//
		enum class Std: uint8_t { A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y };

		enum class Blosum50: uint8_t { FWY, CLVIM, H, AG, ST, DENQ, KR, P };

		enum class HydrophobicityTomii: uint8_t { RKEDQN, GASTPHY, CLVIMFW };

		enum class NormVWTomii: uint8_t { GASTCPD, NVEQIL, MHKFRYW };

		enum class PolarityTomii: uint8_t { LIFWCMVY, PATGS, HQRKNED };

		enum class PolarizabilityTomii: uint8_t { GASDT, CPNVEQIL, KMHFRYW };

		enum class ChargeTomii: uint8_t { KR, ANCQGHILMFPSTWYV, DE };

		enum class SecondStructTomii: uint8_t { EALMQKRH, VIYCWFT, GNPSD };

		enum class SolventAccTomii: uint8_t { ALFCGIVW, RKQEND, MPSTHY };


	}

}

#endif //INPROT_REDUCED_ALPHABETS_H
