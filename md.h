//
// Created by germelcar on 8/31/17.
//

#ifndef INPROT_MD_H
#define INPROT_MD_H


#include <type_traits>
#include <string>
#include <cmath>
#include <iostream>
#include <valarray>
#include "pka_scale.h"
#include "hydrophobicity_scale.h"
#include "hydrophilicity_scale.h"
#include "charge_scale.h"
#include "reduced_alphabets.h"
#include "globals.h"

template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;

namespace md {

	enum MDS: uint8_t {
		LENGTH = 0,                                 // MD # 01
		COMP_STD_F,                                 // MD # 02
		NET_CHRG,                                   // MD # 03
		DIST_NORMVW_TOMII_MHKFRYW_50,               // MD # 04
		COMP_STD_M,                                 // MD # 05
		DIST_NORMVW_TOMII_NVEQIL_75,                // MD # 06
		COMP_STD_Q,                                 // MD # 07
		HMM_EISENBERG,                              // MD # 08
		DIST_POLAR_TOMII_PATGS_0,                   // MD # 09
		DIST_POLAR_TOMII_LIFWCMVY_0,                // MD # 10
		AVG_CHRG_KLEP810101,                        // MD # 11
		DIST_POLAR_TOMII_HQRKNED_25,                // MD # 12
		NET_CHRG_CHAM83108,                         // MD # 13
		AVG_HYDRO_KUHL950101,                       // MD # 14
		AVG_HYDRO_CIDH_CIDH920102,                  // MD # 15
		AVG_HYDRO_CIDH_CIDH920104,                  // MD # 16
		AVG_HYDRO_CIDH_CIDH920105,                  // MD # 17
		AVG_HYDRO_CIDH_MANP780101,                  // MD # 18
		AVG_HYDRO_CIDH_PONP800105,                  // MD # 19
		DIST_POLAR_TOMII_GASDT_75,                  // MD # 20
		AVG_HYDRO_PRAM900101,                       // MD # 21
		AVG_HYDRO_SWER830101,                       // MD # 22
		DIST_POLAR_TOMII_GASDT_100,                 // MD # 23
		AVG_HYDRO_ZIMJ680101,                       // MD # 24
		AVG_HYDRO_WOLR790101,                       // MD # 25
		AVG_HYDRO_CASG920101,                       // MD # 26
		AVG_HYDRO_TOSSI2002,                        // MD # 27
		DIST_SS_TOMII_EALMQKRH_50,                  // MD # 28
		DIST_SS_TOMII_EALMQKRH_100,                 // MD # 29
		COMP_BLOSUM50_CLVIM,                        // MD # 30
		DIST_CHRG_TOMII_DE_0,                       // MD # 31
		COMP_BLOSUM50_FWY,                          // MD # 32
		DIST_CHRG_TOMII_KR_100,                     // MD # 33
		COMP_NORMVW_TOMII_MHKFRYW,                  // MD # 34
		DIST_SOLVENT_TOMII_ALFCGIVW_0,              // MD # 35
		DIST_SOLVENT_TOMII_MPSTHY_0,                // MD # 36
		DIST_SOLVENT_TOMII_RKQEND_0,                // MD # 37
		DIST_SOLVENT_TOMII_RKQEND_25,               // MD # 38
		COMP_POLAR_TOMII_KMHFRYW,                   // MD # 39
		COMP_CHRG_TOMII_DE,                         // MD # 40
		COMP_CHRG_TOMII_KR,                         // MD # 41
		COMP_SS_TOMII_VIYCWFT,                      // MD # 42
		COMP_SA_TOMII_ALFCGIVW,                     // MD # 43
		TRANS_HYDRO_TOMII_CLVIM_RKEDQN,             // MD # 44
		TRIP_HYDRO_TOMII_RKEDQN_CLVIMFW_GASTPHY,    // MD # 45
		TRIP_HYDRO_TOMII_CLVIMFW_CLVIMFW_GASTPHY,   // MD # 46
		TRANS_SA_TOMII_ALFCGIVW_RKQEND,             // MD # 47
		DIST_HYDRO_TOMII_GASTPHY_0,                 // MD # 48
		DIST_HYDRO_TOMII_CLVIMFW_0,                 // MD # 49
		TRIP_HYDRO_TOMII_CLVIMFW_CLVIMFW_CLVIMFW,   // MD # 50
		DIST_HYDRO_TOMII_GASTPHY_75                 // MD # 51

	};


	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr T PI {3.14159265358979323846};

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr T toRadians(T angdeg) { return angdeg / static_cast<T>(180) * PI<T>; };

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr T netCharge(const std::string& seq, int ph = Globals::PH_NET_CHARGE,
	                      const std::unordered_map<char, T>& pKMap = scales::IPC<T>) {
		std::unordered_map<char, T> nj {
				{'Y', 0},
				{'D', 0},
				{'E', 0},
				{'C', 0},
				{'@', 1} // This is equivalent to: COOH (check IPC map)
		};

		std::unordered_map<char, T> ni {
				{'R', 0},
				{'H', 0},
				{'K', 0},
				{'#', 1} // This is equivalent to: NH2 (check IP map)
		};

		std::for_each(seq.cbegin(),seq.cend(), [&] (const auto& c) {
			auto found_nj = nj.find(c);
			auto found_ni = ni.find(c);

			if (found_nj != nj.end())
				found_nj->second++;
			else if (found_ni != ni.end())
				found_ni->second++;
		});

		T pos {0};
		std::for_each(ni.begin(), ni.end(), [&] (auto& pair) {
			auto pKai = pKMap.find(pair.first);
			if (pair.second != 0 && pKai != pKMap.end())
				pos += pair.second * (1 / (1 + std::pow(10, ph - pKai->second)));
		});

		T neg {0};
		std::for_each(nj.begin(), nj.end(), [&] (auto& pair) {
			auto pKai = pKMap.find(pair.first);
			if (pair.second != 0 && pKai != pKMap.end())
				neg += pair.second * (-1 / (1 + std::pow(10, pKai->second - ph)));
		});

		return pos + neg;
	}

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr T hMoment(const std::string& seq, uint angle = Globals::HMM_ANGLE,
	                    uint window = Globals::HMM_WINDOW_SIZE,
	                    const std::unordered_map<char, T>& hydroMap = scales::NormalizedEisenberg<T>) {

		if (window == 0)
			return 0;

		if (window > seq.length())
			return -1;

		T sumHmSin {0};
		T sumHmCos {0};
		T hmMax {std::numeric_limits<T>::lowest()};
		T hM {0};

		for (size_t i = 0, j = seq.length() - window + 1; i < j; ++i) {
			for (size_t k = i, r = (window + i); k < r; ++k) {
				auto hydro_found = hydroMap.find(seq[k]);
				T hv {hydro_found == hydroMap.end() ? 0 : hydro_found->second};
				T rads {toRadians<T>(angle * (k + i + 1))};
				sumHmSin += hv * std::sin(rads);
				sumHmCos += hv * std::cos(rads);
			}

			hM = std::sqrt(std::pow(sumHmSin, 2) + std::pow(sumHmCos, 2)) / window;

			if (hM > hmMax)
				hmMax = hM;

			sumHmSin = {0};
			sumHmCos = {0};
		}

		return hmMax;
	}

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr T sum(const std::string& seq, const std::unordered_map<char, T>& sumMap = scales::ChartonCTDC<T>) {
		T value {0};

		std::for_each(seq.cbegin(), seq.cend(), [&] (const auto& c){
			auto c_found = sumMap.find(c);
			value += (c_found == sumMap.end()) ? 0 : c_found->second;
		});

		return value;
	}

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr T average(const std::string &seq, const std::unordered_map<char, T> &avgMap = scales::Klein<T>) {
		return sum<T>(seq, avgMap) / seq.length();
	}

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr void composition(const std::string& seq, ra::ReducedAlphabet<T>& ra) noexcept {
		auto& counter = ra.getCounter();
		const auto& counterPtr = ra.getCounterPtr();

		std::for_each(counter.begin(), counter.end(), [&](auto& pair) {
			*pair.second += std::count(seq.cbegin(), seq.cend(), pair.first);
		});

		std::for_each(counterPtr.begin(), counterPtr.end(), [&](auto& cntr) {
			*cntr = static_cast<T>((*cntr / seq.length()) * 100);
		});

	}

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr std::valarray<T> compositionReduceAlph(const std::string& seq, ra::ReducedAlphabet<T>& ra) noexcept {
		composition(seq, ra);
		auto& counterPtr = ra.getCounterPtr();
		std::valarray<T> values(ra.getSize());

		size_t i = 0;
		std::for_each(counterPtr.begin(), counterPtr.end(), [&] (const auto& ptr) {
			values[i++] = *ptr;
		});

		return values;
	}

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr std::valarray<T> distributionReduceAlph(const std::string& seq,
	                                                  ra::ReducedAlphabet<T>& ra, T percentage) noexcept {
		composition(seq, ra);
		std::valarray<T> values(ra.getSize());
		auto& counter = ra.getCounter();
		auto& counterPtr = ra.getCounterPtr();


		std::for_each(counterPtr.begin(), counterPtr.end(), [&](auto& ptr){
			auto nRa = 0;
			if (percentage > 0) {
				nRa = static_cast<int>(std::round(((*ptr) * seq.length()) / 100));
				nRa = static_cast<int>(std::round((nRa * percentage) / 100));
			}

			auto aux = (nRa == 0 && percentage == 0) ? -1 : 0;
			if (nRa == 0 && percentage > 0) {
				*ptr = 0;
			} else {
				for (size_t i = 0; i < seq.length(); ++i) {
					auto c_found = counter.find(seq[i]);
					aux += (c_found != counter.end() && c_found->second == ptr) ? 1 : 0;
					if (nRa == aux) {
						*ptr = (static_cast<T>(i + 1) / seq.length()) * 100;
						break;
					}
				}
			}

		});


		for (size_t i = 0; i < counterPtr.size(); ++i)
			values[i] = *counterPtr[i];

		return values;
	}

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr T transitionReduceAlph(const std::string& seq, ra::ReducedAlphabet<T>& ra) noexcept {
		auto& transTuple = ra.getPairTriPep();
		auto& counter = ra.getCounter();
		auto first = transTuple[0];
		auto second = transTuple[1];
		size_t counting = 0;

		for (std::string::size_type i = 1; i < seq.length(); ++i) {

			const auto& c1 = seq[i - 1];
			const auto& c2 = seq[i];
			auto& c1Counter = counter[c1];
			auto& c2Counter = counter[c2];

			if (c1Counter == nullptr || c2Counter == nullptr || c1Counter == c2Counter)
				continue;


			if ((c1Counter == first && c2Counter == second) ||
					(c1Counter == second && c2Counter == first)) {
				++counting;
			}

		}

		return (static_cast<T>(counting) / (seq.length() - 1)) * 100;
	}

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	constexpr T tripeptideReduceAlph(const std::string& seq, ra::ReducedAlphabet<T>& ra) noexcept {
		auto& tri = ra.getPairTriPep();
		auto& counter = ra.getCounter();
		auto first = tri[0];
		auto second = tri[1];
		auto third = tri[2];

		size_t counting = 0;

		for (std::string::size_type i = 2; i < seq.length(); ++i) {

			auto c1 = seq[i - 2];
			auto c2 = seq[i - 1];
			auto c3 = seq[i];
			auto& c1Counter = counter[c1];
			auto& c2Counter = counter[c2];
			auto& c3Counter = counter[c3];

			if ((c1Counter == first) && (c2Counter == second) && (c3Counter == third))
				++counting;

		}

		return (static_cast<T>(counting) / (seq.length() - 2)) * 100;
	}


}

#endif //INPROT_MD_H
