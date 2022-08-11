//
// Created by germelcar on 9/1/17.
//

#ifndef INPROT_RA_H
#define INPROT_RA_H

#include <unordered_map>
#include <map>
#include <memory>
#include <algorithm>
#include <vector>
#include "ras.h"

template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;

namespace md {

	namespace ra {

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		class ReducedAlphabet {

			using alphaPtr = std::shared_ptr<T>;

		public:
			//
			// Constructors
			//
			ReducedAlphabet() = delete;

			ReducedAlphabet(const std::initializer_list<Std>& flags) noexcept {
				std(flags);
			}

			ReducedAlphabet(const std::initializer_list<Blosum50>& flags) noexcept {
				blosum50Murphy(flags);
			}

			ReducedAlphabet(const std::initializer_list<NormVWTomii>& flags) noexcept {
				normVWTomii(flags);
			}

			ReducedAlphabet(const std::initializer_list<PolarityTomii>& flags) noexcept {
				polarityTomii(flags);
			}

			ReducedAlphabet(const std::initializer_list<PolarizabilityTomii>& flags) noexcept {
				polarizabilityTomii(flags);
			}

			ReducedAlphabet(const std::initializer_list<ChargeTomii>& flags) noexcept {
				chargeTomii(flags);
			}

			ReducedAlphabet(const std::initializer_list<SecondStructTomii>& flags) noexcept {
				secondStructTomii(flags);
			}

			ReducedAlphabet(const std::initializer_list<SolventAccTomii>& flags,
			                         const std::tuple<SolventAccTomii, SolventAccTomii>& trans =
			                         std::make_tuple(SolventAccTomii::ALFCGIVW, SolventAccTomii::RKQEND)) noexcept {
				solvAccessTomii(flags, trans);
			}

			ReducedAlphabet(const std::initializer_list<HydrophobicityTomii>& flags) noexcept {
				hydroTomii(flags);
			}


			void resetCounter() const {
				std::for_each(mCounter.begin(), mCounter.end(), [](auto& pair){ *pair.second = 0; });
			}


			//
			// Getters & setters
			//
			std::unordered_map<char, alphaPtr>& getCounter() {
				return mCounter;
			}

			auto getSize() const {
				return mCounterPtr.size();
			}

			const auto& getCounterPtr() const {
				return mCounterPtr;
			}

			auto& getPairTriPep() const {
				return mPairTriPep;
			}

			template <typename PT, EnableIf<std::is_enum<PT>>...>
			void setPairTriPep(std::initializer_list<PT>&& pair, bool alluniques = false) {

				if (pair.size() != 2 && pair.size() != 3)
					return;

				std::vector<PT> vecPair(pair);
				if (alluniques && (std::unique(vecPair.begin(), vecPair.end()) != vecPair.end()))
					return;

				for (const auto& p : vecPair)
					if (!hasPairElement(p))
						return;

				mPairTriPep.clear();
				mPairTriPep.reserve(pair.size());

				for (const auto& p : vecPair) {
					inPairElement(p);
				}

			}


		private:
			//
			// Factory methods
			//
			void std(const std::initializer_list<Std>& flags) noexcept {
				std::unordered_map<char, alphaPtr> counting(3);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());

				for (const auto& f : flags) {

					if (f == Std::A) {
						auto ptrA = std::make_shared<T>(0);
						counting['A'] = ptrA;
						counterPtr.emplace_back(ptrA);
					}

					if (f == Std::C) {
						auto ptrC = std::make_shared<T>(0);
						counting['C'] = ptrC;
						counterPtr.emplace_back(ptrC);
					}

					if (f == Std::D) {
						auto ptrD = std::make_shared<T>(0);
						counting['D'] = ptrD;
						counterPtr.emplace_back(ptrD);
					}

					if (f == Std::E) {
						auto ptrE = std::make_shared<T>(0);
						counting['E'] = ptrE;
						counterPtr.emplace_back(ptrE);
					}

					if (f == Std::F) {
						auto ptrF = std::make_shared<T>(0);
						counting['F'] = ptrF;
						counterPtr.emplace_back(ptrF);
					}

					if (f == Std::G) {
						auto ptrG = std::make_shared<T>(0);
						counting['G'] = ptrG;
						counterPtr.emplace_back(ptrG);
					}

					if (f == Std::H) {
						auto ptrH = std::make_shared<T>(0);
						counting['H'] = ptrH;
						counterPtr.emplace_back(ptrH);
					}

					if (f == Std::I) {
						auto ptrI = std::make_shared<T>(0);
						counting['I'] = ptrI;
						counterPtr.emplace_back(ptrI);
					}

					if (f == Std::K) {
						auto ptrK = std::make_shared<T>(0);
						counting['K'] = ptrK;
						counterPtr.emplace_back(ptrK);
					}

					if (f == Std::L) {
						auto ptrL = std::make_shared<T>(0);
						counting['L'] = ptrL;
						counterPtr.emplace_back(ptrL);
					}

					if (f == Std::M) {
						auto ptrM = std::make_shared<T>(0);
						counting['M'] = ptrM;
						counterPtr.emplace_back(ptrM);
					}

					if (f == Std::N) {
						auto ptrN = std::make_shared<T>(0);
						counting['N'] = ptrN;
						counterPtr.emplace_back(ptrN);
					}

					if (f == Std::P) {
						auto ptrP = std::make_shared<T>(0);
						counting['P'] = ptrP;
						counterPtr.emplace_back(ptrP);
					}

					if (f == Std::Q) {
						auto ptrQ = std::make_shared<T>(0);
						counting['Q'] = ptrQ;
						counterPtr.emplace_back(ptrQ);
					}

					if (f == Std::R) {
						auto ptrR = std::make_shared<T>(0);
						counting['R'] = ptrR;
						counterPtr.emplace_back(ptrR);
					}

					if (f == Std::S) {
						auto ptrS = std::make_shared<T>(0);
						counting['S'] = ptrS;
						counterPtr.emplace_back(ptrS);
					}

					if (f == Std::T) {
						auto ptrT = std::make_shared<T>(0);
						counting['T'] = ptrT;
						counterPtr.emplace_back(ptrT);
					}

					if (f == Std::V) {
						auto ptrV = std::make_shared<T>(0);
						counting['V'] = ptrV;
						counterPtr.emplace_back(ptrV);
					}

					if (f == Std::W) {
						auto ptrW = std::make_shared<T>(0);
						counting['W'] = ptrW;
						counterPtr.emplace_back(ptrW);
					}

					if (f == Std::Y) {
						auto ptrY = std::make_shared<T>(0);
						counting['Y'] = ptrY;
						counterPtr.emplace_back(ptrY);
					}


				} // End of "for"

				mCounter = counting;
				mCounterPtr = counterPtr;
			}

			void blosum50Murphy(const std::initializer_list<Blosum50>& flags) noexcept {
				std::unordered_map<char, alphaPtr> counting(5);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());

				for (const auto& f : flags) {

					if (f == Blosum50::FWY) {
						auto fwy = std::make_shared<T>(0);
						counterPtr.emplace_back(fwy);
						counting['F'] = fwy;
						counting['W'] = fwy;
						counting['Y'] = fwy;
					}

					if (f == Blosum50::CLVIM) {
						auto clvim = std::make_shared<T>(0);
						counterPtr.emplace_back(clvim);
						counting['L'] = clvim;
						counting['M'] = clvim;
						counting['C'] = clvim;
						counting['I'] = clvim;
						counting['V'] = clvim;
					}

					if (f == Blosum50::H) {
						auto h = std::make_shared<T>(0);
						counting['H'] = h;
						counterPtr.emplace_back(h);
					}

					if (f == Blosum50::AG) {
						auto ag = std::make_shared<T>(0);
						counterPtr.emplace_back(ag);
						counting['A'] = ag;
						counting['G'] = ag;
					}

					if (f == Blosum50::ST) {
						auto st = std::make_shared<T>(0);
						counterPtr.emplace_back(st);
						counting['S'] = st;
						counting['T'] = st;
					}

					if (f == Blosum50::DENQ) {
						auto denq = std::make_shared<T>(0);
						counterPtr.emplace_back(denq);
						counting['N'] = denq;
						counting['D'] = denq;
						counting['Q'] = denq;
						counting['E'] = denq;
					}

					if (f == Blosum50::KR) {
						auto kr = std::make_shared<T>(0);
						counterPtr.emplace_back(kr);
						counting['R'] = kr;
						counting['K'] = kr;
					}

					if (f == Blosum50::P) {
						auto p = std::make_shared<T>(0);
						counting['P'] = p;
						counterPtr.emplace_back(p);
					}

				} // End of "for"

				mCounter = counting;
				mCounterPtr = counterPtr;

			}

			void normVWTomii(const std::initializer_list<NormVWTomii>& flags) noexcept {
				std::unordered_map<char, alphaPtr> counting(7);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());

				for (const auto& f : flags) {

					if (f == NormVWTomii::GASTCPD) {
						auto gastcpd = std::make_shared<T>(0);
						counterPtr.emplace_back(gastcpd);
						counting['A'] = gastcpd;
						counting['D'] = gastcpd;
						counting['C'] = gastcpd;
						counting['P'] = gastcpd;
						counting['S'] = gastcpd;
						counting['T'] = gastcpd;
						counting['G'] = gastcpd;
					}

					if (f == NormVWTomii::MHKFRYW) {
						auto mhkfryw = std::make_shared<T>(0);
						counterPtr.emplace_back(mhkfryw);
						counting['R'] = mhkfryw;
						counting['K'] = mhkfryw;
						counting['M'] = mhkfryw;
						counting['F'] = mhkfryw;
						counting['W'] = mhkfryw;
						counting['H'] = mhkfryw;
						counting['Y'] = mhkfryw;
					}

					if (f == NormVWTomii::NVEQIL) {
						auto nveqil = std::make_shared<T>(0);
						counterPtr.emplace_back(nveqil);
						counting['L'] = nveqil;
						counting['N'] = nveqil;
						counting['Q'] = nveqil;
						counting['E'] = nveqil;
						counting['I'] = nveqil;
						counting['V'] = nveqil;
					}


				} // End of "for"

				mCounter = counting;
				mCounterPtr = counterPtr;
			}

			void polarityTomii(const std::initializer_list<PolarityTomii>& flags) noexcept {
				std::unordered_map<char, alphaPtr> counting(8);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());

				for (const auto& f : flags) {

					if (f == PolarityTomii::LIFWCMVY) {
						auto lifwcmvy = std::make_shared<T>(0);
						counterPtr.emplace_back(lifwcmvy);
						counting['L'] = lifwcmvy;
						counting['M'] = lifwcmvy;
						counting['F'] = lifwcmvy;
						counting['C'] = lifwcmvy;
						counting['W'] = lifwcmvy;
						counting['Y'] = lifwcmvy;
						counting['I'] = lifwcmvy;
						counting['V'] = lifwcmvy;
					}

					if (f == PolarityTomii::PATGS) {
						auto patgs = std::make_shared<T>(0);
						counterPtr.emplace_back(patgs);
						counting['A'] = patgs;
						counting['P'] = patgs;
						counting['S'] = patgs;
						counting['T'] = patgs;
						counting['G'] = patgs;
					}

					if (f == PolarityTomii::HQRKNED) {
						auto hqrkned = std::make_shared<T>(0);
						counterPtr.emplace_back(hqrkned);
						counting['R'] = hqrkned;
						counting['K'] = hqrkned;
						counting['N'] = hqrkned;
						counting['D'] = hqrkned;
						counting['Q'] = hqrkned;
						counting['E'] = hqrkned;
						counting['H'] = hqrkned;
					}

				} // End of "for"

				mCounter = counting;
				mCounterPtr = counterPtr;
			}

			void polarizabilityTomii(const std::initializer_list<PolarizabilityTomii>& flags) noexcept {
				std::unordered_map<char, alphaPtr> counting(8);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());

				for (const auto& f : flags) {

					if (f == PolarizabilityTomii::GASDT) {
						auto gasdt = std::make_shared<T>(0);
						counterPtr.emplace_back(gasdt);
						counting['A'] = gasdt;
						counting['D'] = gasdt;
						counting['S'] = gasdt;
						counting['T'] = gasdt;
						counting['G'] = gasdt;
					}

					if (f == PolarizabilityTomii::CPNVEQIL) {
						auto cpnveqil = std::make_shared<T>(0);
						counterPtr.emplace_back(cpnveqil);
						counting['L'] = cpnveqil;
						counting['N'] = cpnveqil;
						counting['C'] = cpnveqil;
						counting['P'] = cpnveqil;
						counting['Q'] = cpnveqil;
						counting['E'] = cpnveqil;
						counting['I'] = cpnveqil;
						counting['V'] = cpnveqil;
					}

					if (f == PolarizabilityTomii::KMHFRYW) {
						auto kmhfryw = std::make_shared<T>(0);
						counterPtr.emplace_back(kmhfryw);
						counting['R'] = kmhfryw;
						counting['K'] = kmhfryw;
						counting['M'] = kmhfryw;
						counting['F'] = kmhfryw;
						counting['W'] = kmhfryw;
						counting['H'] = kmhfryw;
						counting['Y'] = kmhfryw;
					}

				} // End of "for"

				mCounter = counting;
				mCounterPtr = counterPtr;
			}

			void chargeTomii(const std::initializer_list<ChargeTomii>& flags) noexcept {
				std::unordered_map<char, alphaPtr> counting(16);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());

				for (const auto& f : flags) {

					if (f == ChargeTomii::KR) {
						auto kr = std::make_shared<T>(0);
						counterPtr.emplace_back(kr);
						counting['R'] = kr;
						counting['K'] = kr;
					}

					if (f == ChargeTomii::ANCQGHILMFPSTWYV) {
						auto ancqghilmfpstwyv = std::make_shared<T>(0);
						counterPtr.emplace_back(ancqghilmfpstwyv);
						counting['A'] = ancqghilmfpstwyv;
						counting['L'] = ancqghilmfpstwyv;
						counting['N'] = ancqghilmfpstwyv;
						counting['M'] = ancqghilmfpstwyv;
						counting['F'] = ancqghilmfpstwyv;
						counting['C'] = ancqghilmfpstwyv;
						counting['P'] = ancqghilmfpstwyv;
						counting['Q'] = ancqghilmfpstwyv;
						counting['S'] = ancqghilmfpstwyv;
						counting['T'] = ancqghilmfpstwyv;
						counting['G'] = ancqghilmfpstwyv;
						counting['W'] = ancqghilmfpstwyv;
						counting['H'] = ancqghilmfpstwyv;
						counting['Y'] = ancqghilmfpstwyv;
						counting['I'] = ancqghilmfpstwyv;
						counting['V'] = ancqghilmfpstwyv;
					}

					if (f == ChargeTomii::DE) {
						auto de = std::make_shared<T>(0);
						counterPtr.emplace_back(de);
						counting['D'] = de;
						counting['E'] = de;
					}

				} // End of "for"

				mCounter = counting;
				mCounterPtr = counterPtr;
			}

			void secondStructTomii(const std::initializer_list<SecondStructTomii>& flags) noexcept {
				std::unordered_map<char, alphaPtr> counting(8);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());

				for (const auto& f : flags) {

					if (f == SecondStructTomii::EALMQKRH) {
						auto ealmqkrh = std::make_shared<T>(0);
						counterPtr.emplace_back(ealmqkrh);
						counting['A'] = ealmqkrh;
						counting['L'] = ealmqkrh;
						counting['R'] = ealmqkrh;
						counting['K'] = ealmqkrh;
						counting['M'] = ealmqkrh;
						counting['Q'] = ealmqkrh;
						counting['E'] = ealmqkrh;
						counting['H'] = ealmqkrh;
					}

					if (f == SecondStructTomii::VIYCWFT) {
						auto viycwft = std::make_shared<T>(0);
						counterPtr.emplace_back(viycwft);
						counting['F'] = viycwft;
						counting['C'] = viycwft;
						counting['T'] = viycwft;
						counting['W'] = viycwft;
						counting['Y'] = viycwft;
						counting['I'] = viycwft;
						counting['V'] = viycwft;
					}

					if (f == SecondStructTomii::GNPSD) {
						auto gnpsd = std::make_shared<T>(0);
						counterPtr.emplace_back(gnpsd);
						counting['N'] = gnpsd;
						counting['D'] = gnpsd;
						counting['P'] = gnpsd;
						counting['S'] = gnpsd;
						counting['G'] = gnpsd;
					}

				} // End of "for"

				mCounter = counting;
				mCounterPtr = counterPtr;

			}

			void solvAccessTomii(const std::initializer_list<SolventAccTomii>& flags,
			                     const std::tuple<SolventAccTomii, SolventAccTomii >& trans =
			                     std::make_tuple(SolventAccTomii::ALFCGIVW, SolventAccTomii::RKQEND)) noexcept {


				auto trans1 = std::get<0>(trans);
				auto trans2 = std::get<1>(trans);
				auto foundTrans1 = std::find(flags.begin(), flags.end(), trans1);
				auto foundTrans2 = std::find(flags.begin(), flags.end(), trans2);

				if (foundTrans1 == flags.end() ||  foundTrans2 == flags.end())
					return;

				std::unordered_map<char, alphaPtr> counting(8);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());
				std::tuple<alphaPtr, alphaPtr, alphaPtr> transitions;


				for (const auto& f : flags) {

					if (f == SolventAccTomii::ALFCGIVW) {
						auto alfcgivw = std::make_shared<T>(0);
						counterPtr.emplace_back(alfcgivw);
						counting['A'] = alfcgivw;
						counting['L'] = alfcgivw;
						counting['F'] = alfcgivw;
						counting['C'] = alfcgivw;
						counting['G'] = alfcgivw;
						counting['W'] = alfcgivw;
						counting['I'] = alfcgivw;
						counting['V'] = alfcgivw;
					}

					if (f == SolventAccTomii::RKQEND) {
						auto rkqend = std::make_shared<T>(0);
						counterPtr.emplace_back(rkqend);
						counting['R'] = rkqend;
						counting['K'] = rkqend;
						counting['N'] = rkqend;
						counting['D'] = rkqend;
						counting['Q'] = rkqend;
						counting['E'] = rkqend;
					}

					if (f == SolventAccTomii::MPSTHY) {
						auto mpsthy = std::make_shared<T>(0);
						counterPtr.emplace_back(mpsthy);
						counting['M'] = mpsthy;
						counting['P'] = mpsthy;
						counting['S'] = mpsthy;
						counting['T'] = mpsthy;
						counting['H'] = mpsthy;
						counting['Y'] = mpsthy;
					}

				} // End of "for"

				mCounter = counting;
				mCounterPtr = counterPtr;
			}

			void hydroTomii(const std::initializer_list<HydrophobicityTomii>& flags) noexcept {
				std::unordered_map<char, alphaPtr> counting(8);
				std::vector<std::shared_ptr<T>> counterPtr;
				counterPtr.reserve(flags.size());

				for (const auto& f : flags) {

					if (f == HydrophobicityTomii::RKEDQN) {
						auto rkedqn = std::make_shared<T>(0);
						counterPtr.emplace_back(rkedqn);
						counting['R'] = rkedqn;
						counting['K'] = rkedqn;
						counting['N'] = rkedqn;
						counting['D'] = rkedqn;
						counting['Q'] = rkedqn;
						counting['E'] = rkedqn;
					}

					if (f == HydrophobicityTomii::GASTPHY) {
						auto gastphy = std::make_shared<T>(0);
						counterPtr.emplace_back(gastphy);
						counting['A'] = gastphy;
						counting['P'] = gastphy;
						counting['S'] = gastphy;
						counting['T'] = gastphy;
						counting['G'] = gastphy;
						counting['H'] = gastphy;
						counting['Y'] = gastphy;
					}

					if (f == HydrophobicityTomii::CLVIMFW) {
						auto clvimfw = std::make_shared<T>(0);
						counterPtr.emplace_back(clvimfw);
						counting['L'] = clvimfw;
						counting['M'] = clvimfw;
						counting['F'] = clvimfw;
						counting['C'] = clvimfw;
						counting['W'] = clvimfw;
						counting['I'] = clvimfw;
						counting['V'] = clvimfw;
					}
				}

				mCounter = counting;
				mCounterPtr = counterPtr;

			}


			//
			// Private helper methods
			//
			bool hasPairElement(const SolventAccTomii& p) {

				if (p == SolventAccTomii::ALFCGIVW)
					return mCounter['A'] != nullptr;

				if (p == SolventAccTomii::RKQEND)
					return mCounter['R'] != nullptr;

				if (p == SolventAccTomii::MPSTHY)
					return mCounter['M'] != nullptr;

				return false;
			}

			bool hasPairElement(const HydrophobicityTomii& p) {

				if (p == HydrophobicityTomii::RKEDQN)
					return mCounter['R'] != nullptr;

				if (p == HydrophobicityTomii::GASTPHY)
					return mCounter['G'] != nullptr;

				if (p == HydrophobicityTomii::CLVIMFW)
					return mCounter['C'] != nullptr;

				return false;
			}

			void inPairElement(const SolventAccTomii& p) {

				if (p == SolventAccTomii::ALFCGIVW)
					mPairTriPep.emplace_back(mCounter['A']);

				if (p == SolventAccTomii::RKQEND)
					mPairTriPep.emplace_back(mCounter['R']);

				if (p == SolventAccTomii::MPSTHY)
					mPairTriPep.emplace_back(mCounter['M']);
			}

			void inPairElement(const HydrophobicityTomii& p) {

				if (p == HydrophobicityTomii::RKEDQN)
					mPairTriPep.emplace_back(mCounter['R']);

				if (p == HydrophobicityTomii::GASTPHY)
					mPairTriPep.emplace_back(mCounter['G']);

				if (p == HydrophobicityTomii::CLVIMFW)
					mPairTriPep.emplace_back(mCounter['C']);
			}


			//
			// Fields
			//
			std::unordered_map<char, alphaPtr> mCounter;
			std::vector<alphaPtr> mPairTriPep;
			std::vector<alphaPtr> mCounterPtr;

		};


	}

}


#endif //INPROT_RA_H
