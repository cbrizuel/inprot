//
// Created by germelcar on 8/31/17.
//

#ifndef INPROT_KMER_OFFSET_H
#define INPROT_KMER_OFFSET_H

#include "fasta_seq.h"
#include "reduced_alphabets.h"
#include "md.h"
#include <memory>
#include <string>
#include "svm_scaling.h"
#include <tbb/tbb.h>

#ifdef USE_LIBSVM
#include "svm.h"
#else
#include "libsvm.h"
#include "svm_model.h"
using namespace libsvm;
#endif


namespace fasta {

	class KmerOffset {

	public:
		//
		// Constructors & destructors
		//
		explicit KmerOffset(const FastaSeq& fseq, size_t offset, uint k, bool amp = false):
				mFseq(fseq), mOffset(offset), mSize(k), mAmp(amp) { }


		//
		// Operators
		//
		bool operator<(const KmerOffset& rhs) const {
			return compare(rhs) < 0;
		}

		bool operator==(const KmerOffset& rhs) const {
			return compare(rhs) == 0;
		}

		//
		// Methods
		//
		uint size() const {
			return mSize;
		}

		bool isConnected(const std::pair<size_t, size_t>& pair) {
			return (pair.first < mOffset && pair.second == mOffset) || ((pair.second + 1) == mOffset);
		}

		bool intersect(const std::pair<size_t, size_t>& pair) {
			auto end = getEnd();

			if (pair.first <= mOffset && pair.second < end)
				return mOffset <= pair.second && end >= pair.second;

			return pair.first >= mOffset && end <= pair.second;
		}

		bool isInside(const std::pair<size_t, size_t>& pair) {
			return mOffset >= pair.first && getEnd() <= pair.second;
		}

#ifdef USE_LIBSVM
		template<typename T, EnableIf<std::is_floating_point<T>>...>
		void evaluate(const SvmScaling<T>& scaling, const std::shared_ptr<svm_model>& model) {
			std::valarray<T> mds = calculateMD<T>();
			scaling.scale(mds);

			std::unique_ptr<svm_node[]> nodes = std::make_unique<svm_node[]>(mds.size() + 1);

			for (size_t i = 0; i < mds.size(); ++i) {
				nodes[i].index = static_cast<int>((i + 1));
				nodes[i].value = mds[i];
			}

			nodes[mds.size()].index = -1;
			nodes[mds.size()].value = std::numeric_limits<T>::max();

			auto label = svm_predict(model.get(), nodes.get());
			mAmp = static_cast<decltype(Globals::SVM_POSITIVE_LABEL)>(label) == Globals::SVM_POSITIVE_LABEL;
		}

#else
		template<typename T, EnableIf<std::is_floating_point<T>>...>
		void evaluate(SvmScaling<T>& scaling, SvmModel<T>& model) {
			std::valarray<T> mds = calculateMD<T>();
			scaling.scale(mds);
			auto ll = model.predict(mds);

			mAmp = (ll == Globals::SVM_POSITIVE_LABEL);
		}
#endif


		//
		// Getters & setters
		//
		const FastaSeq& getFastaSeq() const {
			return mFseq.get();
		}

		const std::string getKmer() const {
			return mFseq.get().getSeq().substr(mOffset, mSize);
		}

		bool isAMP() const {
			return mAmp;
		}

		size_t getOffset() const {
			return mOffset;
		}

		uint getSize() const {
			return mSize;
		}

		size_t getEnd() const {
			return mOffset + mSize;
		}



	private:

		//
		// Private methods
		//
		int compare(const KmerOffset& rhs) const noexcept {
			if (this == &rhs)
				return 0;

			if (mSize < rhs.mSize)
				return -1;
			else if (mSize > rhs.mSize)
				return 1;

			return std::memcmp(&mFseq.get().getSeq().at(mOffset),
			                   &rhs.mFseq.get().getSeq().at(rhs.mOffset), mSize);
		}

		template<typename T, EnableIf<std::is_floating_point<T>>...>
		std::valarray<T> calculateMD() noexcept {

			std::valarray<T> mds (Globals::NUM_MDS);
			const auto& seq = getKmer();


			//
			// Initialize le neccessary maps (reduced alphabets)
			//
			auto stdMap = md::ra::ReducedAlphabet<T>({md::ra::Std::F,
			                                       md::ra::Std::M,
			                                       md::ra::Std::Q});

			auto normVWMap = md::ra::ReducedAlphabet<T>({md::ra::NormVWTomii::MHKFRYW,
			                                          md::ra::NormVWTomii::NVEQIL});

			auto polarityMap = md::ra::ReducedAlphabet<T>(
					{md::ra::PolarityTomii::PATGS,
					  md::ra::PolarityTomii::LIFWCMVY,
					  md::ra::PolarityTomii::HQRKNED});

			auto polarizaMap = md::ra::ReducedAlphabet<T>({md::ra::PolarizabilityTomii::GASDT,
							 md::ra::PolarizabilityTomii::KMHFRYW});

			auto secStructMap = md::ra::ReducedAlphabet<T>({md::ra::SecondStructTomii::EALMQKRH,
							 md::ra::SecondStructTomii::VIYCWFT});

			auto blosumMap = md::ra::ReducedAlphabet<T>({md::ra::Blosum50::CLVIM, md::ra::Blosum50::FWY});

			auto chrgMap = md::ra::ReducedAlphabet<T>({md::ra::ChargeTomii::DE, md::ra::ChargeTomii::KR});

			auto solvAccMap = md::ra::ReducedAlphabet<T>({md::ra::SolventAccTomii::ALFCGIVW,
			                                              md::ra::SolventAccTomii::MPSTHY,
			                                              md::ra::SolventAccTomii::RKQEND});

			auto hydroMap = md::ra::ReducedAlphabet<T>(
							{md::ra::HydrophobicityTomii::GASTPHY,
							 md::ra::HydrophobicityTomii::CLVIMFW,
							 md::ra::HydrophobicityTomii::RKEDQN});


			// MD #1 ==> Length
			mds[md::LENGTH] = static_cast<T>(seq.length());

			// MD #3 ==> Net charge with Isoelectric Point map
			mds[md::NET_CHRG] = md::netCharge<T>(seq);

			// MD #8 ==> hmm
			mds[md::HMM_EISENBERG] = md::hMoment<T>(seq);

			// MD #11 ==> Average charge with Klein map
			mds[md::AVG_CHRG_KLEP810101] = md::average<T>(seq);

			// MD #13 ==> Total/Sum charge with Charton CTCD map
			mds[md::NET_CHRG_CHAM83108] = md::sum<T>(seq);

			// MD #14 ==> Average charge with Kuhno map
			mds[md::AVG_HYDRO_KUHL950101] = md::average<T>(seq, md::scales::KuhnHydrov<T>);

			// MD #15-19
			mds[md::AVG_HYDRO_CIDH_CIDH920102] = md::average<T>(seq, md::scales::CID2<T>);
			mds[md::AVG_HYDRO_CIDH_CIDH920104] = md::average<T>(seq, md::scales::CID4<T>);
			mds[md::AVG_HYDRO_CIDH_CIDH920105] = md::average<T>(seq, md::scales::CID5<T>);
			mds[md::AVG_HYDRO_CIDH_MANP780101] = md::average<T>(
					seq, md::scales::ManavalanPonnuswamy<T>);
			mds[md::AVG_HYDRO_CIDH_PONP800105] = md::average<T>(
					seq, md::scales::Ponnuswamy5<T>);

			// MD #21
			mds[md::AVG_HYDRO_PRAM900101] = md::average<T>(seq, md::scales::Prabhakaran<T>);

			// MD # 22
			mds[md::AVG_HYDRO_SWER830101] = md::average<T>(seq, md::scales::SweetEisenberg<T>);

			// MD #24-27
			mds[md::AVG_HYDRO_ZIMJ680101] = md::average<T>(seq, md::scales::Zimmerman<T>);
			mds[md::AVG_HYDRO_WOLR790101] = md::average<T>(seq, md::scales::Wolfenden<T>);
			mds[md::AVG_HYDRO_CASG920101] = md::average<T>(seq, md::scales::CasariSippl<T>);
			mds[md::AVG_HYDRO_TOSSI2002] = md::average<T>(seq, md::scales::Tossi<T>);



			// +++
			// MD #2,5,7 ==> Composition of standard aminoacids: F, M, Q
			// +++
			auto stdComps = md::compositionReduceAlph<T>(seq, stdMap);
			mds[md::COMP_STD_F] = stdComps[0];
			mds[md::COMP_STD_M] = stdComps[1];
			mds[md::COMP_STD_Q] = stdComps[2];

			// +++
			// MD #4,6,34 ==> Distribution and composition with NormVWTomii reduced alphabet
			// +++
			// MD #4
			auto normvwValues = md::distributionReduceAlph<T>(seq, normVWMap, 50);
			mds[md::DIST_NORMVW_TOMII_MHKFRYW_50] = normvwValues[0];

			// MD #6
			normVWMap.resetCounter();
			normvwValues = md::distributionReduceAlph<T>(seq, normVWMap, 75);
			mds[md::DIST_NORMVW_TOMII_NVEQIL_75] = normvwValues[1];

			// MD #34
			normVWMap.resetCounter();
			normvwValues = md::compositionReduceAlph<T>(seq, normVWMap);
			mds[md::COMP_NORMVW_TOMII_MHKFRYW] = normvwValues[0];


			// +++
			// MD #9,10,12
			// +++
			auto polarityValues = md::distributionReduceAlph<T>(seq, polarityMap, 0);
			// MD #9,10
			mds[md::DIST_POLAR_TOMII_PATGS_0] = polarityValues[0];
			mds[md::DIST_POLAR_TOMII_LIFWCMVY_0] = polarityValues[1];

			// MD #12
			polarityMap.resetCounter();
			polarityValues = md::distributionReduceAlph<T>(seq, polarityMap, 25);
			mds[md::DIST_POLAR_TOMII_HQRKNED_25] = polarityValues[2];

			// +++
			// MD #20,23,39
			// +++
			auto polarizaValues = md::distributionReduceAlph<T>(seq, polarizaMap, 75);
			// MD #20
			mds[md::DIST_POLAR_TOMII_GASDT_75] = polarizaValues[0];

			// MD #23
			polarizaMap.resetCounter();
			polarizaValues = md::distributionReduceAlph<T>(seq, polarizaMap, 100);
			mds[md::DIST_POLAR_TOMII_GASDT_100] = polarizaValues[0];

			// MD #39
			polarizaMap.resetCounter();
			polarizaValues = md::compositionReduceAlph<T>(seq, polarizaMap);
			mds[md::COMP_POLAR_TOMII_KMHFRYW] = polarizaValues[1];

			// +++
			// MD #28,29,42
			// +++
			auto ssValues = md::distributionReduceAlph<T>(seq, secStructMap, 50);

			// MD #28
			mds[md::DIST_SS_TOMII_EALMQKRH_50] = ssValues[0];

			// MD #29
			secStructMap.resetCounter();
			ssValues = md::distributionReduceAlph<T>(seq, secStructMap, 100);
			mds[md::DIST_SS_TOMII_EALMQKRH_100] = ssValues[0];

			// MD #42
			secStructMap.resetCounter();
			ssValues = md::compositionReduceAlph<T>(seq, secStructMap);
			mds[md::COMP_SS_TOMII_VIYCWFT] = ssValues[1];

			// +++
			// MD #30,32
			// +++
			auto blosumValues = md::compositionReduceAlph<T>(seq, blosumMap);
			mds[md::COMP_BLOSUM50_CLVIM] = blosumValues[0];
			mds[md::COMP_BLOSUM50_FWY] = blosumValues[1];

			// +++
			// MD #31,33,40,41
			// +++
			auto chrgValues = md::distributionReduceAlph<T>(seq, chrgMap, 0);

			// MD # 31
			mds[md::DIST_CHRG_TOMII_DE_0] = chrgValues[0];

			// MD #33
			chrgMap.resetCounter();
			chrgValues = md::distributionReduceAlph<T>(seq, chrgMap, 100);
			mds[md::DIST_CHRG_TOMII_KR_100] = chrgValues[1];

			// MD #40,41
			chrgMap.resetCounter();
			chrgValues = md::compositionReduceAlph<T>(seq, chrgMap);
			mds[md::COMP_CHRG_TOMII_DE] = chrgValues[0];
			mds[md::COMP_CHRG_TOMII_KR] = chrgValues[1];

			// +++
			// MD #35-38,43,47
			// +++
			auto saccValues = md::distributionReduceAlph<T>(seq, solvAccMap, 0);

			// MD #35-37
			mds[md::DIST_SOLVENT_TOMII_ALFCGIVW_0] = saccValues[0];
			mds[md::DIST_SOLVENT_TOMII_MPSTHY_0] = saccValues[1];
			mds[md::DIST_SOLVENT_TOMII_RKQEND_0] = saccValues[2];

			// MD #38
			solvAccMap.resetCounter();
			saccValues = md::distributionReduceAlph<T>(seq, solvAccMap, 25);
			mds[md::DIST_SOLVENT_TOMII_RKQEND_25] = saccValues[2];

			// MD #43
			solvAccMap.resetCounter();
			saccValues = md::compositionReduceAlph<T>(seq, solvAccMap);
			mds[md::COMP_SA_TOMII_ALFCGIVW] = saccValues[0];

			// MD #47
			solvAccMap.resetCounter();
			solvAccMap.setPairTriPep({md::ra::SolventAccTomii::ALFCGIVW, md::ra::SolventAccTomii::RKQEND});
			mds[md::TRANS_SA_TOMII_ALFCGIVW_RKQEND] = md::transitionReduceAlph<T>(seq, solvAccMap);


			// +++
			// MD #44-46,48-51
			// +++
			auto hydroValues = md::distributionReduceAlph<T>(seq, hydroMap, 0);

			// MD #48,49
			mds[md::DIST_HYDRO_TOMII_GASTPHY_0] = hydroValues[0];
			mds[md::DIST_HYDRO_TOMII_CLVIMFW_0] = hydroValues[1];

			// MD #51
			hydroMap.resetCounter();
			hydroValues = md::distributionReduceAlph<T>(seq, hydroMap, 75);
			mds[md::DIST_HYDRO_TOMII_GASTPHY_75] = hydroValues[0];

			// MD #44
			hydroMap.setPairTriPep({md::ra::HydrophobicityTomii::CLVIMFW, md::ra::HydrophobicityTomii::RKEDQN});
			mds[md::TRANS_HYDRO_TOMII_CLVIM_RKEDQN] = md::transitionReduceAlph<T>(seq, hydroMap);

			// MD #45,46,50
			//
			// MD #45
			hydroMap.setPairTriPep({md::ra::HydrophobicityTomii::RKEDQN,
			                         md::ra::HydrophobicityTomii::CLVIMFW,
			                         md::ra::HydrophobicityTomii::GASTPHY});
			mds[md::TRIP_HYDRO_TOMII_RKEDQN_CLVIMFW_GASTPHY] = md::tripeptideReduceAlph<T>(seq, hydroMap);

			// MD #46
			hydroMap.setPairTriPep({md::ra::HydrophobicityTomii::CLVIMFW,
			                         md::ra::HydrophobicityTomii::CLVIMFW,
			                         md::ra::HydrophobicityTomii::GASTPHY});
			mds[md::TRIP_HYDRO_TOMII_CLVIMFW_CLVIMFW_GASTPHY] = md::tripeptideReduceAlph<T>(seq, hydroMap);

			// MD #50
			hydroMap.setPairTriPep({md::ra::HydrophobicityTomii::CLVIMFW,
			                         md::ra::HydrophobicityTomii::CLVIMFW,
			                         md::ra::HydrophobicityTomii::CLVIMFW});
			mds[md::TRIP_HYDRO_TOMII_CLVIMFW_CLVIMFW_CLVIMFW] = md::tripeptideReduceAlph<T>(seq, hydroMap);

			return mds;
		}

		//
		// Fields
		//
		std::reference_wrapper<const FastaSeq> mFseq;
		size_t mOffset;
		uint mSize;
		bool mAmp;

	};

}


#endif //INPROT_KMER_OFFSET_H
