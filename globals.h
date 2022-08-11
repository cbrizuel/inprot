//
// Created by germelcar on 8/30/17.
//

#ifndef INPROT_GLOBALS_H
#define INPROT_GLOBALS_H

#include <thread>

namespace Globals {

	using namespace std::string_literals;

    constexpr size_t     MIN_KMER_SIZE = 10;
    constexpr size_t     MAX_KMER_SIZE = 200;
    const std::string    ALPHABET = "ACDEFGHIKLMNPQRSTVWY"s;
    constexpr size_t     NUM_MDS {51};
	constexpr int        PH_NET_CHARGE {9};
	constexpr uint       HMM_ANGLE {100};
	constexpr uint       HMM_WINDOW_SIZE {10};
	constexpr int        SVM_POSITIVE_LABEL = 1;
	constexpr int        SVM_NEGATIVE_LABEL = -1;

	enum WRITE_PREDICTEDS: uint { WRITE_NONE_PREDS, WRITE_AMPS_PREDS, WRITE_NAMPS_PREDS, WRITE_BOTHS_PREDS };

	// APP - CLI
	const std::string    APP_VERSION = "v0.2.3"s;
	const std::string    APP_NAME = "InProt - In silico Proteolysis " + APP_VERSION;
	const int            MAX_NUM_THREADS = std::thread::hardware_concurrency();


	// Class obtained (and modified) from: https://stackoverflow.com/a/7277333
	class comma_numpunct : public std::numpunct<char> {

	protected:

		virtual char do_thousands_sep() const {
			return ',';
		}

		virtual std::string do_grouping() const {
			return "\03";
		}
	};


}

#endif //INPROT_GLOBALS_H
