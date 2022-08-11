//
// Created by germelcar on 9/13/17.
//

#ifndef INPROT_CLI_H
#define INPROT_CLI_H

#include "CLI11.hpp"
#include "globals.h"
#include "rang.hpp"
#include <memory>
#include <thread>

using namespace rang;

class Cli {

public:

	Cli() {
		generateOptions();
	}

	//
	// Methods
	//
	int parseArgs(int argc, char *argv[], bool printOpts = true) {
		try {

			mApp.parse(argc, argv);

			std::string kmerRangeMsg = "k-mer size should be a value between [" +
			                           std::to_string(Globals::MIN_KMER_SIZE) + "," +
			                           std::to_string(Globals::MAX_KMER_SIZE) + "]";

			if (mLowerKmerSize < Globals::MIN_KMER_SIZE || mLowerKmerSize > Globals::MAX_KMER_SIZE)
				throw CLI::ValidationError("Invalid lower k-mer size. " + kmerRangeMsg);

			if (mUpperKmerSize < Globals::MIN_KMER_SIZE || mUpperKmerSize > Globals::MAX_KMER_SIZE)
				throw CLI::ValidationError("Invalid upper k-mer size. " + kmerRangeMsg);

			if (mLowerKmerSize > mUpperKmerSize)
				throw CLI::ValidationError("Invalid values for lower and upper k-mer size");

			if (mWritePreds > Globals::WRITE_BOTHS_PREDS)
				mWritePreds = Globals::WRITE_NONE_PREDS;


		} catch (CLI::ParseError& e) {
			return mApp.exit(e);
		}

		if (printOpts)
			printOptions();

		return 1; // all OK
	}

	void printOptions() const {
		std::cout << fg::green << style::bold << "-------------------------- " << "CONFIGURATION"
		          << " --------------------------" << style::reset << "\n";

		std::cout << style::bold << fg::green << "Input file: " << style::reset << fg::green
		          << mInputFile << style::reset << "\n";

		std::cout << style::bold << fg::green << "Output's basename: " << style::reset << fg::green
		          << mOutputBaseName << style::reset << "\n";

		std::cout << style::bold << fg::green << "SVM's model file: " << style::reset << fg::green
		          << mModelFile << style::reset << "\n";

		std::cout << style::bold << fg::green << "SVM's scaling file: " << style::reset << fg::green
		          << mScalingFile << style::reset << "\n";

		std::cout << style::bold << fg::green << "Lower k-mer size: " << style::reset << fg::green
		          << std::to_string(mLowerKmerSize) << style::reset << "\n";

		std::cout << style::bold << fg::green << "Upper k-mer size: " << style::reset << fg::green
		          << std::to_string(mUpperKmerSize) << style::reset << "\n";

		std::cout << style::bold << fg::green << "Number of threads: " << style::reset << fg::green
		          << (mNumThreads == -1 ? "automatic" : std::to_string(mNumThreads))
		          << style::reset << "\n";

		std::cout << style::bold << fg::green << "Write predicteds k-mers: " << style::reset << fg::green;

		if (mWritePreds == Globals::WRITE_NONE_PREDS)
			std::cout << "none" << style::reset << "\n";
		else if (mWritePreds == Globals::WRITE_AMPS_PREDS)
			std::cout << "amps" << style::reset << "\n";
		else if (mWritePreds == Globals::WRITE_NAMPS_PREDS)
			std::cout << "namps" << style::reset << "\n";
		else if (mWritePreds == Globals::WRITE_BOTHS_PREDS)
			std::cout << "boths (amps, namps)" << style::reset << "\n";

        std::cout << style::bold << fg::green << "Aware memory mode (low-memory consumption): "
                  << style::reset << fg::green << ((mAware) ? "true" : "false") << style::reset << "\n";

		std::cout << style::bold << fg::green << "Verbose mode (show extra info.): " << style::reset << fg::green
		          << ((mVerbose) ? "true" : "false") << style::reset << "\n";

		std::cout << fg::green << style::bold
		          << "---------------------------------------------------------------------"
		          << style::reset << std::endl;
	}

	//
	// Getters & setters
	//
	const std::string& getInputFile() const {
		return mInputFile;
	}

	const std::string& getOutputFile() const {
		return mOutputBaseName;
	}

	const std::string& getScalingFile() const {
		return mScalingFile;
	}

	const std::string& getModelFile() const {
		return mModelFile;
	}

	uint getLowerKmer() const {
		return mLowerKmerSize;
	}

	uint getUpperKmer() const {
		return mUpperKmerSize;
	}

	int getNumThreads() const {
		return mNumThreads;
	}

	uint getWritePreds() const {
		return mWritePreds;
	}

    bool hasAwareMode() const {
        return mAware;
    }

	bool hasVerboseMode() const {
		return mVerbose;
	}


private:
	void generateOptions() {
		mApp.add_option("-i,--input",
		                mInputFile,
		                "Proteome's file where to extract the k-mers")
				->required()->check(CLI::ExistingFile);

		mApp.add_option("-o,--output",
		                mOutputBaseName,
		                "Output's basename where to put the AMP-predicted k-mers")
				->required()->check([](auto filename) {
					std::ofstream fout(filename, std::ios_base::out | std::ios_base::ate);
					auto allgood = fout.good();
					std::remove(filename.c_str());
					return allgood;
				});

		mApp.add_option("-s,--scaling",
		                mScalingFile,
		                "SVM's scaling factors configuration file")
				->required()->check(CLI::ExistingFile);

		mApp.add_option("-m,--model",
		                mModelFile,
		                "SVM's model file")
				->required()->check(CLI::ExistingFile);

		mApp.add_option("-l,--lower",
		                mLowerKmerSize,
		                "Lower k-mer size (mininum = " + std::to_string(Globals::MIN_KMER_SIZE) + " )")
				->required();

		mApp.add_option("-u,--upper",
		                mUpperKmerSize,
		                "Upper k-mer size (maximum = " + std::to_string(Globals::MAX_KMER_SIZE) + " )")
				->required();

		mApp.add_option("-t,--threads",
		                mNumThreads,
		                "Number of threads (0 = for automatic determination (default). Maximum = " +
				                std::to_string(Globals::MAX_NUM_THREADS) + ").")
				->required()
				->check([&] (std::string numThreads) {

					auto num = std::stoi(numThreads);

					if (num <= 0 || num > Globals::MAX_NUM_THREADS)
						mNumThreads = -1;
					else
						mNumThreads = num;

					return true;
				});

		mApp.add_option("-w,--write", mWritePreds,
		                "Write predicteds k-mers (0 = none, 1 = amps, 2 = amps, 3 = both; default = none");

        mApp.add_flag("-a,--aware", mAware, "Enable aware mode (low-memory consumption; default false)");

		mApp.add_flag("-v,--verbose", mVerbose, "Enable verbose mode (show extra information; default false)");


	}

private:
	CLI::App mApp {Globals::APP_NAME};
	std::string mInputFile;
	std::string mOutputBaseName;
	std::string mScalingFile;
	std::string mModelFile;
	uint mLowerKmerSize;
	uint mUpperKmerSize;
	int mNumThreads = 0;
	uint mWritePreds = Globals::WRITE_NONE_PREDS;
    bool mAware = false;
	bool mVerbose = false;

};


#endif //INPROT_TEST_CLI_H
