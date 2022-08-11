//
// Created by germelcar on 9/16/17.
//

#ifndef INPROT_KMER_MANAGER_H
#define INPROT_KMER_MANAGER_H

#include <vector>
#include <unordered_map>
#include <memory>
#include "fasta_utils.h"
#include "fasta_seq.h"
#include "svm_scaling.h"
#include "group_koff.h"
#include "rang.hpp"
#include <tbb/tbb.h>
#include <tbb/concurrent_unordered_map.h>
#include <mutex>
#include <random>


template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;
using namespace rang;

namespace fasta {


	bool operator==(const FastaSeq& rhs, const FastaSeq& lhs) {
        return &rhs == &lhs;
	}

	class KmersManager {

	public:
		KmersManager(const std::string& inFileName, const std::string& outFileName,
		             uint lowerKSize, uint upperKSize, uint writePreds, bool awareMode, bool verbose):
				mInFileName(inFileName), mOutFileName(outFileName), mLowerKSize(lowerKSize), mUpperKSize(upperKSize),
				mWritePreds(writePreds), mTmpFiles(mUpperKSize - mLowerKSize + 1),
				mKmersMap(mUpperKSize - mLowerKSize + 1), mAwareMode(awareMode), mVerbose(verbose) { }


		template<typename T, typename M, EnableIf<std::is_floating_point<T>>...>
		std::pair<bool, std::string> extract(SvmScaling<T>& scaling, M& model) noexcept {


			if (mFseqs.empty())
				mFseqs = FastaUtils::readFasta(mInFileName);

            if (mVerbose) {
                std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
                          << "A total of " << mFseqs.size() << " sequences" << style::reset << std::endl;
            }

			// Memory mode: AWARE
			// Low memory consumption
			if (mAwareMode) { // ---------- Aware mode ----------

				auto basename = generateTmpBaseName();
                std::string error;
                bool allOK = true;

				for (auto i = mLowerKSize; i <= mUpperKSize; ++i)
                    mTmpFiles[i] = std::to_string(i) + "_" + basename;


				// For each k-mer size (one by one):
				// 1.- Extract the unique k-mers
				// 2.- Evaluate the k-mers (predict against the SVM)
				// 3.- Add if the k-mer is predicted as AMP
				for (auto i = mLowerKSize; i <= mUpperKSize; ++i) {

					if (mVerbose) {
						std::cout << style::bold << fg::blue << "[INFO] " << style::reset
						          << fg::blue << "Extracting unique " << i << "-mers" << std::endl;
					}

					// Extract the uniques k-mers of size "i"
					auto kmers = FastaUtils::uniqKmers(mFseqs, i);

					if (mVerbose) {
						std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
						          << "A total of " << kmers.size() << " unique " << i
						          << "-mers" << std::endl;


						std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
						          << "Evaluating " << kmers.size() << " unique " << i << "-mers" << std::endl;
					}

					// Evaluate k-mers and add to "mKmersMap" those k-mers predicted as AMP
					tbb::parallel_for_each(kmers.begin(), kmers.end(), [&](auto &koff) {
						koff.evaluate(scaling, model);
					});

					// Sort the k-mers by AMP activity
					// The non-AMPs will be at the end
                    tbb::parallel_sort(kmers.begin(), kmers.end(), [] (const auto& ki, const auto& kj) {
						return ki.isAMP() && !kj.isAMP();
					});

					// Write predicteds (AMPs, NAMPs, boths) k-mers in multifasta format
					if (mWritePreds != Globals::WRITE_NONE_PREDS)
						auto okErr = writePreds(kmers);

					// As the k-mers are sorted by first the AMPs and later the Non-AMPs, then, if we find
					// the first non-AMP, we know that at that point, all the k-mers on the left are AMPs
					auto lastNonAMP = std::find_if(kmers.cbegin(), kmers.cend(), [](const auto &k) {
						return k.isAMP() == false;
					});

					// Calculate the total de AMPs
					size_t totalAMPs {0};
					if (lastNonAMP == kmers.cend())
						totalAMPs = kmers.size();
					else
						totalAMPs = static_cast<size_t>(std::distance(kmers.cbegin(), lastNonAMP));

					// Get the k-mer temp filename and write out
					auto fname = mTmpFiles[i];

					if (mVerbose) {

						std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
						          << "Writing " << totalAMPs << " " << i << "-mers predicted as AMPs (of "
						          << kmers.size() << " uniques" << ") to file: " << fname << std::endl;
					}

					std::ofstream fout(fname, std::ios_base::out | std::ios_base::trunc);

					// Check if creating the file was done with success
					// TODO Return if error occurred and information about the error if no exception is thrown
					if (!fout) {
						allOK = false;
                        error = "Failed to open file: " + fname;
                        break;
					}

					// Insert k-mers information as:
					// Fasta sequence index (starting from zero) [space]
					// K-mer offset (start position) [space]
					// k-mer size [end line]
                    for (size_t j {0}; j < totalAMPs; ++j) {

						const auto& koff = kmers[j];
                        const auto fsItr = std::find_if(mFseqs.cbegin(), mFseqs.cend(), [&](const auto &fs) {
                            return fs == koff.getFastaSeq();
                        });

                        // Get the fasta sequence index from the vector of fasta sequences
                        auto fsIdx = std::distance(mFseqs.cbegin(), fsItr);

                        // Write out the fasta sequence index, k-mers offset and k-mer size
                        fout << fsIdx << " "
                             << koff.getOffset() << " "
                             << koff.getSize() << "\n";


                        if (fout.fail() || fout.bad()) {
                            allOK = false;
                            error = "Error while writing k-mer (offset=" + std::to_string(koff.getOffset())  + ";" +
                                    "size=" + std::to_string(koff.getSize()) + ") for fasta sequence " +
                                    mFseqs[fsIdx].getDesc() + " (index=" + std::to_string(fsIdx)
                                    + ") in file: " + fname;
                            break;
                        }

                    }

                    if (!allOK)
                        break;

					fout.flush();

				} // End of for (lower..upper)

                return std::make_pair(allOK, error);

			} else { // ---------- Normal mode ----------

				// For each k-mer size (one by one):
				// 1.- Extract the unique k-mers
				// 2.- Evaluate the k-mers (predict against the SVM)
				// 3.- Add if the k-mer is predicted as AMP
				for (auto i = mLowerKSize; i <= mUpperKSize; ++i) {

					if (mVerbose) {
						std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue <<
						          "Extracting unique " << i << "-mers" << std::endl;
					}

					// Extract the uniques k-mers of size "i"
					auto kmers = FastaUtils::uniqKmers(mFseqs, i);

					if (mVerbose) {

						std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
						          << "A total of " << kmers.size() << " unique " << i << "-mers" << std::endl;

						std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
						          << "Evaluating " << kmers.size() << " unique " << i << "-mers" << std::endl;
					}

					// Evaluate k-mers and add to "mKmersMap" those k-mers predicted as AMP
					tbb::parallel_for_each(kmers.begin(), kmers.end(), [&] (auto& koff) {
						koff.evaluate(scaling, model);
					});

					tbb::parallel_sort(kmers.begin(), kmers.end(), [] (const auto& ki, const auto& kj) {
						return ki.isAMP() && !kj.isAMP();
					});

					auto lastNAMP = std::find_if(kmers.cbegin(), kmers.cend(), [] (const auto& koff) {
						return koff.isAMP() == false;
					});

					// Write predicteds (AMPs, NAMPs, boths) k-mers in multifasta format
					if (mWritePreds != Globals::WRITE_NONE_PREDS)
						auto okErr = writePreds(kmers);

					// Calculate the total de AMPs
					size_t totalAMPs {0};
					if (lastNAMP == kmers.cend())
						totalAMPs = kmers.size();
					else
						totalAMPs = static_cast<size_t>(std::distance(kmers.cbegin(), lastNAMP));


                    // Reserve space for all AMPs and insert them
                    auto& kmersSize = mKmersMap[i];
                    kmersSize.reserve(totalAMPs);

                    tbb::parallel_for_each(kmers.begin(), kmers.end(), [&] (auto& koff) {
                        if (koff.isAMP())
                            kmersSize.emplace_back(koff);
                    });


					// Shrink the vector for reduce memory usage
                    kmersSize.shrink_to_fit();

					if (mVerbose) {

						std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
						          <<  mKmersMap[i].size() << " " << i << "-mers predicted as AMPs (of "
						          << kmers.size() << " uniques)" << std::endl;
					}

				} // End of for (lower..upper)


			} // End of else {


			return std::make_pair(true, std::string()); // ok, errors

		} // End of extract(...)

		bool shrinkProteome() {
			if (mKmersMap.empty() && !mAwareMode)
				return false;

			std::ofstream outFile(mOutFileName, std::ios_base::trunc | std::ios_base::out);

			// If memory save is enabled, then, for each sequence,
			// reduce and write the reduced sequence
			if (mAwareMode) {

				if (mVerbose) {

					size_t numSeq{0};
                    tbb::concurrent_vector<GroupKoff> groups;

					std::for_each(mFseqs.cbegin(), mFseqs.cend(), [&](const auto &fs) {

						tbb::concurrent_vector<std::shared_ptr<KmerOffset>> koffs_fs;

						tbb::parallel_for_each(mTmpFiles.begin(), mTmpFiles.end(), [&](auto &pair) {
							koffFromFsFile(fs, koffs_fs, pair.second);
						});

						reduceKoffs(fs, koffs_fs, groups);
						groups.shrink_to_fit();

                        std::cout << style::bold << fg::blue << "[" << ++numSeq << " / " << mFseqs.size()
                                  << "] " << style::reset << fg::blue << "sequences shrinked\r";

					}); // End of for_each (mFseqs...)

                    std::cout << std::endl << style::bold << fg::blue << "[INFO] " << style::reset
                              << fg::blue << "Writing shrinked sequences" << style::reset << std::endl;

                    // Get the last unique group
                    auto totUniqGroups = reduceGroups(groups);

                    // Write out those overlapped (grouped) k-mers
                    writeGroups(outFile, groups, totUniqGroups);

					std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
					          << "Written a total of " << totUniqGroups << " new sequences "
					          << "to file: " << mOutFileName << std::endl;

				} else {

                    tbb::concurrent_vector<GroupKoff> groups;

					std::for_each(mFseqs.cbegin(), mFseqs.cend(), [&](const auto &fs) {

						tbb::concurrent_vector<std::shared_ptr<KmerOffset>> koffs_fs;

						tbb::parallel_for_each(mTmpFiles.begin(), mTmpFiles.end(), [&](auto &pair) {
							koffFromFsFile(fs, koffs_fs, pair.second);
						});

						koffs_fs.shrink_to_fit();
						reduceKoffs(fs, koffs_fs, groups);
						groups.shrink_to_fit();

					}); // End of for_each (mFseqs...)

                    // Get the last unique group
                    auto totUniqGroups = reduceGroups(groups);

                    // Write out those overlapped (grouped) k-mers
                    writeGroups(outFile, groups, totUniqGroups);

				}

				std::cout << std::endl;

				if (mVerbose) {
					std::cout << style::bold << fg::green << "[STATUS] " << style::reset << fg::green
					          << "Removing temporary files..." << std::endl;
				}

				// Removing temporary files
				std::for_each(mTmpFiles.cbegin(), mTmpFiles.cend(), [] (const auto& pair) {
					std::remove(pair.second.c_str());
				});


			} else { // ---------- Normal mode ----------

				// For each sequence:
				// 1.- Extract the k-mers from that sequence.
				// 2.- Group those overlapped k-mers.
				// 3.- Write out those overlapped k-mers.
				//
				// Perform the three steps above for each sequence, one at a time
				if (mVerbose) {

					size_t numSeq {0};
                    tbb::concurrent_vector<GroupKoff> groups;

					std::for_each(mFseqs.begin(), mFseqs.end(), [&] (const auto& fs) {

						auto koffs = koffFromSeq(fs);
                        reduceKoffs(fs, koffs, groups);

                        std::cout << style::bold << fg::blue << "[" << ++numSeq << " / " << mFseqs.size()
                                  << "] " << style::reset << fg::blue << "sequences shrinked\r";

					});

                    std::cout << std::endl << style::bold << fg::blue << "[INFO] " << style::reset
                              << fg::blue << "Writing shrinked sequences" << style::reset << std::endl;

                    // Get the last unique group
                    auto totUniqGroups = reduceGroups(groups);

                    // Write out those overlapped (grouped) k-mers
                    writeGroups(outFile, groups, totUniqGroups);

					std::cout << style::bold << fg::blue << "[INFO] " << style::reset
					          << fg::blue << "Written a total of " << totUniqGroups << " new sequences "
					          << "to file: " << mOutFileName << std::endl;

				} else {

                    tbb::concurrent_vector<GroupKoff> groups;

					std::for_each(mFseqs.begin(), mFseqs.end(), [&] (const auto& fs) {

						auto koffs = koffFromSeq(fs);
						reduceKoffs(fs, koffs, groups);

					});

                    // Get the last unique group
                    auto totUniqGroups = reduceGroups(groups);

                    // Write out those overlapped (grouped) k-mers
                    writeGroups(outFile, groups, totUniqGroups);

				}

				std::cout << std::endl;

			}

			outFile.close();
			return true; // all OK
		}


	private:
		tbb::concurrent_vector<std::shared_ptr<KmerOffset>> koffFromSeq(const FastaSeq& fs) {

			// Initialize the vector with reference to k-mers for each sequences (fs)
			tbb::concurrent_vector<std::shared_ptr<KmerOffset>> refKmers;
			refKmers.reserve(1'000);

			// For each k-mer, search if the given k-mer, k, was extracted from the sequence "fs"
			// If so, then add to the "refKmers" vector for later add the vector to the
			// "shrinked" unordered_map.
			tbb::parallel_for_each(mKmersMap.begin(), mKmersMap.end(), [&] (const auto& pair) {

				// Get the k-mers of size "pair.first"
				const auto& kmers = pair.second;

				// For earch k-mer, check if it is a reference to the given sequence, fs
				tbb::parallel_for_each(kmers.cbegin(), kmers.cend(), [&] (const auto& koff) {

					// Add those k-mers that were not extracted directly from the seq `fs` but they are present
					// in the sequence (repeateds in other sequences).
					//
					// For example, the k-mers ki could be present in seq1 and seq2, but only one instance
					// will be kept when the duplicates were removed when obtaining the uniques
					// in order to avoid redundant calculation of the molecular descriptors and
					// evaluation of the SVM
					size_t pos {0};
					const auto& kmer = koff.getKmer();

					while ((pos = fs.getSeq().find(kmer, pos)) != std::string::npos) {
						refKmers.emplace_back(std::make_shared<KmerOffset>(fs, pos, koff.getSize(), true));
						++pos;
					}


				});

			});

			refKmers.shrink_to_fit();
			return refKmers;
		}

		void koffFromFsFile(const FastaSeq& fs,
		                    tbb::concurrent_vector<std::shared_ptr<KmerOffset>>& refKmers,
		                    const std::string& filename) {


			std::ifstream inFile(filename);

			if (!inFile) {
				std::cerr << style::bold << fg::yellow << "[WARNING] " << style::reset << fg::yellow
				          << "To run-out from koffFromFsFile for seq: " << fs.getDesc()
				          << "\tand file: " << filename << style::reset << std::endl;

				return;
			}


			unsigned long fsIdx {0};
			size_t koff_offset {0};
			uint koff_size {0};

			while (inFile >> fsIdx >> koff_offset >> koff_size) {

				// Add those k-mers that were not extracted directly from the seq `fs` but they are present
				// in the sequence (repeateds in other sequences).
				//
				// For example, the k-mers ki could be present in seq1 and seq2, but only one instance
				// will be kept when the duplicates were removed when obtaining the uniques
				// in order to avoid redundant calculation of the molecular descriptors and
				// evaluation of the SVM
				size_t pos {0};
				KmerOffset koff(mFseqs[fsIdx], koff_offset, koff_size, true);
				const auto& kmer = koff.getKmer();

				while ((pos = fs.getSeq().find(kmer, pos)) != std::string::npos) {

					refKmers.emplace_back(std::make_shared<KmerOffset>(fs, pos, koff.getSize(), true));
					++pos;
				}

				inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

			} // End of while (inFile >>...)

			inFile.close();

		}

		void reduceKoffs(const FastaSeq& fs, tbb::concurrent_vector<std::shared_ptr<KmerOffset>>& kmers,
		                 tbb::concurrent_vector<GroupKoff>& groups) {

			if (kmers.empty())
				return;

			// Sort k-mers by start position (offset)
			tbb::parallel_sort(kmers.begin(), kmers.end(), [&] (const auto& ki, const auto& kj) {
				return ki->getOffset() < kj->getOffset();
			});


			auto totKmers = kmers.size();
			std::shared_ptr<KmerOffset> kprev = kmers[0];
			std::shared_ptr<KmerOffset> klast = kmers[0];

			// first    => start/offset
			// second   => end
			std::pair<size_t, size_t> kgc = std::make_pair(kprev->getOffset(), kprev->getSize());
			size_t j = 1;

			while (j < totKmers) {
				kgc.second = klast->getEnd();
				std::shared_ptr<KmerOffset> kj = kmers[j];

				if (!kj->isInside(kgc)) {

					if (kj->isConnected(kgc) || kj->intersect(kgc)) {
						klast = kj;
					} else {
						//groups.emplace_back(kgc.first, kgc.second);
						groups.emplace_back(fs, kgc.first, kgc.second);
						kprev = kj;
						klast = kj;
						kgc.first = kprev->getOffset();
					}
				}

				j++;
			}

			if (kprev->getOffset() != klast->getOffset() && kprev->getEnd() != klast->getEnd())
				groups.emplace_back(fs, kprev->getOffset(), klast->getEnd());
			else
				groups.emplace_back(fs, klast->getOffset(), klast->getEnd());

		}

		size_t reduceGroups(tbb::concurrent_vector<GroupKoff>& groups) {

			// Sort by k-mer content (sequence)
			tbb::parallel_sort(groups);

			// Get the iterator ("index") past the end of the last unique group
			auto last = std::unique(groups.begin(), groups.end());

			// Get the distance in order to calculate the total of uniques groups
			return static_cast<size_t>(std::distance(groups.begin(), last));
		}

		void writeGroups(std::ofstream& fout, tbb::concurrent_vector<GroupKoff>& groups, size_t totUniqs) const {

			if (groups.empty())
				return;

			for (size_t i {0}; i < totUniqs; ++i) {

				const auto& g = groups[i];
				const auto& fs = g.getFastaSeq();

				const auto& desc = fs.getDesc();
				const auto& fseq = fs.getSeq();

				const auto begin = g.getBegin();
				const auto end = g.getEnd() > 0 ? g.getEnd() - 1 : g.getEnd();

				const auto len = g.getEnd() - begin;
				const auto header = ">" + desc + "_" + std::to_string(begin) + "_" + std::to_string(end);

				fout << header << "\n" << fseq.substr(begin, len) << "\n";
			}

		}

		std::string generateTmpBaseName() const {
			std::random_device rd;
			std::uniform_int_distribution<int> dist;
			std::uniform_int_distribution<int> dist2 (0, static_cast<int>(mOutFileName.length() - 1));

			return std::to_string(
					((dist(rd) << 3) ^ (std::hash<std::string>()(mOutFileName) >> 5))  ^ (dist2(rd) >> 7));

		}

		std::pair<bool, std::string> writePreds(tbb::concurrent_vector<KmerOffset>& kmers) {

			bool allOK = true;
			std::string error;

			if (kmers.empty())
				return std::make_pair(false, "Error: empty kmers to write");

			auto ksizeStr = std::to_string(kmers[0].size());
			auto lastDotIdx = mOutFileName.find_last_of('.');
			std::string fileBaseName;
			std::string toWriteStr;

			if (lastDotIdx == std::string::npos)
				fileBaseName = mOutFileName;
			else
				fileBaseName = mOutFileName.substr(0, lastDotIdx);


			if (mWritePreds == Globals::WRITE_BOTHS_PREDS) {

				if (mVerbose) {
					std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
					          << "Writing " + ksizeStr + "-mers predicteds as AMPs" << std::endl;
				}

				auto okErrAMPs = writePreds(kmers, ksizeStr, fileBaseName, Globals::WRITE_AMPS_PREDS);

				if (!okErrAMPs.first) {
					std::cout << style::bold << fg::yellow << "[WARNING] " << style::reset << fg::yellow
					          << "Error while writing " + ksizeStr + "-mers predicteds as AMPs" << std::endl;

					return okErrAMPs;
				}

				if (mVerbose) {
					std::cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
					          << "Writing " + ksizeStr + "-mers predicteds as No-AMPs" << std::endl;
				}

				auto okErrNAMPs = writePreds(kmers, ksizeStr, fileBaseName, Globals::WRITE_NAMPS_PREDS);

				if (!okErrNAMPs.first) {
					std::cout << style::bold << fg::yellow << "[WARNING] " << style::reset << fg::yellow
					          << "Error while writing " + ksizeStr + "-mers predicteds as No-AMPs" << std::endl;

					return okErrNAMPs;
				}
			}

			if (mWritePreds == Globals::WRITE_AMPS_PREDS) {

				if (mVerbose) {
					std::cout << std::endl << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
					          << "Writing " + ksizeStr + "-mers predicteds as AMPs" << std::endl;
				}

				auto okErr = writePreds(kmers, ksizeStr, fileBaseName, Globals::WRITE_AMPS_PREDS);

				if (!okErr.first) {
					std::cout << style::bold << fg::yellow << "[WARNING] " << style::reset << fg::yellow
					          << "Error while writing " + ksizeStr + "-mers predicteds as AMPs" << std::endl;

					return okErr;
				}
			}

			if (mWritePreds == Globals::WRITE_NAMPS_PREDS) {

				if (mVerbose) {
					std::cout << std::endl << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
					          << "Writing " + ksizeStr + "-mers predicteds as No-AMPs" << std::endl;
				}

				auto okErr = writePreds(kmers, ksizeStr, fileBaseName, Globals::WRITE_NAMPS_PREDS);

				if (!okErr.first) {
					std::cout << style::bold << fg::yellow << "[WARNING] " << style::reset << fg::yellow
					          << "Error while writing " + ksizeStr + "-mers predicteds as No-AMPs" << std::endl;
					return okErr;
				}
			}

			return std::make_pair(allOK, error);
		}

		std::pair<bool, std::string> writePreds(tbb::concurrent_vector<KmerOffset> &kmers,
		                                        std::string& ksizeStr,
		                                        std::string &fileBaseName, uint toWrite) {

			bool allOK = true;
			std::string error;
			std::string fname;
			std::string msgErr;
			tbb::concurrent_vector<KmerOffset>::const_iterator beginItr;
			tbb::concurrent_vector<KmerOffset>::const_iterator endItr;

			if (toWrite == Globals::WRITE_AMPS_PREDS) {
				fname = fileBaseName + "_" + ksizeStr + "-mers_amps.fasta";
				msgErr = "as AMP to file: ";

				beginItr = kmers.cbegin();
				endItr = std::find_if(kmers.cbegin(), kmers.cend(), [] (const auto& koff) {
					return koff.isAMP() == false;
				});
			}

			if (toWrite == Globals::WRITE_NAMPS_PREDS) {
				fname = fileBaseName + "_" + ksizeStr + "-mers_namps.fasta";
				msgErr = "as No-AMP to file: ";

				beginItr = std::find_if(kmers.cbegin(), kmers.cend(), [] (const auto& koff) {
					return koff.isAMP() == false;
				});

				endItr = kmers.cend();
			}

			std::ofstream fout (fname, std::ios_base::trunc | std::ios_base::out);

			for (auto i = beginItr; i < endItr; ++i) {
				const auto& koff = *i;
				auto header = ">" + koff.getFastaSeq().getDesc() + "_" +
						std::to_string(koff.getOffset()) + "_" + std::to_string(koff.getEnd() - 1) + "\n";

				auto seq = koff.getKmer() + "\n";

				fout << header << seq;

				if (fout.bad() || fout.fail()) {
					allOK = false;
					error = "Unable to write " + ksizeStr + "-mers predicteds " +  msgErr + fname;
					break;
				}
			}

			fout.flush();
			fout.close();

			return std::make_pair(allOK, error);

		}

        // TODO Check data alignment variables declaration
		std::string mInFileName;
		std::string mOutFileName;
		uint mLowerKSize;
		uint mUpperKSize;
		uint mWritePreds;
		tbb::concurrent_vector<FastaSeq> mFseqs;
		tbb::concurrent_unordered_map<uint, std::string> mTmpFiles;
		tbb::concurrent_unordered_map<uint, tbb::concurrent_vector<KmerOffset>> mKmersMap;
		bool mAwareMode;
		bool mVerbose;

	};

}



#endif //INPROT_KMER_MANAGER_H
