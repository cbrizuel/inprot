//
// Created by germelcar on 8/30/17.
//

#ifndef INPROT_FASTA_UTILS_H
#define INPROT_FASTA_UTILS_H


#include <vector>
#include "fasta_seq.h"
#include <fstream>
#include <sstream>
#include "globals.h"
#include "kmer_offset.h"
#include <algorithm>
#include <numeric>
#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>


namespace fasta {

    class FastaUtils {

    public:
        static tbb::concurrent_vector<FastaSeq> readFasta(const std::string& filename) {
            std::ifstream inFile(filename, std::ios_base::in | std::ios_base::binary);

            if (!inFile) {
	            return tbb::concurrent_vector<FastaSeq>();
            }

	        // Create vector of fasta sequences with an initial space for 100,000 fasta sequences
	        tbb::concurrent_vector<FastaSeq> fseqs;
            fseqs.reserve(20'000);
            std::string line;
            std::string seq;
            std::string desc;
            size_t numLine {1};

            while(std::getline(inFile, line)) {
                // If the first character is a space, then, "space_idx" variable will be greater
                // than zero
                auto space_idx = line.find_first_not_of(' ');

                if (space_idx > 0 && line.length() > 0) {
                    throw std::runtime_error("Line " + std::to_string(numLine) + " has invalid characters");
                }

                // If the line begins with ">" character, then, is a description line
                if (line[0] == '>') {
                    space_idx = line.find_first_of(' ');

                    if (desc.length() > 0) {

                        // Shrink to reduce memory usage
                        seq.shrink_to_fit();
                        desc.shrink_to_fit();

                        // Insert new fasta sequence
                        fseqs.emplace_back(seq, desc);
                        desc.clear();
                        seq.clear();
                    }

                    desc += line.substr(1, space_idx == std::string::npos ? line.length() : space_idx);

                } else {

                    // Convert sequence line to upper-case
                    std::transform(line.begin(), line.end(), line.begin(), ::toupper);

                    // Check if only contains valid standard aminoacids
                    if (!isValid(line)) {
                        throw std::runtime_error(
                                "Line " + std::to_string(numLine) + " has invalid characters");
                    }

                    seq += line;
                }

                ++numLine;
            }

            fseqs.emplace_back(seq, desc);
            tbb::parallel_sort(fseqs);

            auto lastUnique = std::unique(fseqs.begin(), fseqs.end(), [] (const auto& fi, const auto& fj) {
                return &fi == &fj || fi.getSeq() == fj.getSeq();
            });

            if (lastUnique != fseqs.end()) {
                auto totUniques = static_cast<size_t>(std::distance(fseqs.begin(), lastUnique));
                auto totDuplicates = fseqs.size() - totUniques;

                std::cout << style::bold << fg::yellow << "[WARNING] " << style::reset << fg::yellow
                          << "A total of " << fseqs.size() << " sequences read with " << totUniques
                          << " uniques (" << totDuplicates << " duplicates)" << style::reset << std::endl;

                tbb::concurrent_vector<FastaSeq> uniqs;
                uniqs.reserve(static_cast<size_t>(totUniques));

                tbb::parallel_for(tbb::blocked_range<size_t>(0, totUniques), [&] (const auto& r) {

                    for (auto i = r.begin(); i != r.end(); ++i)
                        uniqs.emplace_back(fseqs[i]);
                });

                fseqs.clear();
                uniqs.shrink_to_fit();
                return uniqs;

            }

	        fseqs.shrink_to_fit();
            return fseqs;
        }

		static tbb::concurrent_vector<KmerOffset> uniqKmers(
                const tbb::concurrent_vector<FastaSeq>& fseqs, uint32_t ksize) {

			auto total = totalKmers(fseqs, ksize);

			if (total == 0)
				return tbb::concurrent_vector<KmerOffset>();

			tbb::concurrent_vector<KmerOffset> totKmers;
			totKmers.reserve(total);

			tbb::parallel_for_each(fseqs.cbegin(), fseqs.cend(), [&] (const FastaSeq& fs) {

				if (ksize > fs.length())
					return;

				size_t startKmer {0};
				size_t stopKmer {fs.length() - ksize + 1};

				while (startKmer < stopKmer)
					totKmers.emplace_back(fs, startKmer++, ksize);

			});

			tbb::parallel_sort(totKmers);
            auto lastUnique = std::unique(totKmers.begin(), totKmers.end());
            auto totUniques = static_cast<size_t>(std::distance(totKmers.begin(), lastUnique));

            tbb::concurrent_vector<KmerOffset> uniqKmers;
            uniqKmers.reserve(totUniques);

            tbb::parallel_for_each(totKmers.begin(), lastUnique, [&] (auto& koff) {
                uniqKmers.emplace_back(koff);
            });

            totKmers.clear();
            uniqKmers.shrink_to_fit();

			return uniqKmers;
		}

    private:
        static bool isValid(const std::string& seq) {
		    for (const auto& c : seq) {
			    if (Globals::ALPHABET.find_first_of(c) == std::string::npos)
				    return false;
		    }

            return true; // all OK
        }

	    static size_t totalKmers(const tbb::concurrent_vector<FastaSeq>& fseqs, size_t ksize) {
		    size_t total {0};

            std::for_each(fseqs.cbegin(), fseqs.cend(), [&total, &ksize] (const FastaSeq& fs) {
                total += (ksize > fs.length() ? 0 : fs.length() - ksize + 1);
            });

		    return total;
	    }

    };

}



#endif //INPROT_FASTA_UTILS_H
