//
// Created by germelcar on 10/30/17.
//

#ifndef INPROT_GROUP_KOFF_H
#define INPROT_GROUP_KOFF_H

#include "fasta_seq.h"

namespace fasta {

    class GroupKoff {

    public:

        explicit GroupKoff(const FastaSeq& fs, size_t begin, size_t end): mFseq(fs), mBegin(begin), mEnd(end) { }

        bool operator<(const GroupKoff& rhs) const {
            return getKmmer() < rhs.getKmmer();
        }

        bool operator==(const GroupKoff& rhs) const {
            return getKmmer() == rhs.getKmmer();
        }

        std::string getKmmer() const {
            return mFseq.get().getSeq().substr(mBegin, mEnd - mBegin);
        }

        const FastaSeq& getFastaSeq() const {
            return mFseq.get();
        }

        size_t getBegin() const {
            return mBegin;
        }

        size_t getEnd() const {
            return mEnd;
        }

    private:
        std::reference_wrapper<const FastaSeq> mFseq;
        size_t mBegin;
        size_t mEnd;

    };

}

#endif //INPROT_GROUP_KOFF_H
