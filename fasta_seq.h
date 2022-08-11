//
// Created by germelcar on 8/30/17.
//

#ifndef INPROT_FASTASEQ_H
#define INPROT_FASTASEQ_H

#include <string>
#include <utility>
#include <tbb/tbb.h>

namespace fasta {



    class FastaSeq {

    public:
	    //
	    // Constructors & destructors
	    //
        FastaSeq() = default;
        explicit FastaSeq(std::string seq, std::string desc): mSeq(std::move(seq)), mDesc(std::move(desc)) {}

	    //
	    // Methods
	    //
		bool operator<(const FastaSeq& rhs) const {
			return mSeq < rhs.mSeq;
		}

	    size_t length() const {
		    return mSeq.length();
	    }

	    //
	    // Getters & setters
	    //
	    const std::string& getSeq() const {
		    return mSeq;
	    }

	    const std::string& getDesc() const {
		    return mDesc;
	    }

    private:
        std::string mSeq;
        std::string mDesc;

    };

}


#endif //INPROT_FASTASEQ_H
