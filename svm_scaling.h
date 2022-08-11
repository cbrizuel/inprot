//
// Created by germelcar on 9/14/17.
//

#ifndef INPROT_SVM_SCALING_H
#define INPROT_SVM_SCALING_H

#include <valarray>
#include <string>
#include "globals.h"
#include <vector>
#include <utility>

template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;


template<typename T, EnableIf<std::is_floating_point<T>>...>
class SvmScaling {

public:
	SvmScaling() {
		mBounds.reserve(Globals::NUM_MDS);
	}

	bool restore(const std::string& restFile) {
		std::ifstream inFile(restFile, std::ios_base::in | std::ios_base::binary);

		if (!inFile.good())
			return false;

		try {

			if (inFile.peek() == 'y') {
				// We read three lines for the "y" scaling parameters. Example of file with
				// y scaling:
				//
				// y        --> first line of "y" scaling
				// 0 1      --> second line of "y" scaling
				// 0 1      --> third line of "y" scaling
				// x        --> first line of "x" scaling. Read after this "if condition"
				// -1 1     --> second line of "x" scaling. Read after this "if condition"
				// 1 0 297.05
				// 2 -4.555206 581.0731
				// 3 -0.7524385 0.7170606
				// 4 8.157474000000001 180
				//
				inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}

			if (inFile.peek() == 'x') {

				// Skip/ignore the first line of "x" scaling
				inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				T lowX{0};
				T upperX{0};

				if (!(inFile >> lowX >> upperX))
					return false;

				if (lowX != Globals::SVM_NEGATIVE_LABEL || upperX != Globals::SVM_POSITIVE_LABEL)
					return false;

				// Skip until the next line (where the scaling paramters begins)
				inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

				int idx {0};
				int prevIdx {0};
				T minVal {0};
				T maxVal {0};

				while (inFile >> idx >> minVal >> maxVal &&
				       inFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n')) {

					if (std::abs(prevIdx - idx) != 1)
						return false;

					mBounds.emplace_back(minVal, maxVal);
					prevIdx = idx;

				}

			}

			if (!inFile.eof() && inFile.fail())
				return false;

			inFile.close();

		} catch (std::bad_alloc& ex) {
			std::cerr << "Unable to allocate memory" << std::endl;
			return false;
		}


		return true;
	}

	bool scale(std::valarray<T>& mds) const {
		if (mds.size() != mBounds.size())
			return false;

		auto lower = Globals::SVM_NEGATIVE_LABEL;
		auto upper = Globals::SVM_POSITIVE_LABEL;

		for (size_t i = 0; i < mBounds.size(); ++i) {

			const std::pair<T, T>& pair = mBounds[i];
			const auto curr_value = mds[i];

			if (pair.first == pair.second)
				continue;

			if (pair.first == curr_value)
				mds[i] = lower;
			else if (pair.second == curr_value)
				mds[i] = upper;
			else
				mds[i] = lower + (upper - lower) * (curr_value - pair.first) / (pair.second - pair.first);
		}

		return true;
	}


private:
	std::vector<std::pair<T, T>> mBounds;

};



#endif //INPROT_SVM_SCALING_H
