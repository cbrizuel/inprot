//
// Created by germelcar on 9/14/17.
//

#ifndef INPROT_SVM_MODEL_H
#define INPROT_SVM_MODEL_H

#include <fstream>
#include <sstream>
#include <string>
#include <functional>
#include <vector>
#include <algorithm>
#include "libsvm.h"
#include "globals.h"



template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;



namespace libsvm {

	template<typename T, EnableIf<std::is_floating_point<T>>...>
	class SvmModel {

	public:
		SvmModel() = default;

		bool load(const std::string& filename) {
			std::ifstream inFile(filename, std::ios_base::in | std::ios_base::binary);
			readHeader(inFile);

			auto m = static_cast<size_t>(mNrClass - 1);
			auto ml = static_cast<size_t>(mL);

			// Reserve space for SV's coefs and zero-fill them
			mSVCoef.reserve(m);
			for (size_t i = 0; i < m; ++i)
				mSVCoef.emplace_back(0, ml);

			// Reserve space for SVs and zero-fill them
			// Each SV will be an array of size = Globals::NUM_MDS
			// filled with zeroes
			mSV.reserve(ml);
			for (size_t i = 0; i < ml; ++i)
				mSV.emplace_back(0, Globals::NUM_MDS);

			std::string line;
			size_t numLine {0};
			while (numLine < ml) {

				if (!std::getline(inFile, line)) {
					throw std::runtime_error(
							"Unable to read SV's line " + std::to_string(numLine + 1) +
									". Error occurred at line " + std::to_string(__LINE__) +
									" file " + __FILE__);
				}

				// Read first cofficient for the current SVs
				T coef {0};
				std::stringstream ss(line);
				ss >> coef;
				mSVCoef[0][numLine] = coef;

				// Read the rest of the coefficients
				for (size_t i = 1; i < m; ++i) {
					if (!(ss >> coef)) {
						throw std::runtime_error(
								"Unable to read SV's coefficient for SV's line: " + std::to_string(numLine + 1) +
										". Error occured at line " + std::to_string(__LINE__) +
										" file " + __FILE__);
					}

				}

				// Read the SVs
				unsigned idx {0};
				char separator {':'};
				T val {0};
				while (ss >> idx >> separator >> val) {

					if (idx < 1 || idx > Globals::NUM_MDS) {
						throw std::runtime_error(
								"SV's index is outside of range (#SV=" + std::to_string(numLine + 1) +
										"). Current index: " + std::to_string(idx) + ". Valid range: [1," +
										std::to_string(Globals::NUM_MDS) + "]");
					}

					mSV[numLine][--idx] = val;
				}

				numLine++;
			}

			return true; // all OK
		}

		int predict(std::valarray<T>& mds) const noexcept {
			std::vector<T> decValues;
			if (mParams->svm_type == SVM_TYPE::ONE_CLASS ||
					mParams->svm_type == SVM_TYPE::EPSILON_SVR ||
					mParams->svm_type == SVM_TYPE::NU_SVR)
				decValues.resize(1, 0);
			else
				decValues.resize(mNrClass * (mNrClass - 1) / 2, 0);

			if (mParams->svm_type == SVM_TYPE::ONE_CLASS ||
					mParams->svm_type == SVM_TYPE::EPSILON_SVR ||
					mParams->svm_type == SVM_TYPE::NU_SVR) {

				const auto& svCoef = mSVCoef[0];
				T sum {0};
				for (auto i = 0; i < mL; ++i)
					sum += svCoef[i] * 0;

				sum -= mRho[0];
				decValues[0] = sum;

				if (mParams->svm_type == SVM_TYPE::ONE_CLASS)
					return (sum > 0) ? 1 : -1;
				else
					return sum;

			} else {

				std::vector<T> kValues(mL, 0);
				for (auto i = 0; i < mL; ++i)
					kValues[i] = mKernelFnc(mds, mSV[i]);

				std::vector<int> start (mNrClass, 0);
				for (auto i = 1; i < mNrClass; ++i)
					start[i] = start[i - 1] + mNumSV[i - 1];

				std::vector<int> vote (mNrClass, 0);
				size_t p {0};
				for (auto i  = 0; i < mNrClass; ++i) {
					for (auto j  = i + 1; j < mNrClass; ++j) {

						T sum {0};
						auto si = start[i];
						auto sj = start[j];
						auto ci = mNumSV[i];
						auto cj = mNumSV[j];

						const auto& coef1 = mSVCoef[j - 1];
						const auto& coef2 = mSVCoef[i];

						for (size_t k = 0; k < ci; ++k) {
							sum += coef1[si + k] * kValues[si + k];
						}


						for (size_t k = 0; k < cj; ++k) {
							sum += coef2[sj + k] * kValues[sj + k];
						}

						sum -= mRho[p];
						decValues[p] = sum;

						if (decValues[p] > 0)
							++vote[i];
						else
							++vote[j];

						p++;

					} // End for (j..)

				} // End for (i..)

				auto voteMaxIdx = 0;
				for (auto i = 1; i < mNrClass; ++i)
					if (vote[i] > vote[voteMaxIdx])
						voteMaxIdx = i;


				return mLabel[voteMaxIdx];

			} // End else {..

		}

		constexpr T dotKernel(const std::valarray<T>& mds, const std::valarray<T>& y) const noexcept {
            T dot {0};
            for (size_t i = 0; i < mds.size(); ++i)
                dot += (mds[i] * y[i]);

            return dot;
		}

		constexpr T polyKernel(const std::valarray<T>& mds, const std::valarray<T>& y) const noexcept {
			return std::pow(mParams->gamma * dotKernel(mds, y) + mParams->coef0, mParams->degree);
		}

		constexpr T rbfKernel(const std::valarray<T>& mds, const std::valarray<T>& y) const noexcept {
			T rbf {0};
			for (size_t i = 0; i < mds.size(); ++i)
				rbf += std::pow(mds[i] - y[i], 2);

			return rbf;
		}

		constexpr T sigmoidKernel(const std::valarray<T>& mds, const std::valarray<T>& y) const noexcept {
			return std::tanh(mParams->gamma * dotKernel(mds, y) + mParams->coef0);
		}


	private:
		bool readHeader(std::ifstream& modelFile) {
			try {

				mParams = std::make_shared<SvmParameter<T>>();
				modelFile.clear();
				std::string line;

				while (std::getline(modelFile, line, ' ')) {

					if (line == "svm_type") {
						if (!std::getline(modelFile, line))
							return false;

						auto type = svm_type_table.find(line);
						if (type == svm_type_table.end())
							return false;

						mParams->svm_type = type->second;
					}

					if (line == "kernel_type") {
						if (!std::getline(modelFile, line))
							return false;

						auto type = kernel_type_table.find(line);
						if (type == kernel_type_table.end())
							return false;

						mParams->kernel_type = type->second;

						if (type->second == KERNEL_TYPE::LINEAR) {
							mKernelFnc = std::bind(&SvmModel::dotKernel, this,
							                       std::placeholders::_1,
							                       std::placeholders::_2);
						}

						if (type->second == KERNEL_TYPE::RBF) {
							mKernelFnc = std::bind(&SvmModel::rbfKernel, this,
							                      std::placeholders::_1,
							                      std::placeholders::_2);
						}

						if (type->second == KERNEL_TYPE::POLY) {
							mKernelFnc = std::bind(&SvmModel::polyKernel, this,
							                      std::placeholders::_1,
							                      std::placeholders::_2);
						}

						if (type->second == KERNEL_TYPE::SIGMOID) {
							mKernelFnc = std::bind(&SvmModel::sigmoidKernel, this,
							                       std::placeholders::_1,
							                       std::placeholders::_2);
						}

					}

					if (line == "degree") {
						if (!std::getline(modelFile, line))
							return false;

						mParams->degree = std::stoi(line);
					}

					if (line == "gamma") {
						if (!std::getline(modelFile, line))
							return false;

						std::stringstream ss(line);
						ss >> mParams->gamma;

						if (ss.fail())
							return false;
					}

					if (line == "coef0") {

						if (!std::getline(modelFile, line))
							return false;

						std::stringstream ss(line);
						ss >> mParams->coef0;

						if (ss.fail())
							return false;
					}

					if (line == "nr_class") {
						if (!std::getline(modelFile, line))
							return false;

						mNrClass = std::stoi(line);
						if (mNrClass <= 0)
							return false;

					}

					if (line == "total_sv") {
						if (!std::getline(modelFile, line))
							return false;

						mL = std::stoi(line);
					}

					if (line == "rho") {
						auto n = mNrClass * (mNrClass - 1) / 2;
						auto totals = 0;
						mRho.reserve(n);
						T val {0};

						if (!std::getline(modelFile, line))
							return false;

						std::stringstream ss(line);
						while (ss >> val) {
							mRho.emplace_back(val);
							totals++;
						}

						if (totals != n)
							return false;
					}

					if (line == "label") {

						mLabel.reserve(static_cast<size_t>(mNrClass));
						int val {0};
						auto totals = 0;

						if (!std::getline(modelFile, line))
							return false;

						std::stringstream ss(line);
						while (ss >> val) {
							mLabel.emplace_back(val);
							totals++;
						}

						if (totals != mNrClass)
							return false;

						for (const auto& l : mLabel) {
							if (l != Globals::SVM_POSITIVE_LABEL && l != Globals::SVM_NEGATIVE_LABEL) {
								std::cerr << "The labels does not correspond to expected labels" << std::endl;
								return false;
							}
						}
					}

					if (line == "probA") {

						auto n = mNrClass * (mNrClass - 1) / 2;
						mProbA.reserve(n);
						T val {0};
						auto totals = 0;

						if (!std::getline(modelFile, line))
							return false;

						std::stringstream ss(line);
						while (ss >> val) {
							mProbA.emplace_back(val);
							totals++;
						}

						if (totals != n)
							return false;
					}

					if (line == "probB") {

						auto n = mNrClass * (mNrClass - 1) / 2;
						mProbB.reserve(n);
						T val {0};
						auto totals = 0;

						if (!std::getline(modelFile, line))
							return false;

						std::stringstream ss(line);
						while (ss >> val) {
							mProbB.emplace_back(val);
							totals++;
						}

						if (totals != mNrClass)
							return false;
					}

					if (line == "nr_sv") {

						mNumSV.reserve(static_cast<unsigned>(mNrClass));
						int val {0};
						auto totals = 0;

						if (!std::getline(modelFile, line))
							return false;

						std::stringstream ss(line);
						while (ss >> val) {
							mNumSV.emplace_back(val);
							totals++;
						}

						if (totals != mNrClass)
							return false;

						break;
					}
				} // end while(...)

				if (!std::getline(modelFile, line))
					return false;

				if (line != "SV")
					return false;

			} catch (std::bad_alloc& ex) {
				return false;
			} catch (std::invalid_argument& iae) {
				return false;
			} catch (std::out_of_range& oore) {
				return false;
			}

			return true;
		}



		std::shared_ptr<SvmParameter<T>> mParams;
		std::function<T(const std::valarray<T>&, const std::valarray<T>&) noexcept> mKernelFnc;
		std::vector<T> mRho;
		std::vector<T> mProbA;
		std::vector<T> mProbB;
		std::vector<int> mLabel;
		std::vector<std::valarray<T>> mSV; // SVs
		std::vector<unsigned> mNumSV; // Number of SVs for each class
		std::vector<std::valarray<T>> mSVCoef; // Coefficients for SVs
		int mNrClass;
		int mL;
	};

}


#endif //INPROT_SVM_MODEL_H