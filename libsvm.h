//
// Created by germelcar on 9/13/17.
//

#ifndef INPROT_LIBSVM_H
#define INPROT_LIBSVM_H

#include <memory>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <sstream>
#include <vector>



template <typename Condition>
using EnableIf = typename std::enable_if<Condition::value>::type;


namespace libsvm {

	//
	// SVM and KERNEL types
	//
	enum class SVM_TYPE { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };
	enum class KERNEL_TYPE { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED };

	void skipLine(std::ifstream& input) {
		input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}

	std::vector<std::string> split(const std::string& line, const char delim = ' ') {
		std::vector<std::string> tokens;

		std::stringstream ss(line);
		std::string tok;
		while (std::getline(ss, tok, delim))
			tokens.emplace_back(tok);

		return tokens;
	}


	const std::unordered_map<std::string, SVM_TYPE> svm_type_table {
			{"c_svc", SVM_TYPE::C_SVC},
			{"nu_svc", SVM_TYPE::NU_SVC},
			{"one_class", SVM_TYPE::ONE_CLASS},
			{"epsilon_svr", SVM_TYPE::EPSILON_SVR},
			{"nu_svr", SVM_TYPE::NU_SVR}
	};

	const std::unordered_map<std::string, KERNEL_TYPE> kernel_type_table {
			{"linear", KERNEL_TYPE::LINEAR},
			{"polynomial", KERNEL_TYPE::POLY},
			{"rbf", KERNEL_TYPE::RBF},
			{"sigmoid", KERNEL_TYPE::SIGMOID},
			{"precomputed", KERNEL_TYPE::PRECOMPUTED}
	};


	template<typename T, EnableIf<std::is_floating_point<T>>...>
	class SvmParameter {

	public:
		SvmParameter() = default;


		std::string check() const {
			if (svm_type != SVM_TYPE::C_SVC && svm_type != SVM_TYPE::NU_SVC &&
			    svm_type != SVM_TYPE::ONE_CLASS && svm_type != SVM_TYPE::EPSILON_SVR &&
			    svm_type != SVM_TYPE::NU_SVR)
				return std::string("Unknown svm type");

			if (kernel_type != KERNEL_TYPE::LINEAR && kernel_type != KERNEL_TYPE::POLY &&
			    kernel_type != KERNEL_TYPE::RBF && kernel_type != KERNEL_TYPE::SIGMOID &&
			    kernel_type != KERNEL_TYPE::PRECOMPUTED)
				return std::string("Unknown kernel type");

			if (gamma < 0)
				return std::string("gamma < 0");

			if (degree < 0)
				return std::string("Degree of polynomial kernel < 0");

			if (eps <= 0)
				return std::string("eps <= 0");

			if (svm_type == SVM_TYPE::C_SVC || svm_type == SVM_TYPE::EPSILON_SVR ||
			    svm_type == SVM_TYPE::NU_SVR)
				if (C <= 0)
					return std::string("C <= 0");

			if (svm_type == SVM_TYPE::NU_SVC || svm_type == SVM_TYPE::ONE_CLASS ||
			    svm_type == SVM_TYPE::NU_SVR)
				if (nu <= 0 || nu > 1)
					return std::string("nu <= 0 or nu > 1");

			return std::string("");
		}

		std::shared_ptr<size_t> weight_label; // for C_SVC
		std::shared_ptr<T> weight; // for C_SVC
		size_t nr_weight; // for C_SVC

		T gamma; // for poly/rbf/sigmoid kernels
		T coef0; // for poly/sigmoid kernels
		T eps; // stopping criteria
		T C; // for C_SVC, EPSILON_SVR and NU_SVR
		T nu; // for NU_SVC, ONE_CLASS and NU_SVR
		T p; // for EPSILON_SVR

		int degree; // for poly kernel
		SVM_TYPE svm_type;
		KERNEL_TYPE kernel_type;


		bool probability; // do probability estimates
	};

}


#endif //INPROT_LIBSVM_H