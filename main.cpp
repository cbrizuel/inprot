#include <iostream>
#include <locale>
#include "cli.h"
#include "kmers_manager.h"

using namespace fasta;
using namespace std;
using namespace rang;

// Change the data type of the variable for specify the data type used
// during the molecular descriptors calculation
// for example:
//      using MD_T = float;
// if you want to use "float" instead of double
using MD_T = double;

int main(int argc, char *argv[]) {


	try {

        // Start counting elapsed time
        auto start = chrono::steady_clock::now();

		// Create thousand comma separator
		locale comma_locale(locale(), new Globals::comma_numpunct());

		// tell cout to use our new locale.
		cout.imbue(comma_locale);

		Cli cli;
		auto ec = cli.parseArgs(argc, argv);

		if (ec != 1)
			exit(ec);

		SvmScaling<MD_T> scaling;
		auto scalingLoaded = scaling.restore(cli.getScalingFile());

		if (!scalingLoaded) {
			cerr << style::bold << fg::red<< "[ERROR] " << style::reset << fg::red
			     << "Error while loading SVM's scaling file" << endl;
			return 0;
		}

#ifdef USE_LIBSVM

		shared_ptr<svm_model> model (svm_load_model(cli.getModelFile().c_str()));

		if (!model) {
			cerr << style::bold << fg::red << "[ERROR] " << style::reset << fg::red
			<< "Error while loading SVM's model file" << endl;
			return 0;
		}
#else
		SvmModel<MD_T> model;
		auto modelLoaded = model.load(cli.getModelFile());

		if (!modelLoaded) {
			cerr << style::bold << fg::red << "[ERROR] " << style::reset << fg::red
			     << "Error while loading SVM's model file" << endl;
			return 0;
		}
#endif
        if (cli.getNumThreads() == -1) {
            std::cout << style::bold << fg::yellow << "[WARNING] " << style::reset << fg::yellow
                      << "Number of threads set to: automatic" << std::endl;
        }

		// Initialize the task scheduler
		tbb::task_scheduler_init init (cli.getNumThreads());

		KmersManager km(cli.getInputFile(),     // Input file - proteome file
		                cli.getOutputFile(),    // Output file - shrinked proteome
		                cli.getLowerKmer(),     // Lower k-mer size
		                cli.getUpperKmer(),     // Upper k-mer size
		                cli.getWritePreds(),    // Write predicteds k-mers (AMPs, NAMPs, boths, none)
		                cli.hasAwareMode(),     // Aware mode ==> low-memory consumption
						cli.hasVerboseMode());  // Has verbose mode enabled? ==> show extra information

		//
		// Extracting k-mers
		//
		cout << endl << style::bold << fg::green << "[STATUS] " << style::reset << fg::green
		     <<"Extracting k-mers..." << endl;

		auto okErr = km.extract(scaling, model);

		if (!okErr.first) {
            cerr << style::bold << fg::red << "[ERROR] " << style::reset << fg::red << okErr.second << endl;
			return 1;
		}

		//
		// Shrinking proteome
		//
		cout << endl << style::bold << fg::green << "[STATUS] " << style::reset << fg::green
		     <<"Shrinking proteome..." << endl;

		auto shrinked = km.shrinkProteome();

		if (!shrinked) {
			cerr << style::bold << fg::red << "[ERROR] " << style::reset << fg::red
			     << "Error while processing and shrinking proteome" << endl;
			return 0;
		}

        // Calculate elapsed time in milliseconds
        auto stop = chrono::steady_clock::now();
        auto diff = stop - start;
        cout << style::bold << fg::blue << "[INFO] " << style::reset << fg::blue
             << "Elapsed: " << chrono::duration_cast<chrono::milliseconds>(diff).count()
             << " ms" << style::reset << std::endl;

        cout << style::bold << fg::green << "[DONE]" << style::reset << endl;


    } catch (CLI::ValidationError& vee) {
		cerr << vee.what() << endl;
	} catch (std::exception& exc) {
		cerr << exc.what() << endl;
	}

	std::cout << style::reset;

    return 0;
}
