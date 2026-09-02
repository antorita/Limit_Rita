#include <BAT/BCLog.h>
#include <BAT/BCModelManager.h>

#include "bkgModel.h"
#include "sigModel.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<std::string> split(const std::string& value, char delimiter) {
    std::vector<std::string> result;
    std::string item;
    std::stringstream stream(value);
    while (std::getline(stream, item, delimiter)) {
        result.push_back(item);
    }
    return result;
}

int main(int argc, char* argv[]) {
    try {
        if (argc != 7) {
            throw std::invalid_argument(
                "Usage: ./runfit DATA BKG SIG RUN_INDEX REDUCTION_EXPOSURE FLAG\n"
                "Example: ./runfit ../toyMC_txt/1_10/0.txt ../background.txt "
                "../signal.txt 0 1_10 NID");
        }

        const std::string dataFile = argv[1];
        const std::string bkgFile = argv[2];
        const std::string sigFile = argv[3];
        const int runIndex = std::stoi(argv[4]);
        const std::vector<std::string> info = split(argv[5], '_');
        if (info.size() != 2) {
            throw std::invalid_argument("The fifth argument must have the form reduction_exposure, e.g. 1_10.");
        }
        const int reduction = std::stoi(info[0]);
        const int exposure = std::stoi(info[1]);
        const std::string flag = argv[6];

        constexpr int iterations = 100000;
        constexpr int chains = 8;
        const std::string resultDirectory =
            "./results_Merged/run_" + std::string(argv[5]) + "_" +
            std::to_string(runIndex) + "_" +
            std::to_string(iterations * chains) + "/";
        const std::string mkdirCommand = "mkdir -p " + resultDirectory;
        if (std::system(mkdirCommand.c_str()) != 0) {
            throw std::runtime_error("Could not create result directory: " + resultDirectory);
        }

        sigMOD signalBackground("neutrinofit_SB", dataFile, bkgFile, sigFile);
        bkgMOD backgroundOnly("neutrinofit_B", dataFile, bkgFile);

        signalBackground.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        backgroundOnly.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        signalBackground.SetPrecision(BCEngineMCMC::kMedium);
        backgroundOnly.SetPrecision(BCEngineMCMC::kMedium);
        signalBackground.SetNIterationsPreRunMax(100000000);
        backgroundOnly.SetNIterationsPreRunMax(100000000);
        signalBackground.SetNIterationsRun(iterations);
        backgroundOnly.SetNIterationsRun(iterations);

        BCLog::OpenLog(resultDirectory + "neutrinofit_S+B_log.txt",
                       BCLog::detail, BCLog::detail);
        BCLog::OutSummary("Merged S+B and background-only models created");

        BCModelManager manager;
        manager.AddModel(&signalBackground, 0.5);
        manager.AddModel(&backgroundOnly, 0.5);
        manager.SetNChains(chains);
        manager.WriteMarkovChain(resultDirectory + "mcmc_", "RECREATE", true, true);
        manager.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
        manager.SetPrecision(BCEngineMCMC::kHigh);
        manager.SetNIterationsMin(70000);
        manager.SetNIterationsMax(100000);
        manager.SetAbsolutePrecision(1e-8);
        manager.SetRelativePrecision(1e-10);
        manager.MarginalizeAll();

        manager.SetIntegrationMethod(BCIntegrate::BCIntegrationMethod::kIntGrid);
        manager.Integrate();
        const long double bayesFactor = manager.BayesFactor(0, 1);
        manager.PrintModelComparisonSummary();

        std::cout << exposure << "\t" << reduction << "\t" << bayesFactor
                  << " " << std::log10(bayesFactor) << std::endl;
        std::ofstream summary("pippo_newAR_" + std::to_string(exposure) + "_" + flag + ".txt",
                              std::ios::app);
        summary << reduction << "\t" << exposure << "\t" << bayesFactor
                << "\t" << std::log10(bayesFactor) << "\n";

        BCLog::OutSummary("Exiting");
        BCLog::CloseLog();
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << std::endl;
        return 1;
    }
    return 0;
}
