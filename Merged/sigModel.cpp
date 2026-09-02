#include "sigModel.h"

#include <BAT/BCMath.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>

sigMOD::sigMOD(const std::string& name, const std::string& datafile,
               const std::string& bkgfile, const std::string& sigfile)
    : BCModel(name), bkgNorm(0.0), sigNorm(0.0), bkgPriorMean(0.0),
      signalMaximum(0.0), nbins(0) {
    data = readVector(datafile);
    bkg = readVector(bkgfile);
    sig = readVector(sigfile);

    if (data.size() < 3) {
        throw std::runtime_error("The toy data file must contain two metadata values and at least one bin.");
    }
    if (bkg.size() != data.size() - 2) {
        throw std::runtime_error("Background template length does not match the toy data bins.");
    }
    if (sig.size() != data.size() - 2) {
        throw std::runtime_error("Signal template length does not match the toy data bins.");
    }

    for (double value : bkg) {
        bkgNorm += value;
    }
    for (double value : sig) {
        sigNorm += value;
    }
    if (bkgNorm <= 0.0 || sigNorm <= 0.0) {
        throw std::runtime_error("Background and signal templates must have positive sums.");
    }

    bkgPriorMean = data[1];
    if (bkgPriorMean < 0.0 || data[0] < 0.0) {
        throw std::runtime_error("Toy metadata event counts cannot be negative.");
    }

    nbins = static_cast<int>(data.size() - 2);
    double observedTotal = 0.0;
    for (int bin = 0; bin < nbins; ++bin) {
        observedTotal += data[bin + 2];
    }
    signalMaximum = std::max(1.0, observedTotal);

    AddParameter("nb", 0.0, 3.0 * (data[1] + data[0]), "nb", "[counts]");
    AddParameter("ns", 0.0, signalMaximum, "ns", "[counts]");
}

double sigMOD::LogLikelihood(const std::vector<double>& pars) {
    double logLikelihood = 0.0;
    for (int bin = 0; bin < nbins; ++bin) {
        const double expected = pars[0] * bkg[bin] / bkgNorm
                              + pars[1] * sig[bin] / sigNorm;
        logLikelihood += BCMath::LogPoisson(
            static_cast<long long>(data[bin + 2]), expected);
    }
    return logLikelihood;
}

double sigMOD::LogAPrioriProbability(const std::vector<double>& pars) {
    const double backgroundPrior = BCMath::LogPoisson(
        static_cast<long long>(pars[0]), bkgPriorMean);

    // Uniform signal prior on [0, signalMaximum].
    const double signalPrior = -std::log(signalMaximum);
    return backgroundPrior + signalPrior;
}

std::vector<double> sigMOD::readVector(const std::string& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("Could not open input file: " + filename);
    }

    std::vector<double> values;
    double value = 0.0;
    while (input >> value) {
        values.push_back(value);
    }
    if (!input.eof()) {
        throw std::runtime_error("Input file contains a non-numeric value: " + filename);
    }
    return values;
}
