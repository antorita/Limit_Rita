#include "bkgModel.h"

#include <BAT/BCMath.h>

#include <fstream>
#include <stdexcept>

bkgMOD::bkgMOD(const std::string& name, const std::string& datafile,
               const std::string& bkgfile)
    : BCModel(name), bkgNorm(0.0), bkgPriorMean(0.0), nbins(0) {
    data = readVector(datafile);
    bkg = readVector(bkgfile);

    if (data.size() < 3) {
        throw std::runtime_error("The toy data file must contain two metadata values and at least one bin.");
    }
    if (bkg.size() != data.size() - 2) {
        throw std::runtime_error("Background template length does not match the toy data bins.");
    }

    for (double value : bkg) {
        bkgNorm += value;
    }
    if (bkgNorm <= 0.0) {
        throw std::runtime_error("Background template must have a positive sum.");
    }

    bkgPriorMean = data[1];
    if (bkgPriorMean < 0.0) {
        throw std::runtime_error("The background metadata count cannot be negative.");
    }

    nbins = static_cast<int>(data.size() - 2);
    AddParameter("nb", 0.0, 3.0 * (data[1] + data[0]), "nb", "[counts]");
}

double bkgMOD::LogLikelihood(const std::vector<double>& pars) {
    double logLikelihood = 0.0;
    const double expectedBackground = pars[0];

    for (int bin = 0; bin < nbins; ++bin) {
        const double expected = expectedBackground * bkg[bin] / bkgNorm;
        logLikelihood += BCMath::LogPoisson(
            static_cast<long long>(data[bin + 2]), expected);
    }
    return logLikelihood;
}

double bkgMOD::LogAPrioriProbability(const std::vector<double>& pars) {
    // The background prior is Poisson with mean equal to the toy metadata count.
    return BCMath::LogPoisson(static_cast<long long>(pars[0]), bkgPriorMean);
}

std::vector<double> bkgMOD::readVector(const std::string& filename) {
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
