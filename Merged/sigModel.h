#ifndef MERGED_SIGMODEL_H
#define MERGED_SIGMODEL_H

#include <BAT/BCModel.h>

#include <string>
#include <vector>

class sigMOD : public BCModel {
public:
    sigMOD(const std::string& name, const std::string& datafile,
           const std::string& bkgfile, const std::string& sigfile);

    double LogLikelihood(const std::vector<double>& pars) override;
    double LogAPrioriProbability(const std::vector<double>& pars) override;

private:
    static std::vector<double> readVector(const std::string& filename);

    std::vector<double> data;
    std::vector<double> bkg;
    std::vector<double> sig;
    double bkgNorm;
    double sigNorm;
    double bkgPriorMean;
    double signalMaximum;
    int nbins;
};

#endif
