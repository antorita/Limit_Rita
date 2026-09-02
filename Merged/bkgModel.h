#ifndef MERGED_BKGMODEL_H
#define MERGED_BKGMODEL_H

#include <BAT/BCModel.h>

#include <string>
#include <vector>

class bkgMOD : public BCModel {
public:
    bkgMOD(const std::string& name, const std::string& datafile,
           const std::string& bkgfile);

    double LogLikelihood(const std::vector<double>& pars) override;
    double LogAPrioriProbability(const std::vector<double>& pars) override;

private:
    static std::vector<double> readVector(const std::string& filename);

    std::vector<double> data;
    std::vector<double> bkg;
    double bkgNorm;
    double bkgPriorMean;
    int nbins;
};

#endif
