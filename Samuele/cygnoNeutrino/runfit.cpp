//
//  runfit.cpp
//  LIMEPMTfits
//
//  Created by Stefano Piacentini on 23/09/22.
//

#include <BAT/BCLog.h>
#include <BAT/BCModelManager.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <fstream> 
#include <cstring>
#include <ios>

#include "bkgModel.h"
#include "sigModel.h"


std::vector<std::string> split (const std::string &s, char delim);

int main(int argc, char *argv[]) {
    
    // Parsing arguments
    if(argc != 7) {
            throw std::invalid_argument("\nWrong number of arguments. Arguments:\n\t\
                                        1 - Name of data file.\n\t\
                                        2 - Name of bkg  file.\n\t\
                                        3 - Name of sig  file.\n\t\
                                        4 - Index of the run, for the output dir.\n\t\
                                        5 - Additional infos in the name of the output dir. \n\t\
                                        6 - Flag for output file \n\t\
                                        ");
    }

    const std::string datafile_name(argv[1]);
    const std::string  bkgfile_name(argv[2]);
    const std::string  sigfile_name(argv[3]);
    const int run_index  = atoi(argv[4]);
    const std::string addinfo(argv[5]);

    std::vector<std::string> v = split(addinfo, '_');

    int reduction = std::stoi(v[0]);
    int exposure = std::stoi(v[1]);
    
    int NIter = 100000;
    int Nch   = 8;
    
    std::string res_dir = "./results_V4_100toy_ModelComparison/run_" +
      addinfo +
      "_" +
      std::to_string(run_index) +
      "_"+
      std::to_string(Nch*NIter)+"/";
    
    std::system(("mkdir -p "+res_dir).c_str());
    /*int com = std::system(("mkdir "+res_dir).c_str());
    if(com != 0) {
        std::cout<<com<<std::endl;
	}*/
    
    //bkgMOD m("neutrinofit", datafile_name, bkgfile_name);
    sigMOD m("neutrinofit_S+B", datafile_name, bkgfile_name, sigfile_name);
    bkgMOD mbkg("neutrinofit_B",datafile_name,bkgfile_name);
    
    m.WriteMarkovChain(res_dir + m.GetSafeName() + "_mcmc.root", "RECREATE");
    mbkg.WriteMarkovChain(res_dir + mbkg.GetSafeName() + "_mcmc.root", "RECREATE");

    BCLog::OpenLog(res_dir+m.GetSafeName()+"_log.txt", BCLog::detail, BCLog::detail);

    BCLog::OutSummary("S+B model created");
    BCLog::OutSummary("B model created");
    
    m.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNIterationsPreRunMax(100000000);
    m.SetNIterationsRun(NIter);

    mbkg.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    mbkg.SetPrecision(BCEngineMCMC::kMedium);
    mbkg.SetNIterationsPreRunMax(100000000);
    mbkg.SetNIterationsRun(NIter);    

    double modem= m.LogEval(m.FindMode(m.GetBestFitParameters()));
    double modembgk= mbkg.LogEval(mbkg.FindMode(mbkg.GetBestFitParameters()));
    
    m.SetModeEval((modem+modembgk)/2);
    mbkg.SetModeEval((modem+modembgk)/2);
    
    std::cout << m.LogEval(m.FindMode(m.GetBestFitParameters())) << std::endl;
    std::cout << mbkg.LogEval(mbkg.FindMode(mbkg.GetBestFitParameters())) << std::endl;    
    
    BCModelManager mgr;
    mgr.AddModel(&m,0.5);
    mgr.AddModel(&mbkg,0.5);
    
    mgr.SetNChains(Nch);    
    
    // write Markov chains for both models to ROOT files
    mgr.WriteMarkovChain(res_dir + "mcmc_", "RECREATE", true, true);

    mgr.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    mgr.SetPrecision(BCEngineMCMC::kHigh);

    mgr.SetNIterationsMin(70000);
    mgr.SetNIterationsMax(100000);
    mgr.SetAbsolutePrecision(1e-8);
    mgr.SetRelativePrecision(1e-10);
    
    // run MCMC, marginalizing posterior and writing the Markov chains to ROOT files
    mgr.MarginalizeAll();
    
    std::cout << "Int S+B from marg " << mgr.GetModel(0)->GetIntegral() << std::endl;
    std::cout << "Int B from marg " << mgr.GetModel(1)->GetIntegral() << std::endl;
    
    /*
    m.FindMode(m.GetBestFitParameters());
    
    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized(res_dir+m.GetSafeName() + "_plots.pdf");

    // print summary plots
    m.PrintParameterPlot(res_dir+m.GetSafeName() + "_parameters.pdf");
    m.PrintCorrelationPlot(res_dir+m.GetSafeName() + "_correlation.pdf");
    m.PrintCorrelationMatrix(res_dir+m.GetSafeName() + "_correlationMatrix.pdf");
    m.PrintKnowledgeUpdatePlots(res_dir+m.GetSafeName() + "_update.pdf");
    */
    
    /*
    m.PrintAllMarginalized(res_dir+m.GetSafeName() + "_plots.pdf");
    m.PrintKnowledgeUpdatePlots(res_dir+m.GetSafeName() + "_update.pdf");  
    mbkg.PrintAllMarginalized(res_dir+mbkg.GetSafeName() + "_plots.pdf");  
    */
    
    //m.SetCubaIntegrationMethod(BCIntegrate::BCCubaMethod::kCubaDefault);
    //mbkg.SetCubaIntegrationMethod(BCIntegrate::BCCubaMethod::kCubaDefault);

    m.SetIntegrationMethod(BCIntegrate::BCIntegrationMethod::kIntGrid);
    mbkg.SetIntegrationMethod(BCIntegrate::BCIntegrationMethod::kIntGrid);
    
    mgr.Integrate();    

    long double B = mgr.BayesFactor(0, 1);
    mgr.PrintModelComparisonSummary();

    std::cout<< exposure <<"\t" << reduction << "\t" << B << " " << log10(B) << std::endl;
    
    std::ofstream myfile(Form("pippo_newAR_%i_%s.txt",exposure,argv[6]),std::ios_base::app);
    myfile<<reduction << "\t" << exposure <<"\t" << B << "\t" << log10(B) <<"\n";

    //sm.PrintResults(res_dir+m.GetSafeName() + "_results.txt");
    // print results of the analysis into a text file
    //mgr.PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();
    
    return 0;
}


//Borexino has a signal to noise ratio of 3.8e-5 3e-5





///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
////////////////////////// utility functions //////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

std::vector<std::string> split (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}






