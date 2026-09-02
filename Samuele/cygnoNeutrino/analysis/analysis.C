#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TMath.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <iterator>
#include <algorithm>

Double_t ConfLowInt(TH1D* h){
  Double_t integral=h->Integral();
  Double_t sum=0;
  Int_t counter = 1;
  while(sum < 0.16*integral){
    sum+= h->Integral(counter,counter);
    counter++;
  };
  return h->GetBinCenter(counter);
}


Double_t Median(TH1D* h){
  Double_t integral=h->Integral();
  Double_t sum=0;
  Int_t counter = 1;
  while(sum < 0.5*integral){
    sum+= h->Integral(counter,counter);
    counter++;
  };
  return h->GetBinCenter(counter);
}



void analysis(){
  
  //std::vector<Int_t> BkgLv={10, 100, 500, 1000,5000, 10000, 100000, 1000000};
  //std::vector<Int_t> BkgLv={10,100};
  
  std::vector<Int_t> SigLv={1, 5, 10, 20, 40, 60, 100, 200, 400, 600, 1000};
  //std::vector<Int_t> SigLv={10, 20, 100, 600, 1000};

  std::vector<Int_t> BkgLv={10};
  //std::vector<Int_t> SigLv={1};

  TFile* f;
  TTree* tree;
  TH1D* hSignal;
  //TGraphErrors* g[BkgLv.size()];
  std::vector<TGraphErrors*> g(BkgLv.size(), nullptr);

  std::vector<Double_t> nSigma;
  std::vector<Double_t> snSigma;  
  std::vector<Double_t> nSignal;

  std::vector<Double_t> Ratio;
  std::vector<Double_t> sRatio;

  Double_t sum, mean, smean, rms, srms;
    
  Double_t ns,nb;
  Int_t phase;
  
  for(int b=0;b<BkgLv.size();b++){

    for(int s=0;s<SigLv.size();s++){
      ///for(int file=0;file<100;file++){
      for(int file=0;file<4;file++){

	std::cout << BkgLv[b] << "\t" << SigLv[s] << "\t" << file << std::endl;
	
  // Updated path to match current results directory naming
  ///f = TFile::Open(Form("../results_V4_100toy_ModelComparison/run_%i_%i_%i_800000/mcmc_neutrinofit_SB.root",SigLv[s],BkgLv[b],file));
  f = TFile::Open(Form("results_V4_100toy_ModelComparison/run_%i_%i_%i_800000/mcmc_neutrinofit_SB.root",SigLv[s],BkgLv[b],file));
	if (!f || f->IsZombie()) {
	  std::cout << "Missing file for run_" << SigLv[s] << "_" << BkgLv[b] << "_" << file << "\n";
	  if (f) f->Close();
	  continue;
	}
	
	tree = (TTree*)f->Get("neutrinofit_SB_mcmc");
	if (!tree) {
	  std::cout << "Missing tree in file for run_" << SigLv[s] << "_" << BkgLv[b] << "_" << file << "\n";
	  f->Close();
	  continue;
	}
	
	tree->SetBranchAddress("nb",&nb);
	tree->SetBranchAddress("ns",&ns);
	tree->SetBranchAddress("Phase",&phase);
	
	///hSignal = new TH1D("hSignal","hSignal",500,0,1200);
  hSignal = new TH1D("hSignal","hSignal",500,0,1200);
  //hSignal->SetDirectory(0);
	
	for(int i=0;i<tree->GetEntries();i++){
	  tree->GetEntry(i);
	  if(phase==1){
	    hSignal->Fill(ns);
	  }//chiudo if fase
	}//chiudo for entries
	
	mean = hSignal->GetMean();
	smean = hSignal->GetMeanError();
	rms = hSignal->GetRMS();
	srms = hSignal->GetRMSError();
	
	if (rms != 0 && mean != 0) {
	  Ratio.push_back(mean/rms);
	  sRatio.push_back( mean/rms * sqrt( smean/mean * smean/mean + srms/rms * srms/rms ) );
	} else {
	  std::cout << "Skipped invalid statistics for run_" << SigLv[s] << "_" << BkgLv[b] << "_" << file << "\n";
	}
	f->Close();
      }//chiudo for on file index

      if (Ratio.empty() || sRatio.empty()) {
        std::cout << "No valid toys for combo " << SigLv[s] << "_" << BkgLv[b] << ", skipping" << std::endl;
      } else {
        nSignal.push_back(SigLv[s]);
        sum = std::accumulate(std::begin(Ratio), std::end(Ratio), 0.0);      
        nSigma.push_back(sum/Ratio.size());

        sum = std::accumulate(std::begin(sRatio), std::end(sRatio), 0.0);  
        snSigma.push_back(sum/sRatio.size());
      }

      Ratio.clear();
      sRatio.clear();
      
    }//chiudo for s

    if (!nSignal.empty()) {
      g[b]= new TGraphErrors((Int_t)nSignal.size(), &nSignal[0], &nSigma[0], 0, &snSigma[0]);
      g[b]->SetTitle(Form("BkgLv_%i",BkgLv[b]));
      g[b]->SetName(g[b]->GetTitle());
    } else {
      g[b] = nullptr;
    }

    nSignal.clear();
    nSigma.clear();
    snSigma.clear();

    
  }//chiudo for bkg

  
    
  gSystem->mkdir("Out", true);
  TFile* pippo = new TFile("Out/Plot_NSigmaVsSignal.root","recreate");
  pippo->cd();

  for(int b=0;b<BkgLv.size();b++){
    if (!g[b]) continue;
    g[b]->SetMarkerStyle(8);
    g[b]->SetMarkerColor(40+b*2);
    g[b]->Write();
  }

  pippo->Save();
  pippo->Close();
  
}//chiudo funzione
