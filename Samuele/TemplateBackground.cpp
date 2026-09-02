#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <THStack.h>
#include <Math/RootFinder.h>

using namespace ROOT::Math;

std::vector<Double_t> GenerateBinningE(){

  Double_t Thr=10;
  Double_t p0=10.249379;
  Double_t p1=21.230281;
  
  TF1* pippo = new TF1("pippo","x - ( [0] + [1] / sqrt(x) ) / 100 * x - [2] ");

  pippo->SetParameter(0,p0);
  pippo->SetParameter(1,p1);
		      
  std::vector<Double_t> BinE;

  BinE.push_back(Thr);

  int i=0;
  Double_t c=0;
    
  do{

    pippo->SetParameter(2,BinE[i]);

    RootFinder *k = new RootFinder();
    ROOT::Math::Functor1D f(*pippo);
    k->SetFunction(f, 1, 1600);
    k->Solve();
    c = k->Root();

    BinE.push_back( c + (c-BinE[i]));

    std::cout << BinE[i] << std::endl;
    i++;
  }while(c<1300);


  return BinE;
  
}

Double_t SmearE(Double_t Eele, TRandom* r){
  
  Double_t p0=10.249379;
  Double_t p1=21.230281;

  TF1* pippo = new TF1("pippo","[0] + [1] / sqrt(x)");

  pippo->SetParameter(0,p0);
  pippo->SetParameter(1,p1);

  Double_t ESmeared = r->Gaus(Eele,pippo->Eval(Eele)/100*Eele);
  
  return ESmeared;
}

double GetAngReso(TF1*AngReso, double Eele){

  if(Eele<215){ //215 NID   290 ED
    return AngReso->Eval(Eele);
  } else {
    return 3; //3 degrees NID     5 ED
  }
  
}


Double_t SmearTheta(Double_t Theta, Double_t E, TRandom* r){

  Double_t p0 = -4.37; //NID
  Double_t p1 = 108.5;

  //Double_t p0 = -2.97;
  //Double_t p1 = 135.9;

  TF1* AngReso = new TF1("AngReso","[0]+[1]/sqrt(x)");
  AngReso->SetParameter(0,p0);
  AngReso->SetParameter(1,p1);

  //std::cout << "energy: "<< E <<" theta "<< Theta  <<" reso " << GetAngReso(AngReso,E) <<" smeared " << r->Gaus(Theta,GetAngReso(AngReso,E)) << std::endl;    
  
  Double_t ThetaSmeared = r->Gaus(Theta,GetAngReso(AngReso,E));

  return ThetaSmeared;
}


std::vector<std::string> ReadData(std::string FileName){
  
  std::vector<std::string> FileList;
  std::ifstream file(FileName);
  std::string str; 
    
  if(!file.is_open()){std::cout << "file non trovato" << std::endl;}
  
  while (std::getline(file, str))
    {
      FileList.push_back(str);
    }
  
  return FileList;
  
}



double GetElectronEnergy(double Eneut, double cosTheta){

  double Ee=2*Eneut*Eneut*0.51099895*cosTheta*cosTheta/((Eneut+0.51099895)*(Eneut+0.51099895)-Eneut*Eneut*cosTheta*cosTheta);

  return Ee;
}


























  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////// Code for signal sim
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






void GenSignal(){


  std::vector<Double_t> BinHistoEReso;

  BinHistoEReso = GenerateBinningE();

  TH1D *HistoE = new TH1D("HistoE","HistoE",BinHistoEReso.size()-2,&BinHistoEReso[0]);

  TH1D* HistoCosThetaTot = new TH1D("HistoCosThetaTot","HistoCosThetaTot",150,-1,1);
  TH1D* HistoERecoTot = new TH1D("HistoERecoTot","HistoERecoTot",100,0,500);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////// Building angular resolution histogram. The number of histogram will be BinHistoEReso.size()-2
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  TF1* AngReso = new TF1("AngReso","[0]+[1]/sqrt(x)");
  AngReso->SetParameter(0,-10.11);
  AngReso->SetParameter(1,184.3);
  
  
  TH2D* AngHistogram[BinHistoEReso.size()-2];
  TH2D* AngHistogram_test[BinHistoEReso.size()-2]; 
  TH2D* AngHistogramCos[BinHistoEReso.size()-2];
  TH2D* AngHistogramCos_test[BinHistoEReso.size()-2];
  std::vector<Double_t> CosBin;
  
  Double_t ECenter;
  Int_t NbinsX;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////// Building Phi Vs Theta Histo  ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int i=0;i<BinHistoEReso.size()-2;i++){
    
    ECenter = (BinHistoEReso[i]+BinHistoEReso[i+1])/2;
    NbinsX=(int)( 180/(2*GetAngReso(AngReso,ECenter)) );
    if(NbinsX<4) NbinsX=4;
    AngHistogram[i] = new TH2D(Form("Histo%.0f",ECenter),Form("Histo%.0f",ECenter),NbinsX,0,180,(int)( 360/(2*GetAngReso(AngReso,ECenter)) ),-180,180);
    AngHistogram[i]->GetXaxis()->SetTitle("Theta [deg]");
    AngHistogram[i]->GetYaxis()->SetTitle("Phi [deg]");

    AngHistogram_test[i] = new TH2D(Form("Histo%.0f_test",ECenter),Form("Histo%.0f_test",ECenter),100,0,180,100,-180,180);
    AngHistogram_test[i]->GetXaxis()->SetTitle("Theta [deg]");
    AngHistogram_test[i]->GetYaxis()->SetTitle("Phi [deg]");
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////// Building Phi Vs CosTheta Histo  ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i=0;i<BinHistoEReso.size()-2;i++){
    
    ECenter = (BinHistoEReso[i]+BinHistoEReso[i+1])/2;
    NbinsX=(int)( 360/(2*GetAngReso(AngReso,ECenter)) );
    CosBin.clear();
    
    for(int j=AngHistogram[i]->GetXaxis()->GetNbins();j>=0;j--){
      CosBin.push_back(TMath::Cos(  (AngHistogram[i]->GetXaxis()->GetBinLowEdge(j+1))/180*TMath::Pi()  ) );
      std::cout << i <<"  "<<AngHistogram[i]->GetXaxis()->GetBinLowEdge(j+1)<< "  "<<(AngHistogram[i]->GetXaxis()->GetBinLowEdge(j+1))/180*TMath::Pi() << "  " << TMath::Cos((AngHistogram[i]->GetXaxis()->GetBinLowEdge(j+1))/180*TMath::Pi()) << std::endl;
    }

    AngHistogramCos[i]= new TH2D(Form("HistoCos%.0f",ECenter),Form("HistoCos%.0f",ECenter),CosBin.size()-1,&CosBin[0],(int)( 360/(2*GetAngReso(AngReso,ECenter)) ),-180,180);
    AngHistogramCos[i]->GetXaxis()->SetTitle("Theta [deg]");
    AngHistogramCos[i]->GetYaxis()->SetTitle("Phi [deg]");

    AngHistogramCos_test[i]= new TH2D(Form("HistoCos%.0f_test",ECenter),Form("HistoCos%.0f_test",ECenter),100,-1,1,100,-180,180);
    AngHistogramCos_test[i]->GetXaxis()->SetTitle("CosTheta [deg]");
    AngHistogramCos_test[i]->GetYaxis()->SetTitle("Phi [deg]");
  }

  std::cout << "here" <<std::endl;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////// Reading neutrino Spectrum And building histogram ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<std::string> FileValue;
  std::vector<double> ENu;
  std::vector<double> PE;

  int size= 20; 
  
  FileValue=ReadData("spectrum2.txt");

  ///////////////////////////////// Get Spectrum /////////////////////////////////
  
  for(int i=0; i<(int)FileValue.size();i++){
    std::string arr[size];
    int k = 0;
    std::stringstream ssin(FileValue.at(i));

    while (ssin.good() && k < size){
      ssin >> arr[k];
      //std::cout << arr[k]; 
      ++k;
    }//while
    ENu.push_back(std::stod(arr[0]));
    PE.push_back(std::stod(arr[2]));
    
  }//for
  
  /*
  for(int i=0;i<(int)FileValue.size();i++){
    std::cout << ENu[i] << "  " << PE[i] << std::endl;
  }
  */

  Double_t step = ENu[0];
  TH1D* HistoSolarSpec = new TH1D("HistoSolarSpec","HistoSolarSpec",ENu.size()-1,0+step+step/2,ENu.back()+step/2);
  HistoSolarSpec->GetXaxis()->SetTitle("EnuSun");

  for(int i=0;i<ENu.size();i++){
    HistoSolarSpec->SetBinContent(i,PE[i]);
    //  std::cout<< ENu[i] << "\t"<<PE[i] <<"\t"<< HistoSolarSpec->GetBinCenter(i)<< "\t" << HistoSolarSpec->GetBinContent(i) <<std::endl;
  }
  
  HistoSolarSpec->Draw("");  

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////// Interaction Simulation ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  double gL=-0.5+0.23142;
  double gR=0.23142;
  
  
  TF1* DiffCS = new TF1("pippo","2*(1.1663787e-5)^2*[0]*[1]*[1]*x/[2]* ([1]*[1]*[0]+[0]*[0]*[0]+2*[0]*[0]*[1]) / (([1]+[0])^2 -[1]*[1]*x^2)^2 *( (2+[3]+[4])^2+([3]+[4])^2*(1-(2*[1]*[0]*x^2)/(([1]+[0])^2-[1]*[1]*x^2) )^2  - ([3]-[4])*([3]+[4]+2)*(2*[0]*[0]*x^2/(([1]+[0])^2 -[1]*[1]*x^2)))*0.2*0.2*10^(-26)",0,1);
  
  
  DiffCS->SetParameter(0,0.51099895e-3);
  DiffCS->SetParameter(2,3.14159265358);
  DiffCS->SetParameter(3,gL+gR);
  DiffCS->SetParameter(4,gL-gR);

  TRandom* r = new TRandom();
  
  int NTries =  50000;

  Double_t Eneut;
  Double_t cosTheta;

  Double_t Eele;
  Double_t Theta;
  Double_t Phi;
  
  for(int i=0;i<NTries;i++){
    if(i%1000==0) std::cout << i << " / " << NTries<< std::endl;
    
    //HistoSolarSpec->GetXaxis()->SetRangeUser(400,600);
    Eneut=HistoSolarSpec->GetRandom();
    DiffCS->SetParameter(1,Eneut*0.001); //Energia in GeV

    cosTheta=DiffCS->GetRandom(0,1);

    Eele=GetElectronEnergy(Eneut,cosTheta)*1000; 

    //std::cout << "Eele= " << Eele << std::endl;       
    
    Eele=SmearE(Eele,r);
    
    //std::cout << "EeleSm= " << Eele << std::endl;     
    
    if(Eele>10){
    
      //if(i%10000==0)std::cout << Eele << std::endl;
      HistoE->Fill(Eele);
      HistoERecoTot->Fill(Eele);
      
      Theta = TMath::ACos(cosTheta)/TMath::Pi()*180;//* (arc4random() % 2 ? 1 : -1);
      //Qui theta è in gradi
      //std::cout << "Theta= " << Theta << std::endl;
      //std::cout << "Thetarad= " << Theta/180*TMath::Pi() << std::endl;     
      //std::cout << "CosTheta= " << TMath::Cos(Theta/180*TMath::Pi()) << std::endl;
      Theta= SmearTheta(Theta,Eele,r);
      //std::cout << "ThetaSm= " << Theta << std::endl;     
      //std::cout << "ThetaSmRad= " << Theta/180*TMath::Pi() << std::endl;
      //std::cout << "CosThetaSmRad= " << TMath::Cos(Theta/180*TMath::Pi()) << std::endl;     
      //std::cout << "cosTh=" <<cosTheta << "  Th" << Theta << std::endl; 
      Phi = (float)rand()/RAND_MAX*360-180;
      
      for(int j=0;j<BinHistoEReso.size()-2;j++){
	if(Eele>BinHistoEReso[j] && Eele<BinHistoEReso[j+1]){
	  
	  AngHistogram[j]->Fill(Theta+90,Phi);
	  AngHistogram_test[j]->Fill(Theta+90,Phi);
	  AngHistogramCos[j]->Fill(TMath::Cos(Theta/180.*TMath::Pi()),Phi);
	  AngHistogramCos_test[j]->Fill(TMath::Cos(Theta/180.*TMath::Pi()),Phi);
	  HistoCosThetaTot->Fill(TMath::Cos(Theta/180.*TMath::Pi()));
	}//chiudo if
      }//chiudo for
    } 
    //std::cout << std::endl;
  }//chiudo for
  
  HistoE->Draw();

  TCanvas* c = new TCanvas();
  c->Divide(4,3);

  TCanvas* cCos = new TCanvas();
  cCos->Divide(4,3);
  
  for(int i=0;i<BinHistoEReso.size()-2;i++){
    c->cd(i+1);
    AngHistogram[i]->Draw("COLZ");
    cCos->cd(i+1);
    AngHistogramCos[i]->Draw("COLZ");
  }
  
  TFile* fOutTemplate = new TFile("OutTemplate_Sign.root","recreate");
  fOutTemplate->cd();

  HistoE->Write();
   
  for(int i=0;i<BinHistoEReso.size()-2;i++){
    AngHistogram[i]->Write();
    AngHistogramCos[i]->Write();
  }
  
  fOutTemplate->Save();
  fOutTemplate->Close();


  TFile* fOutTemplate_test = new TFile("OutTemplate_Sign_test.root","recreate");
  fOutTemplate_test->cd();
  HistoERecoTot->Write();
  HistoCosThetaTot->Write();
  
  for(int i=0;i<BinHistoEReso.size()-2;i++){
    AngHistogram_test[i]->Write();
    AngHistogramCos_test[i]->Write();
  }
  
}




































  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////// Code for background sim
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






int main(int argc, char** argv){


  std::vector<Double_t> BinHistoEReso;

  BinHistoEReso = GenerateBinningE();

  TH1D *HistoE = new TH1D("HistoE","HistoE",BinHistoEReso.size()-2,&BinHistoEReso[0]);

  TH1D* HistoCosThetaTot = new TH1D("HistoCosThetaTot","HistoCosThetaTot",150,-1,1);
  TH1D* HistoERecoTot = new TH1D("HistoERecoTot","HistoERecoTot",100,0,350);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////// Building angular resolution histogram. The number of histogram will be BinHistoEReso.size()-2
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  TF1* AngReso = new TF1("AngReso","[0]+[1]/sqrt(x)");
  AngReso->SetParameter(0,-4.37); //NID
  AngReso->SetParameter(1,108.5);

  //AngReso->SetParameter(0,-2.97); //ED
  //AngReso->SetParameter(1,135.9);
  
  
  TH2D* AngHistogram[BinHistoEReso.size()-2];
  TH2D* AngHistogram_test[BinHistoEReso.size()-2]; 
  TH2D* AngHistogramCos[BinHistoEReso.size()-2];
  TH2D* AngHistogramCos_test[BinHistoEReso.size()-2];
  std::vector<Double_t> CosBin;
  
  Double_t ECenter;
  Int_t NbinsX;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////// Building Phi Vs Theta Histo  ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  for(int i=0;i<BinHistoEReso.size()-2;i++){
    
    ECenter = (BinHistoEReso[i]+BinHistoEReso[i+1])/2;
    NbinsX=(int)( 180/(2*GetAngReso(AngReso,ECenter)) );
    if(NbinsX<4) NbinsX=4;
    AngHistogram[i] = new TH2D(Form("Histo%.0f",ECenter),Form("Histo%.0f",ECenter),NbinsX,0,180,(int)( 360/(2*GetAngReso(AngReso,ECenter)) ),-180,180);
    AngHistogram[i]->GetXaxis()->SetTitle("Theta [deg]");
    AngHistogram[i]->GetYaxis()->SetTitle("Phi [deg]");

    AngHistogram_test[i] = new TH2D(Form("Histo%.0f_test",ECenter),Form("Histo%.0f_test",ECenter),100,0,180,100,-180,180);
    AngHistogram_test[i]->GetXaxis()->SetTitle("Theta [deg]");
    AngHistogram_test[i]->GetYaxis()->SetTitle("Phi [deg]");
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////// Building Phi Vs CosTheta Histo  ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i=0;i<BinHistoEReso.size()-2;i++){
    
    ECenter = (BinHistoEReso[i]+BinHistoEReso[i+1])/2;
    NbinsX=(int)( 360/(2*GetAngReso(AngReso,ECenter)) );
    CosBin.clear();
    
    for(int j=AngHistogram[i]->GetXaxis()->GetNbins();j>=0;j--){
      CosBin.push_back(TMath::Cos(  (AngHistogram[i]->GetXaxis()->GetBinLowEdge(j+1))/180*TMath::Pi()  ) );
      std::cout << i <<"  "<<AngHistogram[i]->GetXaxis()->GetBinLowEdge(j+1)<< "  "<<(AngHistogram[i]->GetXaxis()->GetBinLowEdge(j+1))/180*TMath::Pi() << "  " << TMath::Cos((AngHistogram[i]->GetXaxis()->GetBinLowEdge(j+1))/180*TMath::Pi()) << std::endl;
    }

    AngHistogramCos[i]= new TH2D(Form("HistoCos%.0f",ECenter),Form("HistoCos%.0f",ECenter),CosBin.size()-1,&CosBin[0],(int)( 360/(2*GetAngReso(AngReso,ECenter)) ),-180,180);
    AngHistogramCos[i]->GetXaxis()->SetTitle("Theta [deg]");
    AngHistogramCos[i]->GetYaxis()->SetTitle("Phi [deg]");

    AngHistogramCos_test[i]= new TH2D(Form("HistoCos%.0f_test",ECenter),Form("HistoCos%.0f_test",ECenter),100,-1,1,100,-180,180);
    AngHistogramCos_test[i]->GetXaxis()->SetTitle("CosTheta [deg]");
    AngHistogramCos_test[i]->GetYaxis()->SetTitle("Phi [deg]");
  }

  std::cout << "here" <<std::endl;

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////// Reading neutrino Spectrum And building histogram ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TRandom* r = new TRandom();
  r->SetSeed(0);
  
  int NTries =  1000000;

  Double_t Eele;
  Double_t cosTheta;

  Double_t Theta;
  Double_t CosTheta;
  Double_t Phi;

  TFile* fBkg= TFile::Open("NormalizedHistoEFCu.root");
  THStack* stack;
  stack = (THStack*)fBkg->Get("Hstack_cut");

  TH1D* SpectrumSum = (TH1D*)stack->GetStack()->Last();
    
  for(int i=0;i<NTries;i++){
    if(i%1000==0) std::cout << i << " / " << NTries<< std::endl;

    Eele = SpectrumSum->GetRandom();
    
    Eele=SmearE(Eele,r);

    
    if(Eele>10){
    
      //if(i%10000==0)std::cout << Eele << std::endl;
      HistoE->Fill(Eele);
      HistoERecoTot->Fill(Eele);

      CosTheta= 1-2*(float)rand()/RAND_MAX;
      
      if(rand()>0.5){
	Theta = TMath::ACos(CosTheta);
      } else {
	Theta = -TMath::ACos(CosTheta);
      }
      //Qui theta è in gradi
      //std::cout << "Theta= " << Theta << std::endl;
      //std::cout << "Thetarad= " << Theta/180*TMath::Pi() << std::endl;     
      //std::cout << "CosTheta= " << TMath::Cos(Theta/180*TMath::Pi()) << std::endl;
      //Theta= SmearTheta(Theta,r);

      //if(Theta<0) Theta = TMath::Abs(Theta);
      //if(Theta>180) Theta = 180-(Theta-180);
      
      //std::cout << "ThetaSm= " << Theta << std::endl;     
      //std::cout << "ThetaSmRad= " << Theta/180*TMath::Pi() << std::endl;
      //std::cout << "CosThetaSmRad= " << TMath::Cos(Theta/180*TMath::Pi()) << std::endl;     
      //std::cout << "cosTh=" <<cosTheta << "  Th" << Theta << std::endl; 
      Phi = (float)rand()/RAND_MAX*360-180;
      
      for(int j=0;j<BinHistoEReso.size()-2;j++){
	if(Eele>BinHistoEReso[j] && Eele<BinHistoEReso[j+1]){
	  
	  AngHistogram[j]->Fill(Theta,Phi);
	  AngHistogram_test[j]->Fill(Theta,Phi);
	  AngHistogramCos[j]->Fill(CosTheta,Phi);
	  AngHistogramCos_test[j]->Fill(CosTheta,Phi);
	  HistoCosThetaTot->Fill(CosTheta);
	}//chiudo if
      }//chiudo for
    } 
    //std::cout << std::endl;
  }//chiudo for
  
  HistoE->Draw();

  TCanvas* c = new TCanvas();
  c->Divide(4,3);

  TCanvas* cCos = new TCanvas();
  cCos->Divide(4,3);
  
  for(int i=0;i<BinHistoEReso.size()-2;i++){
    c->cd(i+1);
    AngHistogram[i]->Draw("COLZ");
    cCos->cd(i+1);
    AngHistogramCos[i]->Draw("COLZ");
  }
  
  TFile* fOutTemplate = new TFile(Form("Out_V8_CosThetaFlat_HSBkg_NID/OutTemplate_Bkg_%s.root",argv[1]),"recreate");
  fOutTemplate->cd();

  HistoE->Write();
   
  for(int i=0;i<BinHistoEReso.size()-2;i++){
    AngHistogram[i]->Write();
    AngHistogramCos[i]->Write();
  }

  fOutTemplate->Save();
  fOutTemplate->Close();


  TFile* fOutTemplate_test = new TFile(Form("Out_V8_CosThetaFlat_HSBkg_NID/OutTemplate_Bkg_%s_test.root",argv[1]),"recreate");
  fOutTemplate_test->cd();
  HistoERecoTot->Write();
  HistoCosThetaTot->Write();
  
  for(int i=0;i<BinHistoEReso.size()-2;i++){
    AngHistogram_test[i]->Write();
    AngHistogramCos_test[i]->Write();
  }
  
}


