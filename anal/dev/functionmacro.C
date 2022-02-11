/*
  functionmacro.C : provide QCD weights at drawing level
  modified from : https://github.com/CMS-HTT/QCDModelingEMu

  4 April 2017 Y.T
 */


#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TROOT.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TFile.h"
#include "TLorentzVector.h"


//RooWorkspace *w;
//TLorentzVector vec;
//TLorentzVector vec_check;

TFile *f;
TH1F *mc;
TH1F *weight;
//TH1F *weight;

void ReadFile(){

  std::cout << "Read file";

  f = new TFile("datacard/inclusive/b_mass_sf.root");

  mc = (TH1F*) f->Get("inclusive/bg_ul");
  //  mc->Scale(1./mc->GetSumOfWeights());

  weight = (TH1F*) f->Get("inclusive/data_obs");  
  //  weight->Scale(1./weight->GetSumOfWeights());  
  
  weight->Divide(mc);

  //  f->Close();

  std::cout << "... end" << std::endl;

}


Float_t getWeight(float val){

  //  std::cout << "getWeight for " << val << std::endl;

  Int_t binid = weight->GetXaxis()->FindBin(val);
  
  //  std::cout <<"binid=" << binid << std::endl;
  
  Float_t w = weight->GetBinContent(binid);
  //  Float_t w = 1.;

  //  std::cout << "weight: b mass = "<<  val << ", w = "  << w << std::endl;
  return w;

}


Float_t deltaPhi(Float_t p1, Float_t p2){

  Float_t res = p1 - p2;
  while(res > TMath::Pi()){
    res -= 2*TMath::Pi();
  }
  while(res < -TMath::Pi()){
    res += 2*TMath::Pi();
  }
  
  return res;
}









void functionmacro(){
  std::cout << std::endl;
  std::cout << "Initialize functionmacro.C ..." << std::endl;
  std::cout << std::endl;
  ReadFile();
}
