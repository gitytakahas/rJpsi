/*
  functionmacro.C : provide QCD weights at drawing level
  modified from : https://github.com/CMS-HTT/QCDModelingEMu

  4 April 2017 Y.T
 */


#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
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


TFile *tauf;
TH1F *taubr_up;
TH1F *taubr_down;
TH2F *taubr_ref;
TH1F *Bcweight1D_up;
TH1F *Bcweight1D_down;
TFile *fmap;

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

float getCorrection(int year, float xgbs){
  float p0 = -1;
  float p1 = -1;
  float p2 = -1;

  if(year==2016){
    p0 = 0.922312;
    p1 = 0.0165410;
    p2 = 0.0166504;
  }else if (year==2017){
    p0 = 0.897591;
    p1 = 0.0247323;
    p2 = 0.0115164;
  }else if (year==2018){
    p0 = 0.986485;
    p1 = 0.0545809;
    p2 = 0.00711602;
  }


  float weight = p0 + p1*xgbs + p2*xgbs*xgbs;

  //  std::cout << year << " " << weight << std::endl;

  return weight;


}


void ReadFileTau(){

  std::cout << "Read tau BR file";

  tauf = new TFile("tauola/correction.root");

  //  mc = (TH1F*) f->Get("inclusive/bg_ul");
  //  mc->Scale(1./mc->GetSumOfWeights());

  taubr_up = (TH1F*) tauf->Get("envelope_up");  
  taubr_down = (TH1F*) tauf->Get("envelope_down");  
  taubr_ref = (TH2F*) tauf->Get("evt1_set1_2d");
  //  weight->Scale(1./weight->GetSumOfWeights());  
  
  //  weight->Divide(mc);

  //  f->Close();
 

  std::cout << ".... end" << std::endl;

}

void ReadFileBcWeight(){

  std::cout << "Read Bc weight file";

  fmap = new TFile("correction_complete_BcPt.root");
  Bcweight1D_up = (TH1F*) fmap->Get("mc_weight_up");   
  Bcweight1D_down = (TH1F*) fmap->Get("mc_weight_down");

  std::cout << ".... end" << std::endl;
}

Float_t getBcWeight(float pt, float direction){
  Float_t Bcw = 1; 
  Float_t binweight=1; 
  if (direction > 0.5 ) binweight = Bcweight1D_up->GetBinContent(Bcweight1D_up->FindBin(pt));  
  if (direction < -0.5 ) binweight = Bcweight1D_down->GetBinContent(Bcweight1D_down->FindBin(pt));  
  Bcw=binweight; 
  return Bcw;				       
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


Float_t getTauBrWeight_up(float val1, float val2){

  //  std::cout << "getWeight for " << val << std::endl;

  //  Int_t binid = taubr_up->GetXaxis()->FindBin(val);
  
  if(val1==-1 || val2==-1){
    //    std::cout << "negative value!" << std::endl;
    return 1;
  }

  Int_t binid_x = taubr_ref->GetXaxis()->FindBin(val1);
  Int_t binid_y = taubr_ref->GetYaxis()->FindBin(val2);
  Int_t binid = taubr_ref->GetXaxis()->GetNbins()*(binid_y - 1) + binid_x; 
  //  std::cout <<val1 << " " << val2 << ", binid_x =" << binid_x << ", binid_y= " << binid_y << ", binid=" << binid << std::endl;
  
  Float_t w = taubr_up->GetBinContent(binid);
  //  Float_t w = 1.;

  //  std::cout << "weight: b mass = "<<  val << ", w = "  << w << std::endl;
  return w;

}

Float_t getTauBrWeight_down(float val1, float val2){

  if(val1==-1 || val2==-1){
    //    std::cout << "negative value!" << std::endl;
    return 1;
  }

  //  std::cout << "getWeight for " << val << std::endl;

  //  Int_t binid = taubr_down->GetXaxis()->FindBin(val); 
  Int_t binid_x = taubr_ref->GetXaxis()->FindBin(val1);
  Int_t binid_y = taubr_ref->GetYaxis()->FindBin(val2);
  Int_t binid = taubr_ref->GetXaxis()->GetNbins()*(binid_y - 1) + binid_x; 
  //  Int_t binid = taubr_ref->FindBin(val1, val2);
 
  //  std::cout <<"binid=" << binid << std::endl;
  
  Float_t w = taubr_down->GetBinContent(binid);
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
  //  ReadFile();
  ReadFileTau();
  ReadFileBcWeight();
}
