#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <stdio.h>
#include <time.h>


using namespace TMVA;

inline void InitRand()
{
  //  srand((unsigned int)time(NULL));
  // seeding set on 24 Jan. 
  srand(20220124);
}




void pick(TString filename, TString output, TString output_training, TString frac){

  InitRand();
  
  //  std::cout << "input file = " <<  filename << std::endl;

  float fr = std::stof((string) frac);

  cout << "Processing " << filename << " -> output = " << output << " with fraction = " << fr << endl;

  //  TFile *wFile = new TFile("weights/weight.root");
  //  TH2F *hist = (TH2F*) gROOT->FindObject("ratio");

  TFile *lFile = new TFile();
  lFile = TFile::Open(filename);
  TTree *lTree = (TTree*) lFile->FindObjectAny("tree");

  Int_t nentries = lTree->GetEntriesFast();
  std::cout << "number of events = " << nentries << std::endl;


  //  Double_t tau_mass;

  //  lTree->SetBranchAddress( "tau_mass", &tau_mass);
  //  lTree->SetBranchAddress( "b_eta", &b_eta);


  TFile *lOFile = new TFile(output,"RECREATE");
  TTree *lOTree = lTree->CloneTree(0);


  //  delete lFile_analysis; 
  //  delete lTree_analysis;


  //  TFile *lFile_training = new TFile(filename);
  //  TTree *lTree_training = (TTree*) lFile_training->FindObjectAny("tree");

  TFile *lOFile_training = new TFile(output_training,"RECREATE");
  TTree *lOTree_training = lTree->CloneTree(0);
   

  //  Float_t lweight;

  //  lOTree->Branch("weight", &lweight, "weight/F");
  
  Int_t counter = 0;
  Int_t counter_training = 0;

  for (Long64_t i0=0; i0<lTree->GetEntries(); i0++) {
    lTree->GetEntry(i0);

    
    //    std::cout << num << std::endl;


    if(i0%10000==0) std::cout << i0 << " " << lTree->GetEntries() << " (" << double(i0)/double(lTree->GetEntries()) << ") has been processed ..." << std::endl;

    //    Float_t histw = hist->GetBinContent(hist->GetXaxis()->FindBin(TMath::Abs(b_eta)), hist->GetYaxis()->FindBin(b_pt));
    //    std::cout << b_pt << " " << b_eta << " "  << histw  << std::endl;
    //    lweight = histw;
    Float_t num = (double)rand()/RAND_MAX;

    if(num > fr){
    //    if(0.75 < tau_mass && tau_mass < 1.5 && !(0.8 < tau_mass && tau_mass < 1.45)){
      lOTree->Fill();    
      counter+=1;
    }else{
      lOTree_training->Fill();    
      counter_training+=1;
    }



//    }else{
//      lOTree_analysis->Fill();
//      counter_analysis+=1;
//    }
  }

  std::cout << counter << ", " << lOTree->GetEntriesFast() << "/" << nentries << "  (" << (Float_t) counter/nentries << ") is used for the analysis" << std::endl;
  std::cout << counter_training << ", " << lOTree_training->GetEntriesFast() << "/" << nentries << "  (" << (Float_t) counter_training/nentries << ") is used for the training" << std::endl;
  //  std::cout << counter_analysis << ", " << lOTree_analysis->GetEntries()  << "/" << nentries << "  (" << (Float_t) counter_analysis/nentries << ") is used for the analysis" << std::endl;

  lOFile->cd();
  lOTree->Write();
  lOFile->Close();

  lOFile_training->cd();
  lOTree_training->Write();
  lOFile_training->Close();


  delete lFile;

  lFile = 0; lTree = 0; lOTree = 0;

}

