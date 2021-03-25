#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace TMVA;


void addvariable(TString filename, TString output, TString wfile, TString tdir, TString idstr){

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "input file = " <<  filename << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "output file = " << output << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "weight file = " << wfile << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "transient dir = " << tdir << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "id dir = " << idstr << std::endl;
  std::cout << "------------------------------------------" << std::endl;

  int id = std::stoi((string) idstr);
  //  std::cout << "ID = " << idstr << " " << id << std::endl;

  bool isData = false;
  bool isBG = false;
  bool isSignal_truth = false;
  bool isSignal_all = false;

  if(string(filename).find("Data")!=std::string::npos){
    std::cout << "This is data!" << std::endl;
    isData = true;
  }
  else if(string(filename).find("BcJpsiX")!=std::string::npos){
    std::cout << "This is BG!" << std::endl;
    isBG = true;
  }
  else if(string(filename).find("all")!=std::string::npos){
    std::cout << "This is signal all!" << std::endl;
    isSignal_all = true;
  }
  else if(string(filename).find("truth")!=std::string::npos){
    std::cout << "This is signal truth!" << std::endl;
    isSignal_truth = true;
  }else{
    std::cout << "This is invalid!" << std::endl;
  }

  TFile *wFile = new TFile(wfile);
  TH2F *hist = (TH2F*) gROOT->FindObject("ratio");

  TFile *lFile = new TFile(filename);
  TTree *lTree = (TTree*) lFile->FindObjectAny("tree");

  Double_t b_pt, b_eta, b_puweight, b_hammer;
  Double_t b_toy[3000];
  Bool_t b_isveto;
  //  Bool_t b_tau_isRight;
  Double_t b_rhomass1, b_rhomass2; 
  Int_t b_tau_index;

  TString output_;
  output_ = output;

  if(!isData){
    lTree->SetBranchAddress( "puweight", &b_puweight);
  }
  if(isSignal_truth || isSignal_all){
    lTree->SetBranchAddress( "hammer_ebe", &b_hammer);
  }

  if(isSignal_truth){
    lTree->SetBranchAddress( "tau_index", &b_tau_index);
  }
  
  if(isSignal_all){
    lTree->SetBranchAddress( "hammer_ebe_toy", &b_toy);
    output_ = tdir + "/transient_large_" + idstr + ".root";
  }

  if(isBG || isData){
    lTree->SetBranchAddress( "b_pt", &b_pt);
    lTree->SetBranchAddress( "b_eta", &b_eta);
    lTree->SetBranchAddress( "isveto", &b_isveto);
    lTree->SetBranchAddress( "tau_rhomass1", &b_rhomass1);
    lTree->SetBranchAddress( "tau_rhomass2", &b_rhomass2);
    //    lTree->SetBranchAddress( "tau_isRight", &b_tau_isRight);
  }

  TFile *lOFile = new TFile(output_,"RECREATE");
  TTree *lOTree = lTree->CloneTree(0);
   
  Float_t lweight;

  lOTree->Branch("weight", &lweight, "weight/F");
  

  for (Long64_t i0=0; i0<lTree->GetEntries(); i0++) {
    lTree->GetEntry(i0);

    if(i0%10000==0) std::cout << i0 << " " << lTree->GetEntries() << " (" << double(i0)/double(lTree->GetEntries()) << ") has been processed ..." << std::endl;
    //    if(i0==20000) break;

    
    if(isBG){
      if(b_isveto){
	std::cout << "This event has signal !!! so removed." <<std::endl;
	continue;
      }

      //      if(b_rhomass1 > 1.2 || b_rhomass2 > 1.2) continue;

      Float_t histw = hist->GetBinContent(hist->GetXaxis()->FindBin(TMath::Abs(b_eta)), hist->GetYaxis()->FindBin(b_pt));
      lweight = histw*b_puweight;

    }else if(isSignal_truth){
      //      std::cout << "hammer = " << b_hammer << std::endl;
//      if(b_hammer==-1){
//	std::cout <<"This cannot be possible!" << std::endl;
//	continue;
//      }

      //      if(b_tau_index!=0) continue;

      lweight =b_puweight;
      //      std::cout << b_toy[3000] << std::endl;
    }else if(isSignal_all){

      //      std::cout << b_toy[0] << " " << b_toy[1] << " " <<b_toy[2999] << std::endl;
      
//      if(b_hammer==-1){
//	std::cout <<"This cannot be possible!" << std::endl;
//	continue;
//      }

      if(id==-1) lweight = b_puweight;
      else{
	lweight = b_toy[id]*b_puweight;
      }

    }else if(isData){

      Float_t histw = hist->GetBinContent(hist->GetXaxis()->FindBin(TMath::Abs(b_eta)), hist->GetYaxis()->FindBin(b_pt));
      lweight = histw;

      //      lweight = 1;
    }else{
      std::cout << "Invalid data type" << std::endl;
    }

    lOTree->Fill();    
  }

  lOTree->Write();
  lOFile->Close();

  delete lFile;

  lFile = 0; lTree = 0; lOTree = 0;


  if(isSignal_all){

    

    std::cout << "cleaning up the branch" << std::endl;
    TFile *input = TFile::Open(output_);
    TTree *inputtree; input->GetObject("tree",inputtree);
    TFile *foutput = TFile::Open(tdir + "/transient_" + idstr + ".root","RECREATE");
    //    //    TFile *output = TFile::Open(tdir + "/transient_" + idstr + ".root","RECREATE");
    //    TFile *output = TFile::Open("transient_" + idstr + ".root","RECREATE");
    inputtree->SetBranchStatus("hammer_ebe_toy",0);
    TTree *outputtree = inputtree->CloneTree(-1,"fast"); 
    foutput->Write(); 
    delete input; 
    delete foutput;

    //    char passbuf_mv[] = "mv " + tdir + "/transient_" + idstr + ".root " + output;
    TString passbuf_mv = "mv -f " + tdir + "/transient_" + idstr + ".root " + output;

    system(passbuf_mv);

    //char passbuf_rm[] = "rm " + output_;
    TString passbuf_rm = "rm -f " + output_;

    system(passbuf_rm);

//    std::cout << "cleaning up the branch" << std::endl;
//    TFile *input = TFile::Open(output, "update");
//    TTree *inputtree; input->GetObject("tree",inputtree);
//    //    TFile *output = TFile::Open("test.root","RECREATE");
//    inputtree->SetBranchStatus("hammer_ebe_toy",0);
//    //    TTree *outputtree = inputtree->CloneTree(-1,"fast"); 
//    //    output->Write(); 
//    inputtree->Write();
//    delete input; 
//    delete output;

  }


  
//  TFile f(output, "update");
//  TTree *T = (TTree*)f.Get("tree");
//  TBranch *b = T->GetBranch("hammer_ebe_toy");
//  T->GetListOfBranches()->Remove(b);
//  T->Write();
//  
//  f.Close();
}

