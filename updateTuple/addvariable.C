#include <iostream>
#include <string>
#include <algorithm> 

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace TMVA;


//void addvariable(TString filename, TString wfile, TString tdir, TString idstr){
void addvariable(TString filename, TString wfile){

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "input filename = " <<  filename << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "weight file = " << wfile << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  //  std::cout << "transient dir = " << tdir << std::endl;
  //  std::cout << "------------------------------------------" << std::endl;
  //  std::cout << "id dir = " << idstr << std::endl;
  //  std::cout << "------------------------------------------" << std::endl;

  //  int id = std::stoi((string) idstr);

  bool isBG = false;
  bool isSignal = false;
  bool isData = false;

//  if(string(filename).find("BcJpsiX")!=std::string::npos){
//    std::cout << "This is BG!" << std::endl;
//    isBG = true;
//  }else 
  if(string(filename).find("BcJpsiTau")!=std::string::npos){
    std::cout << "This is signal !!" << std::endl;
    isSignal = true;
  }else if(string(filename).find("Data")!=std::string::npos){
    std::cout << "This is Data !!" << std::endl;
    isData = true;
  }else if(string(filename).find("BJpsiX")!=std::string::npos){
    std::cout << "This is BG !!" << std::endl;
    isBG = true;
  }else{
    std::cout << "This is invalid!" << std::endl;
  }

  std::cout << "reading weight file" << std::endl;

  TFile *wFile = new TFile(wfile);
  TH2F *hist = (TH2F*) gROOT->FindObject("ratio");

  std::cout << "input file: ";

  TFile *lFile = new TFile(filename);
  TTree *lTree = (TTree*) lFile->FindObjectAny("tree");

  std::cout << "... done" << std::endl;

  Double_t b_pt, b_eta;

  //  TString _output;
  string output = (string)filename;

  //  std::cout << "output file = " <<  output << std::endl;
  string tbr = ".root";
  auto pos = output.find(tbr);
  auto len = tbr.length();

  if(pos != std::string::npos){
    output.replace(pos, len, "_weightAdded.root");
  }

  std::cout << "output file = " <<  output << std::endl;

    
//  if(isSignal){
//    lTree->SetBranchAddress( "puweight", &b_puweight);
//    lTree->SetBranchAddress( "hammer_ebe", &b_hammer);
//    //    lTree->SetBranchAddress( "hammer_ebe_toy", &b_toy);
//    //    _output = tdir + "/transient_large_" + idstr + ".root";
//  }

  if(isBG){
    lTree->SetBranchAddress( "b_pt", &b_pt);
    lTree->SetBranchAddress( "b_eta", &b_eta);
  }

//    lTree->SetBranchAddress( "genWeightBkgB", &b_bkgBweight);
//    _output = (TString)output;
//  }

  TFile *lOFile = new TFile((TString)output,"RECREATE");
  TTree *lOTree = lTree->CloneTree(0);
   
  Float_t lweight;

  lOTree->Branch("weight", &lweight, "weight/F");
  

  for (Long64_t i0=0; i0<lTree->GetEntries(); i0++) {
    lTree->GetEntry(i0);

    if(i0%10000==0) std::cout << i0 << " " << lTree->GetEntries() << " (" << double(i0)/double(lTree->GetEntries()) << ") has been processed ..." << std::endl;
    //    if(i0==20000) break;

    
//    if(isBG){
//

//      lweight = histw*b_puweight*b_bkgBweight;
//
//    }else 
//    if(isSignal){
//
//      //if(id==-1)
//      lweight = b_hammer*b_puweight/0.55;
////      else{
////	lweight = b_toy[id]*b_puweight;
////      }
//
//    }else if(isData){

    if(isBG){
      Float_t histw = hist->GetBinContent(hist->GetXaxis()->FindBin(TMath::Abs(b_eta)), hist->GetYaxis()->FindBin(b_pt));      
      lweight = histw;
    }
      //      std::cout << "Invalid data type" << std::endl;
      //    }

    lOTree->Fill();    
  }

  lOTree->Write();
  lOFile->Close();

  delete lFile;

  lFile = 0; lTree = 0; lOTree = 0;


////  if(isSignal){
////
////   
////    std::cout << "cleaning up the branch" << std::endl;
////    TFile *input = TFile::Open(_output);
////    TTree *inputtree; input->GetObject("tree",inputtree);
////    TFile *foutput = TFile::Open(tdir + "/transient_" + idstr + ".root","RECREATE");
////    inputtree->SetBranchStatus("hammer_ebe_toy",0);
////    TTree *outputtree = inputtree->CloneTree(-1,"fast"); 
////    foutput->Write(); 
////    delete input; 
////    delete foutput;
////
////    //    TString passbuf_mv = "mv -f " + tdir + "/transient_" + idstr + ".root " + output;
////    TString passbuf_mv = "xrdcp -f " + tdir + "/transient_" + idstr + ".root root://t3dcachedb03.psi.ch/" + output;
////
//////    if(idstr!="-1"){
//////      string output_tbr = (string)output;
//////      output_tbr.replace(output_tbr.begin(), output_tbr.end(), ".root", "_" + idstr + ".root");
//////
////////      auto pos = output_tbr.find(".root");
////////      auto len = t.length();
////////      if (pos != std::string::npos) {
////////	.replace(pos, len, "|"); // s == "a|b"
////////      }
//////
//////
//////      //      TString output_wid = output;
//////      //      string output_wid = std::regex_replace(output, regex(".root"), "_" + idstr + ".root");
//////      std::cout << "check!!! " << output_tbr << std::endl;
//////      passbuf_mv = "mv -f " + tdir + "/transient_" + idstr + ".root " + output_tbr;
//////    }
////
////    system(passbuf_mv);
////
////    //char passbuf_rm[] = "rm " + output_;
////    TString passbuf_rm = "rm -f " + _output;
////
////    system(passbuf_rm);
////
////    passbuf_rm = "rm -f " + tdir + "/transient_" + idstr + ".root";
////
////    system(passbuf_rm);
////
////  }


  
}

