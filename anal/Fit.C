#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooAddPdf.h"
#include "RooAbsReal.h"
#include "RooBinning.h"
#include "RooCBShape.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooConstVar.h"
#include "RooDLLSignificanceMCSModule.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooDecay.h"
#include "RooErrorVar.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "RooNLLVar.h"
#include "RooNumConvPdf.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooProfileLL.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "TColor.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/NumberCountingPdfFactory.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/SPlot.h"
#include "RooVoigtian.h"
#include "RooWorkspace.h"
#include "RooArgusBG.h"
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH3.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPaveLabel.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "TSystem.h"

using namespace RooFit ;
using namespace RooStats ;

void fit();
void Fit() { fit(); }
void fit() {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  //gSystem->Load("common/functionmacro_C.so");
  //  gROOT->Macro("common/functionmacro.C+");

  double lower_limit = 2.0, upper_limit =7.0;
  RooRealVar xgbs ("xgbs", "xgbs", lower_limit, upper_limit);

  RooRealVar puweight("puweight", "puweight", 0., 100.); 
  RooRealVar mu1_SFReco ("mu1_SFReco", "mu1_SFReco", 0., 100);
  RooRealVar mu2_SFReco ("mu2_SFReco", "mu2_SFReco", 0, 100);
  RooRealVar mu1_SFID   ("mu1_SFID", "mu1_SFID", 0., 100);
  RooRealVar mu2_SFID   ("mu2_SFID", "mu2_SFID", 0., 100);
  RooRealVar genWeightBkgB ("genWeightBkgB", "genWeightBkgB", 0, 5000);
  RooRealVar b_pt       ("b_pt", "b_pt", 0, 10000);
  RooRealVar b_mass       ("b_mass", "b_mass", 4.,6.);
  RooRealVar tau_pt       ("tau_pt", "tau_pt", 2, 10000);
  RooRealVar mu1_isLoose ("mu1_isLoose", "mu1_isLoose", 0.5, 2);
  RooRealVar mu2_isLoose ("mu2_isLoose", "mu2_isLoose", 0.5, 2);
 
 
  int Bc_bins = (int)(upper_limit - lower_limit)*10;

  ////1) Yuta: select the file you want
  TFile *ntuple_data = new TFile("/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_pt_2018/BJpsiX/bkg.root");
  //TFile *ntuple_data = new TFile("/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi/job_inv_pt_2018/Data/data.root"); 
  //=====
  
  TTree* tree_data   = (TTree*) ntuple_data->Get("tree");

  RooArgSet Variables(xgbs, tau_pt, mu1_isLoose, mu2_isLoose);
  Variables.add(mu1_SFID); Variables.add(mu2_SFID); Variables.add(genWeightBkgB); Variables.add(puweight); Variables.add(mu1_SFReco);  Variables.add(mu2_SFReco);
  Variables.add(b_mass);


  RooDataSet *data          = new RooDataSet("data", "data", tree_data, Variables);

  TCut SelectionCut = "tau_pt>3.&&mu1_isLoose==1&&mu2_isLoose==1";

  RooDataSet *cut_data0      = (RooDataSet*)data     ->reduce(SelectionCut);

  ////2) Yuta: select if it is data or MC
  //For MC:
  ////RooFormulaVar n1Func("n1Func","@0*@1*@2*@3*@4*@5*getWeight(@6)",RooArgList(puweight,mu1_SFReco,mu2_SFReco,mu1_SFID,mu2_SFID,genWeightBkgB,b_mass));     
  RooFormulaVar n1Func("n1Func","@0*@1*@2*@3*@4*@5",RooArgList(puweight,mu1_SFReco,mu2_SFReco,mu1_SFID,mu2_SFID,genWeightBkgB));
  RooRealVar* w = (RooRealVar*) cut_data0->addColumn(n1Func) ;
  RooDataSet cut_data(cut_data0->GetName(),cut_data0->GetTitle(),cut_data0,*cut_data0->get(),0,w->GetName()) ;

  //For data:
  //RooDataSet cut_data(cut_data0->GetName(),cut_data0->GetTitle(),cut_data0,*cut_data0->get()) ;  

  //=====
  TCut SelectionCut_numerator   = "xgbs>=4.0"; RooDataSet *cut_data_numerator = (RooDataSet*)cut_data.reduce(SelectionCut_numerator);
  TCut SelectionCut_denominator = "xgbs>2.5 && xgbs<3.5"; RooDataSet *cut_data_denominator = (RooDataSet*)cut_data.reduce(SelectionCut_denominator);

  //RooFormulaVar n1Func("n1Func","@0*@1*@2*@3*@4*@5*@6",RooArgList(puweight,mu1_SFReco,mu2_SFReco,mu1_SFID,mu2_SFID,genWeightBkgB,getWeight(b_mass)));
  std::cout<<"after cut num den definition "<<std::endl;
  // * * * * * * 
  // exponential
  RooRealVar        c1("c1", "c1", -1.5, -20., 0.);
  RooExponential mBkg1("mBkg1", "exponential", xgbs, c1);
  RooRealVar N_mBkg1 ("N_mBkg1", "N_mBkg1", 10., 0., 10000000.);
  RooExtendPdf e_mBkg1 ("e_mBkg1", "e_mBkg1", mBkg1, N_mBkg1);
  float  Range_fit_low =2.5, Range_fit_high=3.5;
  RooFitResult *fr = e_mBkg1.fitTo(cut_data, Range(Range_fit_low, Range_fit_high), NumCPU(4, kTRUE), Save(), Extended());

  std::cout<<"after fit exponential "<<std::endl;  
//  // * * * * * *
//  // polunomial p0
//  RooRealVar a00("a00", "a00", -10000, 10000);
//  RooPolynomial p0("p0","p0",xgbs,RooArgList(a00),0) ;
//  RooRealVar N_p0 ("N_p0", "N_p0", 10., 0., 1000000.);
//  RooExtendPdf e_p0 ("e_p0", "e_p0", p0, N_p0);
//  RooFitResult *frp0 = e_p0.fitTo(*cut_data, NumCPU(4, kTRUE), Save(), Extended());
//  // * * * * * *
//  // polunomial p5
  RooRealVar a0("a0", "a0",1, -100000, 100000);
  RooRealVar a31("a31", "a31",1, -10000, 10000);
  RooRealVar a32("a32", "a32",1, -10000, 10000);
  RooRealVar a33("a33", "a33",1, -10000, 10000);
  RooRealVar a34("a34", "a34",1, -10000, 10000);
  RooChebychev pol5("pol5","pol5",xgbs,RooArgList(a0,a31,a32,a33,a34)) ;
  RooRealVar N_p5 ("N_p5", "N_p5", 10., 0., 1000000.);
  RooExtendPdf e_p5 ("e_p5", "e_p5", pol5, N_p5);
  //  RooFitResult *frp5 = e_p5.fitTo(cut_data, NumCPU(4, kTRUE), Save(), Extended());

  RooKeysPdf kern01 ("kern01", "kern01", xgbs, cut_data, RooKeysPdf::NoMirror,0.6);
  RooKeysPdf kern02 ("kern02", "kern02", xgbs, cut_data, RooKeysPdf::MirrorLeft,0.6);
  RooKeysPdf kern03 ("kern03", "kern03", xgbs, cut_data, RooKeysPdf::MirrorRight,0.6);
  RooKeysPdf kern04 ("kern04", "kern04", xgbs, cut_data, RooKeysPdf::MirrorBoth,0.6);
  RooKeysPdf kern05 ("kern05", "kern05", xgbs, cut_data, RooKeysPdf::MirrorAsymLeft,0.6);
  RooKeysPdf kern06 ("kern06", "kern06", xgbs, cut_data, RooKeysPdf::MirrorAsymLeftRight,0.6);
  RooKeysPdf kern07 ("kern07", "kern07", xgbs, cut_data, RooKeysPdf::MirrorAsymRight, 0.6);
  RooKeysPdf kern08 ("kern08", "kern08", xgbs, cut_data, RooKeysPdf::MirrorLeftAsymRight,0.6);
  RooKeysPdf kern09 ("kern09", "kern09", xgbs, cut_data, RooKeysPdf::MirrorAsymBoth,0.6);
  std::cout<<"after RooKeysPdf kern09 "<<std::endl;
  RooKeysPdf kern11 ("kern11", "kern11", xgbs, cut_data, RooKeysPdf::NoMirror,0.4);
  RooKeysPdf kern12 ("kern12", "kern12", xgbs, cut_data, RooKeysPdf::MirrorLeft,0.4);
  RooKeysPdf kern13 ("kern13", "kern13", xgbs, cut_data, RooKeysPdf::MirrorRight,0.4);
  RooKeysPdf kern14 ("kern14", "kern14", xgbs, cut_data, RooKeysPdf::MirrorBoth,0.4);
  RooKeysPdf kern15 ("kern15", "kern15", xgbs, cut_data, RooKeysPdf::MirrorAsymLeft,0.4);
  RooKeysPdf kern16 ("kern16", "kern16", xgbs, cut_data, RooKeysPdf::MirrorAsymLeftRight,0.4);
  RooKeysPdf kern17 ("kern17", "kern17", xgbs, cut_data, RooKeysPdf::MirrorAsymRight,0.4);
  RooKeysPdf kern18 ("kern18", "kern18", xgbs, cut_data, RooKeysPdf::MirrorLeftAsymRight,0.4);
  RooKeysPdf kern19 ("kern19", "kern19", xgbs, cut_data, RooKeysPdf::MirrorAsymBoth,0.4);
  std::cout<<"after RooKeysPdf kern19 "<<std::endl;

  RooRealVar N_kern01 ("N_kern01", "N_kern01", 1000, 0, 10000000.); 
  RooRealVar N_kern02 ("N_kern02", "N_kern02", 1000, 0, 10000000.); 
  RooRealVar N_kern03 ("N_kern03", "N_kern03", 1000, 0, 10000000.); 
  RooRealVar N_kern04 ("N_kern04", "N_kern04", 1000, 0, 10000000.); 
  RooRealVar N_kern05 ("N_kern05", "N_kern05", 1000, 0, 10000000.); 
  RooRealVar N_kern06 ("N_kern06", "N_kern06", 1000, 0, 10000000.); 
  RooRealVar N_kern07 ("N_kern07", "N_kern07", 1000, 0, 10000000.); 
  RooRealVar N_kern08 ("N_kern08", "N_kern08", 1000, 0, 10000000.); 
  RooRealVar N_kern09 ("N_kern09", "N_kern09", 1000, 0, 10000000.); 

  RooRealVar N_kern11 ("N_kern11", "N_kern11", 1000, 0, 10000000.); 
  RooRealVar N_kern12 ("N_kern12", "N_kern12", 1000, 0, 10000000.); 
  RooRealVar N_kern13 ("N_kern13", "N_kern13", 1000, 0, 10000000.); 
  RooRealVar N_kern14 ("N_kern14", "N_kern14", 1000, 0, 10000000.); 
  RooRealVar N_kern15 ("N_kern15", "N_kern15", 1000, 0, 10000000.); 
  RooRealVar N_kern16 ("N_kern16", "N_kern16", 1000, 0, 10000000.); 
  RooRealVar N_kern17 ("N_kern17", "N_kern17", 1000, 0, 10000000.); 
  RooRealVar N_kern18 ("N_kern18", "N_kern18", 1000, 0, 10000000.); 
  RooRealVar N_kern19 ("N_kern19", "N_kern19", 1000, 0, 10000000.); 
  std::cout<<"after RooRealVar kern19 "<<std::endl;

  RooExtendPdf e_kern01 ("e_kern01", "e_kern01", kern01, N_kern01); 
  RooExtendPdf e_kern02 ("e_kern02", "e_kern02", kern02, N_kern02); 
  RooExtendPdf e_kern03 ("e_kern03", "e_kern03", kern03, N_kern03); 
  RooExtendPdf e_kern04 ("e_kern04", "e_kern04", kern04, N_kern04); 
  RooExtendPdf e_kern05 ("e_kern05", "e_kern05", kern05, N_kern05); 
  RooExtendPdf e_kern06 ("e_kern06", "e_kern06", kern06, N_kern06); 
  RooExtendPdf e_kern07 ("e_kern07", "e_kern07", kern07, N_kern07); 
  RooExtendPdf e_kern08 ("e_kern08", "e_kern08", kern08, N_kern08); 
  RooExtendPdf e_kern09 ("e_kern09", "e_kern09", kern09, N_kern09); 
  std::cout<<"after  RooExtendPdf kern09 "<<std::endl;
  RooExtendPdf e_kern11 ("e_kern11", "e_kern11", kern11, N_kern11); 
  RooExtendPdf e_kern12 ("e_kern12", "e_kern12", kern12, N_kern12); 
  RooExtendPdf e_kern13 ("e_kern13", "e_kern13", kern13, N_kern13); 
  RooExtendPdf e_kern14 ("e_kern14", "e_kern14", kern14, N_kern14); 
  RooExtendPdf e_kern15 ("e_kern15", "e_kern15", kern15, N_kern15); 
  RooExtendPdf e_kern16 ("e_kern16", "e_kern16", kern16, N_kern16); 
  RooExtendPdf e_kern17 ("e_kern17", "e_kern17", kern17, N_kern17); 
  RooExtendPdf e_kern18 ("e_kern18", "e_kern18", kern18, N_kern18); 
  RooExtendPdf e_kern19 ("e_kern19", "e_kern19", kern19, N_kern19); 
  std::cout <<" finishing definitio kerns " <<std::endl;
  // plot
  RooPlot *frame_main_fit1 = xgbs.frame(Title("xgbs 1 fit"), Bins(Bc_bins));
  cut_data.plotOn(frame_main_fit1, XErrorSize(0), Name("plotdata"));
  e_mBkg1.plotOn(frame_main_fit1, LineColor(kBlack), Name("totalpdf"));
  cut_data.plotOn(frame_main_fit1, XErrorSize(0));

  std::cout <<" finishing plotting fit 1 " <<std::endl;     

  RooPlot *frame_p5_fit = xgbs.frame(Title("xgbs plo5 fit"), Bins(Bc_bins));
  cut_data.plotOn(frame_p5_fit, XErrorSize(0), Name("plotdata"));
  cut_data.plotOn(frame_p5_fit, XErrorSize(0));

  RooPlot *frame_kern01_fit = xgbs.frame(Title("frame_kern01_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern01_fit, XErrorSize(0)); e_kern01.plotOn(frame_kern01_fit, LineColor(kBlack)); 
  RooPlot *frame_kern02_fit = xgbs.frame(Title("frame_kern02_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern02_fit, XErrorSize(0)); e_kern02.plotOn(frame_kern02_fit, LineColor(kBlack));
  RooPlot *frame_kern03_fit = xgbs.frame(Title("frame_kern03_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern03_fit, XErrorSize(0)); e_kern03.plotOn(frame_kern03_fit, LineColor(kBlack));
  RooPlot *frame_kern04_fit = xgbs.frame(Title("frame_kern04_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern04_fit, XErrorSize(0)); e_kern04.plotOn(frame_kern04_fit, LineColor(kBlack));
  RooPlot *frame_kern05_fit = xgbs.frame(Title("frame_kern05_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern05_fit, XErrorSize(0)); e_kern05.plotOn(frame_kern05_fit, LineColor(kBlack));
  RooPlot *frame_kern06_fit = xgbs.frame(Title("frame_kern06_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern06_fit, XErrorSize(0)); e_kern06.plotOn(frame_kern06_fit, LineColor(kBlack));
  RooPlot *frame_kern07_fit = xgbs.frame(Title("frame_kern07_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern07_fit, XErrorSize(0)); e_kern07.plotOn(frame_kern07_fit, LineColor(kBlack));
  RooPlot *frame_kern08_fit = xgbs.frame(Title("frame_kern08_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern08_fit, XErrorSize(0)); e_kern08.plotOn(frame_kern08_fit, LineColor(kBlack));
  RooPlot *frame_kern09_fit = xgbs.frame(Title("frame_kern09_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern09_fit, XErrorSize(0)); e_kern09.plotOn(frame_kern09_fit, LineColor(kBlack));

  RooPlot *frame_kern11_fit = xgbs.frame(Title("frame_kern11_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern11_fit, XErrorSize(0)); e_kern11.plotOn(frame_kern11_fit, LineColor(kBlack)); 
  RooPlot *frame_kern12_fit = xgbs.frame(Title("frame_kern12_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern12_fit, XErrorSize(0)); e_kern12.plotOn(frame_kern12_fit, LineColor(kBlack));
  RooPlot *frame_kern13_fit = xgbs.frame(Title("frame_kern13_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern13_fit, XErrorSize(0)); e_kern13.plotOn(frame_kern13_fit, LineColor(kBlack));
  RooPlot *frame_kern14_fit = xgbs.frame(Title("frame_kern14_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern14_fit, XErrorSize(0)); e_kern14.plotOn(frame_kern14_fit, LineColor(kBlack));
  RooPlot *frame_kern15_fit = xgbs.frame(Title("frame_kern15_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern15_fit, XErrorSize(0)); e_kern15.plotOn(frame_kern15_fit, LineColor(kBlack));
  RooPlot *frame_kern16_fit = xgbs.frame(Title("frame_kern16_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern16_fit, XErrorSize(0)); e_kern16.plotOn(frame_kern16_fit, LineColor(kBlack));
  RooPlot *frame_kern17_fit = xgbs.frame(Title("frame_kern17_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern17_fit, XErrorSize(0)); e_kern17.plotOn(frame_kern17_fit, LineColor(kBlack));
  RooPlot *frame_kern18_fit = xgbs.frame(Title("frame_kern18_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern18_fit, XErrorSize(0)); e_kern18.plotOn(frame_kern18_fit, LineColor(kBlack));
  RooPlot *frame_kern19_fit = xgbs.frame(Title("frame_kern19_fit"), Bins(Bc_bins)); cut_data.plotOn(frame_kern19_fit, XErrorSize(0)); e_kern19.plotOn(frame_kern19_fit, LineColor(kBlack));
  std::cout <<" finishing plotting kerns " <<std::endl; 

  TCanvas *c_xgbs_1 = new TCanvas("c_xgbs_1", "c_xgbs_1", 1200, 1200); c_xgbs_1->Divide(1,2);
  c_xgbs_1->cd(1);
  frame_main_fit1->GetXaxis()->SetNdivisions(504);
  frame_main_fit1->Draw(); frame_main_fit1->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_1->cd(1)->SetLogy(1);
  c_xgbs_1->cd(2);
  frame_p5_fit->GetXaxis()->SetNdivisions(504);
  frame_p5_fit->Draw(); frame_p5_fit->GetYaxis()->SetRangeUser(0.01, 5000);  c_xgbs_1 ->cd(2)->SetLogy(1);
  c_xgbs_1->SaveAs("plots/canvas_01.pdf");
  
  std::cout <<" finishing canvas 01 " <<std::endl; 
  TCanvas *c_xgbs_2 = new TCanvas("c_xgbs_2", "c_xgbs_2", 1200, 1200); c_xgbs_2->Divide(3,3);
  c_xgbs_2->cd(1); frame_kern01_fit->Draw(); frame_kern01_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(1)->SetLogy(1); 
  c_xgbs_2->cd(2); frame_kern02_fit->Draw(); frame_kern02_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(2)->SetLogy(1);
  c_xgbs_2->cd(3); frame_kern03_fit->Draw(); frame_kern03_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(3)->SetLogy(1);
  c_xgbs_2->cd(4); frame_kern04_fit->Draw(); frame_kern04_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(4)->SetLogy(1);
  c_xgbs_2->cd(5); frame_kern05_fit->Draw(); frame_kern05_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(5)->SetLogy(1);
  c_xgbs_2->cd(6); frame_kern06_fit->Draw(); frame_kern06_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(6)->SetLogy(1);
  c_xgbs_2->cd(7); frame_kern07_fit->Draw(); frame_kern07_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(7)->SetLogy(1);
  c_xgbs_2->cd(8); frame_kern08_fit->Draw(); frame_kern08_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(8)->SetLogy(1);
  c_xgbs_2->cd(9); frame_kern09_fit->Draw(); frame_kern09_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_2->cd(9)->SetLogy(1);
  c_xgbs_2->SaveAs("plots/canvas_02.pdf");

  std::cout <<" finishing canvas 02 " <<std::endl; 
  TCanvas *c_xgbs_3 = new TCanvas("c_xgbs_3", "c_xgbs_3", 1200, 1200); c_xgbs_3->Divide(3,3);
  c_xgbs_3->cd(1); frame_kern11_fit->Draw(); frame_kern11_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(1)->SetLogy(1); 
  c_xgbs_3->cd(2); frame_kern12_fit->Draw(); frame_kern12_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(2)->SetLogy(1);
  c_xgbs_3->cd(3); frame_kern13_fit->Draw(); frame_kern13_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(3)->SetLogy(1);
  c_xgbs_3->cd(4); frame_kern14_fit->Draw(); frame_kern14_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(4)->SetLogy(1);
  c_xgbs_3->cd(5); frame_kern15_fit->Draw(); frame_kern15_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(5)->SetLogy(1);
  c_xgbs_3->cd(6); frame_kern16_fit->Draw(); frame_kern16_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(6)->SetLogy(1);
  c_xgbs_3->cd(7); frame_kern17_fit->Draw(); frame_kern17_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(7)->SetLogy(1);
  c_xgbs_3->cd(8); frame_kern18_fit->Draw(); frame_kern18_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(8)->SetLogy(1);
  c_xgbs_3->cd(9); frame_kern19_fit->Draw(); frame_kern19_fit->GetYaxis()->SetRangeUser(0.01, 5000); c_xgbs_3->cd(9)->SetLogy(1);
  c_xgbs_3->SaveAs("plots/canvas_03.pdf");

////3) Yuta: select the regions boundaries.
  std::cout <<" reporting values " <<std::endl;
  xgbs.setRange("numerator", 4.0, 7.);
  xgbs.setRange("denominator", 2.5, 3.5);
  //=====
  cout << "Histogram: " << cut_data_numerator->sumEntries() << " / " << cut_data_denominator->sumEntries() << " = " << cut_data_numerator->sumEntries()/cut_data_denominator->sumEntries() << endl;

  cout << " - - - - - \n";
  RooAbsReal* a = e_mBkg1.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b = e_mBkg1.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "Exponential: " << a->getVal()*N_mBkg1.getVal() << " / " << b->getVal()*N_mBkg1.getVal() << " = " << (a->getVal() / b->getVal()) << endl;


  cout << " - - - - - \n";
  RooAbsReal* a_kern01 = e_kern01.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern01 = e_kern01.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern01: " << a_kern01->getVal()*N_kern01.getVal() << " / " << b_kern01->getVal()*N_kern01.getVal() << " = " << (a_kern01->getVal() / b_kern01->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern02 = e_kern02.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern02 = e_kern02.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern02: " << a_kern02->getVal()*N_kern02.getVal() << " / " << b_kern02->getVal()*N_kern02.getVal() << " = " << (a_kern02->getVal() / b_kern02->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern03 = e_kern03.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern03 = e_kern03.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern03: " << a_kern03->getVal()*N_kern03.getVal() << " / " << b_kern03->getVal()*N_kern03.getVal() << " = " << (a_kern03->getVal() / b_kern03->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern04 = e_kern04.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern04 = e_kern04.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern04: " << a_kern04->getVal()*N_kern04.getVal() << " / " << b_kern04->getVal()*N_kern04.getVal() << " = " << (a_kern04->getVal() / b_kern04->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern05 = e_kern05.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern05 = e_kern05.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern05: " << a_kern05->getVal()*N_kern05.getVal() << " / " << b_kern05->getVal()*N_kern05.getVal() << " = " << (a_kern05->getVal() / b_kern05->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern06 = e_kern06.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern06 = e_kern06.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern06: " << a_kern06->getVal()*N_kern06.getVal() << " / " << b_kern06->getVal()*N_kern06.getVal() << " = " << (a_kern06->getVal() / b_kern06->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern07 = e_kern07.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern07 = e_kern07.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern07: " << a_kern07->getVal()*N_kern07.getVal() << " / " << b_kern07->getVal()*N_kern07.getVal() << " = " << (a_kern07->getVal() / b_kern07->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern08 = e_kern08.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern08 = e_kern08.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern08: " << a_kern08->getVal()*N_kern08.getVal() << " / " << b_kern08->getVal()*N_kern08.getVal() << " = " << (a_kern08->getVal() / b_kern08->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern09 = e_kern09.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern09 = e_kern09.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern09: " << a_kern09->getVal()*N_kern09.getVal() << " / " << b_kern09->getVal()*N_kern09.getVal() << " = " << (a_kern09->getVal() / b_kern09->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern11 = e_kern11.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern11 = e_kern11.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern11: " << a_kern11->getVal()*N_kern11.getVal() << " / " << b_kern11->getVal()*N_kern11.getVal() << " = " << (a_kern11->getVal() / b_kern11->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern12 = e_kern12.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern12 = e_kern12.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern12: " << a_kern12->getVal()*N_kern12.getVal() << " / " << b_kern12->getVal()*N_kern12.getVal() << " = " << (a_kern12->getVal() / b_kern12->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern13 = e_kern13.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern13 = e_kern13.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern13: " << a_kern13->getVal()*N_kern13.getVal() << " / " << b_kern13->getVal()*N_kern13.getVal() << " = " << (a_kern13->getVal() / b_kern13->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern14 = e_kern14.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern14 = e_kern14.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern14: " << a_kern14->getVal()*N_kern14.getVal() << " / " << b_kern14->getVal()*N_kern14.getVal() << " = " << (a_kern14->getVal() / b_kern14->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern15 = e_kern15.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern15 = e_kern15.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern15: " << a_kern15->getVal()*N_kern15.getVal() << " / " << b_kern15->getVal()*N_kern15.getVal() << " = " << (a_kern15->getVal() / b_kern15->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern16 = e_kern16.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern16 = e_kern16.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern16: " << a_kern16->getVal()*N_kern16.getVal() << " / " << b_kern16->getVal()*N_kern16.getVal() << " = " << (a_kern16->getVal() / b_kern16->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern17 = e_kern17.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern17 = e_kern17.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern17: " << a_kern17->getVal()*N_kern17.getVal() << " / " << b_kern17->getVal()*N_kern17.getVal() << " = " << (a_kern17->getVal() / b_kern17->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern18 = e_kern18.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern18 = e_kern18.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern18: " << a_kern18->getVal()*N_kern18.getVal() << " / " << b_kern18->getVal()*N_kern18.getVal() << " = " << (a_kern18->getVal() / b_kern18->getVal()) << endl;

  cout << " - - - - - \n";
  RooAbsReal* a_kern19 = e_kern19.createIntegral(xgbs, NormSet(xgbs), Range("numerator"));
  RooAbsReal* b_kern19 = e_kern19.createIntegral(xgbs, NormSet(xgbs), Range("denominator"));
  cout << "kern19: " << a_kern19->getVal()*N_kern19.getVal() << " / " << b_kern19->getVal()*N_kern19.getVal() << " = " << (a_kern19->getVal() / b_kern19->getVal()) << endl;




// fit with polX
// fit with gaussian kernel 




}
