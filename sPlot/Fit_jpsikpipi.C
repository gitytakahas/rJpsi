// signal histo ===> "bg_res"
//
// please fit data_obs with parameterized signal histo + whatever function you want

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

//#define lowermcstats
//#define DrawResiduals

using namespace RooFit ;
using namespace RooStats ;

void fit_jpsikpipi();
void Fit_jpsikpipi() { fit_jpsikpipi(); }
void fit_jpsikpipi() {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  RooRealVar jpsi_kpipi ("jpsi_kpipi", "jpsi_kpipi", 5.0, 5.5);
  RooRealVar tau_pt ("tau_pt", "tau_pt", 3.0, 10000);
  RooRealVar mu1_isLoose ("mu1_isLoose", "mu1_isLoose", 0.5, 1.5);
  RooRealVar mu2_isLoose ("mu2_isLoose", "mu2_isLoose", 0.5, 1.5);
  RooRealVar jpsi_isrestos ("jpsi_isrestos", "jpsi_isrestos", 0.5, 1.5);
  RooRealVar procid ("procid", "procid", -100, 600);
  RooRealVar puweight ("puweight", "puweight", -10000, 100000);
  //  RooRealVar genWeightBkgB ("genWeightBkgB", "genWeightBkgB", -10000, 100000);
  RooRealVar xgbs ("xgbs", "xgbs", -15, 8);

  int Bc_bins = 50;

  TFile *ntuple_mc = new TFile("/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi_Validation/job_pt_BcJpsiTauNu/BcJpsiX_ul_2018/bkg.root");

  TTree* tree_mc   = (TTree*) ntuple_mc->Get("tree");

  RooArgSet Variables(jpsi_kpipi, tau_pt, mu1_isLoose, mu2_isLoose, jpsi_isrestos); //, jpsi_kpipi_psi2smass);
  Variables.add(procid); Variables.add(xgbs); 

#ifdef lowermcstats
  RooRealVar evt ("evt", "evt", 76000107, 5700010700);
  Variables.add(evt);
#endif

  //  RooDataSet *mc          = new RooDataSet("mc", "mc", tree_mc, Variables); //, puweight.GetName());
  RooDataSet *mc          = new RooDataSet("mc", "mc", tree_mc, Variables, puweight.GetName());

  TCut SelectionCut = "(procid==2 || procid==3 || procid==4 || procid==5 || procid==15 || procid==16 || procid==28 || procid==30 || procid==32 || procid==37 || procid==54 || procid==55)";

  RooDataSet *cut_mc    = (RooDataSet*)mc   ->reduce(SelectionCut);

  // exponential
  RooRealVar        c0("c0", "c0", -2.77542e+00); //-.5, -20., 0.);
  RooExponential mBkg0("mBkg0", "exponential", jpsi_kpipi, c0);
  RooRealVar N_mBkg0 ("N_mBkg0", "N_mBkg0", 10000., 0., 100000.);
  RooExtendPdf e_mBkg0 ("e_mBkg0", "e_mBkg0", mBkg0, N_mBkg0);

  //Error Function
  RooRealVar ErfSlope ("ErfSlope", "Erf Slope", 6.47636e+00); //5., 2., 10.);
  RooRealVar meanErf ("meanErf", "mean of the Erf gaussian", 5.30000e+00); //5.224, 5.100, 5.300);  //5300
  RooRealVar sigmaErf ("sigmaErf", "width of the Erf gaussian", 4.00000e-02); //0.0386, 0.010, 0.040); // 200
  RooRealVar ErfOffset ("ErfOffset", "Offset of Erf exponential", 5.07743e+00); //5.100, 5.000, 5.200); //4930.
  RooGenericPdf Erf("Erf", "Error Function", "TMath::Exp(TMath::Abs(ErfSlope)*(jpsi_kpipi-ErfOffset))*TMath::Erfc((jpsi_kpipi-meanErf)/sigmaErf)", RooArgSet(jpsi_kpipi, meanErf, sigmaErf,ErfSlope,ErfOffset));
//  RooGenericPdf Erf("Erf", "Error Function", "TMath::Erfc((jpsi_kpipi-meanErf)/sigmaErf)", RooArgSet(jpsi_kpipi, meanErf, sigmaErf, ErfOffset));

  RooRealVar N_Erf ("N_Erf", "N_Erf", 30000, 0, 100000);
  RooExtendPdf e_Erf ("e_Erf", "e_Erf", Erf, N_Erf);

  // signal
  RooRealVar mean_m  ("mean_m","mean of gaussian", 5.2749, 5.26, 5.3);
  RooRealVar sigma_m ("sigma_m","Scale Factor 1",  2.45912e-02); //0.034, 0.01, 0.04);
  RooGaussian mSig1   ("mSig1","signal p.d.f.", jpsi_kpipi, mean_m, sigma_m);
  RooRealVar sigma_m2 ("sigma_m2","Scale Factor 1",  0.08, 0.0, 0.2);
  RooGaussian mSig2   ("mSig2","signal p.d.f.", jpsi_kpipi, mean_m, sigma_m2);
  RooRealVar frac ("frac", "frac", 6.10512e-01); //0.8, 0., 1);

  RooRealVar sigma_cb ("sigma_cb","Scale Factor 1",  5.75354e-02); //0.034, 0.02, 0.075);
  RooRealVar alpha("alpha", "Alpha", 1.2); //5,0.1,10);
  RooRealVar n("n", "Order", 10);//6,0.1,10);
  RooCBShape  mSig3   ("mSig3","signa p.d.f.", jpsi_kpipi, mean_m, sigma_cb, alpha,n);

  RooAddPdf mSig ("mSig", "mSig", RooArgList(mSig1, mSig3), frac);
  RooRealVar N_mSig  ("N_mSig", "N_mSig", 80000., 0., 1000000.);
  RooExtendPdf e_mSig ("e_mSig", "e_mSig", mSig, N_mSig);

  RooRealVar frac1 ("frac1", "frac1", 1.80389e-01); //0.1, 0., 1.);
  RooRealVar frac2 ("frac2", "frac2", 4.66709e-01); //0.1, 0., 1.);
//  RooAddPdf total ("total", "total", RooArgSet(e_mBkg0, e_mSig, e_Erf), RooArgList(frac1, frac2));

  RooRealVar signalYield ("signalYield", "signalYield", 100000, 0, 100000);
  RooAddPdf totalp ("totalp", "totalp", RooArgSet(mBkg0, mSig, Erf), RooArgList(frac1, frac2));
  RooExtendPdf total ("total", "total", totalp, signalYield);

  RooFitResult *fr = total.fitTo(*cut_mc, NumCPU(4, kTRUE), Save(), Extended());

  // plot
  RooPlot *frame_main_fit1 = jpsi_kpipi.frame(Title("jpsi_kpipi 1 fit"), Bins(Bc_bins));
  cut_mc->plotOn(frame_main_fit1, XErrorSize(0), Name("plotmc"));
  total.plotOn(frame_main_fit1, LineColor(kBlack), Name("totalpdf"));
  total.plotOn(frame_main_fit1, Components("e_mBkg0"), LineColor(kAzure+2));
  total.plotOn(frame_main_fit1, Components("e_Erf"), LineColor(kAzure+2));
  total.plotOn(frame_main_fit1, Components(e_mSig), LineColor(kRed));
  cut_mc->plotOn(frame_main_fit1, XErrorSize(0));

  TCanvas *c_jpsi_kpipi_1 = new TCanvas("c_jpsi_kpipi_1", "c_jpsi_kpipi_1", 1200, 1200); //c_jpsi_kpipi_1->Divide(2,2);
  frame_main_fit1->GetXaxis()->SetNdivisions(504);
  frame_main_fit1->Draw(); //c_jpsi_kpipi_1->cd(1)->SetLogy(1);
  c_jpsi_kpipi_1->SaveAs("plots/canvas_jpsi_kpipi.pdf");

  // ** ** ** ** 
  // start the data fit

  TFile *ntuple_data = new TFile("/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi_Validation/job_pt_BcJpsiTauNu/Data_2018/data4Stefanos.root");
  //TFile *ntuple_data = new TFile("/pnfs/psi.ch/cms/trivcat/store/user/ytakahas/RJpsi_Validation/job_pt_BcJpsiTauNu/Data_2018/data_training_v2.root");
  TTree* tree_data   = (TTree*) ntuple_data->Get("tree");

  RooArgSet Variablesd(jpsi_kpipi, tau_pt, mu1_isLoose, mu2_isLoose, jpsi_isrestos, xgbs); 

  RooDataSet *data = new RooDataSet("data", "data", tree_data, Variablesd); //, puweight.GetName());

  // signal template!
  frac1.setConstant();
  frac2.setConstant();

  c0.setConstant();

  ErfSlope .setConstant(); 
  meanErf  .setConstant();
  sigmaErf .setConstant();
  ErfOffset.setConstant();

//  RooRealVar mean_m2  ("mean_m2","mean of gaussian", 5.2749, mean_m.getVal() - 2.*mean_m.getError(), mean_m.getVal() + 2.*mean_m.getError());
//  mean_m  .setConstant();
  sigma_m .setConstant();
  sigma_cb.setConstant(); 
  alpha   .setConstant(); 
  n       .setConstant();
  frac    .setConstant();

  RooRealVar signalYield_data ("signalYield_data", "signalYield_data", 100000, 0, 100000);
  RooAddPdf total_signal_data ("total_signal_data", "total_signal_data", RooArgSet(mBkg0, mSig, Erf), RooArgList(frac1, frac2));
  RooExtendPdf e_total_signal_data ("e_total_signal_data", "e_total_signal_data", total_signal_data, signalYield_data);

  // background
  RooRealVar        c1("c1", "c1", -.5, -20., 10.);
  RooExponential mBkg1("mBkg1", "exponential", jpsi_kpipi, c1);
  RooRealVar N_mBkg1 ("N_mBkg1", "N_mBkg1", 10000., 0., 100000.);
  RooExtendPdf e_mBkg1 ("e_mBkg1", "e_mBkg1", mBkg1, N_mBkg1);

  //Error Function
  RooRealVar    data_ErfSlope   ("data_ErfSlope", "Erf Slope", 5., 2., 10.);
  RooRealVar    data_meanErf    ("data_meanErf", "mean of the Erf gaussian", 5.224, 5.100, 5.300);  //5300
  RooRealVar    data_sigmaErf   ("data_sigmaErf", "width of the Erf gaussian", 0.0386, 0.010, 0.040); // 200
  RooRealVar    data_ErfOffset  ("data_ErfOffset", "Offset of Erf exponential", 5.100, 5.000, 5.200); //4930.
  RooGenericPdf data_Erf        ("data_Erf", "Error Function", "TMath::Exp(TMath::Abs(data_ErfSlope)*(jpsi_kpipi-data_ErfOffset))*TMath::Erfc((jpsi_kpipi-data_meanErf)/data_sigmaErf)", RooArgSet(jpsi_kpipi, data_meanErf, data_sigmaErf, data_ErfSlope, data_ErfOffset));
//  RooGenericPdf Erf("Erf", "Error Function", "TMath::Erfc((jpsi_kpipi-meanErf)/sigmaErf)", RooArgSet(jpsi_kpipi, meanErf, sigmaErf, ErfOffset));
  RooRealVar N_data_Erf ("N_data_Erf", "N_data_Erf", 1000, 0, 1000000);
  RooExtendPdf e_data_Erf ("e_data_Erf", "e_data_Erf", data_Erf, N_data_Erf);

  RooAddPdf total_data ("total_data", "total_data", RooArgSet(e_total_signal_data, e_mBkg1, e_data_Erf));

  RooFitResult *fr_data = total_data.fitTo(*data, NumCPU(4, kTRUE), Save(), Extended());

  // plot
  RooPlot *frame_data_fit = jpsi_kpipi.frame(Title("jpsi_kpipi 1 fit"), Bins(Bc_bins));
  data->plotOn(frame_data_fit, XErrorSize(0), Name("plotmc"));
  total_data.plotOn(frame_data_fit, LineColor(kAzure+1), Name("total_data"));
  total_data.plotOn(frame_data_fit, Components("e_total_signal_data"), LineColor(kAzure+2));
  total_data.plotOn(frame_data_fit, Components("e_mBkg1"), LineColor(kAzure+3));
  total_data.plotOn(frame_data_fit, Components("e_data_Erf"), LineColor(kAzure+4));
//  total_data.plotOn(frame_data_fit, Components("e_data_Erf"), LineColor(kAzure+2));
  data->plotOn(frame_data_fit, XErrorSize(0));

  TCanvas *c_data = new TCanvas("c_data", "c_data", 1200, 1200); //c_data->Divide(2,2);
  c_data->cd(1);
  frame_data_fit->GetXaxis()->SetNdivisions(504);
  frame_data_fit->Draw(); //c_data->cd(1)->SetLogy(1);
  c_data->SaveAs("plots/canvas_jpsi_kpipi.pdf");

  // ******
  // splot

  RooWorkspace* ws = new RooWorkspace("myWS");
  ws->import(total_data);
  ws->import(*data, Rename("sPlotdata"));

  RooAbsPdf* sPlottotal_data = ws->pdf("total_data");
  RooRealVar* SignalPromptYield  = ws->var("signalYield_data");
  RooRealVar* Signal2PromptYield = ws->var("N_mBkg1");
  RooRealVar* Signal3PromptYield = ws->var("N_data_Erf");

  SignalPromptYield->setConstant();
  Signal2PromptYield->setConstant();
  Signal3PromptYield->setConstant();

  RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, sPlottotal_data, RooArgList(*SignalPromptYield, *Signal2PromptYield, *Signal3PromptYield));

  RooDataSet * data_weighted = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "signalYield_data_sw");

//  double weight = signalYield_data.getVal()/signalYield.getVal();
  double weight = signalYield_data.getVal()/cut_mc->sumEntries();
//  RooPlot *frame_xgbs_gc = xgbs.frame(Title("xgbs"), Bins(26));

  RooPlot *frame_xgbs_weighted = xgbs.frame(Title("xgbs"), Bins(26));
  data_weighted->plotOn(frame_xgbs_weighted, XErrorSize(0));
  cut_mc->plotOn(frame_xgbs_weighted, XErrorSize(0), Rescale(weight), MarkerColor(kAzure+1));

  TCanvas *c = new TCanvas("c", "c", 900, 900); 
  c->cd(); //frame_xgbs_weighted->GetXaxis()->SetRangeUser(0.1, 2.5);  frame_xgbs_weighted->GetYaxis()->SetRangeUser(0., 10.);  

//#ifdef DrawResiduals
//  TH1* tmp1 = data_weighted->createHistogram("xgbs",24) ;
//  TH1* tmp2 = cut_mc->createHistogram("xgbs",24) ; tmp2->Scale(weight);
//  TH1F *diff = new TH1F("diff", "diff", 24, -20, 20);
//  diff->Divide(tmp1, tmp2, 1., 1.);
//
//  TPad *pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.05,1.0,1.0,21);
//  TPad *pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.1,22);
//  pad1->SetFillColor(0); pad2->SetFillColor(0);
//  pad1->Draw();
//  pad2->Draw();
//  pad2->cd(); diff->Draw(); diff->GetYaxis()->SetRangeUser(0, 2);
//  TF1 *f_straighline = new TF1("f_straighline", "1", -20, 20);
//  f_straighline->SetLineColor(kBlack);
//  f_straighline->Draw("same");
//  pad1->cd();
//#endif

  frame_xgbs_weighted->Draw();
  frame_xgbs_weighted->GetYaxis()->SetRange(-50, 1200);
  c->SaveAs("~/Desktop/xgbs.pdf");

  jpsi_kpipi.setRange("signal", 0, 10) ;
  jpsi_kpipi.setRange("signal2",5.,5.5) ;
  jpsi_kpipi.setRange("signal3", 0, 100) ;
  jpsi_kpipi.setRange("signal4", 4.5, 6.0) ;

//  RooAbsReal* fracInt = e_total_signal_data.createIntegral(jpsi_kpipi, Range("signal")) ;  cout << "first integral: " << fracInt->getVal() << endl;
  RooAbsReal* fracIn3 = e_total_signal_data.createIntegral(jpsi_kpipi, Range("signal3")) ; cout << "first integral: " << fracIn3->getVal()*data->sumEntries() << endl;
//  RooAbsReal* fracIn4 = e_total_signal_data.createIntegral(jpsi_kpipi, Range("signal4")) ; cout << "first integral: " << fracIn4->getVal() << endl;
//  RooAbsReal* fracIn3 = e_total_signal_data.createIntegral(jpsi_kpipi, NormSet(jpsi_kpipi), Range("signal2")) ; cout << "first integral: " << fracIn3->getVal() << endl;
//  cout << " - - - - - - - \n";
//  RooAbsReal* fracIn2  = e_total_signal_data.createIntegral(jpsi_kpipi, Range("signal2")) ; cout << "first integral: " << fracIn2->getVal() << endl;
//  RooAbsReal* fracIn20 = e_mBkg1.createIntegral(jpsi_kpipi, Range("signal2")) ; cout << "first integral: " << fracIn20->getVal() << endl;
//  RooAbsReal* fracIn21 = e_data_Erf.createIntegral(jpsi_kpipi, Range("signal2")) ; cout << "first integral: " << fracIn21->getVal() << endl;
}
