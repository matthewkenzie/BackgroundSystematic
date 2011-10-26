#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "Rtypes.h"

#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"
#include "RooFit.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"

using namespace std;
using namespace RooFit;

//string fcnNames[12] = {"exp1","exp3","exp5","pow1","pow3","pow5","lau1","lau3","lau5","pol1","pol3","pol5"};

// in and out files
TFile *inFile = new TFile("CMS-HGG_1665pb.root");
TFile *outFile = new TFile("BkgSystOut.root","RECREATE");

// get workspace and mass data
RooWorkspace *dataWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
RooRealVar *mass = (RooRealVar*)dataWS->var("CMS_hgg_mass");

// -------------- Global fit parameters ----------------------
  // dummy varaible for 1
  RooRealVar dumY("dumY","dumY",1.0,1.0,1.0);
  // single exp
  RooRealVar exp_p0("exp_p0","exp_p0",-0.02,-10.,10.);
  //RooGenericPdf sin_exp("sin_exp","sin_exp","exp(@1*@0)",RooArgSet(*mass,exp_p0));
  RooExponential temp_sin_exp("temp_sin_exp","temp_sin_exp",*mass,exp_p0);
  // double exp
  RooRealVar exp_p1("exp_p1","exp_p1",-0.02,-10.,10.); //find starting vals
  RooRealVar exp_p2("exp_p2","exp_p2",-0.02,-10.,10.); //find starting vals
  //RooGenericPdf dbl_exp("dbl_exp","dbl_exp","exp(@1*@0)+@2*exp(@3*@0)",RooArgSet(*mass,exp_p0,exp_p1,exp_p2));
  // triple exp
  RooRealVar exp_p3("exp_p3","exp_p3",-0.02,-10.,10.); //find starting vals
  RooRealVar exp_p4("exp_p4","exp_p4",-0.02,-10.,10.); //find starting vals
  RooGenericPdf trip_exp("trip_exp","trip_exp","exp(@1*@0)+@2*exp(@3*@0)+@4*exp(@5*@0)",RooArgSet(*mass,exp_p0,exp_p1,exp_p2,exp_p3,exp_p4));
  // single pow
  RooRealVar pl_p0("pl_p0","pl_p0",-0.01,-10.,10.);
  RooGenericPdf temp_sin_pow("temp_sin_pow","temp_sin_pow","pow(@0,@1)",RooArgSet(*mass,pl_p0));
  // double pow
  RooRealVar pl_p1("pl_p1","pl_p1",-0.01,-10.,10.);
  RooRealVar pl_p2("pl_p2","pl_p2",-0.01,-10.,10.);
  RooGenericPdf dbl_pow("dbl_pow","dbl_pow","pow(@0,@1)+@2*pow(@0,@3)",RooArgSet(*mass,pl_p0,pl_p1,pl_p2));
  // triple pow
  RooRealVar pl_p3("pl_p3","pl_p3",-0.01,-10.,10.);
  RooRealVar pl_p4("pl_p4","pl_p4",-0.01,-10.,10.);
  RooGenericPdf trip_pow("trip_pow","trip_pow","pow(@0,@1)+@2*pow(@0,@3)+@4*pow(@0,5)",RooArgSet(*mass,pl_p0,pl_p1,pl_p2,pl_p3,pl_p4));
  // two term laurent
  RooRealVar lau_p0("lau_p0","lau_p0",-0.01,-10.,0.);
  RooGenericPdf two_lau("two_lau","two_lau","pow(@0,-4.0)+@1*pow(@0,-5.0)",RooArgSet(*mass,lau_p0));
  // four term laurent
  RooRealVar lau_p1("lau_p1","lau_p1",-0.01,-10.,10.);
  RooRealVar lau_p2("lau_p2","lau_p2",-0.01,-10.,10.);
  RooGenericPdf four_lau("four_lau","four_lau","pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)",RooArgSet(*mass,lau_p0,lau_p1,lau_p2));
  // six term laurent
  RooRealVar lau_p3("lau_p3","lau_p3",-0.01,-10.,10.);
  RooRealVar lau_p4("lau_p4","lau_p4",-0.01,-10.,10.);
  RooGenericPdf six_lau("six_lau","six_lau","pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)+@4*pow(@0,-2.0)+@5*pow(@0,-7.0)",RooArgSet(*mass,lau_p0,lau_p1,lau_p2,lau_p3,lau_p4));
  // single poly
  RooRealVar pol_p0("pol_p0","pol_p0",-0.01,-10.,0.);
  RooGenericPdf sin_pol("sin_pol","sin_pol","1.0+(@1*@0)",RooArgSet(*mass,pol_p0));
  // double poly
  RooRealVar pol_p1("pol_p1","pol_p1",-0.01,-10.,10.);
  RooRealVar pol_p2("pol_p2","pol_p2",-0.01,-10.,10.);
  RooGenericPdf dbl_pol("dbl_pol","dbl_pol","1.0+(@1*@0)+(@2*pow(@0,2.0))+(@3*pow(@0,3.0))",RooArgSet(*mass,pol_p0,pol_p1,pol_p2));
  // triple poly
  RooRealVar pol_p3("pol_p3","pol_p3",-0.01,-10.,10.);
  RooRealVar pol_p4("pol_p4","pol_p4",-0.01,-10.,10.);
  RooGenericPdf trip_pol("trip_pol","trip_pol","1.0+(@1*@0)+(@2*pow(@0,2.0))+(@3*pow(@0,3.0))+(@4*pow(@0,4.0))+(@5*pow(@0,5.0))",RooArgSet(*mass,pol_p0,pol_p1,pol_p2,pol_p3,pol_p4));
// ----------------------------------------------------------

// ---------- functions which get Pdf for each function ------

RooAddPdf getExp(int nFs=0){
  if (nFs<0 || nFs>2) {
    cout << "Wrong number of params" << endl;
    cout << "Returning default: single exponential" << endl;
    nFs=0;
  }
  // single exp
  RooAddPdf sin_exp("sin_exp","sin_exp",RooArgList(temp_sin_exp),RooArgList(RooConst(1)));
  //RooGenericPdf sin_exp = (RooGenericPdf)temp_sin_exp;
  if (nFs==0) return sin_exp;
  // double exp
  //RooAddPdf dbl_exp("dbl_exp","dbl_exp",RooArgList(sin_exp,temp_exp_1),RooArgList(dumY,exp_p1));
  //if (nFs==1) return dbl_exp;
  // triple exp 
  //RooAddPdf trip_exp("dbl_exp","dbl_exp",RooArgList(dbl_exp,temp_exp_2),RooArgList(dumY,exp_p3));
  //if (nFs==2) return trip_exp;
}

RooAddPdf getPow(int nFs=0){
  if (nFs<0 || nFs>2) {
    cout << "Wrong number of params" << endl;
    cout << "Returning default: single power law" << endl;
    nFs=0;
  }
  // single pow
  RooAddPdf sin_pow("sin_pow","sin_pow",RooArgList(temp_sin_pow),RooArgList(dumY));
  if (nFs==0) return sin_pow;
  // double pow
  //RooAddPdf dbl_pow("dbl_pow","dbl_pow",RooArgList(sin_pow,temp_pl_1),RooArgList(dumY,pl_p1));
  //if (nFs==1) return dbl_pow;
  // triple pow
  //RooAddPdf trip_pow("trip_pow","trip_pow",RooArgList(sin_pow,temp_pl_2),RooArgList(dumY,pl_p3));
  //if (nFs==2) return trip_pow;
}

RooGenericPdf getLau(int nFs=0){
  if (nFs<0 || nFs>2) {
    cout << "Wrong number of params" << endl;
    cout << "Returning default: two term laurent" << endl;
  nFs=0;
  }
  // two term laurent
  //RooAddPdf two_lau("two_lau","two_lau",RooArgList(temp_lau_0),RooArgList(dumY));
  if (nFs==0) return two_lau;
  // four term laurent
  //RooAddPdf four_lau("four_lau","four_lau",RooArgList(two_lau,temp_lau_1),RooArgList(dumY,dumY));
  if (nFs==1) return four_lau;
  // six term laurent
  //RooAddPdf six_lau("six_lau","six_lau",RooArgList(four_lau,temp_lau_2),RooArgList(dumY,dumY));
  if (nFs==2) return six_lau;
}

RooGenericPdf getPol(int nFs=0){
  if (nFs<0 || nFs>2) {
    cout << "Wrong number of params" << endl;
    cout << "Returning default: single polynomial" << endl;
    nFs=0;
  }
  // single polynomial
  //RooAddPdf sin_pol("sin_pol","sin_pol",RooArgList(temp_pol_0),RooArgList(dumY));
  if (nFs==0) return sin_pol;
  // double polynomial
  //RooAddPdf dbl_pol("dbl_pol","dbl_pol",RooArgList(sin_pol,temp_pol_1),RooArgList(dumY,dumY));
  if (nFs==1) return dbl_pol;
  // double polynomial
  //RooAddPdf trip_pol("trip_pol","trip_pol",RooArgList(dbl_pol,temp_pol_2),RooArgList(dumY,dumY));
  if (nFs==2) return trip_pol;
}
// ----------------------------------------------------------

// ----------- plotting functions ---------------------------
void Plot(RooRealVar *mass, RooDataSet* hMassData, RooAddPdf dataFitFcn, RooDataSet* genDat, RooAddPdf sigAndBkg, RooAddPdf genFitFcn, int MCmass, int cat, int toy){
 
  TCanvas *canv = new TCanvas();

  RooPlot *mFrame = mass->frame(Title(Form("Fit of %s to data cat %d",dataFitFcn.GetName(),cat)));
  hMassData->plotOn(mFrame,DataError(RooDataSet::SumW2));
  dataFitFcn.plotOn(mFrame);
  mFrame->Draw();
  canv->Print(Form("plots/Fit_%s_cat%d.pdf",dataFitFcn.GetName(),cat),"pdf");
  canv->Clear();

  RooPlot *gFrame = mass->frame(Title(Form("Gen data from %s fitted with %s plus signal at M %d in cat %d toy %d",dataFitFcn.GetName(),genFitFcn.GetName(),MCmass,cat,toy)));
  genDat->plotOn(gFrame,DataError(RooDataSet::SumW2));
  sigAndBkg.plotOn(gFrame);
  gFrame->Draw();
  canv->Print(Form("plots/Gen_%s_fit_%s_m%d_cat%d_toy%d.pdf",dataFitFcn.GetName(),genFitFcn.GetName(),MCmass,cat,toy),"pdf");

  delete canv;
}

void histPlot(TH1F *hist){

  TCanvas *canv = new TCanvas();
  hist->Draw();
  canv->Print(Form("plots/%s.pdf",hist->GetName()),"pdf");
  delete canv;
}
void histPlot(TH2F *hist){

  gStyle->SetPalette(1);
  TCanvas *canv = new TCanvas();
  hist->Draw("colz");
  canv->Print(Form("plots/%s.pdf",hist->GetName()),"pdf");
  delete canv;
}

// ----------------------------------------------------------

// ---- utility function which requests correct function ----

RooAddPdf getFunction(string name){
  
  // exponentials
  if (name=="sin_exp") { 
    RooAddPdf returnFcn = getExp(0);
    return returnFcn;
  }
  /*
  else if (name=="dbl_exp") { 
    RooGenericPdf returnFcn = getExp(1);
    return returnFcn;
  }
  else if (name=="trip_exp") { 
    RooGenericPdf returnFcn = getExp(2);
    return returnFcn;
  }
  */
  // power laws
  else if (name=="sin_pow") { 
    RooAddPdf returnFcn = getPow(0);
    return returnFcn;
  }
  /*
  else if (name=="dbl_pow") { 
    RooGenericPdf returnFcn = getPow(1);
    return returnFcn;
  }
  else if (name=="trip_pow") { 
    RooGenericPdf returnFcn = getPow(2);
    return returnFcn;
  }
  // laurent series
  else if (name=="two_lau") { 
    RooGenericPdf returnFcn = getLau(0);
    return returnFcn;
  }
  else if (name=="four_lau") { 
    RooGenericPdf returnFcn = getLau(1);
    return returnFcn;
  }
  else if (name=="six_lau") { 
    RooGenericPdf returnFcn = getLau(2);
    return returnFcn;
  }
  // polynomials
  else if (name=="sin_pol") { 
    RooGenericPdf returnFcn = getPol(0);
    return returnFcn;
  }
  else if (name=="dbl_pol") { 
    RooGenericPdf returnFcn = getPol(1);
    return returnFcn;
  }
  else if (name=="trip_pol") { 
    RooGenericPdf returnFcn = getPol(2);
    return returnFcn;
  }
  */
  else {
    cout << "ERROR: this function name is not recognised: " << name << endl;
    cout << "Existing program" << endl;
    exit(1);
  }
}

void checkInput(string name){
  if (name!="sin_exp" && name!="dbl_exp" && name!="trip_exp" && name!="sin_pow" && name!="dbl_pow" && name!="trip_pow" && name!="two_lau" && name!="four_lau" && name!="six_lau" && name!="sin_pol" && name!="dbl_pol" && name!="trip_pol") 
    cout << name << " is not a valid function" << endl;
  assert(name=="sin_exp" || name=="dbl_exp" || name=="trip_exp" || name=="sin_pow" || name=="dbl_pow" || name=="trip_pow" || name=="two_lau" || name=="four_lau" || name=="six_lau" || name=="sin_pol" || name=="dbl_pol" || name=="trip_pol");  
}

int main(int argc, char* argv[]){
  
  if (argc!=7) cout << "Not enough arguments" << endl;
  assert(argc==7);
  string fitFunc(argv[1]);
  string genFunc(argv[2]);
  checkInput(fitFunc);
  checkInput(genFunc);

  ofstream diagFile("diags.txt");
  const int nToys=atoi(argv[3]);
  const int nCats=atoi(argv[4]) ;
  const int nMasses=9;
  const int mRlow=atoi(argv[5]);
  const int mRhigh=atoi(argv[6]);
  cerr << "Running " << nToys << " toys over " << nCats << " catgeories around " << fitFunc << " fit to data and fitting this with " << genFunc << endl;
  cerr << "Signal mass range = [" << mRlow << "-" << mRhigh << "]" << endl;
  TH1F *sigNormHist[nMasses][nCats];
  TH1F *sigYieldHist[nMasses][nCats];
  TH1F *bkgNormHist[nMasses][nCats];
  TH2F *bkgSigNormCorr[nMasses][nCats];

  for (int cat=0; cat<nCats; cat++){
    cerr << "Category " << cat << endl;
    int mIt=0;
    for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
      if (mMC==145) continue;
      // ------- make histos for bias study -------
      sigNormHist[mIt][cat] = new TH1F(Form("sigNorm_m%d_cat%d",mMC,cat),Form("sigNorm_m%d_cat%d",mMC,cat),200,-50,50);
      sigYieldHist[mIt][cat] = new TH1F(Form("sigYield_m%d_cat%d",mMC,cat),Form("sigYield_m%d_cat%d",mMC,cat),200,-50,50);
      bkgNormHist[mIt][cat] = new TH1F(Form("bkgNorm_m%d_cat%d",mMC,cat),Form("bkgNorm_m%d_cat%d",mMC,cat),200,-50,50);
      bkgSigNormCorr[mIt][cat] = new TH2F(Form("bkgSigCorr_m%d_cat%d",mMC,cat),Form("bkgSigCorr_m%d_cat%d",mMC,cat),20,-50,50,20,-50,50);
      bkgSigNormCorr[mIt][cat]->GetXaxis()->SetTitle("sigNorm");
      bkgSigNormCorr[mIt][cat]->GetYaxis()->SetTitle("bkgNorm");
      // ------------------------------------------
      mIt++;
    }
  }

  diagFile << setw(6) << "Mass" << setw(6) << "Cat" << setw(10) << "datNevt" << setw(10) << "datFitInt" << setw(10) << "genNevt" << setw(10) << "genFitInt" << setw(10) << "sigNorm" << setw(10) << "sigYield" << setw(10) << "bkgNormFit" << setw(10) << "bkgNormGen" << setw(10) << "bkgDiff" << setw(10) << "s+bNorm" << endl;
    for (int cat=0; cat<nCats; cat++){

      // Get data for each cat and fit with function
      RooDataSet *hMassData = (RooDataSet*)dataWS->data(Form("data_mass_cat%d",cat));
      RooAddPdf dataFitFcn = getFunction(fitFunc); //getExp(0);
      dataFitFcn.fitTo(*hMassData);

      // generate toy data from this
      for (int itToy=0; itToy<nToys; itToy++){
        if (itToy%200==0) cerr << Form("%2f %% of toys thown",double(itToy)/double(nToys)) << endl;
        RooDataSet *genDat = dataFitFcn.generate(*mass,hMassData->numEntries(),Extended());
        // fit to sig + bkg
        int mIt=0;
        for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
          if (mMC==145) continue;
          // define sig region
          double lowBand = 0.98*double(mMC);
          double highBand = 1.02*double(mMC);
          RooRealVar intRange(*mass);
          intRange.setRange("sigWindow",lowBand,highBand);
          RooRealVar wholeRange(*mass);
          wholeRange.setRange("wholeRange",100,160);
          // find norm of bkg in sig region
          RooAbsReal* int_dataFitFcn = dataFitFcn.createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
          double bkgFitInt = int_dataFitFcn->getVal()*hMassData->numEntries();
          RooAbsReal* wholeInt_dataFitFcn = dataFitFcn.createIntegral(*mass,NormSet(*mass),Range("wholeRange"));
          double bkgIntegral = wholeInt_dataFitFcn->getVal()*hMassData->numEntries();

          RooDataHist *sigMCHist = (RooDataHist*)dataWS->data(Form("roohist_sig_mass_m%d_cat%d",mMC,cat));
          RooAddPdf genFitFcn = getFunction(genFunc); //getPow(0);
          RooRealVar bkgYield("nBkg","nBkg",1001.,1000.,10000.);
          RooRealVar sigYield("nSig","nSig",0.01,-50.,50.);
          RooHistPdf sigMC("sigMC","sigMC",*mass,*sigMCHist);
          RooAddPdf sigAndBkg(Form("sigAndBkg%d",mMC),Form("sigAndBkg%d",mMC),RooArgList(genFitFcn,sigMC),RooArgList(bkgYield,sigYield));
          sigAndBkg.fitTo(*genDat);

          RooArgSet n(*mass);
          RooAbsReal* int_genFitFcn = genFitFcn.createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
          double bkgGenInt = (int_genFitFcn->getVal()*bkgYield.getVal());
          RooAbsReal* int_sigMC = sigMC.createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
          double sigInt = (int_sigMC->getVal()*sigYield.getVal());
          RooAbsReal* wholeInt_sigAndBkg = sigAndBkg.createIntegral(*mass,NormSet(*mass),Range("wholeRange"));
          double sandbIntegral = wholeInt_sigAndBkg->getVal()*genDat->numEntries();
          RooAbsReal* windowInt_sigAndBkg = sigAndBkg.createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
          double sandbWinIntegral = windowInt_sigAndBkg->getVal()*genDat->numEntries();
          
          // --- print to diagnostic file -------
          diagFile << setw(6) << mMC << setw(6) << cat << setw(10) << hMassData->numEntries() << setw(10) << bkgIntegral << setw(10) << genDat->numEntries() << setw(10) << sandbIntegral << setw(10) << sigInt << setw(10) << sigYield.getVal() << setw(10) << bkgFitInt << setw(10) << bkgGenInt << setw(10) << bkgFitInt-bkgGenInt << setw(10) << sandbWinIntegral << endl;
          // -----------------------------------

          sigNormHist[mIt][cat]->Fill(sigInt);
          sigYieldHist[mIt][cat]->Fill(sigYield.getVal());
          bkgNormHist[mIt][cat]->Fill(bkgFitInt-bkgGenInt);
          bkgSigNormCorr[mIt][cat]->Fill(sigInt,bkgFitInt-bkgGenInt);

          if (itToy==nToys-1) Plot(mass,hMassData,dataFitFcn,genDat,sigAndBkg,genFitFcn,mMC,cat,itToy);
        mIt++;
        }
      }
      
      int mIt=0;
      for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
        if (mMC==145) continue;
        histPlot(sigNormHist[mIt][cat]);
        histPlot(sigYieldHist[mIt][cat]);
        histPlot(bkgNormHist[mIt][cat]);
        histPlot(bkgSigNormCorr[mIt][cat]);
        outFile->cd();
        sigNormHist[mIt][cat]->Write();
        sigYieldHist[mIt][cat]->Write();
        bkgNormHist[mIt][cat]->Write();
        bkgSigNormCorr[mIt][cat]->Write();
        mIt++;
      }
    }
    
    inFile->Close();
    outFile->Close();
    diagFile.close();

    return 0;
  }

  /*
  void fitSingExp(RooDataSet *hMassData, RooDataHist *sigMCHist, RooRealVar *mass, int cat){
    
    // single exponential
    RooRealVar exp_p1("exp_p1","exp_p1",-0.02,-1.,0.);
  RooExponential exp_1("exp1","exp1",*mass,exp_p1);
  RooRealVar expY("expY","expY",1.0,1.0,1.0);
  RooAddPdf exp_pdf("exp_pdf","exp_pdf",RooArgList(exp_1),RooArgList(expY));

  // single power law
  //RooRealVar m(*mass);
  RooRealVar pl_p1("pl_p1","pl_p1",-0.01,-10.,10.);
  RooGenericPdf pl_1("pl1","pl1","pow(@0,@1)",RooArgSet(*mass,pl_p1));
  
  // fit exp and plot
  exp_pdf.fitTo(*hMassData);
  RooPlot *Mframe = mass->frame(Title(Form("Single exp fit to actual data cat %d",cat)));
  hMassData->plotOn(Mframe,DataError(RooDataSet::SumW2));
  exp_pdf.plotOn(Mframe);

  TCanvas *c1 = new TCanvas();
  Mframe->Draw();
  c1->Print(Form("plots/fit_cat%d.pdf",cat),"pdf");

  // gen new dataset, fit and plot
  RooDataSet *genDat = exp_1.generate(*mass,hMassData->numEntries(),Extended());

  RooRealVar bkgYield("nBkg","nBkg",0.01,0.,10000.);
  RooRealVar sigYield("nSig","nSig",0.01,0.,50.);
  RooHistPdf sigMC("sigMC","sigMC",*mass,*sigMCHist);
  RooAddPdf sigAndBkg("sigAndBkg","sigAndBkg",RooArgList(pl_1,sigMC),RooArgList(bkgYield,sigYield));
  sigAndBkg.fitTo(*genDat);
  //pl_p1.fitTo(*genDat);
  
  RooPlot *Gframe = mass->frame(Title(Form("Generated data from single exp cat %d",cat)));
  genDat->plotOn(Gframe,DataError(RooDataSet::SumW2));
  sigAndBkg.plotOn(Gframe);
  c1->Clear();
  Gframe->Draw();
  c1->Print(Form("plots/gen_cat%d.pdf",cat),"pdf");

}
*/
// basic plan
/*
RooExponential getExp(RooDataSet *hMassData, int nFs=0){

  if (nPars<0 || nPars>2) {
    cout << "Wrong number of params" << endl;
    cout << "Returning default: single exponential" << endl;
  }
  // single exp
  if (nFs==0 || nFs==1 || nFs==2){
    RooRealVar exp_p0("exp_p0","exp_p0",-0.02,-1.,0.);
    RooExponential sin_exp("exp0","exp0",*mass,exp_p0);
    if (nFs==0) return sin_exp;
  }
  if (nFs==1 || nFs==2){
    RooRealVar exp_p1("exp_p1","exp_p1",-0.02,-1.,0.); //find starting vals
    RooRealVar exp_p2("exp_p2","exp_p2",-0.02,-1.,0.); //find starting vals
    RooExponential temp_exp_1("temp_exp_1","temp_exp_1",*mass,exp_p2);
    RooAddPdf dbl_exp("dbl_exp","dbl_exp",RooArgList(sin_exp,temp_exp_1),RooArgList(1,exp_p1));
    if (nFs==1) return dbl_exp;
  }
    
  if (nFs==2){
    RooRealVar exp_p3("exp_p3","exp_p3",-0.02,-1.,0.); //find starting vals
    RooRealVar exp_p4("exp_p4","exp_p4",-0.02,-1.,0.); //find starting vals
    RooExponential temp_exp_2("temp_exp_2","temp_exp_2",*mass,exp_p4);
    RooAddPdf trip_exp("dbl_exp","dbl_exp",RooArgList(dbl_exp,exp_1),RooArgList(1,exp_p3));
    if (nFs==2) return trip_exp;
  } 
}

// same for power laws, laurent and poly

void fitToData(RooDataSet*, RooAddPdf){}
void genNewData(){} // prob don't need function for this

int main(){

  for (int categories loop){
    for (int nFits loop){ // do job split here
      fitToData(actualData);
      for (int nFits loop){
        for (int nToys){
          genNewData();
          for (int nMasses (5 GeV steps)){
            fitToData(genData);
            getNormOfSignal();
            fillHisto with norm of signal
          }
        }
      }
    }
  }

}
*/
