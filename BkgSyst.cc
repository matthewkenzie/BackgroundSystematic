#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
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
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooNLLVar.h"

using namespace std;
using namespace RooFit;

const int nFits=12;
string fcnNames[nFits] = {"sin_exp","dbl_exp","trip_exp","sin_pow","dbl_pow","trip_pow","two_lau","four_lau","six_lau","sin_pol","dbl_pol","trip_pol"};
// in file

TFile *inFile = new TFile("CMS-HGG_1665pb.root");

// get workspace and mass data
RooWorkspace *dataWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
RooRealVar *mass = (RooRealVar*)dataWS->var("CMS_hgg_mass");

// -------------- Global fit parameters ----------------------
  // dummy varaible for 1
  RooRealVar dumY("dumY","dumY",1.0,1.0,1.0);
  // single exp
  RooRealVar exp0_p0("exp0_p0","exp0_p0",-0.01,-10.,0.);
  RooGenericPdf sin_exp("sin_exp","sin_exp","exp(@1*@0)",RooArgSet(*mass,exp0_p0));
  //RooExponential sin_exp("sin_exp","sin_exp",*mass,exp_p0);
  // double exp
  RooRealVar exp1_p0("exp1_p0","exp1_p0",-0.,-1.,0.);
  RooRealVar exp1_p1("exp1_p1","exp1_p1",0.,-0.5,0.5); //find starting vals
  RooRealVar exp1_p2("exp1_p2","exp1_p2",-0.,-1.,0.); //find starting vals
 // RooExponential temp_dbl_exp("temp_dbl_exp","temp_dbl_exp",*mass,exp_p2);
  RooGenericPdf dbl_exp("dbl_exp","dbl_exp","exp(@1*@0)+@2*exp(@3*@0)",RooArgSet(*mass,exp1_p0,exp1_p1,exp1_p2));
  // triple exp
  RooRealVar exp2_p0("exp2_p0","exp2_p0",-0.,-1.,0.);
  RooRealVar exp2_p1("exp2_p1","exp2_p1",0.,-1.,1.); //find starting vals
  RooRealVar exp2_p2("exp2_p2","exp2_p2",-0.,-1.,0.); //find starting vals
  RooRealVar exp2_p3("exp2_p3","exp2_p3",0.,-0.5,0.5); //find starting vals
  RooRealVar exp2_p4("exp2_p4","exp2_p4",-0.,-1.,0.); //find starting vals
  RooGenericPdf trip_exp("trip_exp","trip_exp","exp(@1*@0)+@2*exp(@3*@0)+@4*exp(@5*@0)",RooArgSet(*mass,exp2_p0,exp2_p1,exp2_p2,exp2_p3,exp2_p4));
  // single pow
  RooRealVar pl0_p0("pl0_p0","pl0_p0",-0.,-10.,0.);
  RooGenericPdf sin_pow("sin_pow","sin_pow","pow(@0,@1)",RooArgSet(*mass,pl0_p0));
  // double pow
  RooRealVar pl1_p0("pl1_p0","pl1_p0",-0.,-10.,0.);
  RooRealVar pl1_p1("pl1_p1","pl1_p1",0.,-2.,2.);
  RooRealVar pl1_p2("pl1_p2","pl1_p2",.0,-10.,0.);
  RooGenericPdf dbl_pow("dbl_pow","dbl_pow","pow(@0,@1)+@2*pow(@0,@3)",RooArgSet(*mass,pl1_p0,pl1_p1,pl1_p2));
  // triple pow
  RooRealVar pl2_p0("pl2_p0","pl2_p0",-0.,-10.,0.);
  RooRealVar pl2_p1("pl2_p1","pl2_p1",-0.,-2.,2.);
  RooRealVar pl2_p2("pl2_p2","pl2_p2",0.,-10.,10.);
  RooRealVar pl2_p3("pl2_p3","pl2_p3",0.,-0.1,0.1);
  RooRealVar pl2_p4("pl2_p4","pl2_p4",0.,-10.,10.);
  RooGenericPdf trip_pow("trip_pow","trip_pow","pow(@0,@1)+@2*pow(@0,@3)+@4*pow(@0,5)",RooArgSet(*mass,pl2_p0,pl2_p1,pl2_p2,pl2_p3,pl2_p4));
  // two term laurent
  RooRealVar lau0_p0("lau0_p0","lau0_p0",0.104,-10.,10.);
  RooGenericPdf two_lau("two_lau","two_lau","pow(@0,-4.0)+@1*pow(@0,-5.0)",RooArgSet(*mass,lau0_p0));
  // four term laurent
  RooRealVar lau1_p0("lau1_p0","lau1_p0",0.159,-10.,10.);
  RooRealVar lau1_p1("lau1_p1","lau1_p1",0.263,-10.,10.);
  RooRealVar lau1_p2("lau1_p2","lau1_p2",0.03,-10.,10.);
  RooGenericPdf four_lau("four_lau","four_lau","pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)",RooArgSet(*mass,lau1_p0,lau1_p1,lau1_p2));
  // six term laurent
  RooRealVar lau2_p0("lau2_p0","lau2_p0",0.,-10.,10.);
  RooRealVar lau2_p1("lau2_p1","lau2_p1",0.,-10.,10.);
  RooRealVar lau2_p2("lau2_p2","lau2_p2",0.,-10.,10.);
  RooRealVar lau2_p3("lau2_p3","lau2_p3",0.,-10.,10.);
  RooRealVar lau2_p4("lau2_p4","lau2_p4",0.,-10.,10.);
  RooGenericPdf six_lau("six_lau","six_lau","pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)+@4*pow(@0,-2.0)+@5*pow(@0,-7.0)",RooArgSet(*mass,lau2_p0,lau2_p1,lau2_p2,lau2_p3,lau2_p4));
  // single poly
  RooRealVar pol0_p0("pol0_p0","pol0_p0",-0.001,-0.05,0.);
  RooGenericPdf sin_pol("sin_pol","sin_pol","1.0+(@1*@0)",RooArgSet(*mass,pol0_p0));
  // double poly
  RooRealVar pol1_p0("pol1_p0","pol1_p0",-0.,-0.05,0.);
  RooRealVar pol1_p1("pol1_p1","pol1_p1",-0.,-0.001,0.001);
  RooRealVar pol1_p2("pol1_p2","pol1_p2",0.,-0.01,0.01);
  RooGenericPdf dbl_pol("dbl_pol","dbl_pol","1.0+(@1*@0)+(@2*pow(@0,2.0))+(@3*pow(@0,3.0))",RooArgSet(*mass,pol1_p0,pol1_p1,pol1_p2));
  // triple poly
  RooRealVar pol2_p0("pol2_p0","pol2_p0",-0.001,-10.,10.);
  RooRealVar pol2_p1("pol2_p1","pol2_p1",-0.0,-1.,1.);
  RooRealVar pol2_p2("pol2_p2","pol2_p2",0.,-1.,1.);
  RooRealVar pol2_p3("pol2_p3","pol2_p3",-0.,-0.1,0.1);
  RooRealVar pol2_p4("pol2_p4","pol2_p4",0.,-0.003,0.003);
  RooGenericPdf trip_pol("trip_pol","trip_pol","1.0+(@1*@0)+(@2*pow(@0,2.0))+(@3*pow(@0,3.0))+(@4*pow(@0,4.0))+(@5*pow(@0,5.0))",RooArgSet(*mass,pol2_p0,pol2_p1,pol2_p2,pol2_p3,pol2_p4));
// ----------------------------------------------------------

// ---------- functions which get Pdf for each function ------
RooGenericPdf getFunction(string name){
  
  // exponentials
  if (name=="sin_exp") { 
    return sin_exp;
    //RooGenericPdf returnFcn = getExp(0);
    //return returnFcn;
  }
  else if (name=="dbl_exp") { 
    
    return dbl_exp;
    //RooGenericPdf returnFcn = getExp(1);
    //return returnFcn;
  }
  else if (name=="trip_exp") { 
    return trip_exp;
    //RooGenericPdf returnFcn = getExp(2);
    //return returnFcn;
  }
  // power laws
  else if (name=="sin_pow") { 
    return sin_pow;
    //RooGenericPdf returnFcn = getPow(0);
    //return returnFcn;
  }
  else if (name=="dbl_pow") { 
    return dbl_pow;
    //RooGenericPdf returnFcn = getPow(1);
    //return returnFcn;
  }
  else if (name=="trip_pow") { 
    return trip_pow;
    //RooGenericPdf returnFcn = getPow(2);
    //return returnFcn;
  }
  // laurent series
  else if (name=="two_lau") { 
    return two_lau;
  //  RooGenericPdf returnFcn = getLau(0);
  //  return returnFcn;
  }
  else if (name=="four_lau") { 
    return four_lau;
    //RooGenericPdf returnFcn = getLau(1);
    //return returnFcn;
  }
  else if (name=="six_lau") { 
    return six_lau;
    //RooGenericPdf returnFcn = getLau(2);
    //return returnFcn;
  }
  // polynomials
  else if (name=="sin_pol") { 
    return sin_pol;
    //RooGenericPdf returnFcn = getPol(0);
    //return returnFcn;
  }
  else if (name=="dbl_pol") { 
    return dbl_pol;
    //RooGenericPdf returnFcn = getPol(1);
    //return returnFcn;
  }
  else if (name=="trip_pol") { 
    return trip_pol;
    //RooGenericPdf returnFcn = getPol(2);
    //return returnFcn;
  }
  else {
    cout << "ERROR: this function name is not recognised: " << name << endl;
    cout << "Existing program" << endl;
    exit(1);
  }
}

// ----------- plotting functions ---------------------------
void Plot(RooRealVar *mass, RooDataSet* hMassData, RooGenericPdf *dataFitFcn, int cat){
 
  system("mkdir -p plots/Fits");
  TCanvas *canv = new TCanvas();

  RooPlot *mFrame = mass->frame(Title(Form("Fit of %s to data cat %d",dataFitFcn->GetName(),cat)));
  hMassData->plotOn(mFrame,DataError(RooDataSet::SumW2));
  dataFitFcn->plotOn(mFrame);
  mFrame->Draw();
  canv->Print(Form("plots/Fits/Fit_%s_cat%d.pdf",dataFitFcn->GetName(),cat),"pdf");
  canv->Clear();
  delete canv;
}

void Plot(RooRealVar *mass, RooDataSet* genDat, RooAddPdf sigAndBkg, RooGenericPdf genFitFcn, int MCmass, int cat, int toy, string fitFuncName, string genFuncName,int sigOption){
 
  string temp[2]={"","pos"};
  TCanvas *canv = new TCanvas();

  RooPlot *gFrame = mass->frame(Title(Form("Gen data from %s fitted with %s plus %s signal at M %d in cat %d toy %d",fitFuncName.c_str(),genFuncName.c_str(),temp[sigOption].c_str(),MCmass,cat,toy)));
  genDat->plotOn(gFrame,DataError(RooDataSet::SumW2));
  sigAndBkg.plotOn(gFrame);
  gFrame->Draw();
  canv->Print(Form("plots/%s/%s/%sGen_%s_fit_%s_m%d_cat%d_toy%d.pdf",fitFuncName.c_str(),genFuncName.c_str(),temp[sigOption].c_str(),fitFuncName.c_str(),genFuncName.c_str(),MCmass,cat,toy),"pdf");
  canv->Clear();

  delete canv;
}

void Plot(RooRealVar *mass, RooDataSet* genDat, RooGenericPdf genFitFcn, int MCmass, int cat, int toy, string fitFuncName, string genFuncName){
 
  // set up sidebands
  double lowBand = 0.93*double(MCmass);
  double highBand = 1.07*double(MCmass);
  RooRealVar wholeRange(*mass);
  wholeRange.setRange("wholeRange",100,160);
  RooRealVar lowRange(*mass);
  lowRange.setRange("lowRange",100,lowBand);
  RooRealVar highRange(*mass);
  highRange.setRange("highRange",highBand,160);
  
  TCanvas *canv = new TCanvas();

  RooPlot *gFrame = mass->frame(Title(Form("Gen data from %s fitted with %s outside window at M %d in cat %d toy %d",fitFuncName.c_str(),genFuncName.c_str(),MCmass,cat,toy)));
  genDat->plotOn(gFrame,DataError(RooDataSet::SumW2));
  genFitFcn.plotOn(gFrame,Range("lowband","highband"),NormRange("wholeRange"));
  gFrame->Draw();
  canv->Print(Form("plots/%s/%s/bdtGen_%s_fit_%s_m%d_cat%d_toy%d.pdf",fitFuncName.c_str(),genFuncName.c_str(),fitFuncName.c_str(),genFuncName.c_str(),MCmass,cat,toy),"pdf");
  canv->Clear();

  delete canv;
}

void histPlot(TH1F *hist, string fitFuncName, string genFuncName){

  TCanvas *canv = new TCanvas();
  hist->Draw();
  canv->Print(Form("plots/%s/%s/%s.pdf",fitFuncName.c_str(),genFuncName.c_str(),hist->GetName()),"pdf");
  delete canv;
}
void histPlot(TH2F *hist, string fitFuncName, string genFuncName){

  gStyle->SetPalette(1);
  TCanvas *canv = new TCanvas();
  hist->Draw("colz");
  canv->Print(Form("plots/%s/%s/%s.pdf",fitFuncName.c_str(),genFuncName.c_str(),hist->GetName()),"pdf");
  delete canv;
}

// ----------------------------------------------------------

// ---- utility function which requests correct function ----

void checkInput(string name){
  if (name!="sin_exp" && name!="dbl_exp" && name!="trip_exp" && name!="sin_pow" && name!="dbl_pow" && name!="trip_pow" && name!="two_lau" && name!="four_lau" && name!="six_lau" && name!="sin_pol" && name!="dbl_pol" && name!="trip_pol") 
    cout << name << " is not a valid function" << endl;
  assert(name=="sin_exp" || name=="dbl_exp" || name=="trip_exp" || name=="sin_pow" || name=="dbl_pow" || name=="trip_pow" || name=="two_lau" || name=="four_lau" || name=="six_lau" || name=="sin_pol" || name=="dbl_pol" || name=="trip_pol");  
}

int main(int argc, char* argv[]){

// --------------- Instructions ---------------------------
// Can run with two main options (initial fit on or off)
//   Initial fit on run with ./BkgSyst.exe 1 $1 $2 $3
//      where $1=nToys
//            $2=cat (0-3)
//            $3=fit index (0-11)
//  Secondary fit on run with ./BkgSyst.exe 0 $1 $2 $3 $4 $5
//      where $1=fitFunc name
//            $2=postGenFunc name
//            $3=nToys
//            $4=cat (0-3)
//            $5=mLow
//            $6=mHigh
//            $7=sigOption (0 - for signal as pos or neg
//                          1 - for signal pos only
//                          2 - for bdt method of window +-7%
// -----------------------------------------------------------

  for (int n1=0; n1<12; n1++){
    for (int n2=0; n2<12; n2++){
      system(("mkdir -p plots/"+fcnNames[n1]+"/"+fcnNames[n2]).c_str());
      system(("mkdir -p fitFiles/"+fcnNames[n1]+"/"+fcnNames[n2]).c_str());
    }
  }
  // ---------------- setup ------------------------
  bool doFit;
  if (argv[1]=="true" || atoi(argv[1])==1){
    doFit=true;
    cout << "Fitting method selected" << endl;
    cout << argv[2] << " toys being written" << endl;
  }
  else {
    doFit=false;
    cout << "Generating method selected" << endl;
  }
  // ---------------------------------------------
  

  if (doFit){
    const int nCats=4;
    const int nToys = atoi(argv[2]);
    const int cat = atoi(argv[3]);
    const int itFits = atoi(argv[4]);
    TFile *fitFile = new TFile(Form("fitFiles/FitAndToys_%s_cat%d.root",fcnNames[itFits].c_str(),cat),"RECREATE");
    RooWorkspace *fitWS = new RooWorkspace("cms-hgg-fits");
    fitWS->import(*mass);

    //for (int cat=0; cat<nCats; cat++){
      RooDataSet *hMassData = (RooDataSet*)dataWS->data(Form("data_mass_cat%d",cat));
      cout << "Category " << cat+1 << " of " << nCats << endl;
      // Get data for each cat and fit with function
      fitWS->import(*hMassData,PrintLevel(-1));
     // for (int itFits=0; itFits<nFits; itFits++){
        cout << "fit " << itFits+1 << " of " << nFits << endl;
        RooGenericPdf dataFitFcn = getFunction(fcnNames[itFits]);
        dataFitFcn.fitTo(*hMassData,PrintLevel(-1));
        dataFitFcn.SetName(Form("%s_cat%d",fcnNames[itFits].c_str(),cat));
        fitWS->import(dataFitFcn,PrintLevel(-1));
        Plot(mass,hMassData,&dataFitFcn,cat);
        for (int itToy=0; itToy<nToys; itToy++){
          if (nToys>1 && itToy%(nToys/10)==0) cout << Form("%2.0f %% of toys thown",100*double(itToy)/double(nToys)) << endl;
          RooDataSet *genDat = dataFitFcn.generate(*mass,hMassData->numEntries(),Extended());
          genDat->SetName(Form("%s_cat%d_toy%d",fcnNames[itFits].c_str(),cat,itToy));
          fitWS->import(*genDat,PrintLevel(-1));
        }
     // }
   // }
    fitFile->cd();
    fitWS->Write();
    fitFile->Close();
  }
  else {
  // ---- setup stuff ------
    string fitFunc(argv[2]);
    string genFunc(argv[3]);
    checkInput(fitFunc);
    checkInput(genFunc);
    ofstream diagFile("diags.txt");
    ofstream fitResFile("fitResults.txt");
    const int nToys=atoi(argv[4]);
    const int cat=atoi(argv[5]) ;
    const int nMasses=9;
    const int mRlow=atoi(argv[6]);
    const int mRhigh=atoi(argv[7]);
    const int sigOption=atoi(argv[8]);
    TFile *outFile = new TFile(Form("fitFiles/BkgSystOut_%s_%s_cat%d.root",fitFunc.c_str(),genFunc.c_str(),cat),"RECREATE");
    cout << "Running " << nToys << " toys over in category" << cat+1 << " around " << fitFunc << " fit to data and fitting this with " << genFunc << endl;
    cout << "Signal mass range = [" << mRlow << "-" << mRhigh << "]" << endl;
    diagFile << setw(6) << "Mass" << setw(6) << "Cat" << setw(10) << "datNevt" << setw(10) << "datFitInt" << setw(10) << "genNevt" << setw(10) << "sandBFitInt" << setw(10) << "winBkgDat" << setw(10) << "winSandBGen" << setw(10) << "winBkgGen" << setw(10) << "winSigGen" << setw(10) << "bkgDiff" << endl;
 
    // --- declate out workspace
    RooWorkspace *compWS = new RooWorkspace("cms_hgg_bkgsys");
    compWS->import(*mass,PrintLevel(-1));
    // --- declare histograms--------------
    TH1F *sigNormHist[nMasses];
    TH1F *bkgNormHist[nMasses];
    TH2F *bkgSigNormCorr[nMasses];
    TH1F *posSigNormHist[nMasses];
    TH1F *posBkgNormHist[nMasses];
    TH2F *posBkgSigNormCorr[nMasses];
    TH1F *bdtBkgNormHist[nMasses];

    std::vector<std::pair<string,double> > chi2fitResults;
    std::vector<std::pair<string,double> > NLLfitResults;

    //for (int cat=0; cat<nCats; cat++){
      int mIt=0;
      for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
        if (mMC==145) continue;
        // ------- make histos for bias study -------
        if (sigOption==0){
          sigNormHist[mIt] = new TH1F(Form("sigNorm_m%d_cat%d",mMC,cat),Form("sigNorm_m%d_cat%d",mMC,cat),200,-50,50);
          bkgNormHist[mIt] = new TH1F(Form("bkgNorm_m%d_cat%d",mMC,cat),Form("bkgNorm_m%d_cat%d",mMC,cat),200,-50,50);
          bkgSigNormCorr[mIt] = new TH2F(Form("bkgSigCorr_m%d_cat%d",mMC,cat),Form("bkgSigCorr_m%d_cat%d",mMC,cat),20,-50,50,20,-50,50);
          bkgSigNormCorr[mIt]->GetXaxis()->SetTitle("sigNorm");
          bkgSigNormCorr[mIt]->GetYaxis()->SetTitle("bkgNorm");
        }
        if (sigOption==1){
          posSigNormHist[mIt] = new TH1F(Form("posSigNorm_m%d_cat%d",mMC,cat),Form("posSigNorm_m%d_cat%d",mMC,cat),200,-50,50);
          posBkgNormHist[mIt] = new TH1F(Form("posBkgNorm_m%d_cat%d",mMC,cat),Form("posBkgNorm_m%d_cat%d",mMC,cat),200,-50,50);
          posBkgSigNormCorr[mIt] = new TH2F(Form("posBkgSigCorr_m%d_cat%d",mMC,cat),Form("posBkgSigCorr_m%d_cat%d",mMC,cat),20,-50,50,20,-50,50);
          posBkgSigNormCorr[mIt]->GetXaxis()->SetTitle("sigNorm");
          posBkgSigNormCorr[mIt]->GetYaxis()->SetTitle("bkgNorm");
        }
        if (sigOption==2){
          bdtBkgNormHist[mIt] = new TH1F(Form("bdtBkgNorm_m%d_cat%d",mMC,cat),Form("bdtBkgNorm_m%d_cat%d",mMC,cat),200,-50,50);
        }
        // ------------------------------------------
        mIt++;
      }
   // }

    // ----- get fits -----
    //for (int cat=0; cat<nCats; cat++){
      TFile *fitFile = new TFile(Form("data/fitFiles/FitAndToys_%s_cat%d.root",fitFunc.c_str(),cat));
      RooWorkspace *fitWS = (RooWorkspace*)fitFile->Get("cms-hgg-fits");
      RooDataSet *hMassData = (RooDataSet*)dataWS->data(Form("data_mass_cat%d",cat));
      RooGenericPdf *dataFitFcn = (RooGenericPdf*)fitWS->pdf(Form("%s_cat%d",fitFunc.c_str(),cat));
      compWS->import(*hMassData,PrintLevel(-1));
      compWS->import(*dataFitFcn,PrintLevel(-1));
      // get toy data from this
      for (int itToy=0; itToy<nToys; itToy++){
        if (nToys>1 && itToy%(nToys/10)==0) cout << Form("%2.0f %% of toys thown",100*double(itToy)/double(nToys)) << endl;
        //RooDataSet *genDat = dataFitFcn->generate(*mass,hMassData->numEntries(),Extended());
        RooDataSet *genDat = (RooDataSet*)fitWS->data(Form("%s_cat%d_toy%d",fitFunc.c_str(),cat,itToy));
        compWS->import(*genDat,PrintLevel(-1));
        int mIt=0; 
        for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
          if (mMC==145) continue;
          RooRealVar wholeRange(*mass);
          wholeRange.setRange("wholeRange",100,160);
          RooAbsReal* wholeInt_dataFitFcn = dataFitFcn->createIntegral(*mass,NormSet(*mass),Range("wholeRange"));
          double bkgIntegral = wholeInt_dataFitFcn->getVal()*hMassData->numEntries();
          
          // let signal go pos and neg
          double sigYieldMin;
          if (sigOption==0) sigYieldMin==-50;
          if (sigOption==1) sigYieldMin==0;

          if (sigOption==0 || sigOption==1){
            // get signal MC histo
            RooDataHist *sigMCHist = (RooDataHist*)dataWS->data(Form("roohist_sig_mass_m%d_cat%d",mMC,cat));
            RooHistPdf sigMC("sigMC","sigMC",*mass,*sigMCHist);
            // get background func
            RooGenericPdf genFitFcn = getFunction(genFunc); //getPow(0);
            // --- contruct s+b model and fit allowing signal to go negative
            RooRealVar bkgYield("nBkg","nBkg",2000.,1000.,2500.);
            RooRealVar sigYield("nSig","nSig",0.,sigYieldMin,50.);
            RooAddPdf sigAndBkg(Form("sigAndBkg%d",mMC),Form("sigAndBkg%d",mMC),RooArgList(genFitFcn,sigMC),RooArgList(bkgYield,sigYield));
            sigAndBkg.fitTo(*genDat,PrintLevel(-1),PrintEvalErrors(-1));
            compWS->import(sigAndBkg);
            RooDataHist *dh = genDat->binnedClone();
            RooChi2Var chi2_sigBkg("chi2_sigBkg","chi2_sigBkg",sigAndBkg,*dh);
            RooNLLVar nll_sigBkg("nll_sigBkg","nll_sigBkg",sigAndBkg,*dh);
            string temp1(sigAndBkg.GetName());
            double temp2(chi2_sigBkg.getVal());
            double temp4(nll_sigBkg.getVal());
            pair<string,double> temp3(temp1,temp2);
            pair<string,double> temp5(temp1,temp4);
            chi2fitResults.push_back(temp3);
            NLLfitResults.push_back(temp5);
            
            // --- calc integrals in window and across whole range for checking.
            double lowBand = 0.95*double(mMC);
            double highBand = 1.05*double(mMC);
            RooRealVar intRange(*mass);
            intRange.setRange("sigWindow",lowBand,highBand);
            // ---- background in sig region before gen
            RooAbsReal* int_dataFitFcn = dataFitFcn->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
            double bkgFitInt = int_dataFitFcn->getVal()*hMassData->numEntries();
            // ---- background in sig region after gen
            RooAbsReal* int_genFitFcn = genFitFcn.createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
            double bkgGenInt = (int_genFitFcn->getVal()*bkgYield.getVal());
            // ---- signal in sig region after gen
            RooAbsReal* int_sigMC = sigMC.createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
            double sigInt = (int_sigMC->getVal()*sigYield.getVal());
            // ---- s and b in sig region
            RooAbsReal* windowInt_sigAndBkg = sigAndBkg.createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
            double sandbWinIntegral = windowInt_sigAndBkg->getVal()*genDat->numEntries();
            // ---- s and b across whole range 
            RooAbsReal* wholeInt_sigAndBkg = sigAndBkg.createIntegral(*mass,NormSet(*mass),Range("wholeRange"));
            double sandbIntegral = wholeInt_sigAndBkg->getVal()*genDat->numEntries();

            // bias between before and after
            double bkgDiff = bkgGenInt-bkgFitInt;
            
            if (sigOption==0){
              sigNormHist[mIt]->Fill(sigInt);
              bkgNormHist[mIt]->Fill(bkgDiff);
              bkgSigNormCorr[mIt]->Fill(sigInt,bkgDiff);
            }
            else {
              posSigNormHist[mIt]->Fill(sigInt);
              posBkgNormHist[mIt]->Fill(bkgDiff);
              posBkgSigNormCorr[mIt]->Fill(sigInt,bkgDiff);
            }
            if (itToy==nToys-1) Plot(mass,genDat,sigAndBkg,genFitFcn,mMC,cat,itToy,fitFunc,genFunc,sigOption);
            // --- print to diagnostic file -------
            diagFile << setw(6) << mMC << setw(6) << cat << setw(10) << hMassData->numEntries() << setw(10) << bkgIntegral << setw(10) << genDat->numEntries() << setw(10) << sandbIntegral << setw(10) << bkgFitInt << setw(10) << sandbWinIntegral << setw(10) << bkgGenInt << setw(10) << sigInt << setw(10) << bkgDiff << endl;
            // -----------------------------------
          }
          // use BDT sliding window method
          if (sigOption==2){
            // get background func
            RooGenericPdf genFitFcn = getFunction(genFunc); 
            // set up sidebands
            double lowBand = 0.93*double(mMC);
            double highBand = 1.07*double(mMC);
            RooRealVar intRange(*mass);
            intRange.setRange("sigWindow",lowBand,highBand);
            RooRealVar lowRange(*mass);
            lowRange.setRange("lowRange",100,lowBand);
            RooRealVar highRange(*mass);
            highRange.setRange("highRange",highBand,160);
            // fit outside mass window
            genFitFcn.fitTo(*genDat,Range("lowRange","highRange"),PrintLevel(-1));
            compWS->import(genFitFcn,PrintLevel(-1));

            // ---- whole background after gen
            RooAbsReal* int_wholeGenFitFcn = genFitFcn.createIntegral(*mass,NormSet(*mass),Range("wholeRange"));
            double wholeBkgGen = int_wholeGenFitFcn->getVal()*genDat->numEntries();
            // ---- background in sig region before gen
            RooAbsReal* int_dataFitFcn = dataFitFcn->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
            double bkgFitInt = int_dataFitFcn->getVal()*hMassData->numEntries();
            // ---- background in sig region after gen
            RooAbsReal* int_genFitFcn = genFitFcn.createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
            double bkgGenInt = (int_genFitFcn->getVal()*genDat->numEntries());
            
            // bias between before and after
            double bkgDiff = bkgGenInt-bkgFitInt;

            bdtBkgNormHist[mIt]->Fill(bkgDiff);
            
            if (itToy==nToys-1) Plot(mass,genDat,genFitFcn,mMC,cat,itToy,fitFunc,genFunc);
            // --- print to diagnostic file -------
            diagFile << setw(6) << mMC << setw(6) << cat << setw(10) << hMassData->numEntries() << setw(10) << bkgIntegral << setw(10) << genDat->numEntries() << setw(10) << wholeBkgGen << setw(10) << bkgFitInt << setw(10) << " " << setw(10) << bkgGenInt << setw(10) << " " << setw(10) << bkgDiff << endl;
            // -----------------------------------
          }
          
          mIt++;
        }
      }
          
      mIt=0;
      for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
        if (mMC==145) continue;
          
        outFile->cd();
        if (sigOption==0){
          histPlot(sigNormHist[mIt],fitFunc,genFunc);
          histPlot(bkgNormHist[mIt],fitFunc,genFunc);
          histPlot(bkgSigNormCorr[mIt],fitFunc,genFunc);
          sigNormHist[mIt]->Write();
          bkgNormHist[mIt]->Write();
          bkgSigNormCorr[mIt]->Write();
        }
        if (sigOption==1){
          histPlot(posSigNormHist[mIt],fitFunc,genFunc);
          histPlot(posBkgNormHist[mIt],fitFunc,genFunc);
          histPlot(posBkgSigNormCorr[mIt],fitFunc,genFunc);
          posSigNormHist[mIt]->Write();
          posBkgNormHist[mIt]->Write();
          posBkgSigNormCorr[mIt]->Write();
        }
        if (sigOption==2){
          histPlot(bdtBkgNormHist[mIt],fitFunc,genFunc);
          bdtBkgNormHist[mIt]->Write();
        }
        mIt++;
      }
    fitFile->Close();
    outFile->cd();
    compWS->Write();
   // }
    outFile->Close();
    diagFile.close();
    for (int i=0; i<chi2fitResults.size(); i++){
      fitResFile << chi2fitResults.at(i).first << " " << chi2fitResults.at(i).second << " " << NLLfitResults.at(i).second << endl;
    }
  }
  inFile->Close();

  return 0;
}
