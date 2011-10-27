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

using namespace std;
using namespace RooFit;

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
void Plot(RooRealVar *mass, RooDataSet* hMassData, RooGenericPdf dataFitFcn, int MCmass, int cat, int toy, string fitFuncName, string genFuncName){
 
  TCanvas *canv = new TCanvas();

  RooPlot *mFrame = mass->frame(Title(Form("Fit of %s to data cat %d",dataFitFcn.GetName(),cat)));
  hMassData->plotOn(mFrame,DataError(RooDataSet::SumW2));
  dataFitFcn.plotOn(mFrame);
  mFrame->Draw();
  canv->Print(Form("plots/test/Fit_%s_cat%d.pdf",dataFitFcn.GetName(),cat),"pdf");
  //canv->Print(Form("plots/%s/%s/Fit_%s_cat%d.pdf",fitFuncName.c_str(),genFuncName.c_str(),dataFitFcn.GetName(),cat),"pdf");
  canv->Clear();
  delete canv;
}

void Plot(RooRealVar *mass, RooDataSet* hMassData, RooExponential dataFitFcn, int MCmass, int cat, int toy, string fitFuncName, string genFuncName){
 
  TCanvas *canv = new TCanvas();

  RooPlot *mFrame = mass->frame(Title(Form("Fit of %s to data cat %d",dataFitFcn.GetName(),cat)));
  hMassData->plotOn(mFrame,DataError(RooDataSet::SumW2));
  dataFitFcn.plotOn(mFrame);
  mFrame->Draw();
  canv->Print(Form("plots/test/Fit_%s_cat%d.pdf",dataFitFcn.GetName(),cat),"pdf");
  //canv->Print(Form("plots/%s/%s/Fit_%s_cat%d.pdf",fitFuncName.c_str(),genFuncName.c_str(),dataFitFcn.GetName(),cat),"pdf");
  canv->Clear();

  /*
  RooPlot *gFrame = mass->frame(Title(Form("Gen data from %s fitted with %s plus signal at M %d in cat %d toy %d",dataFitFcn.GetName(),genFitFcn.GetName(),MCmass,cat,toy)));
  genDat->plotOn(gFrame,DataError(RooDataSet::SumW2));
  sigAndBkg.plotOn(gFrame);
  gFrame->Draw();
  canv->Print(Form("plots/%s/%s/Gen_%s_fit_%s_m%d_cat%d_toy%d.pdf",fitFuncName.c_str(),genFuncName.c_str(),dataFitFcn.GetName(),genFitFcn.GetName(),MCmass,cat,toy),"pdf");
  canv->Clear();

  RooPlot *pFrame = mass->frame(Title(Form("Gen data from %s fitted with %s plus positive only signal at M %d in cat %d toy %d",dataFitFcn.GetName(),genFitFcn.GetName(),MCmass,cat,toy)));
  genDat->plotOn(pFrame,DataError(RooDataSet::SumW2));
  posSigAndBkg.plotOn(pFrame);
  pFrame->Draw();
  canv->Print(Form("plots/%s/%s/pos_Gen_%s_fit_%s_m%d_cat%d_toy%d.pdf",fitFuncName.c_str(),genFuncName.c_str(),dataFitFcn.GetName(),genFitFcn.GetName(),MCmass,cat,toy),"pdf");
  */

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
 
  // ---------------- setup ------------------------
  if (argc!=7) cout << "Not enough arguments" << endl;
  assert(argc==7);
  string fitFunc(argv[1]);
  string genFunc(argv[2]);
  checkInput(fitFunc);
  checkInput(genFunc);

  string fcnNames[12] = {"sin_exp","dbl_exp","trip_exp","sin_pow","dbl_pow","trip_pow","two_lau","four_lau","six_lau","sin_pol","dbl_pol","trip_pol"};
  for (int n1=0; n1<12; n1++){
    for (int n2=0; n2<12; n2++){
      system(("mkdir -p plots/"+fcnNames[n1]+"/"+fcnNames[n2]).c_str());
    }
  }
  system("mkdir -p plots/test");

  ofstream diagFile("diags.txt");
  const int nToys=atoi(argv[3]);
  const int nCats=atoi(argv[4]) ;
  const int nMasses=9;
  const int mRlow=atoi(argv[5]);
  const int mRhigh=atoi(argv[6]);
  cerr << "Running " << nToys << " toys over " << nCats << " catgeories around " << fitFunc << " fit to data and fitting this with " << genFunc << endl;
  cerr << "Signal mass range = [" << mRlow << "-" << mRhigh << "]" << endl;
  // ---------------------------------------------

  // --- declare histograms--------------
  TH1F *sigNormHist[nMasses][nCats];
  TH1F *bkgNormHist[nMasses][nCats];
  TH2F *bkgSigNormCorr[nMasses][nCats];
  TH1F *posSigNormHist[nMasses][nCats];
  TH1F *posBkgNormHist[nMasses][nCats];
  TH2F *posBkgSigNormCorr[nMasses][nCats];

  for (int cat=0; cat<nCats; cat++){
    int mIt=0;
    for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
      if (mMC==145) continue;
      // ------- make histos for bias study -------
      sigNormHist[mIt][cat] = new TH1F(Form("sigNorm_m%d_cat%d",mMC,cat),Form("sigNorm_m%d_cat%d",mMC,cat),200,-50,50);
      bkgNormHist[mIt][cat] = new TH1F(Form("bkgNorm_m%d_cat%d",mMC,cat),Form("bkgNorm_m%d_cat%d",mMC,cat),200,-50,50);
      bkgSigNormCorr[mIt][cat] = new TH2F(Form("bkgSigCorr_m%d_cat%d",mMC,cat),Form("bkgSigCorr_m%d_cat%d",mMC,cat),20,-50,50,20,-50,50);
      bkgSigNormCorr[mIt][cat]->GetXaxis()->SetTitle("sigNorm");
      bkgSigNormCorr[mIt][cat]->GetYaxis()->SetTitle("bkgNorm");
      posSigNormHist[mIt][cat] = new TH1F(Form("posSigNorm_m%d_cat%d",mMC,cat),Form("posSigNorm_m%d_cat%d",mMC,cat),200,-50,50);
      posBkgNormHist[mIt][cat] = new TH1F(Form("posBkgNorm_m%d_cat%d",mMC,cat),Form("posBkgNorm_m%d_cat%d",mMC,cat),200,-50,50);
      posBkgSigNormCorr[mIt][cat] = new TH2F(Form("posBkgSigCorr_m%d_cat%d",mMC,cat),Form("posBkgSigCorr_m%d_cat%d",mMC,cat),20,-50,50,20,-50,50);
      posBkgSigNormCorr[mIt][cat]->GetXaxis()->SetTitle("sigNorm");
      posBkgSigNormCorr[mIt][cat]->GetYaxis()->SetTitle("bkgNorm");
      // ------------------------------------------
      mIt++;
    }
  }

  diagFile << setw(6) << "Mass" << setw(6) << "Cat" << setw(10) << "datNevt" << setw(10) << "datFitInt" << setw(10) << "genNevt" << setw(10) << "genFitInt" << setw(10) << "sigNorm" << setw(10) << "sigYield" << setw(10) << "bkgNmFit" << setw(10) << "bkgNmGen" << setw(10) << "bkgDiff" << setw(10) << "s+bNorm" << setw(10) << "+sigN" << setw(10) << "+bkgN" << endl;
    for (int cat=0; cat<nCats; cat++){
      cerr << "Category " << cat << " of " << nCats << endl;

      int mMC=120;
      int itToy=0;
      int mIt=0;
      // Get data for each cat and fit with function
      RooDataSet *hMassData = (RooDataSet*)dataWS->data(Form("data_mass_cat%d",cat));
      //RooAddPdf dataFitFcn = getFunction(fitFunc); //getExp(0);
      //dataFitFcn.fitTo(*hMassData);
      /*
      sin_exp.fitTo(*hMassData);
      Plot(mass,hMassData,sin_exp,mMC,cat,itToy,fitFunc,genFunc);
      double e0p0 = exp0_p0.getVal();
      dbl_exp.fitTo(*hMassData);
      Plot(mass,hMassData,dbl_exp,mMC,cat,itToy,fitFunc,genFunc);
      double e1p0 = exp1_p0.getVal();
      double e1p1 = exp1_p1.getVal();
      double e1p2 = exp1_p2.getVal();
      trip_exp.fitTo(*hMassData);
      Plot(mass,hMassData,trip_exp,mMC,cat,itToy,fitFunc,genFunc);
      double e2p0 = exp2_p0.getVal();
      double e2p1 = exp2_p1.getVal();
      double e2p2 = exp2_p2.getVal();
      double e2p3 = exp2_p3.getVal();
      double e2p4 = exp2_p4.getVal();
      */

      for (int i=0; i<12; i++){
        RooGenericPdf datFitFunc = getFunction(fcnNames[i]);
        datFitFunc.fitTo(*hMassData,PrintEvalErrors(-1));
        /*
        RooAbsReal *resultNLL = datFitFunc.createNLL(*hMassData);
        cout << resultNLL->getVal() << endl;
        RooAbsReal *resultChi2 = datFitFunc.createChi2(*hMassData);
        cout << resultChi2->getVal() << endl;
        */
        Plot(mass,hMassData,datFitFunc,mMC,cat,itToy,fitFunc,genFunc);
      }

      
      //RooAbsReal *chi2result = sin_exp.createChi2(*mass);
      //cerr << chi2result->getVal() << endl;
     
      /*
      cout << setw(15) << "Fit" << setw(15) << "p0" << setw(15) << "p1" << setw(15) << "p2" << setw(15) << "p3" << setw(15) << "p4" << endl;
      cout << setw(15) << "e0" << setw(15) << e0p0 << endl;
      cout << setw(15) << "e1" << setw(15) << e1p0 << setw(15) << e1p1 << setw(15) << e1p2 << endl;
      cout << setw(15) << "e2" << setw(15) << e2p0 << setw(15) << e2p1 << setw(15) << e2p2 << setw(15) << e2p3 << setw(15) << e2p4 << endl;
      */
      
      for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
        if (mMC==145) continue;
        histPlot(sigNormHist[mIt][cat],fitFunc,genFunc);
        histPlot(bkgNormHist[mIt][cat],fitFunc,genFunc);
        histPlot(bkgSigNormCorr[mIt][cat],fitFunc,genFunc);
        histPlot(posSigNormHist[mIt][cat],fitFunc,genFunc);
        histPlot(posBkgNormHist[mIt][cat],fitFunc,genFunc);
        histPlot(posBkgSigNormCorr[mIt][cat],fitFunc,genFunc);
        outFile->cd();
        sigNormHist[mIt][cat]->Write();
        bkgNormHist[mIt][cat]->Write();
        bkgSigNormCorr[mIt][cat]->Write();
        posSigNormHist[mIt][cat]->Write();
        posBkgNormHist[mIt][cat]->Write();
        posBkgSigNormCorr[mIt][cat]->Write();
        mIt++;
      }
    }
    
    inFile->Close();
    outFile->Close();
    diagFile.close();

    return 0;
  }
