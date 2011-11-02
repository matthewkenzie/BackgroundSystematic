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

TFile *inFile = new TFile("CMS-HGG_1665pb.root");

// get workspace and mass data
RooWorkspace *dataWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
RooRealVar *mass = (RooRealVar*)dataWS->var("CMS_hgg_mass");

// -------------- Global fit parameters ----------------------
  // single exp
  RooRealVar exp0_p0("exp0_p0","exp0_p0",-0.03,-0.1,0.);
  RooGenericPdf sin_exp("sin_exp","sin_exp","exp(@1*@0)",RooArgSet(*mass,exp0_p0));
  //RooExponential sin_exp("sin_exp","sin_exp",*mass,exp_p0);
  // double exp
  RooRealVar exp1_p0("exp1_p0","exp1_p0",-0.,-25.,0.);
  RooRealVar exp1_p1("exp1_p1","exp1_p1",0.,0.,0.6); //find starting vals
  RooRealVar exp1_p2("exp1_p2","exp1_p2",-0.03,-0.1,0.); //find starting vals
 // RooExponential temp_dbl_exp("temp_dbl_exp","temp_dbl_exp",*mass,exp_p2);
  RooGenericPdf dbl_exp("dbl_exp","dbl_exp","exp(@1*@0)+@2*exp(@3*@0)",RooArgSet(*mass,exp1_p0,exp1_p1,exp1_p2));
  // triple exp
  RooRealVar exp2_p0("exp2_p0","exp2_p0",-0.08,-0.1,0.);
  RooRealVar exp2_p1("exp2_p1","exp2_p1",0.,-0.2,0.2); //find starting vals
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

// ---- utility function which requests correct function ----

void checkInput(string name){
  if (name!="sin_exp" && name!="dbl_exp" && name!="trip_exp" && name!="sin_pow" && name!="dbl_pow" && name!="trip_pow" && name!="two_lau" && name!="four_lau" && name!="six_lau" && name!="sin_pol" && name!="dbl_pol" && name!="trip_pol") 
    cout << name << " is not a valid function" << endl;
  assert(name=="sin_exp" || name=="dbl_exp" || name=="trip_exp" || name=="sin_pow" || name=="dbl_pow" || name=="trip_pow" || name=="two_lau" || name=="four_lau" || name=="six_lau" || name=="sin_pol" || name=="dbl_pol" || name=="trip_pol");  
}

int main(int argc, char* argv[]){

  string fitFunc(argv[2]);
  string genFunc(argv[3]);
  checkInput(fitFunc);
  checkInput(genFunc);
  const int nToys=atoi(argv[4]);
  const int cat=atoi(argv[5]) ;
  const int nMasses=9;
  const int mRlow=atoi(argv[6]);
  const int mRhigh=atoi(argv[7]);

  // make outputs directories file
  for (int i=0; i<nFits; i++){
    for (int j=0; j<nFits; j++){
      system(Form("mkdir -p BkgSystResults/%s/%s/ToyPlots",fcnNames[i].c_str(),fcnNames[j].c_str()));
      system(Form("mkdir -p BkgSystResults/%s/%s/HistoPlots",fcnNames[i].c_str(),fcnNames[j].c_str()));
    }
  }

  
  TFile *outFile = new TFile(Form("BkgSystResults/%s/%s/HistoFile_cat%d.root",fitFunc.c_str(),genFunc.c_str(),cat),"RECREATE");
  TFile *fitFile = new TFile(Form("data/fitFiles/FitAndToys_%s_cat%d.root",fitFunc.c_str(),cat));
  RooWorkspace *fitWS = (RooWorkspace*)fitFile->Get("cms-hgg-fits");
  
  for (int mMC=mRlow; mMC<=mRhigh; mMC+=5){
    if (mMC==145) continue;
    // ---- declare histos -----
    TH1F *sigNorm = new TH1F(Form("sigNorm_%s_%s_m%d_cat%d",fitFunc.c_str(),genFunc.c_str(),mMC,cat),Form("sigNorm_%s_%s_m%d_cat%d",fitFunc.c_str(),genFunc.c_str(),mMC,cat),100,-100,100);
    sigNorm->GetXaxis()->SetTitle("No. of fitted signal events");
    // -------------------------
    for (int itToy=0; itToy<nToys; itToy++){
      RooDataSet *genDat = (RooDataSet*)fitWS->data(Form("%s_cat%d_toy%d",fitFunc.c_str(),cat,itToy));
      

      RooGenericPdf genFitFcn = getFunction(genFunc);
      RooDataHist *sigMCHist = (RooDataHist*)dataWS->data(Form("roohist_sig_mass_m%d_cat%d",mMC,cat));
      RooHistPdf sigMC("sigMC","sigMC",*mass,*sigMCHist);
      RooRealVar bkgYield("nBkg","nBkg",1000.,0.,5000.);
      RooRealVar sigYield("nSig","nSig",0.,-100,100.);
      RooAddPdf sigAndBkg(Form("sigAndBkg%d_toy%d",mMC,itToy),Form("sigAndBkg%d_toy%d",mMC,itToy),RooArgList(genFitFcn,sigMC),RooArgList(bkgYield,sigYield));
      sigAndBkg.fitTo(*genDat/*,PrintLevel(-1),PrintEvalErrors(-1)*/);
      sigNorm->Fill(sigYield.getVal());

      // ------ plot every 100th Toy ------  
      if (nToys>1 && itToy%(nToys/10)==0) {
        cout << Form("%2.0f %% of toys thown",100*double(itToy)/double(nToys)) << endl;
        TCanvas *canv = new TCanvas();
        RooPlot *gFrame = mass->frame(Title(Form("Gen data from %s fitted with %s outside window at M %d in cat %d toy %d",fitFunc.c_str(),genFunc.c_str(),mMC,cat,itToy)));
        genDat->plotOn(gFrame,DataError(RooDataSet::SumW2));
        sigAndBkg.plotOn(gFrame);
        //sigMC.plotOn(gFrame,LineColor(kRed));
        gFrame->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV/c^{2})");
        gFrame->Draw();
        canv->Print(Form("BkgSystResults/%s/%s/ToyPlots/Gen_%s_fit_%s_m%d_cat%d_toy%d.png",fitFunc.c_str(),genFunc.c_str(),fitFunc.c_str(),genFunc.c_str(),mMC,cat,itToy),"png");
      // -----------------------------------
      }
    }
    TCanvas *canv = new TCanvas();
    sigNorm->Draw();
    canv->Print(Form("BkgSystResults/%s/%s/HistoPlots/sigNorm_%s_%s_m%d_cat%d.png",fitFunc.c_str(),genFunc.c_str(),fitFunc.c_str(),genFunc.c_str(),mMC,cat),"png");
    outFile->cd();
    sigNorm->Write();
  }
  
  inFile->Close();
  fitFile->Close();
  outFile->Close();
 
 return 0;
}
