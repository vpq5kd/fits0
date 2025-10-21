#include "TRandom2.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TGClient.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include <iostream>
using namespace std;

using TMath::Log;

//parms
const double xmin=1;
const double xmax=20;
const int npoints=12;
const double sigma=0.2;

double f(double x){
  const double a=0.5;
  const double b=1.3;
  const double c=0.5;
  return a+b*Log(x)+c*Log(x)*Log(x);
}

void getX(double *x){
  double step=(xmax-xmin)/npoints;
  for (int i=0; i<npoints; i++){
    x[i]=xmin+i*step;
  }
}

void getY(const double *x, double *y, double *ey){
  static TRandom2 tr(0);
  for (int i=0; i<npoints; i++){
    y[i]=f(x[i])+tr.Gaus(0,sigma);
    ey[i]=sigma;
  }
}


void leastsq(){
  double x[npoints];
  double y[npoints];
  double ey[npoints];
  getX(x);
  getY(x,y,ey);
  auto tg = new TGraphErrors(npoints,x,y,0,ey);
  tg->Draw("alp");
}

int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  // ******************************************************************************

  gStyle->SetOptStat(0); // turn off histogram stats box


  TCanvas *tc = new TCanvas("c1","Sample dataset",dw,dh);

  double lx[npoints];
  double ly[npoints];
  double ley[npoints];

  getX(lx);
  getY(lx,ly,ley);
  auto tgl = new TGraphErrors(npoints,lx,ly,0,ley);
  tgl->SetTitle("Pseudoexperiment;x;y");
  
  // An example of one pseudo experiment
  tgl->Draw("alp");
  tc->Draw();


  
  // *** modify and add your code here ***
  TMatrixD SolveLSQ(const TMatrixD &A, const TMatrixD &y){
  	TMatrixD At = (A);
	At.T();
	TMatrixD ATAi(AT,TMatrixD::kMult,A);
	ATAi.invert();
	TMatrixD Adag(ATAi,TMatrixD::kMult,AT);
	TMatrixD theta(Adag, TMatrixD::kMult,y);
	return theta;
  }

  TVectorD bestfit(double *x, double *y, double *err, int npoints){
  	TMatrixD A(npoints, 3);
	TMatrixD Y(npoints, 1);

	for (int i = 0; i<npoints; i++){
		double lx = Log(x[i]);
		A(i, 0) = 1.0;
		A(i, 1) = lx;
	        A(i, 2) = lx*lx;
		Y(i, 0) = y[i];	
	}
	TMatrixD theta = SolveLSQ(A,Y);

	double a = theta(0,0);
	double b = theta(1,0);
	double c = theta(2,0);

	double y_fit = a +b*lx+c*lx*lx;
	double chi2 = pow((y[i]-y_fit)/err[i],2); //developed with chatGPT
	double chi2_reduced = chi2 / (npoints -3 );
	
	TVectorD results(4);
	results[0] = a;
	results[1] = b;
	results[2] = c;
	results[3] = chi2_reduced;

	return results; 


  }  

  

  TH2F *h1 = new TH2F("h1","Parameter b vs a;a;b",100,0,1,100,0,1);
  TH2F *h2 = new TH2F("h2","Parameter c vs a;a;c",100,0,1,100,0,1);
  TH2F *h3 = new TH2F("h3","Parameter c vs b;b;c",100,0,1,100,0,1);
  TH1F *h4 = new TH1F("h4","reduced chi^2;;frequency",100,0,1);
	
  double x[npoints], y[npoints], err[npoints];
  for (int i = 0; i < nexperiments; i++){
  	getX(x);
	getY(x,y,err);
	TVectorD results = bestfit(x, y, err, npoints);
	h1->Fill(results[0],results[1]);
	h2->Fill(results[0],results[2]);
	h3->Fill(results[1],results[2]);
	h4->Fill(results[3]);

  }
  // perform many least squares fits on different pseudo experiments here
  // fill histograms w/ required data
  
  TCanvas *tc2 = new TCanvas("c2","my study results",200,200,dw,dh);
  tc2->Divide(2,2);
  tc2->cd(1); h1->Draw("colz");
  tc2->cd(2); h2->Draw("colz");
  tc2->cd(3); h3->Draw("colz");
  tc2->cd(4); h4->Draw();
  
  tc2->Draw();

  // **************************************
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}
