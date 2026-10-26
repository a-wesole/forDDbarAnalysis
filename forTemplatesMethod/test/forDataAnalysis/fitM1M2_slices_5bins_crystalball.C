// how to run:  root 'fitM1M2_slices_5bins_crystalball.C(6,3)' or ./run_for_all_big_production_fits.sh 
#include <iostream>
#include "TRandom.h"
#include "TFile.h"
#include <random>
#include "TNtuple.h"
#include "TTree.h"
#include <TBranch.h>
#include <vector>
#include <TRandom3.h>
#include "TROOT.h"
#include "TH1D.h"
#include "TH1.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"
#include "Math/Functions.h"

double histo_entries = 1.0e+05;
// const float fit_range_low = 1.74, fit_range_high = 1.99, D0_mass = 1.8648;
const float fit_range_low = 1.86484-0.125, fit_range_high = 1.86484 + 0.125, D0_mass = 1.8648;
// const float fit_range_low = 1.75, fit_range_high = 1.86484 + 0.125, D0_mass = 1.8648;
const int number_of_slices = 50;
const int number_of_dphi_bins = 10;

using namespace std;

void GetInclusiveComponents(TFile *f, int numBins, TH2F *&inclusiveSignal, TH2F *&inclusiveSwap, TH2F *&inclusiveBkg, bool projY)
{
    inclusiveSignal = nullptr;
    inclusiveSwap   = nullptr;
    inclusiveBkg    = nullptr;

    for (int i = 1; i < numBins+1; ++i) 
    {
        TH2F *S1S2   = (TH2F *)f->Get(Form("S1S2Mass_bin%d", i));
        TH2F *SW1S2  = (TH2F *)f->Get(Form("SW1S2Mass_bin%d", i));
        TH2F *B1S2   = (TH2F *)f->Get(Form("B1S2Mass_bin%d", i));
        TH2F *S1SW2  = (TH2F *)f->Get(Form("S1SW2Mass_bin%d", i));
        TH2F *SW1SW2 = (TH2F *)f->Get(Form("SW1SW2Mass_bin%d", i));
        TH2F *B1SW2  = (TH2F *)f->Get(Form("B1SW2Mass_bin%d", i));
        TH2F *S1B2   = (TH2F *)f->Get(Form("S1B2Mass_bin%d", i));
        TH2F *SW1B2  = (TH2F *)f->Get(Form("SW1B2Mass_bin%d", i));
        TH2F *B1B2   = (TH2F *)f->Get(Form("B1B2Mass_bin%d", i));

        if (!S1S2) {
            cout << "Warning: Missing components for bin " << i << endl;
            continue;
        }

        TH2F *binSig = nullptr;
        TH2F *binSwap = nullptr;
        TH2F *binBkg = nullptr;

        if (projY) {
            binSig  = (TH2F *)S1S2->Clone("temp_sig");
            binSig->Add(SW1S2);
            binSig->Add(B1S2);

            binSwap = (TH2F *)S1SW2->Clone("temp_swap");
            binSwap->Add(SW1SW2);
            binSwap->Add(B1SW2);

            binBkg  = (TH2F *)S1B2->Clone("temp_bkg");
            binBkg->Add(SW1B2);
            binBkg->Add(B1B2);
        } 
        else {
            binSig  = (TH2F *)S1S2->Clone("temp_sig");
            binSig->Add(S1SW2);
            binSig->Add(S1B2);

            binSwap = (TH2F *)SW1S2->Clone("temp_swap");
            binSwap->Add(SW1SW2); 
            binSwap->Add(SW1B2);

            binBkg  = (TH2F *)B1B2->Clone("temp_bkg");
            binBkg->Add(B1SW2);
            binBkg->Add(B1S2);
        }

        if (i == 1 || inclusiveSignal == nullptr) {
            inclusiveSignal = (TH2F *)binSig->Clone("inclusive_signal_th2");
            inclusiveSwap   = (TH2F *)binSwap->Clone("inclusive_swap_th2");
            inclusiveBkg    = (TH2F *)binBkg->Clone("inclusive_bkg_th2");
            inclusiveSignal->SetDirectory(0);
            inclusiveSwap->SetDirectory(0);
            inclusiveBkg->SetDirectory(0);
            
        } 
        else {
            inclusiveSignal->Add(binSig);
            inclusiveSwap->Add(binSwap);
            inclusiveBkg->Add(binBkg);
        }

        // Clean up temporary bin-level clones
        delete binSig;
        delete binSwap;
        delete binBkg;
    }
}

void GetComponentSlices(TFile *f, int dphi_bin, TH2F *&signalBin, TH2F *&swapBin, TH2F *&bkgBin, bool projY)
{

  cout << "made it to GetCom" << endl;
  TH2F *S1S2Mass = (TH2F *)f->Get(Form("S1S2Mass_bin%d", dphi_bin));
  if (!S1S2Mass)
  {
    cout << "Error: S1S2Mass histogram not found for dphi_bin = " << dphi_bin << endl;
    cout << "Looking for histogram name: " << Form("S1S2Mass_%d", dphi_bin) << endl;
    return;
  }
  TH2F *SW1S2Mass = (TH2F *)f->Get(Form("SW1S2Mass_bin%d", dphi_bin));
  TH2F *B1S2Mass = (TH2F *)f->Get(Form("B1S2Mass_bin%d", dphi_bin));
  TH2F *S1SW2Mass = (TH2F *)f->Get(Form("S1SW2Mass_bin%d", dphi_bin));
  TH2F *SW1SW2Mass = (TH2F *)f->Get(Form("SW1SW2Mass_bin%d", dphi_bin));
  TH2F *B1SW2Mass = (TH2F *)f->Get(Form("B1SW2Mass_bin%d", dphi_bin));
  TH2F *S1B2Mass = (TH2F *)f->Get(Form("S1B2Mass_bin%d", dphi_bin));
  TH2F *SW1B2Mass = (TH2F *)f->Get(Form("SW1B2Mass_bin%d", dphi_bin));
  TH2F *B1B2Mass = (TH2F *)f->Get(Form("B1B2Mass_bin%d", dphi_bin));

  if (projY){
  S1S2Mass->Add(SW1S2Mass);
  S1S2Mass->Add(B1S2Mass);
  signalBin = S1S2Mass;
  S1SW2Mass->Add(SW1SW2Mass);
  S1SW2Mass->Add(B1SW2Mass);
  swapBin = S1SW2Mass;
  S1B2Mass->Add(SW1B2Mass);
  S1B2Mass->Add(B1B2Mass);
  bkgBin = S1B2Mass;
  }

  else if (!projY){
  S1S2Mass->Add(S1SW2Mass);
  S1S2Mass->Add(S1B2Mass);
  signalBin = S1S2Mass;
  SW1SW2Mass->Add(SW1S2Mass);
  SW1SW2Mass->Add(SW1B2Mass);
  swapBin = SW1SW2Mass;
  B1B2Mass->Add(B1SW2Mass);
  B1B2Mass->Add(B1S2Mass);
  bkgBin = B1B2Mass;

  delete signalBin;
  delete swapBin;
  delete bkgBin;

  }

}

void FitMassPlotsTF1(TF1 *fitFunction, TH1D *MassPlot)
{
  for (int i = 0; i < 20; i++)
  {
    MassPlot->Fit(fitFunction, "Lq");
  }
  MassPlot->Fit(fitFunction, "q R");
}

void Inclusive_dphi_fitting(TF1 *f, TH1D *SignalMass, TH1D *SignalAndSwapMass, TH1D *MassPlot, TF1 *swap_fit, TH1D *swap_mass, TCanvas *canvas, TF1 *&signal, TF1 *&swap1, TF1 *&background, TH1D *inclusiveSignal_TH1F, TH1D *inclusiveSwap_TH1F, TH1D *inclusiveBkg_TH1F)
{

  cout << "inclusiveSignal_TH1F entries = " << inclusiveSignal_TH1F->GetEntries() << endl;
  cout << "inclusiveSwap_TH1F entries = " << inclusiveSwap_TH1F->GetEntries() << endl;
  cout << "inclusiveBkg_TH1F entries = " << inclusiveBkg_TH1F->GetEntries() << endl;

  canvas->cd();
  canvas->Clear();
  swap_fit->SetLineColor(kGreen);
  canvas->Divide(2, 2);
  canvas->cd(1);
  f->SetLineColor(2);
  f->SetLineWidth(2);
  // //abby
  // f->SetParNames("normalization", "mean mass", "sigma 1", "sigma 2", "gaus ratio", "signalfrac", "smearing", "sigma swap", "x0", "x1", "x2");
  // f->SetParName(11, "x3");
  // f->SetParName(12, "raio guas 2/3");
  // f->SetParName(13, "swap gaus ratio");
  // f->SetParName(14, "sigma 3");
  // f->SetParName(15, "sigma2 swap");
  // f->SetParName(16, "swap mean");
  // f->SetParName(17, "CB scaling");
  // // f->SetParName(20, "CB left ");
  // f->SetParName(18, "CB r sigma");
  // f->SetParName(19, "CB l sigma");
  // //abby
  f->SetParameter(0, 100);     // Normalization
  f->SetParameter(1, D0_mass); // mean mass
  f->SetParameter(2, 0.01);    // sigma 1
  f->SetParameter(3, 0.005);   // sigma 2
  f->SetParameter(14, 0.03);   // sigma 3 signal
  f->SetParameter(4, 0.5);     // ratio between first and second guassian
  f->SetParameter(12, 0.5);    // ratio gaus 2&3

  f->FixParameter(5, 1.0); // fraction that is signal
  f->FixParameter(6, 0);   // scaling factor----------always 0 in MC
  f->FixParameter(7, 1.0); // sigma of swap-------------- does not really mater here as yield is fix to 0
  // f->SetParameter(8,D0_mass);//mean of swap
  f->FixParameter(8, 0);        // polynomial
  f->FixParameter(9, 0);        // polynomial
  f->FixParameter(10, 0);       // polynomail
  f->FixParameter(11, 0);       // poly
  f->FixParameter(13, 0);       // ratio btw gaussian swap
  f->FixParameter(15, 1);       // sigma2 swap
  f->FixParameter(16, D0_mass); // swap mean
  f->FixParameter(17, 0);

  f->SetParLimits(1, 1.85, 1.88); // sigma 1
  f->SetParLimits(4, 0, 1);  // amplitude/ratio
  f->SetParLimits(12, 0, 1); // amplitude/ratio

  f->SetParameter(1, 1.8648); // mean mass

  // for CB

  SignalMass->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
  cout << "---------------- fitting overall signal " << endl;
  for (int i = 0; i < 29; i++)
  {
    SignalMass->Fit(f, "q", "", fit_range_low, fit_range_high);
  }
  SignalMass->Fit(f, "M", "", fit_range_low, fit_range_high);
  SignalMass->Draw("ep");
  canvas->cd(2);
  swap_mass->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
  swap_mass->Draw();
  swap_fit->Draw("same");
  canvas->cd(3);

  f->FixParameter(1, f->GetParameter(1));   // mean mass
  f->FixParameter(2, f->GetParameter(2));   // sigma 1
  f->FixParameter(3, f->GetParameter(3));   // sigma 2
  f->FixParameter(4, f->GetParameter(4));   // ratio/amplitude
  f->FixParameter(12, f->GetParameter(12)); // ratio guaes 2/3
  f->FixParameter(14, f->GetParameter(14)); // sigma 3

  f->SetParameter(7, swap_fit->GetParameter(2));  // sigma 1
  f->SetParameter(15, swap_fit->GetParameter(3)); // sigma 2
  f->SetParameter(13, swap_fit->GetParameter(4)); // raio
  f->SetParameter(16, swap_fit->GetParameter(1)); // swap mean

  f->ReleaseParameter(5);
  f->ReleaseParameter(7);
  f->ReleaseParameter(15);
  f->ReleaseParameter(13);
  f->ReleaseParameter(16);
  f->SetParLimits(5, 0, 1);
  f->SetParLimits(7, 0.1, 1.3 * swap_fit->GetParameter(2));
  f->SetParLimits(15, 0.0, 1.3 * swap_fit->GetParameter(3));
  f->SetParLimits(13, 0, 1); // ratio

  SignalAndSwapMass->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
  cout << "---------------- fitting overall signal + swap" << endl;
  for (int i = 0; i < 29; i++)
  {
    SignalAndSwapMass->Fit(f, "q", "", fit_range_low, fit_range_high);
  }
  SignalAndSwapMass->Fit(f, "M", "", fit_range_low, fit_range_high);

  TF1 *swap = new TF1("swap", "[0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
  swap->FixParameter(0, f->GetParameter(0));  // norm
  swap->FixParameter(1, f->GetParameter(16)); // mean mass
  swap->FixParameter(2, f->GetParameter(7));  // sigma1
  swap->FixParameter(3, f->GetParameter(15)); // sigma2
  swap->FixParameter(4, f->GetParameter(13)); // ratio
  swap->FixParameter(5, f->GetParameter(5));  // signal fraction
  swap->FixParameter(6, f->GetParameter(6));  // smearing
  swap->SetLineColor(8);

  SignalAndSwapMass->Draw("ep");
  swap->Draw("same");
  canvas->Update();
  canvas->cd(4);

  f->FixParameter(0, f->GetParameter(0));   // fraction of signal
  f->FixParameter(5, f->GetParameter(5));   // fraction of signal
  f->FixParameter(15, f->GetParameter(15)); // fraction of signal
  f->FixParameter(13, f->GetParameter(13)); // fraction of signal
  f->FixParameter(16, f->GetParameter(16)); // swap mean
  f->FixParameter(7, f->GetParameter(7));   // sigma swap

  f->ReleaseParameter(8);  // poly +C
  f->ReleaseParameter(9);  // poly  *x
  f->ReleaseParameter(10); // poly *x^2
  f->ReleaseParameter(11); // 12 *x^3
  f->ReleaseParameter(0); // Normalization
  MassPlot->Fit(f, " L q", "", fit_range_low, fit_range_high);
  MassPlot->Fit(f, "L q", "", fit_range_low, fit_range_high);

  f->ReleaseParameter(1);    // mean mass ------allow data to have different mass peak mean than MC
  f->SetParLimits(1, 1.84, 1.90);
  f->ReleaseParameter(6);    // smearing factor
  f->SetParameter(6, 0);     // smearing factor
  f->SetParLimits(6, -1, 1); // ssmearing factor
  f->ReleaseParameter(17); //ratio of kk,pipi peaks 
  f->SetParLimits(17, 0, f->GetParameter(5));

  cout << "fitting allowing to scale mc to data" << endl;
  for (int i = 0; i < 10; i++)
  {
    MassPlot->Fit(f, "L q", "", fit_range_low, fit_range_high);
    //for m2 no L
  }
  //before the last round of fitting we need to release all of the parameters so that the error is not underestimated.
  f->ReleaseParameter(1);   // mean mass
  double m1 = f->GetParameter(1);
  f->SetParLimits(1, m1 - 0.03, m1 + 0.03);

  f->ReleaseParameter(2);   // sigma 1
  double s1 = f->GetParameter(2);
  f->SetParLimits(2, 0.8 * s1, 1.2 * s1);

  f->ReleaseParameter(3);   // sigma 2
  double s2 = f->GetParameter(3);
  f->SetParLimits(3, 0.8 * s2, 1.2 * s2);

  f->ReleaseParameter(14);  // sigma 3
  double s3 = f->GetParameter(14);
  f->SetParLimits(14, 0.8 * s3, 1.2 * s3);

  f->ReleaseParameter(4);   // ratio/amplitude 1
  f->SetParLimits(4, 0.0, 1.0); // Physical probability limit

  f->ReleaseParameter(12); // ratio gaus 2/3
  f->SetParLimits(12, 0.0, 1.0);


  // --- 2. RELEASE SWAP & FRACTION PARAMETERS SAFELY ---
  f->ReleaseParameter(5);   // fraction of signal
  double init_s = f->GetParameter(5);
  double low_s = std::max(0.0, init_s* 0.9);  //we add these lines to stay close to the original fit vlaue but not go  outside (0,1) as well
  double hi_s = std::min(1.0, init_s* 1.1);
  f->SetParLimits(5, low_s, hi_s);

  f->ReleaseParameter(16);  // swap mean
  double sm = f->GetParameter(16);
  f->SetParLimits(16, sm - 0.02, sm + 0.02);

  f->ReleaseParameter(7);   // sigma swap 
  double sw1 = f->GetParameter(7);
  f->SetParLimits(7, 0.9 * sw1, 1.1 * sw1);

  f->ReleaseParameter(15);  // sigma2 swap
  double sw2 = f->GetParameter(15);
  f->SetParLimits(15, 0.9 * sw2, 1.1 * sw2);

  f->ReleaseParameter(13);  // swap ratio
  double init = f->GetParameter(13);
  double low_l = std::max(0.0, init* 0.98);  //we add these lines to stay close to the original fit vlaue but not go  outside (0,1) as well
  double hi_l = std::min(1.0, init* 1.02);
  f->SetParLimits(13, low_l, hi_l);

  f->ReleaseParameter(17);  // kk,pipi ratio

  MassPlot->Fit(f, "L m N 0", "", fit_range_low, fit_range_high); ///
  // MassPlot->Fit(f, "L m N 0", "", fit_range_low, fit_range_high);
  //for m2 no L 

  signal = new TF1("signal_updated", "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6]))))", fit_range_low, fit_range_high);

  signal->FixParameter(0, f->GetParameter(0));
  signal->FixParameter(1, f->GetParameter(1));
  signal->FixParameter(2, f->GetParameter(2));
  signal->FixParameter(3, f->GetParameter(3));
  signal->FixParameter(4, f->GetParameter(4));
  signal->FixParameter(5, f->GetParameter(5));
  signal->FixParameter(6, f->GetParameter(6));
  signal->FixParameter(12, f->GetParameter(12));
  signal->FixParameter(14, f->GetParameter(14));

  swap1 = new TF1("swap1", "[0]*(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6]))))", fit_range_low, fit_range_high);

  swap1->FixParameter(0, f->GetParameter(0));
  swap1->FixParameter(5, f->GetParameter(5));
  swap1->FixParameter(6, f->GetParameter(6));
  swap1->FixParameter(7, f->GetParameter(7));
  swap1->FixParameter(13, f->GetParameter(13));
  swap1->FixParameter(15, f->GetParameter(15));
  swap1->FixParameter(16, f->GetParameter(16));

  swap1->SetLineColor(8);

  // TF1 *background1 = new TF1("raw background", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", fit_range_low, fit_range_high);
  TF1 *background1 = new TF1("raw background", "cheb3(0)", fit_range_low, fit_range_high);
  background1->FixParameter(0, f->GetParameter(8));
  background1->FixParameter(1, f->GetParameter(9));
  background1->FixParameter(2, f->GetParameter(10));
  background1->FixParameter(3, f->GetParameter(11));

  TF1 *kkpipi = new TF1("kk_pipi", "[0]*[17]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[1]), 1.96*(1+[2])) + 1.26*[0]*[17]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
  kkpipi->SetLineColor(kViolet - 4);
  kkpipi->SetLineWidth(1);
  kkpipi->SetLineStyle(1);
  kkpipi->SetFillColorAlpha(kViolet - 4, 0.3);
  kkpipi->SetFillStyle(1001);
  kkpipi->FixParameter(0, f->GetParameter(0));
  kkpipi->FixParameter(17, f->GetParameter(17));
  kkpipi->FixParameter(1, f->GetParameter(6));
  kkpipi->FixParameter(2, f->GetParameter(18));
  kkpipi->FixParameter(3, f->GetParameter(19));

  // background = new TF1("background", "[4] + [5]*x + [6]*x*x + [7]*x*x*x + [0]*[17]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[1]), 1.97*(1+[2])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
  background = new TF1("background", 
    // "cheb3(4) + [0]*[17]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[1]), 1.97*(1+[2])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[1]), 1.7734*(1+[3])))", 
    "cheb3(8) + [0]*[17]*ROOT::Math::crystalball_function(x, [20], [21], 0.025*(1+[6]), 1.95*(1+[18])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, [22], [23], 0.0146*(1+[6]), 1.7734*(1+[19])))",
    fit_range_low, fit_range_high);
  background->FixParameter(0, f->GetParameter(0));
  background->FixParameter(17, f->GetParameter(17));
  background->FixParameter(24, f->GetParameter(24));
  background->FixParameter(6, f->GetParameter(6));
  background->FixParameter(18, f->GetParameter(18));
  background->FixParameter(19, f->GetParameter(19));
  background->FixParameter(20, f->GetParameter(20));
  background->FixParameter(21, f->GetParameter(21));
  background->FixParameter(22, f->GetParameter(22));
  background->FixParameter(23, f->GetParameter(23));

  background->FixParameter(8, f->GetParameter(8));
  background->FixParameter(9, f->GetParameter(9));
  background->FixParameter(10, f->GetParameter(10));
  background->FixParameter(11, f->GetParameter(11));



  MassPlot->SetMinimum(0.0);
  // MassPlot->SetMaximum(1500);
  MassPlot->Draw("EP");
  swap1->Draw("same");
  signal->SetLineColor(4);
  signal->Draw("same");
  background->SetLineColor(kBlack);
  background->Draw("same");
  kkpipi->Draw("same");

  inclusiveSignal_TH1F->SetMarkerStyle(4); inclusiveSignal_TH1F->SetMarkerColor(kBlue+2);
  inclusiveSwap_TH1F->SetMarkerStyle(26); inclusiveSwap_TH1F->SetMarkerColor(kGreen+2);
  inclusiveBkg_TH1F->SetMarkerStyle(25); inclusiveBkg_TH1F->SetMarkerColor(kGray+2);
  inclusiveSignal_TH1F->Draw("same P");
  inclusiveSwap_TH1F->Draw("same P");
  inclusiveBkg_TH1F->Draw("same P");

  canvas->Update();
  canvas->SaveAs("abby.png");
  canvas->SaveAs("abby.pdf");

}

void CB_fitting(TF1 *f, TH1D *kk_mass)
{
  // TF1 *kkpipi = new TF1("kk_pipi", "[0]*ROOT::Math::crystalball_function(x, [4], [5], 0.025*(1+[1]), 1.95*(1+[2])) + [8]*[0]*(ROOT::Math::crystalball_function(x, [6], [7], 0.0146*(1+[9]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
  cout << "-------------------" << endl;
  cout << "fitting CB" << endl;
  cout << "-------------------" << endl;

  f->SetLineColor(kRed);
  f->SetLineWidth(2);

  // f->SetParNames("norm", "rsigma", "rmean", "lmean", "r alpha", "r n", "l apla", "l n", "l scal", "lsig");
  //              0       1          2         3         4         5       6       7       8        9

  f->ReleaseParameter(0);
  f->ReleaseParameter(2);
  f->ReleaseParameter(3);
  f->SetParLimits(4, -3.5, -1.5);
  f->SetParameter(4, -2.2);
  f->SetParLimits(5, 3, 9);
  f->SetParameter(5, 6);
  f->SetParLimits(6, 0.0, 1.0);
  f->SetParameter(6, 0.34);
  f->SetParLimits(7, 5, 10);
  f->SetParameter(7, 8);
  f->SetParameter(8, 1.26);

  // f->FixParameter(4, 2.2);
  // f->FixParameter(5, 17);
  // f->FixParameter(6, 0.34);
  // f->FixParameter(7, 5);

  for (int i = 0; i < 20; i++)
  {
    kk_mass->Fit(f, "q R N 0", "", fit_range_low, fit_range_high);
  }

  kk_mass->Fit(f, "L M R q N 0", "", fit_range_low, fit_range_high);
  // kk_mass->Fit(f, "L M ", "", fit_range_low, fit_range_high);

  // TCanvas *cab = new TCanvas("cab", "cab", 800, 600);
  // cab->cd();
  // kk_mass->Draw();
  // cab->SaveAs("abby.pdf");
  // return;

}

void Fit_FBin(TF1 *f, TH1D *SignalMass, TH1D *SignalAndSwapMass, TH1D *MassPlot, TF1 *swap_fit, TH1D *swap_mass, TCanvas *canvas, TF1 *&signal, TF1 *&swap1, TF1 *&background)
{

  canvas->cd();
  canvas->Clear();
  swap_fit->SetLineColor(kGreen);
  // canvas->Divide(2, 2);
  canvas->cd(1);
  f->SetLineColor(2);
  f->SetLineWidth(2);
  // f->SetParNames("normalization", "mean mass", "sigma 1", "sigma 2", "gaus ratio", "signalfrac", "smearing", "sigma swap", "const", "pol1", "pol2");
  // f->SetParName(11, "pol3");
  // f->SetParName(12, "raio guas 2/3");
  // f->SetParName(13, "swap gaus ratio");
  // f->SetParName(14, "sigma 3");
  // f->SetParName(15, "sigma2 swap");
  // f->SetParName(16, "swap mean");
  f->FixParameter(11,0);
  f->ReleaseParameter(0); // Normalization
  cout << "---------------- fitting background" << endl;

  MassPlot->Fit(f, "q", "", fit_range_low, fit_range_high);
  MassPlot->Fit(f, "q", "", fit_range_low, fit_range_high);
  MassPlot->Fit(f, "q", "", fit_range_low, fit_range_high);
  MassPlot->Fit(f, "q", "", fit_range_low, fit_range_high);

  f->ReleaseParameter(8);  // poly +C
  f->ReleaseParameter(9);  // poly  *x
  f->ReleaseParameter(10); // poly *x^2
  // f->ReleaseParameter(11); // poly *x^3
  for (int i = 0; i < 10; i++)
  {
    MassPlot->Fit(f, "L q", "", fit_range_low, fit_range_high);
  }
  cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAbby look" << endl;
    MassPlot->Fit(f, "L ", "", fit_range_low, fit_range_high);
  // MassPlot->Fit(f, "L m N 0 q", "", fit_range_low, fit_range_high); ///

  signal = new TF1("signal_updated", "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[7]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[7])*TMath::Gaus(x,[1],[8]*(1.0 +[6]))/(sqrt(2*3.14159)*[8]*(1.0 +[6]))))", fit_range_low, fit_range_high);
  signal->FixParameter(0, f->GetParameter(0));  // normailzation
  signal->FixParameter(1, f->GetParameter(1));  // mean mass
  signal->FixParameter(2, f->GetParameter(2));  // sigma 1
  signal->FixParameter(3, f->GetParameter(3));  // sigma 2
  signal->FixParameter(4, f->GetParameter(4));  // ratio btw gaus 1 & 2
  signal->FixParameter(5, f->GetParameter(5));  // signal fraction
  signal->FixParameter(6, f->GetParameter(6));  // smearing
  signal->FixParameter(7, f->GetParameter(12)); // ratio btw gaus 2 & 3
  signal->FixParameter(8, f->GetParameter(14)); // sigma 3

  /// the below signal is wronng!! i will come back and find out why exactly
  swap1 = new TF1("swap", "[0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
  swap1->FixParameter(0, f->GetParameter(0));  // norm
  swap1->FixParameter(1, f->GetParameter(16)); // mean mass
  swap1->FixParameter(2, f->GetParameter(7));  // sigma1
  swap1->FixParameter(3, f->GetParameter(15)); // sigma2
  swap1->FixParameter(4, f->GetParameter(13)); // ratio
  swap1->FixParameter(5, f->GetParameter(5));  // signal fraction
  swap1->FixParameter(6, f->GetParameter(6));  // smearing
  swap1->SetLineColor(8);

  TF1 *background_only = new TF1("background_only", "cheb3(0)", fit_range_low, fit_range_high);
  background_only->FixParameter(0, f->GetParameter(8));
  background_only->FixParameter(1, f->GetParameter(9));
  background_only->FixParameter(2, f->GetParameter(10));
  background_only->FixParameter(3, f->GetParameter(11));

  TF1 *kkpipi = new TF1("kk_pipi", "[0]*[17]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[1]), 1.97*(1+[2])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
  kkpipi->SetLineColor(kViolet - 4);
  kkpipi->SetLineWidth(1);
  kkpipi->SetLineStyle(1);
  kkpipi->SetFillColorAlpha(kViolet - 4, 0.3);
  kkpipi->SetFillStyle(1001);
  kkpipi->FixParameter(0, f->GetParameter(0));
  kkpipi->FixParameter(17, f->GetParameter(17));
  kkpipi->FixParameter(1, f->GetParameter(6));
  kkpipi->FixParameter(2, f->GetParameter(18));
  kkpipi->FixParameter(3, f->GetParameter(19));

  background = new TF1("background", "cheb3(8) + [0]*[17]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[1]), 1.97*(1+[2])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
  background->FixParameter(0, f->GetParameter(0));
  background->FixParameter(17, f->GetParameter(17));
  background->FixParameter(1, f->GetParameter(6));
  background->FixParameter(2, f->GetParameter(18));
  background->FixParameter(3, f->GetParameter(19));
  background->FixParameter(8, f->GetParameter(8));
  background->FixParameter(9, f->GetParameter(9));
  background->FixParameter(10, f->GetParameter(10));
  background->FixParameter(11, f->GetParameter(11));

  MassPlot->SetMinimum(0.0);
  // MassPlot->SetMaximum(15000);
  MassPlot->Draw("EP");
  swap1->Draw("same");
  signal->SetLineColor(4);
  signal->Draw("same");
  background_only->SetLineColor(kGray + 1);
  background_only->Draw("same");
  background->SetLineColor(kBlack);
  background->Draw("same");
  kkpipi->Draw("same");
  canvas->Update();

  canvas->SaveAs("TF1_outputs/M2_check.pdf");
  std::cin.get();
}

void ProcessSlicesForAxis(int dphi_bin,
                          TFile *inf1,
                          TFile *results,
                          TF1 *F_bin_ref,               // F1_bin or F2_bin (already created & fixed)
                          TF1 *fitted_signal_ref,      // fitted_signal from overall fit
                          TF1 *fitted_swap_ref,        // fitted_swap from overall fit
                          TF1 *fitted_background_ref,
                          double sigSwapratio,  // fitted_background from overall fit
                          TCanvas *canvas,
                          bool projectY,               // false => use ProjectionX (M1); true => ProjectionY (M2)
                          const char *pdfPattern,      // e.g. "M1_fits_bin%d.pdf"
                          const char *prefixF,         // e.g. "F1" or "F2"
                          const char *prefixS,         // e.g. "S1" or "S2"
                          const char *prefixSW,        // e.g. "SW1" or "SW2"
                          const char *prefixB,        // e.g. "B1" or "B2"
                          int number_of_slices_local,  // pass number_of_slices
                          // TGraphErrors *gr,
                          TH2F *signalMC = nullptr, TH2F *swapMC = nullptr, TH2F *bkgMC = nullptr)  // graph
{
  // open multipage PDF
  TString pdfName = Form(pdfPattern, dphi_bin);
  canvas->cd();
  canvas->Print(Form("%s[", pdfName.Data()));
  float maxX = 0, maxY = 0;

  for (int i = 0; i < number_of_slices_local; ++i)
  // for (int i = 25; i < 26; ++i)
  {
    TH2D *sliceA = (TH2D *)inf1->Get(Form("M1M2Mass_bin%d", dphi_bin));
    TH2D *sliceB = (TH2D *)inf1->Get(Form("M1M2Mass_bin%d", dphi_bin));

    TH1D *hX = nullptr;
    TH1D *hY = nullptr;

    if (number_of_slices == 50){
    if (sliceA) hX = sliceA->ProjectionX(Form("%s_M1M2hx_bin%d_slice%d", prefixF, dphi_bin, i), i + 1, i + 1);
    if (sliceB) hY = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), i + 1, i + 1);
    if (i==0) {
      TH1D *hx_25 = sliceA->ProjectionX(Form("%s_M1M2hx_bin%d_slice25", prefixF, dphi_bin), 25, 25);
      maxX = hx_25->GetMaximum();
      TH1D *hy_25 = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice25", prefixF, dphi_bin), 25, 25);
      maxY = hy_25->GetMaximum();
    }
  }
  
    if (number_of_slices == 25){

    if (sliceA) hX = sliceA->ProjectionX(Form("%s_M1M2hx_bin%d_slice%d", prefixF, dphi_bin, i), i + 1, i + 1);
    if (sliceB) hY = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), i + 1, i + 1);
    // hX->RebinX(2); hY->RebinX(2);
    if (i==0) {
      TH1D *hx_25 = sliceA->ProjectionX(Form("%s_M1M2hx_bin%d_slice25", prefixF, dphi_bin), 13, 13);
      maxX = 1.2*hx_25->GetMaximum();
      TH1D *hy_25 = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice25", prefixF, dphi_bin), 13, 13);
      maxY = 1.2*hy_25->GetMaximum();
    }
    }

    // pick the histogram to fit/draw
    TH1D *M_slice = projectY ? hY : hX;
    if (!M_slice) {
      cout << "Error: Could not retrieve histogram for bin " << dphi_bin << ", slice " << i << endl;
      cout << "and projectY = " << projectY << endl;
      continue;
    }

    cout << "M_slice histogram: " << M_slice->GetName() << ", nbins: " << M_slice->GetNbinsX() << endl;
    M_slice->SetTitle(Form("%s mass, bin %d, slice %d", projectY ? "M2" : "M1", dphi_bin, i + 1));

    // clone the F_bin reference for this slice (protected)
    TF1 *f_ref_clone = (TF1 *)F_bin_ref->Clone(Form("f31_protected_%s_bin%d_slice%d", prefixF, dphi_bin, i));

    double ref_par8 = f_ref_clone->GetParameter(8); 
    double ref_par9 = f_ref_clone->GetParameter(9); 
    double ref_par10 = f_ref_clone->GetParameter(10); 
    double ref_par11 = f_ref_clone->GetParameter(11); 
    double ref_par17 = f_ref_clone->GetParameter(17); 

    // fitted_foreground_slice: p0 scale, p1 ratio (same lambda structure as original)
    TF1 *fitted_foreground_slice = new TF1(Form("%s_bin%d_slice%d", prefixF, dphi_bin, i),
      // [f_ref_clone](double *x, double *p) {
      [f_ref_clone, ref_par8, ref_par9, ref_par10, ref_par11, ref_par17](double *x, double *p) {
        f_ref_clone->SetParameter(8, ref_par8*p[0]); //const 
        f_ref_clone->SetParameter(9, ref_par9*p[0]); //(x)
        f_ref_clone->SetParameter(10, ref_par10*p[0]); //(x^2)
        f_ref_clone->SetParameter(11, ref_par11*p[0]); //(x^3)
        f_ref_clone->SetParameter(17, ref_par17*p[0]); //(multiplies the kk, pipi peaks)

        // f_ref_clone->SetParameter(8, p[0]);
        f_ref_clone->SetParameter(5, p[1]); //signal/swap ratio 
        f_ref_clone->SetParameter(0, p[2]);  //yield of s+sw
        // return p[0] * f_ref_clone->Eval(x[0]);
        return f_ref_clone->Eval(x[0]);
      }, fit_range_low, fit_range_high, 3);

    fitted_foreground_slice->SetParameter(0, (M_slice->GetBinContent(1))/(f_ref_clone->Eval(fit_range_low)));

    fitted_foreground_slice->FixParameter(2, f_ref_clone->GetParameter(0) * (M_slice->GetBinContent(1))/(f_ref_clone->Eval(fit_range_low)));

    fitted_foreground_slice->FixParameter(1, sigSwapratio);
    // fitted_foreground_slice->SetParLimits(1, 0, 0.99999);
    

    // perform the same three quiet fits as in your code
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    fitted_foreground_slice->ReleaseParameter(1);
    fitted_foreground_slice->SetParLimits(1, 0, 0.99999);
    fitted_foreground_slice->SetParameter(1, sigSwapratio);
    if (number_of_slices==50) if (i < 14 ||i >34) fitted_foreground_slice->FixParameter(1, sigSwapratio);
    if (number_of_slices==25) if (i < 7 ||i >17) fitted_foreground_slice->FixParameter(1, sigSwapratio);
    // fitted_foreground_slice->SetParameter(2, f_ref_clone->GetParameter(0)); 
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    fitted_foreground_slice->FixParameter(1, fitted_foreground_slice->GetParameter(1));
    fitted_foreground_slice->FixParameter(0, fitted_foreground_slice->GetParameter(0));
    fitted_foreground_slice->ReleaseParameter(2);
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LM");
    fitted_foreground_slice->ReleaseParameter(0);
    // fitted_foreground_slice->ReleaseParameter(1);
    M_slice->Fit(fitted_foreground_slice, "LM");

    float yield = fitted_foreground_slice->GetParameter(2);
    float ratio = fitted_foreground_slice->GetParameter(1);
    float bkg_scaling = fitted_foreground_slice->GetParameter(0);

    cout << "yield = " << yield << ", ratio = " << ratio << ", bkg_scaling = " << bkg_scaling << endl;

    if (ratio > 0.999) {
      cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
      cout << "ratio is too large, ratio = " << ratio << endl;
      cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    }

    //   int n = gr->GetN();
    //   gr->SetPoint(n, i+1, bkg_scaling*ref_par9);
    // // gr->SetPointError(n, 0, fabs(bkg_scaling * ref_par9) * sqrt(pow(fitted_foreground_slice->GetParError(0) / fitted_foreground_slice->GetParameter(0), 2) + pow(f_ref_clone->GetParError(9) / f_ref_clone->GetParameter(9), 2)));
    // cout << "000000000000000000000000000000" << endl;
    //   cout << "SetPoint = " << i+1 << ", " << bkg_scaling*ref_par9 << endl;
    //   cout << "bkg_scaling = " << bkg_scaling << " and ref_par9 = " << ref_par9 << endl;
    // cout << "000000000000000000000000000000" << endl;


    // create fitted_signal_slice (clone via lambda, same as original)
    TF1 *fitted_signal_slice = new TF1(Form("fitted_signal_slice_%s_%d", prefixS, i),
      [fitted_signal_ref](double *x, double *p) {
        fitted_signal_ref->SetParameter(5, p[0]);
        fitted_signal_ref->SetParameter(0, p[1]);
        // return p[0] * fitted_signal_ref->Eval(x[0]);
        return fitted_signal_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 2);
    fitted_signal_slice->FixParameter(0, ratio);
    fitted_signal_slice->FixParameter(1, yield);


    // create fitted_swap_slice
    TF1 *fitted_swap_slice = new TF1(Form("fitted_swap_slice_%s_%d", prefixSW, i),
      [fitted_swap_ref](double *x, double *p) {
        fitted_swap_ref->SetParameter(5, p[0]);
        fitted_swap_ref->SetParameter(0, p[1]);
        // return p[0] * fitted_swap_ref->Eval(x[0]);
        return  fitted_swap_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 2);
    fitted_swap_slice->FixParameter(0, ratio);
    fitted_swap_slice->FixParameter(1, yield);

    // create fitted_background_slice
    TF1 *fitted_background_slice = new TF1(Form("fitted_background_slice_%s_%d", prefixB, i),
      [fitted_background_ref, yield](double *x, double *p) {
        fitted_background_ref->SetParameter(0, yield);
        return p[0]* fitted_background_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 1);
    fitted_background_slice->FixParameter(0, bkg_scaling);

    cout << "PPPPP -- bin 1 content: F= " << fitted_foreground_slice->Eval(1.75) << ", S=" << fitted_signal_slice->Eval(1.75)
         << ", SW=" << fitted_swap_slice->Eval(1.75) << ", B=" << fitted_background_slice->Eval(1.75) << endl;

    gStyle->SetOptStat(0);
    // gStyle->SetOptFit(0);
    // draw exactly as before
    M_slice->SetLineColor(kBlack); M_slice->SetMarkerColor(kBlack);
    M_slice->SetMarkerStyle(20); M_slice->SetMinimum(0);

    if (projectY) {
      M_slice->SetMaximum(3000);
    } else {
      M_slice->SetMaximum(maxX * 1.2);
    }
    M_slice->Draw("E");
    fitted_foreground_slice->SetLineColor(kRed); fitted_foreground_slice->Draw("same");
    fitted_signal_slice->SetLineColor(kBlue); fitted_signal_slice->Draw("same");
    fitted_swap_slice->SetLineColor(kGreen); fitted_swap_slice->Draw("same");
    fitted_background_slice->SetLineColor(kBlack); fitted_background_slice->Draw("same");

    if (signalMC && swapMC && bkgMC) {
    signalMC->SetMarkerStyle(4); signalMC->SetMarkerColor(kBlue+2);
    if (!projectY)signalMC->ProjectionX("a",i+1, i+1)->Draw("EP same");
    if (projectY)signalMC->ProjectionY("a",i+1, i+1)->Draw("EP same");
    swapMC->SetMarkerStyle(22); swapMC->SetMarkerColor(kGreen+2);
    if (!projectY)swapMC->ProjectionX("b",i+1, i+1)->Draw("EP same");
    if (projectY)swapMC->ProjectionY("b",i+1, i+1)->Draw("EP same");
    bkgMC->SetMarkerStyle(20);  bkgMC->SetMarkerColor(kGray+2);
    cout << endl;
    cout << "ABBY _____________________________ projectY = " << projectY << endl;
    cout << endl;
    if (!projectY) bkgMC->ProjectionX("c",i+1, i+1)->Draw("EP same");
    if (projectY) bkgMC->ProjectionY("c",i+1, i+1)->Draw("EP same");
    }

    
    TLegend *leg = new TLegend(0.58, 0.55, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0); // Transparent background so it looks clean
    leg->SetTextSize(0.032);

    leg->AddEntry(M_slice, "Foreground", "lep");
    leg->AddEntry(fitted_foreground_slice, "Foreground Fit", "l");
    leg->AddEntry(fitted_signal_slice, "Signal Fit", "l");
    leg->AddEntry(fitted_swap_slice, "Swap Fit", "l");
    leg->AddEntry(fitted_background_slice, "Background Fit", "l");

    if (signalMC && swapMC && bkgMC) {
        leg->AddEntry(signalMC, "Signal MC", "p");
        leg->AddEntry(swapMC, "Swap MC", "p");
        leg->AddEntry(bkgMC, "Bkg MC", "p");
    }

    leg->Draw();

    // print page to PDF
    canvas->Update();
    canvas->Print(pdfName);

    // write TF1s and histograms to results file (same names pattern as original)
    results->cd();
    fitted_foreground_slice->Write(Form("%s_bin%d_slice%d", prefixF, dphi_bin, i));
    fitted_signal_slice->Write(Form("%s_bin%d_slice%d", prefixS, dphi_bin, i));
    fitted_swap_slice->Write(Form("%s_bin%d_slice%d", prefixSW, dphi_bin, i));
    fitted_background_slice->Write(Form("%s_bin%d_slice%d", prefixB, dphi_bin, i));

    fitted_foreground_slice->SetNpx(number_of_slices);
    fitted_signal_slice->SetNpx(number_of_slices);
    fitted_swap_slice->SetNpx(number_of_slices);
    fitted_background_slice->SetNpx(number_of_slices);
    TH1 *hF = (TH1 *)fitted_foreground_slice->GetHistogram()->Clone(Form("%shisto_bin%d_slice%d", prefixF, dphi_bin, i));
    TH1 *hS = (TH1 *)fitted_signal_slice->GetHistogram()->Clone(Form("%shisto_bin%d_slice%d", prefixS, dphi_bin, i));
    TH1 *hSW = (TH1 *)fitted_swap_slice->GetHistogram()->Clone(Form("%shisto_bin%d_slice%d", prefixSW, dphi_bin, i));
    TH1 *hB = (TH1 *)fitted_background_slice->GetHistogram()->Clone(Form("%shisto_bin%d_slice%d", prefixB, dphi_bin, i));

    hF->Write(); hS->Write(); hSW->Write(); hB->Write();

    // optional cleanup to avoid memory growth (safe): delete clones created per slice
    delete fitted_foreground_slice;
    delete fitted_signal_slice;
    delete fitted_swap_slice;
    delete fitted_background_slice;
    delete f_ref_clone;
    delete hF;
    delete hS;
    delete hSW;
    delete hB;
    if (hX) { delete hX; hX = nullptr; }
    if (hY) { delete hY; hY = nullptr; }
  } // end slice loop

  // close PDF
  canvas->cd();
  canvas->Print(Form("%s]", pdfName.Data()));
}
double getSigSwapRatio(int dphi_bin,
                          TFile *inf1,
                          TF1 *F_bin_ref, TF1 *fitted_signal_ref, TF1 *fitted_swap_ref, TF1 *fitted_background_ref,       // F1_bin or F2_bin (already created & fixed)
                          bool projectY,               // false => use ProjectionX (M1); true => ProjectionY (M2)
                          const char *prefixF,         // e.g. "F1" or "F2"
                          const char *prefixS,         // e.g. "S1" or "S2"
                          const char *prefixSW,        // e.g. "SW1" or "SW2"
                          const char *prefixB,        // e.g. "B1" or "B2"
                          int number_of_slices_local, 
                          double pT_leading, double pT_subleading)  // pass number_of_slices
{
  // open multipage PDF
  float maxX = 0, maxY = 0;

    TH2D *sliceA = (TH2D *)inf1->Get(Form("M1M2Mass_bin%d", dphi_bin));
    TH2D *sliceB = (TH2D *)inf1->Get(Form("M1M2Mass_bin%d", dphi_bin));

    TH1D *hX = nullptr;
    TH1D *hY = nullptr;

    int i=0;
    if (number_of_slices == 50){
    if (sliceA) hX = sliceA->ProjectionX(Form("%s_M1M2hx_bin%d_slice%d", prefixF, dphi_bin, i), i + 1, i + 1);
    if (sliceB) hY = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), i + 1, i + 1);
    if (i==0) {
      TH1D *hx_25 = sliceA->ProjectionX(Form("%s_M1M2hx_bin%d_slice25", prefixF, dphi_bin), 25, 25);
      maxX = hx_25->GetMaximum();
      TH1D *hy_25 = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice25", prefixF, dphi_bin), 25, 25);
      maxY = hy_25->GetMaximum();
    }
  }
  
    if (number_of_slices==25){
    if (projectY) {
    hX = sliceA->ProjectionY(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), 1, 7);
    hY = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), 18, 25);
    hY->Add(hX);
    }
    if (!projectY) {
    hX = sliceA->ProjectionX(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), 1, 7);
    hY = sliceB->ProjectionX(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), 18, 25);
    hY->Add(hX);
    }
  }

    else {
    if (projectY) {
    hX = sliceA->ProjectionY(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), 1, 14);
    hY = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), 36, 50);
    hY->Add(hX);
    }
    if (!projectY) {
    hX = sliceA->ProjectionX(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), 1, 14);
    hY = sliceB->ProjectionX(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i), 36, 50);
    hY->Add(hX);
    }
  }
    // if (number_of_slices==25) {hX->RebinX(2); hY->RebinX(2);}
    

    // pick the histogram to fit/draw
    TH1D *M_slice = projectY ? hY : hY; 
    if (!M_slice) {
      cout << "Error: Could not retrieve histogram for bin " << dphi_bin << ", slice " << i << endl;
      cout << "and projectY = " << projectY << endl;
      return -999.9;
    }

    // clone the F_bin reference for this slice (protected)
    TF1 *f_ref_clone = (TF1 *)F_bin_ref->Clone(Form("f31_protected_%s_bin%d_slice%d", prefixF, dphi_bin, i));
    cout << "--------------------------" << endl;
    cout << "the original ratio = " << f_ref_clone->GetParameter(5) << " +/- " << f_ref_clone->GetParError(5) << endl;

    double ref_par8 = f_ref_clone->GetParameter(8);
    double ref_par9 = f_ref_clone->GetParameter(9);
    double ref_par10 = f_ref_clone->GetParameter(10);
    double ref_par11 = f_ref_clone->GetParameter(11);
    double ref_par17 = f_ref_clone->GetParameter(17);



    // fitted_foreground_slice: p0 scale, p1 ratio (same lambda structure as original)
    TF1 *fitted_foreground_slice = new TF1(Form("%s_bin%d_slice%d", prefixF, dphi_bin, i),
            [f_ref_clone, ref_par8, ref_par9, ref_par10, ref_par11, ref_par17](double *x, double *p) {

        f_ref_clone->SetParameter(8, ref_par8*p[0]); //const 
        f_ref_clone->SetParameter(9, ref_par9*p[0]); //(x)
        f_ref_clone->SetParameter(10, ref_par10*p[0]); //(x^2)
        f_ref_clone->SetParameter(11, ref_par11*p[0]); //(x^3)
        f_ref_clone->SetParameter(17, ref_par17*p[0]); //(multiplies the kk, pipi peaks)

      // [f_ref_clone](double *x, double *p) {
        f_ref_clone->SetParameter(5, p[1]);
        f_ref_clone->SetParameter(0, p[2]);
        // return p[0] * f_ref_clone->Eval(x[0]);
        return f_ref_clone->Eval(x[0]);
      }, fit_range_low, fit_range_high, 3);

    fitted_foreground_slice->SetParameter(0, (M_slice->GetBinContent(1))/(f_ref_clone->Eval(fit_range_low))); //bkg scaling
    fitted_foreground_slice->FixParameter(2, f_ref_clone->GetParameter(0) * (M_slice->GetBinContent(1))/(f_ref_clone->Eval(fit_range_low))); //signal+swap
    fitted_foreground_slice->FixParameter(1, f_ref_clone->GetParameter(5)); //sig swap ratio


    // perform the same three quiet fits as in your code
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    // fitted_foreground_slice->FixParameter(0, fitted_foreground_slice->GetParameter(0));
    if (i >=15 && i <= 35) {
    fitted_foreground_slice->ReleaseParameter(1);
    fitted_foreground_slice->SetParLimits(1, 0, 0.99999);
    fitted_foreground_slice->SetParameter(1, f_ref_clone->GetParameter(5));
    // fitted_foreground_slice->SetParameter(2, f_ref_clone->GetParameter(0)); 
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    fitted_foreground_slice->FixParameter(1, fitted_foreground_slice->GetParameter(1));
    }
    else {
      if (number_of_slices==50) fitted_foreground_slice->FixParameter(1, f_ref_clone->GetParameter(5));
    if (number_of_slices==25) fitted_foreground_slice->FixParameter(1, f_ref_clone->GetParameter(5));
  }
    fitted_foreground_slice->FixParameter(0, fitted_foreground_slice->GetParameter(0));
    fitted_foreground_slice->ReleaseParameter(2);
    M_slice->Fit(fitted_foreground_slice, "LQ0");
    M_slice->Fit(fitted_foreground_slice, "L");
    fitted_foreground_slice->ReleaseParameter(0);
    fitted_foreground_slice->ReleaseParameter(1);
    M_slice->Fit(fitted_foreground_slice, "L");

    float yield = fitted_foreground_slice->GetParameter(2);
    float ratio = fitted_foreground_slice->GetParameter(1);
    float bkg_scaling = fitted_foreground_slice->GetParameter(0);

        TF1 *fitted_signal_slice = new TF1(Form("fitted_signal_slice_%s_%d", prefixS, i),
      [fitted_signal_ref](double *x, double *p) {
        fitted_signal_ref->SetParameter(5, p[0]);
        fitted_signal_ref->SetParameter(0, p[1]);
        // return p[0] * fitted_signal_ref->Eval(x[0]);
        return fitted_signal_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 2);
    fitted_signal_slice->FixParameter(0, ratio);
    fitted_signal_slice->FixParameter(1, yield);

    // create fitted_swap_slice
    TF1 *fitted_swap_slice = new TF1(Form("fitted_swap_slice_%s_%d", prefixSW, i),
      [fitted_swap_ref](double *x, double *p) {
        fitted_swap_ref->SetParameter(5, p[0]);
        fitted_swap_ref->SetParameter(0, p[1]);
        // return p[0] * fitted_swap_ref->Eval(x[0]);
        return fitted_swap_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 2);
    fitted_swap_slice->FixParameter(0, ratio);
    fitted_swap_slice->FixParameter(1,  yield);

    // create fitted_background_slice
    TF1 *fitted_background_slice = new TF1(Form("fitted_background_slice_%s_%d", prefixB, i),
      [fitted_background_ref, yield](double *x, double *p) {
        fitted_background_ref->SetParameter(0, yield);
        // fitted_background_ref->SetParameter(8, p[0]);
        return p[0] * fitted_background_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 1);
    fitted_background_slice->FixParameter(0, bkg_scaling);

  TCanvas *abi = new TCanvas("abi", "abi", 1000, 1000);
  abi->cd();
    M_slice->SetMinimum(0);
    if (projectY) M_slice->SetTitle(Form("M2 outer slices merged, dphi bin %d", dphi_bin));
    if (!projectY) M_slice->SetTitle(Form("M1 outer slices merged, dphi bin %d", dphi_bin));
    M_slice->Draw("EP");
    fitted_foreground_slice->SetLineColor(kRed);fitted_foreground_slice->Draw("same");
    fitted_signal_slice->SetLineColor(kBlue); fitted_signal_slice->Draw("same");
    fitted_swap_slice->SetLineColor(kGreen); fitted_swap_slice->Draw("same");
    fitted_background_slice->SetLineColor(kBlack); fitted_background_slice->Draw("same");

  //  if (!projectY) abi->SaveAs(Form("big_production_fits_%ibins_%islices/pT_%d_%d/m1_bin%d.pdf", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading, dphi_bin));
  // if (projectY) abi->SaveAs(Form("big_production_fits_%ibins_%islices/pT_%d_%d/m2_bin%d.pdf", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading, dphi_bin)); 

   if (!projectY) abi->SaveAs(Form("for_Off_MC/big_production_fits_%ibins_%islices/pT_%d_%d/m1_bin%d.pdf", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading, dphi_bin));
  if (projectY) abi->SaveAs(Form("for_Off_MC/big_production_fits_%ibins_%islices/pT_%d_%d/m2_bin%d.pdf", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading, dphi_bin)); 

    return fitted_foreground_slice->GetParameter(1);

    if (hX) { delete hX; hX = nullptr; }
    if (hY) { delete hY; hY = nullptr; }

}

void fitM1M2_slices_5bins_crystalball(double pT_leading, double pT_subleading)
{

  TCanvas *M1_canvas = new TCanvas("M1_canvas", "M1_canvas", 1000, 1000);
  TCanvas *M2_canvas = new TCanvas("M2_canvas", "M2_canvas", 1000, 1000);

  // TString outfile = TString("fited_M1M2_with5bins_updatedRidge.root");
  // TString outfile = TString("fited_M1M2_with5bins_pT6_3_30slices_updated14May_fixing_signalswapratio.root");
  // TString outfile = TString("fited_M1M2_with5bins_pT6_3_30slices_updated18May_added_kkpipi_peaks.root");
  // TString outfile = TString("fited_M1M2_with5bins_pT3_3_30slices_forEvtMixing_incl_kkpipi_test_updated_for_simple_method_good.root");
  // Use %g to print decimals nicely (e.g., 6.5) or %.1f to force a specific precision

  TString outfile = TString(Form("for_Off_MC/big_production_fitted_outputs/pT_%d_%d_with_%i_bins_and_%i_slicesOct13.root", (int)pT_leading, (int)pT_subleading, number_of_dphi_bins, number_of_slices));
  // TString outfile = TString(Form("big_production_fitted_outputs/pT_%d_%d_with_%i_bins_and_%i_slices_v2.root", (int)pT_leading, (int)pT_subleading, number_of_dphi_bins, number_of_slices));
  // TString outfile = TString(Form("big_production_EM_Off_MC_fitted_outputs/pT_%d_%d_with_%i_bins_and_%i_slices.root", (int)pT_leading, (int)pT_subleading, number_of_dphi_bins, number_of_slices));
  // TString outfile = TString("test.root");
  cout << "-------------------------------------------------------Output file: " << outfile << endl;
  cout << "---------------------------------------------------------------" << endl;
  // TString outfile = TString("oo.root");
  TFile *results = new TFile(outfile, "recreate");


  double F1F2Entries, F1SB2Entries, SB1F2Entries, SB1SB2Entries;
  // TString file_containing_M1M2_slices = "/scratch/negishi/awesole/DDbar_Jan2026/ROOT/M1M2.root";
  // TString file_containing_M1M2_slices = "/scratch/negishi/awesole/DDbar_analysis_17Oct/pT_6_pT_3/ROOT/wBins/M1M2.root";
  // TString file_containing_M1M2_slices = "/scratch/negishi/awesole/DDbar_Jan2026/pt_10_5/ROOT/M1M2_hadd.root";
  // TString file_containing_M1M2_slices = "/scratch/negishi/awesole/DDbar_Jan2026/ROOT/M1M2_30slices_everything_we_need.root";
  // TString file_containing_M1M2_slices = "/scratch/negishi/awesole/DDbar_Jan2026/ROOT/M1M2_includingSharedDtrk.root";
  // TString file_containing_M1M2_slices = "/scratch/negishi/awesole/DDbar_July_2026_pT_3_3/ROOT/M1M2_skippingSharedTracks.root"; 
  // TString file_containing_M1M2_slices = "/scratch/negishi/awesole/DDbar_July_2026_pT_6_3/ROOT/M1M2.root"; 

  // TString file_containing_M1M2_slices = TString(Form("/scratch/negishi/awesole/DDbar_Sept_2026_xT_pT_%d_%d/ROOT/M1M2.root", (int)pT_leading, (int)pT_subleading));
  // TString file_containing_M1M2_slices = TString(Form("/scratch/negishi/awesole/DDbar_July_2026_pT_%d_%d/ROOT/M1M2.root", (int)pT_leading, (int)pT_subleading));
  // TString file_containing_M1M2_slices = TString(Form("/scratch/negishi/awesole/DDbar_eventmixing_v3_matchesMilan/pT_%d_%d/ROOT/M1M2.root", (int)pT_leading, (int)pT_subleading));
  // TString file_containing_M1M2_slices = TString(Form("/scratch/negishi/awesole/DDbar_eventmixing_official_MC/pT_%d_%d/ROOT/M1M2_size20_updated.root", (int)pT_leading, (int)pT_subleading));
  // TString file_containing_M1M2_slices = TString(Form("/scratch/negishi/awesole/DDbar_eventmixing_official_MC/pT_%d_%d/ROOT/M1M2_Sept28v2.root", (int)pT_leading, (int)pT_subleading));
  TString file_containing_M1M2_slices = TString(Form("/scratch/negishi/awesole/DDbar_eventmixing_official_MC/pT_%d_%d/ROOT/M1M2_Oct13.root", (int)pT_leading, (int)pT_subleading));


  // TString file_containing_M1M2_slices = "/scratch/negishi/awesole/DDbar_June_2026_pT_3_3/ROOT/M1M2.root";
  TFile *inf1 = TFile::Open(file_containing_M1M2_slices);

  // TString file_containing_M1_signal_swap = "/scratch/negishi/awesole/DDbar_analysis/signal_swap_histos_150bins/ROOT/massPlots_pt_" + TString::Format("%d", (int)pT_leading) + "_1000_cent_0_90_y_0_1.root";
  TString file_containing_M1_signal_swap = "/scratch/negishi/awesole/DDbar_analysis/signal_swap_histos_150bins/ROOT/Oct13_massPlots_pt_" + TString::Format("%d", (int)pT_leading) + "_1000_cent_0_90_y_0_1.root";
  TFile *inf2 = TFile::Open(file_containing_M1_signal_swap);

  // TString file_containing_M2_signal_swap = "/scratch/negishi/awesole/DDbar_analysis/signal_swap_histos_150bins/ROOT/massPlots_pt_" + TString::Format("%d", (int)pT_subleading) + "_1000_cent_0_90_y_0_1.root";
  TString file_containing_M2_signal_swap = "/scratch/negishi/awesole/DDbar_analysis/signal_swap_histos_150bins/ROOT/Oct13_massPlots_pt_" + TString::Format("%d", (int)pT_subleading) + "_1000_cent_0_90_y_0_1.root";
  TFile *inf3 = TFile::Open(file_containing_M2_signal_swap);

  TString file_containing_kk = "/home/awesole/forDDbarAnalysis/forTemplatesMethod/dataFiles/d0ana_mc_D0toKK_total.root";
  TFile *inf4 = TFile::Open(file_containing_kk);
  TTree *kk_tree = (TTree *)inf4->Get("d0ana_mc_genmatch/VertexCompositeNtuple");
  TH1D *kk_mass = new TH1D("kk_mass", "kk_mass", 100, fit_range_low, fit_range_high); // Adjust bins/range as needed

  Float_t mass[50000]; // assuming max candSize ≤ 100
  Int_t candSize;
  kk_tree->SetBranchAddress("mass", mass);
  kk_tree->SetBranchAddress("candSize", &candSize);

  Long64_t nentries = kk_tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    kk_tree->GetEntry(i);
    for (Int_t j = 0; j < candSize; ++j)
    {
      kk_mass->Fill(mass[j]);
    }
  }

  TString file_containing_pipi = "/home/awesole/forDDbarAnalysis/forTemplatesMethod/dataFiles/d0ana_mc_D0toPiPi_total.root";
  TFile *inf5 = TFile::Open(file_containing_pipi);
  TTree *pipi_tree = (TTree *)inf5->Get("d0ana_mc_genmatch/VertexCompositeNtuple");

  pipi_tree->SetBranchAddress("mass", mass);
  pipi_tree->SetBranchAddress("candSize", &candSize);

  nentries = pipi_tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i)
  {
    pipi_tree->GetEntry(i);
    for (Int_t j = 0; j < candSize; ++j)
    {
      kk_mass->Fill(mass[j]);
    }
  }

  int a = 2;
  // TH1D *Signal_M1 = (TH1D *)inf2->Get("signalMass_pt_6_1000_cent_0_90_y_0_1");
  // TH1D *Signal_M2 = (TH1D *)inf3->Get("signalMass_pt_3_1000_cent_0_90_y_0_1");
  TH1D *Signal_M1 = (TH1D *)inf2->Get("signalMass_pt_" + TString::Format("%d", (int)pT_leading) + "_1000_cent_0_90_y_0_1");
  TH1D *Signal_M2 = (TH1D *)inf3->Get("signalMass_pt_" + TString::Format("%d", (int)pT_subleading) + "_1000_cent_0_90_y_0_1");

  TH1D *Swap_M1 = (TH1D *)inf2->Get("swapMass_pt_" + TString::Format("%d", (int)pT_leading) + "_1000_cent_0_90_y_0_1");
  TH1D *Swap_M2 = (TH1D *)inf3->Get("swapMass_pt_" + TString::Format("%d", (int)pT_subleading) + "_1000_cent_0_90_y_0_1");

  TH1D *SignalSwap_M1 = (TH1D *)inf2->Get("signalSwapMass_pt_" + TString::Format("%d", (int)pT_leading) + "_1000_cent_0_90_y_0_1");
  TH1D *SignalSwap_M2 = (TH1D *)inf3->Get("signalSwapMass_pt_" + TString::Format("%d", (int)pT_subleading) + "_1000_cent_0_90_y_0_1");

  // TH2D *M1M2Mass = (TH2D *)inf1->Get("M1M2"); // Full Range of F1F2 Mass

  TH2D *M1M2Mass = (TH2D *)inf1->Get("M1M2"); // Full Range of F1F2 Mass
  TH1D *M1M2hx = M1M2Mass->ProjectionX("M1M2hx");
  M1M2hx->SetTitle("M1M2 total projX");
  TH1D *M1M2hy = M1M2Mass->ProjectionY("M1M2hy");
  M1M2hy->SetTitle("M1M2 total projY");


  TF1 *SwapSW1 = new TF1("SwapSW1", "[0]*([4]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + (1-[4])*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]))", fit_range_low, fit_range_high);
  TF1 *SwapSW2 = new TF1("SwapSW2", "[0]*([4]*TMath::Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2]) + (1-[4])*TMath::Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]))", fit_range_low, fit_range_high);
  SwapSW1->SetLineColor(2);
  SwapSW1->SetParameter(0, 100);      // normalization
  SwapSW1->SetParameter(1, D0_mass);  // mean
  SwapSW1->SetParameter(2, 0.1);      // sigma1
  SwapSW1->SetParameter(3, 0.01);     // sigma 2
  SwapSW1->SetParameter(4, 0.7);      // ratio
  SwapSW1->SetParLimits(1, 1.8, 1.9); // mass
  SwapSW1->SetParLimits(2, 0.5, 1.9); // sigma1
  SwapSW1->SetParLimits(3, 0.0, 0.5); // sigma2
  SwapSW1->SetParLimits(4, 0.0, 1.0); // ratio
  FitMassPlotsTF1(SwapSW1, Swap_M1);
  // FitMassPlotsTF1(SwapSW1, SWhx);

  SwapSW2->SetLineColor(2);
  SwapSW2->SetParameter(0, 100);      // normalization
  SwapSW2->SetParameter(1, D0_mass);  // mean
  SwapSW2->SetParameter(2, 0.1);      // sigma1
  SwapSW2->SetParameter(3, 0.01);     // sigma 2
  SwapSW2->SetParameter(4, 0.7);      // ratio btw the twogg/Ca gauses
  SwapSW2->SetParLimits(1, 1.8, 1.9); // mass
  SwapSW2->SetParLimits(2, 0.5, 1.9); // sigma 1
  SwapSW2->SetParLimits(3, 0.0, 0.5); // sigma 2
  SwapSW2->SetParLimits(4, 0.0, 1.0); // ratio
  FitMassPlotsTF1(SwapSW2, Swap_M2);  // first fix swap according to the MC
  // FitMassPlotsTF1(SwapSW2, SWhy);  // first fix swap according to the MC


  // TF1 *kkpipi = new TF1("kk_pipi", "[0]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[1]), 1.97*(1+[2])) + 1.26*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
  TF1 *kkpipi = new TF1("kk_pipi", "[0]*ROOT::Math::crystalball_function(x, [4], [5], 0.025*(1+[1]), 1.95*(1+[2])) + [8]*[0]*(ROOT::Math::crystalball_function(x, [6], [7], 0.0146*(1+[9]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
  CB_fitting(kkpipi, kk_mass);
  // return;

  TF1 *fit_signal = nullptr;
  TF1 *fit_swap = nullptr;
  TF1 *fit_background = nullptr;

  TF1 *F1_total = new TF1("F1_total", 
    "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6])) + (1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6]))) + (1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6]))))) + cheb3(8) + [17]*[0]*ROOT::Math::crystalball_function(x, [20], [21], 0.025*(1+[6]), 1.95*(1+[18])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, [22], [23], 0.0146*(1+[6]), 1.7734*(1+[19])))", 
    fit_range_low, fit_range_high);

TF1 *F2_total = new TF1("F2_total", 
    "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6])) + (1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6]))) + (1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6]))))) + cheb3(8) + [17]*[0]*ROOT::Math::crystalball_function(x, [20], [21], 0.025*(1+[6]), 1.95*(1+[18])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, [22], [23], 0.0146*(1+[6]), 1.7734*(1+[19])))", 
    fit_range_low, fit_range_high);

  // TF1 *F1_total = new TF1("F1_total", "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x + [17]*[0]*ROOT::Math::crystalball_function(x, [20], [21], 0.025*(1+[6]), 1.95*(1+[18])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, [22], [23], 0.0146*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);
                                                 // 4     5                                    8                                                  6      7               9 


  F1_total->FixParameter(18, kkpipi->GetParameter(2));
  F1_total->FixParameter(19, kkpipi->GetParameter(3));
  F1_total->FixParameter(20, kkpipi->GetParameter(4));
  F1_total->FixParameter(21, kkpipi->GetParameter(5));
  F1_total->FixParameter(22, kkpipi->GetParameter(6));
  F1_total->FixParameter(23, kkpipi->GetParameter(7));
  F1_total->FixParameter(24, kkpipi->GetParameter(9));
  F1_total->FixParameter(25, kkpipi->GetParameter(8));

  F1_total->FixParameter(18,0.00);
  F1_total->FixParameter(19,0.00);
  F1_total->FixParameter(24, 3.4);

  // TF1 *F2_total = new TF1("F2_total", "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x + [17]*[0]*ROOT::Math::crystalball_function(x, [20], [21], 0.025*(1+[6]), 1.95*(1+[18])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, [22], [23], 0.0146*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);
  F2_total->FixParameter(18, kkpipi->GetParameter(2));
  F2_total->FixParameter(19, kkpipi->GetParameter(3));

  F2_total->FixParameter(20, kkpipi->GetParameter(4));
  F2_total->FixParameter(21, kkpipi->GetParameter(5));
  F2_total->FixParameter(22, kkpipi->GetParameter(6));
  F2_total->FixParameter(23, kkpipi->GetParameter(7));
  F2_total->FixParameter(24, kkpipi->GetParameter(9));
  F2_total->FixParameter(25, kkpipi->GetParameter(8));

  F2_total->FixParameter(18,0.00);
  F2_total->FixParameter(19,0.00);
  F2_total->FixParameter(24, 3.4);

  TH2F *inclusiveSignal = nullptr; TH2F *inclusiveSwap = nullptr; TH2F *inclusiveBkg = nullptr;
  GetInclusiveComponents(inf1, 10, inclusiveSignal, inclusiveSwap, inclusiveBkg, false);
  TH1D *inclusiveSignal_TH1F = inclusiveSignal->ProjectionX("inclusiveSignal_TH1F");
  TH1D *inclusiveSwap_TH1F = inclusiveSwap->ProjectionX("inclusiveSwap_TH1F");
  TH1D *inclusiveBkg_TH1F = inclusiveBkg->ProjectionX("inclusiveBkg_TH1F");
  cout << "inclusiveSignal_TH1F entries = " << inclusiveSignal_TH1F->GetEntries() << endl;
  Inclusive_dphi_fitting(F1_total, Signal_M1, SignalSwap_M1, M1M2hx, SwapSW1, Swap_M1, M1_canvas, fit_signal, fit_swap, fit_background, inclusiveSignal_TH1F, inclusiveSwap_TH1F, inclusiveBkg_TH1F);
  M1_canvas->SaveAs(Form("for_Off_MC/big_production_fits_%ibins_%islices/pT_%d_%d/m1TOTAL.pdf", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading));
  M1_canvas->SaveAs(Form("for_Off_MC/big_production_fits_%ibins_%islices/pT_%d_%d/m1TOTAL.root", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading));
  results->cd();
  F1_total->Write("F1_total");
  fit_signal->Write("S1_total");
  fit_swap->Write("SW1_total");
  fit_background->Write("B1_total");
  delete inclusiveSignal; delete inclusiveSwap; delete inclusiveBkg;
  delete inclusiveSignal_TH1F; delete inclusiveSwap_TH1F; delete inclusiveBkg_TH1F;

  inclusiveSignal = nullptr;
  inclusiveSwap = nullptr;
  inclusiveBkg = nullptr;
  GetInclusiveComponents(inf1, 10, inclusiveSignal, inclusiveSwap, inclusiveBkg, true);
  inclusiveSignal_TH1F = inclusiveSignal->ProjectionY("inclusiveSignal_TH1F");
  inclusiveSwap_TH1F = inclusiveSwap->ProjectionY("inclusiveSwap_TH1F");
  inclusiveBkg_TH1F = inclusiveBkg->ProjectionY("inclusiveBkg_TH1F");
  Inclusive_dphi_fitting(F2_total, Signal_M2, SignalSwap_M2, M1M2hy, SwapSW2, Swap_M2, M2_canvas, fit_signal, fit_swap, fit_background, inclusiveSignal_TH1F, inclusiveSwap_TH1F, inclusiveBkg_TH1F);
  results->cd();
  F2_total->Write("F2_total");
  fit_signal->Write("S2_total");
  fit_swap->Write("SW2_total");
  fit_background->Write("B2_total");

  // M1_canvas->SaveAs(Form("big_production_fits_%ibins_%islices/pT_%d_%d/m1TOTAL.pdf", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading));
  // M2_canvas->SaveAs(Form("big_production_fits_%ibins_%islices/pT_%d_%d/m2TOTAL.pdf", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading));
  // M1_canvas->SaveAs("m1.root");
  M2_canvas->SaveAs(Form("for_Off_MC/big_production_fits_%ibins_%islices/pT_%d_%d/m2TOTAL.pdf", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading));
  M2_canvas->SaveAs(Form("for_Off_MC/big_production_fits_%ibins_%islices/pT_%d_%d/m2TOTAL.root", number_of_dphi_bins, number_of_slices, (int)pT_leading, (int)pT_subleading));
  // M1_canvas->SaveAs("kkpipi_test_m1.pdf");
  // M2_canvas->SaveAs("kkpipi_test_m2.pdf");
  // return;
  // std::cin.get();

  // std::cin.get();

  // return;
  //   std::vector<TGraphErrors *> ratioGraphs;

  // // for (int dphi_bin = 1; dphi_bin < number_of_dphi_bins +1; dphi_bin++)
  // for (int dphi_bin = 4; dphi_bin < 5; dphi_bin++)
  // {
  //   TH2F *signalMC_x, *swapMC_x, *bkgMC_x;
  //   TH2F *signalMC_y, *swapMC_y, *bkgMC_y;
  //   GetComponentSlices(inf1, dphi_bin, signalMC_x, swapMC_x, bkgMC_x, false);
  //   GetComponentSlices(inf1, dphi_bin, signalMC_y, swapMC_y, bkgMC_y, true);
  //   signalMC_x->SetTitle("signalMC_x"); swapMC_x->SetTitle("swapMC_x"); bkgMC_x->SetTitle("bkgMC_x");
  //   signalMC_y->SetTitle("signalMC_y"); swapMC_y->SetTitle("swapMC_y"); bkgMC_y->SetTitle("bkgMC_y");
  //   // cout << "signal MC entries = " << signalMC->GetEntries() << endl;

  //   TF1 *fitted_signal = nullptr;
  //   TF1 *fitted_swap = nullptr;
  //   TF1 *fitted_background = nullptr;
  //   TF1 *fitted_foreground_reference = nullptr;


  //   float scale1 = 0.0;
  //   float ratio = 0.0;

  //   TGraphErrors *gr_ratio = new TGraphErrors();
  //   gr_ratio->SetName(Form("slope_par9_bin%d", dphi_bin));
  //   gr_ratio->SetTitle(Form("Bkg slope (par9) vs Slice - bin %d;Slice;Parameter 9", dphi_bin));
  //   gr_ratio->SetMarkerStyle(20);
  //   gr_ratio->SetMarkerSize(0.9);

  //   int colors[] = {kRed, kBlue, kGreen, kMagenta, kOrange, kCyan, kBlack, kGray, kViolet, kRed-1};
  //   int color = colors[(dphi_bin - 1) % (sizeof(colors) / sizeof(colors[0]))];
  //   gr_ratio->SetMarkerColor(color);
  //   gr_ratio->SetLineColor(color);


  //   TH2D *M1M2Mass_bin = (TH2D *)inf1->Get(Form("M1M2Mass_bin%d", dphi_bin)); 
  //   if (number_of_slices==25) M1M2Mass_bin->RebinX(2);
  //   if (number_of_slices==25) M1M2Mass_bin->RebinY(2);
  //   TH1D *M1M2hx_bin = M1M2Mass_bin->ProjectionX("M1M2hx_bin");
  //   M1M2hx_bin->SetTitle(Form("M1M2 projX bin %d", dphi_bin));
  //   TH1D *M1M2hy_bin = M1M2Mass_bin->ProjectionY("M1M2hy_bin");
  //   M1M2hy_bin->SetTitle(Form("M1M2 projY bin %d", dphi_bin));
    

  //   ///// fit M1 slices first
  //   TString F1_bin_name = Form("F1_bin%d", dphi_bin);
  //   // TF1 *F1_bin = new TF1(F1_bin_name, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*(x) + [10]*(x)*(x) + [11]*(x)*(x)*(x) + [17]*[0]*ROOT::Math::crystalball_function(x, [20], [21], 0.025*(1+[6]), 1.95*(1+[18])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, [22], [23], 0.0146*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);
  //   TF1 *F1_bin = new TF1(F1_bin_name, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ cheb3(8) + [17]*[0]*ROOT::Math::crystalball_function(x, [20], [21], 0.025*(1+[6]), 1.95*(1+[18])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, [22], [23], 0.0146*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);
  //   //abby fix
  //   for (a = 1; a < 25; a++)
  //   {
  //     if (a==8 || a==9 || a==10 || a==11) continue;
  //     F1_bin->FixParameter(a, F1_total->GetParameter(a));
  //   }
  //   // Fit_FBin(F1_bin, Signal_M1, SignalSwap_M1, M1M2hx_bin, SwapSW1, Swap_M1, M1_canvas, fitted_signal, fitted_swap, fitted_background);

  //   // double sigSwapRatioForM1 = getSigSwapRatio(dphi_bin, inf1, F1_bin, fitted_signal, fitted_swap, fitted_background, 
  //   //                      false, // projectY = true -> use ProjectionY (M2)
  //   //                      "F1", "S1", "SW1", "B1",
  //   //                      number_of_slices, pT_leading, pT_subleading);
  //   //  ProcessSlicesForAxis(dphi_bin, inf1, results, F1_bin, fitted_signal, fitted_swap, fitted_background, sigSwapRatioForM1,
  //   //                       M1_canvas,
  //   //                       false, // projectY = false -> use ProjectionX (M1)
  //   //                       "for_Off_MC/big_production_fits_" + TString::Format("%i", number_of_dphi_bins) + "bins_" + TString::Format("%i", number_of_slices) +"slices/pT_" + TString::Format("%d", (int)pT_leading) + "_" + TString::Format("%d", (int)pT_subleading) + "/M1_fits_bin%d.pdf",
  //   //                       // "big_production_fits_" + TString::Format("%i", number_of_dphi_bins) + "bins_" + TString::Format("%i", number_of_slices) +"slices/pT_" + TString::Format("%d", (int)pT_leading) + "_" + TString::Format("%d", (int)pT_subleading) + "/M1_fits_bin%d.pdf",
  //   //                       "F1", "S1", "SW1", "B1",
  //   //                     //  number_of_slices);
  //   //                       number_of_slices,
  //   //                       // 2,
  //   //                     //  gr_ratio,
  //   //                      signalMC_x, swapMC_x, bkgMC_x);
  //   //                     //  ratioGraphs.push_back(gr_ratio);


  
  //   // // // now fit m2 slices 
  //   // TString F2_bin_name = Form("F2_bin%d", dphi_bin);
  //   // TF1 *F2_bin = new TF1(F2_bin_name, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ cheb3(8) + [17]*[0]*ROOT::Math::crystalball_function(x, [20], [21], 0.025*(1+[6]), 1.95*(1+[18])) + [24]*[17]*[0]*(ROOT::Math::crystalball_function(x, [22], [23], 0.0146*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);
  //   // for (a = 1; a < 25; a++)
  //   // {
  //   //   if (a==8 || a==9 || a==10 || a==11) continue;
  //   //   F2_bin->FixParameter(a, F2_total->GetParameter(a));
  //   // }
  //   // Fit_FBin(F2_bin, Signal_M2, SignalSwap_M2, M1M2hy_bin, SwapSW2, Swap_M2, M2_canvas, fitted_signal, fitted_swap, fitted_background);
  //   // double sigSwapRatioForM2 = getSigSwapRatio(dphi_bin, inf1, F2_bin, fitted_signal, fitted_swap, fitted_background, 
  //   //                      true, // projectY = true -> use ProjectionY (M2)
  //   //                      "F2", "S2", "SW2", "B2",
  //   //                      number_of_slices, 
  //   //                     pT_leading, pT_subleading);
  //   // ProcessSlicesForAxis(dphi_bin, inf1, results, F2_bin, fitted_signal, fitted_swap, fitted_background, sigSwapRatioForM2,
  //   //                      M2_canvas,
  //   //                      true, // projectY = true -> use ProjectionY (M2)
  //   //                      "for_Off_MC/big_production_fits_" + TString::Format("%i", number_of_dphi_bins) + "bins_" + TString::Format("%i", number_of_slices) +"slices/pT_" + TString::Format("%d", (int)pT_leading) + "_" + TString::Format("%d", (int)pT_subleading) + "/M2_fits_bin%d.pdf",
  //   //                     //  "big_production_fits_" + TString::Format("%i", number_of_dphi_bins) + "bins_" + TString::Format("%i", number_of_slices) +"slices/pT_" + TString::Format("%d", (int)pT_leading) + "_" + TString::Format("%d", (int)pT_subleading) + "/M2_fits_bin%d.pdf",
  //   //                      "F2", "S2", "SW2", "B2",
  //   //                     //  number_of_slices);
  //   //                     //  1);
  //   //                      1,
  //   //                      signalMC_y, swapMC_y, bkgMC_y);

  // } // loop over each dphi bin

  //  gStyle->SetOptStat(0);
  //  TCanvas *c = new TCanvas("c", "Param5 vs Slice", 900, 700);

  //  TH1F *frame = new TH1F("frame",
  //                         "background slope p9 vs Slice - M1;Slice;Parameter 9",
  //                         number_of_slices, 1, number_of_slices);
  //  frame->SetMinimum(-10000);
  //  frame->SetMaximum(4000);
  //  frame->Draw();

  //  TLegend *leg = new TLegend(0.15, 0.1, 0.35, 0.3);
  //  int markers[10] = {20, 21, 22, 23, 33, 20, 21, 22, 23, 33};

  //  for (int b = 0; b < number_of_dphi_bins; b++)
  // //  for (int b = 0; b < 1; b++)
  //  {
  //    ratioGraphs[b]->SetMarkerStyle(markers[b]);
  //    ratioGraphs[b]->Draw("p SAME");
  //    leg->AddEntry(ratioGraphs[b], Form("Bin %d", b + 1), "lp");
  //  }
  // //  TLine *line = new TLine(1, 1.0, number_of_slices, 1.0);line->SetLineColor(kRed); line->Draw("same");

  //  leg->Draw();

  //  c->SaveAs("checking_slope_per_slice_bkg_par9.pdf");

  // //  c->SaveAs("big_production_EM_fits_" + TString::Format("%i", number_of_dphi_bins) + "bins_" + TString::Format("%i", number_of_slices) +"slices/pT_" + TString::Format("%d", (int)pT_leading) + "_" + TString::Format("%d", (int)pT_subleading) + "/sigSwapRatioPerSlice.pdf");
  //  c->SaveAs("for_Off_MC/big_production_fits_" + TString::Format("%i", number_of_dphi_bins) + "bins_" + TString::Format("%i", number_of_slices) +"slices/pT_" + TString::Format("%d", (int)pT_leading) + "_" + TString::Format("%d", (int)pT_subleading) + "/sigSwapRatioPerSlice.pdf");

  results->Close();
}
