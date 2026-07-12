// how to run: root -l -b -q fitM1M2.C > TF1_outputs/fitResults_pt_5_centrality_0_10_y_0_1_updatedSignal_25August_150bins.txt 2>&1
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
// const float fit_range_low = 1.75, fit_range_high = 1.86484 + 0.125, D0_mass = 1.8648;
const float fit_range_low = 1.55, fit_range_high = 2.2, D0_mass = 1.8648;
const int number_of_slices = 30;

using namespace std;

void FitMassPlotsTF1(TF1 *fitFunction, TH1D *MassPlot)
{
  for (int i = 0; i < 20; i++)
  {
    MassPlot->Fit(fitFunction, "Lq");
  }
  MassPlot->Fit(fitFunction, "q R");
}

void FitMassPlots_total(TF1 *f, TH1D *SignalMass, TH1D *SignalAndSwapMass, TH1D *MassPlot, TF1 *swap_fit, TH1D *swap_mass, TCanvas *canvas, TF1 *&signal, TF1 *&swap1, TF1 *&background)
{

  canvas->cd();
  canvas->Clear();
  swap_fit->SetLineColor(kGreen);
  canvas->Divide(2, 2);
  canvas->cd(1);
  f->SetLineColor(2);
  f->SetLineWidth(2);
  f->SetParNames("normalization", "mean mass", "sigma 1", "sigma 2", "gaus ratio", "signalfrac", "smearing", "sigma swap", "pol0", "pol1", "pol2");
  f->SetParName(11, "pol3");
  f->SetParName(12, "raio guas 2/3");
  f->SetParName(13, "swap gaus ratio");
  f->SetParName(14, "sigma 3");
  f->SetParName(15, "sigma2 swap");
  f->SetParName(16, "swap mean");
  f->SetParName(17, "CB scaling");
  // f->SetParName(20, "CB left ");
  f->SetParName(18, "CB r sigma");
  f->SetParName(19, "CB l sigma");
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

  // f->SetParLimits(2, 0.01, 0.1);   // sigma 1
  f->SetParLimits(1, 1.85, 1.88); // sigma 1
  // f->SetParLimits(3, 0.001, 0.05); // sigma 2
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
  SignalMass->Fit(f, "qM", "", fit_range_low, fit_range_high);
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
  f->SetParLimits(7, 0.05, 0.2); // sigma1

  // f->SetParLimits(15, 0.5, 2.0);
  SignalAndSwapMass->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);
  cout << "---------------- fitting overall signal + swap" << endl;
  for (int i = 0; i < 29; i++)
  {
    SignalAndSwapMass->Fit(f, "q", "", fit_range_low, fit_range_high);
  }
  SignalAndSwapMass->Fit(f, "qM", "", fit_range_low, fit_range_high);

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
  // swap->Draw("same");
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
  // f->ReleaseParameter(11); // 12 *x^3
  f->ReleaseParameter(0); // Normalization
  cout << "---------------- fitting overall background" << endl;

  MassPlot->Fit(f, "q", "", fit_range_low, fit_range_high);
  MassPlot->Fit(f, "q", "", fit_range_low, fit_range_high);

  f->ReleaseParameter(1);    // mean mass ------allow data to have different mass peak mean than MC
  f->ReleaseParameter(6);    // ratio ------- allow data to have different peak width than MC
  f->SetParameter(6, 0);     // scaling factor
  f->SetParLimits(6, -1, 1); // scaling factor
  f->ReleaseParameter(17);
  // f->SetParLimits(17, 0, f->GetParameter(0) * f->GetParameter(5));
  f->SetParLimits(17, 0, f->GetParameter(5));
  // f->SetParLimits(20, 0, f->GetParameter(0) * f->GetParameter(5));

  cout << "fitting allowing to scale mc to data" << endl;
  for (int i = 0; i < 10; i++)
  {
    MassPlot->Fit(f, "L q", "", fit_range_low, fit_range_high);
  }
  MassPlot->Fit(f, "L m N 0", "", fit_range_low, fit_range_high); ///

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

  background = new TF1("background", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", fit_range_low, fit_range_high);
  background->FixParameter(0, f->GetParameter(8));
  background->FixParameter(1, f->GetParameter(9));
  background->FixParameter(2, f->GetParameter(10));
  background->FixParameter(3, f->GetParameter(11));

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

  MassPlot->SetMinimum(0.0);
  // MassPlot->SetMaximum(1500);
  MassPlot->Draw("EP");
  swap1->Draw("same");
  signal->SetLineColor(4);
  signal->Draw("same");
  background->SetLineColor(kBlack);
  background->Draw("same");
  kkpipi->Draw("same");
  canvas->Update();

}

void CB_fitting(TF1 *f, TH1D *kk_mass)
{
  cout << "-------------------" << endl;
  cout << "fitting CB" << endl;
  cout << "-------------------" << endl;

  f->SetLineColor(kRed);
  f->SetLineWidth(2);

  f->FixParameter(1, 0);
  f->ReleaseParameter(3);
  f->ReleaseParameter(0);
  f->ReleaseParameter(2);


  for (int i = 0; i < 20; i++)
  {
    kk_mass->Fit(f, "q R N 0", "", fit_range_low, fit_range_high);
  }
  kk_mass->Fit(f, "L M R q N 0", "", fit_range_low, fit_range_high);

}

void GetComponentSlices(TFile *f, int dphi_bin, TH2D *&signalBin, TH2D *&swapBin, TH2D *&bkgBin)
{

  cout << "made it to GetCom" << endl;
  TH2D *S1S2Mass = (TH2D *)f->Get(Form("S1S2Mass_%d", dphi_bin));
  TH2D *SW1S2Mass = (TH2D *)f->Get(Form("SW1S2Mass_%d", dphi_bin));
  TH2D *B1S2Mass = (TH2D *)f->Get(Form("B1S2Mass_%d", dphi_bin));
  TH2D *S1SW2Mass = (TH2D *)f->Get(Form("S1SW2Mass_%d", dphi_bin));
  TH2D *SW1SW2Mass = (TH2D *)f->Get(Form("SW1SW2Mass_%d", dphi_bin));
  TH2D *B1SW2Mass = (TH2D *)f->Get(Form("B1SW2Mass_%d", dphi_bin));
  TH2D *S1B2Mass = (TH2D *)f->Get(Form("S1B2Mass_%d", dphi_bin));
  TH2D *SW1B2Mass = (TH2D *)f->Get(Form("SW1B2Mass_%d", dphi_bin));
  TH2D *B1B2Mass = (TH2D *)f->Get(Form("B1B2Mass_%d", dphi_bin));

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

void FitMassPlots(TF1 *f, TH1D *SignalMass, TH1D *SignalAndSwapMass, TH1D *MassPlot, TF1 *swap_fit, TH1D *swap_mass, TCanvas *canvas, TF1 *&signal, TF1 *&swap1, TF1 *&background, TH1D *signalMC, TH1D *swapMC, TH1D *bkgMC)
{

  canvas->cd();
  canvas->Clear();
  swap_fit->SetLineColor(kGreen);
  // canvas->Divide(2, 2);
  canvas->cd(1);
  f->SetLineColor(2);
  f->SetLineWidth(2);
  f->SetParNames("normalization", "mean mass", "sigma 1", "sigma 2", "gaus ratio", "signalfrac", "smearing", "sigma swap", "pol0", "pol1", "pol2");
  f->SetParName(11, "pol3");
  f->SetParName(12, "ratio gaus 2/3");
  f->SetParName(13, "swap gaus ratio");
  f->SetParName(14, "sigma 3");
  f->SetParName(15, "sigma2 swap");
  f->SetParName(16, "swap mean");
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
  cout << "fitting allowing to scale mc to data" << endl;
  for (int i = 0; i < 10; i++)
  {
    MassPlot->Fit(f, "L q", "", fit_range_low, fit_range_high);
  }
  MassPlot->Fit(f, "L m N 0", "", fit_range_low, fit_range_high); ///

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

  background = new TF1("background", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", fit_range_low, fit_range_high);
  background->FixParameter(0, f->GetParameter(8));
  background->FixParameter(1, f->GetParameter(9));
  background->FixParameter(2, f->GetParameter(10));
  background->FixParameter(3, f->GetParameter(11));

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

  MassPlot->SetMinimum(0.0);
  // MassPlot->SetMaximum(15000);
  MassPlot->SetTitle("Mass Plot with Fit Components");
  MassPlot->Draw("EP");
  signalMC->SetMarkerStyle(4); signalMC->SetMarkerColor(kBlack);
  signalMC->RebinX(10);
  signalMC->Draw("EP same");
  swapMC->RebinX(10);
  swapMC->SetMarkerStyle(20); swapMC->SetMarkerColor(kBlack);
  swapMC->Draw("EP same");
  bkgMC->SetMarkerStyle(20);  bkgMC->SetMarkerColor(kBlack);
  bkgMC->RebinX(10);
  bkgMC->Draw("EP same");
  swap1->Draw("same");
  signal->SetLineColor(4);
  signal->Draw("same");
  background->SetLineColor(kBlack);
  background->Draw("same");
  kkpipi->Draw("same");
  // swap_mass->SetLineColor(kGreen + 2);
  // swap_mass->Draw("same");
  canvas->Update();

  // c1->SaveAs("TF1_outputs/M2_check.pdf");
}

// Process slices for one axis (M1 or M2). Behavior matches your existing per-slice code.
void ProcessSlicesForAxis(int dphi_bin,
                          TFile *inf1,
                          TFile *results,
                          TF1 *F_bin_ref,               // F1_bin or F2_bin (already created & fixed)
                          TF1 *fitted_signal_ref,      // fitted_signal from overall fit
                          TF1 *fitted_swap_ref,        // fitted_swap from overall fit
                          TF1 *fitted_background_ref,  // fitted_background from overall fit
                          TCanvas *canvas,
                          const char *sliceNameA,      // e.g. "M1_bin%d_slice%d"
                          const char *sliceNameB,      // e.g. "M2_bin%d_slice%d"
                          bool projectY,               // false => use ProjectionX (M1); true => ProjectionY (M2)
                          const char *pdfPattern,      // e.g. "M1_fits_bin%d.pdf"
                          const char *prefixF,         // e.g. "F1" or "F2"
                          const char *prefixS,         // e.g. "S1" or "S2"
                          const char *prefixSW,        // e.g. "SW1" or "SW2"
                          const char *prefixB,        // e.g. "B1" or "B2"
                          int number_of_slices_local,  // pass number_of_slices
                          TGraphErrors *gr = nullptr)  // graph
{
  // open multipage PDF
  TString pdfName = Form(pdfPattern, dphi_bin);
  canvas->cd();
  canvas->Print(Form("%s[", pdfName.Data()));

  for (int i = 0; i < number_of_slices_local; ++i)
  {
    TH2D *sliceA = (TH2D *)inf1->Get(Form(sliceNameA, dphi_bin, i + 1));
    TH2D *sliceB = (TH2D *)inf1->Get(Form(sliceNameB, dphi_bin, i + 1));

    TH1D *hX = nullptr;
    TH1D *hY = nullptr;
    if (sliceA) hX = sliceA->ProjectionX(Form("%s_M1M2hx_bin%d_slice%d", prefixF, dphi_bin, i));
    if (sliceB) hY = sliceB->ProjectionY(Form("%s_M1M2hy_bin%d_slice%d", prefixF, dphi_bin, i));

    // pick the histogram to fit/draw
    TH1D *toFit = projectY ? hY : hX;
    if (!toFit) continue; // safety: if missing, skip

    cout << "toFit histogram: " << toFit->GetName() << ", nbins: " << toFit->GetNbinsX() << endl;
    toFit->Rebin(10);
    toFit->SetTitle(Form("%s mass, bin %d, slice %d", projectY ? "M2" : "M1", dphi_bin, i + 1));

    // clone the F_bin reference for this slice (protected)
    TF1 *f_ref_clone = (TF1 *)F_bin_ref->Clone(Form("f31_protected_%s_bin%d_slice%d", prefixF, dphi_bin, i));

    // fitted_foreground_slice: p0 scale, p1 ratio (same lambda structure as original)
    TF1 *fitted_foreground_slice = new TF1(Form("%s_bin%d_slice%d", prefixF, dphi_bin, i),
      [f_ref_clone](double *x, double *p) {
        f_ref_clone->SetParameter(5, p[1]);
        return p[0] * f_ref_clone->Eval(x[0]);
      }, fit_range_low, fit_range_high, 2);

    fitted_foreground_slice->SetParameter(0, 1.0);
    fitted_foreground_slice->SetParameter(1, f_ref_clone->GetParameter(5));
    fitted_foreground_slice->SetParLimits(1, 0.0, 1.0);

    // perform the same three quiet fits as in your code
    toFit->Fit(fitted_foreground_slice, "RQ0");
    toFit->Fit(fitted_foreground_slice, "RQ0");
    toFit->Fit(fitted_foreground_slice, "RQ0");

    float scale = fitted_foreground_slice->GetParameter(0);
    float ratio = fitted_foreground_slice->GetParameter(1);

    if(gr) {
      int n = gr->GetN();
      gr->SetPoint(n, i, ratio);
      gr->SetPointError(n, 0, fitted_foreground_slice->GetParError(1));
    }

    // create fitted_signal_slice (clone via lambda, same as original)
    TF1 *fitted_signal_slice = new TF1(Form("fitted_signal_slice_%s_%d", prefixS, i),
      [fitted_signal_ref](double *x, double *p) {
        fitted_signal_ref->SetParameter(5, p[1]);
        return p[0] * fitted_signal_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 2);
    fitted_signal_slice->FixParameter(0, scale);
    fitted_signal_slice->FixParameter(1, ratio);

    // create fitted_swap_slice
    TF1 *fitted_swap_slice = new TF1(Form("fitted_swap_slice_%s_%d", prefixSW, i),
      [fitted_swap_ref](double *x, double *p) {
        fitted_swap_ref->SetParameter(5, p[1]);
        return p[0] * fitted_swap_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 2);
    fitted_swap_slice->FixParameter(0, scale);
    fitted_swap_slice->FixParameter(1, ratio);

    // create fitted_background_slice
    TF1 *fitted_background_slice = new TF1(Form("fitted_background_slice_%s_%d", prefixB, i),
      [fitted_background_ref](double *x, double *p) {
        return p[0] * fitted_background_ref->Eval(x[0]);
      }, fit_range_low, fit_range_high, 1);
    fitted_background_slice->FixParameter(0, scale);

    // draw exactly as before
    toFit->SetLineColor(kGray); toFit->SetMarkerColor(kGray);
    toFit->SetMarkerStyle(20); toFit->SetMinimum(0);

    toFit->Draw("E");
    fitted_foreground_slice->SetLineColor(kRed); fitted_foreground_slice->Draw("same");
    fitted_signal_slice->SetLineColor(kBlue); fitted_signal_slice->Draw("same");
    fitted_swap_slice->SetLineColor(kGreen); fitted_swap_slice->Draw("same");
    fitted_background_slice->SetLineColor(kBlack); fitted_background_slice->Draw("same");

    // print page to PDF
    canvas->Update();
    canvas->Print(pdfName);

    // write TF1s and histograms to results file (same names pattern as original)
    results->cd();
    fitted_foreground_slice->Write(Form("%s_bin%d_slice%d", prefixF, dphi_bin, i));
    fitted_signal_slice->Write(Form("%s_bin%d_slice%d", prefixS, dphi_bin, i));
    fitted_swap_slice->Write(Form("%s_bin%d_slice%d", prefixSW, dphi_bin, i));
    fitted_background_slice->Write(Form("%s_bin%d_slice%d", prefixB, dphi_bin, i));

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


void fitM1M2_slices_5bins_crystalball_pythiaMC()
{

  TCanvas *M1_canvas = new TCanvas("M1_canvas", "M1_canvas", 1000, 1000);
  TCanvas *M2_canvas = new TCanvas("M2_canvas", "M2_canvas", 1000, 1000);

  // TString outfile = TString("fited_M1M2_MC_with5bins_includingSharedDtrks.root");
  TString outfile = TString("fited_M1M2_MC_test.root");
  TFile *results = new TFile(outfile, "recreate");

  TFile *f_mc_comp = TFile::Open("updatedM1M2_with30Slices.root");

  double F1F2Entries, F1SB2Entries, SB1F2Entries, SB1SB2Entries;
  // TString file_containing_M1M2_slices = "updatedM1M2_withSlices.root";
  TString file_containing_M1M2_slices = "updatedM1M2_with30Slices.root"; //skipping uncorrelated d0s
  TFile *inf1 = TFile::Open(file_containing_M1M2_slices);

  TString file_containing_M1_signal_swap = "/scratch/negishi/awesole/DDbar_analysis/signal_swap_histos_150bins/ROOT/massPlots_pt_6_1000_cent_0_90_y_0_1.root";
  TFile *inf2 = TFile::Open(file_containing_M1_signal_swap);

  TString file_containing_M2_signal_swap = "/scratch/negishi/awesole/DDbar_analysis/signal_swap_histos_150bins/ROOT/massPlots_pt_3_1000_cent_0_90_y_0_1.root";
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

  int a = 6; // total sum of all 5 bins
  TH1D *Signal_M1 = (TH1D *)inf2->Get("signalMass_pt_6_1000_cent_0_90_y_0_1");
  TH1D *Signal_M2 = (TH1D *)inf3->Get("signalMass_pt_3_1000_cent_0_90_y_0_1");

  TH1D *Swap_M1 = (TH1D *)inf2->Get("swapMass_pt_6_1000_cent_0_90_y_0_1");
  TH1D *Swap_M2 = (TH1D *)inf3->Get("swapMass_pt_3_1000_cent_0_90_y_0_1");

  TH1D *SignalSwap_M1 = (TH1D *)inf2->Get("signalSwapMass_pt_6_1000_cent_0_90_y_0_1");
  TH1D *SignalSwap_M2 = (TH1D *)inf3->Get("signalSwapMass_pt_3_1000_cent_0_90_y_0_1");


  TH2D *B1SW2Mass = (TH2D *)inf1->Get(Form("B1SW2Mass_%d",a));
  TH2D *B1S2Mass = (TH2D *)inf1->Get(Form("B1S2Mass_%d",a));
  TH2D *SW1S2Mass = (TH2D *)inf1->Get(Form("SW1S2Mass_%d",a)); // Full Range of F1F2 Mass
  SW1S2Mass->SetMinimum(0.5);
  SW1S2Mass->RebinX(1);
  SW1S2Mass->RebinY(1);
  TH2D *SW1B2Mass = (TH2D *)inf1->Get(Form("SW1B2Mass_%d",a));
  TH2D *S1B2Mass = (TH2D *)inf1->Get(Form("S1B2Mass_%d",a));
  TH2D *S1SW2Mass = (TH2D *)inf1->Get(Form("S1SW2Mass_%d",a)); // Full Range of F1F2 Mass
  S1SW2Mass->SetMinimum(0.5);
  S1SW2Mass->RebinX(1);
  S1SW2Mass->RebinY(1);

  TH2D *M1M2Mass = (TH2D *)inf1->Get(Form("M1M2Mass_%d",a));
  TH1D *M1M2hx = M1M2Mass->ProjectionX("M1M2hx");
  cout << " --------------------------------------M1M2hx nbins: " << M1M2hx->GetNbinsX() << endl;
  TH1D *M1M2hy = M1M2Mass->ProjectionY("M1M2hy");

  TH2D *S1S2Mass = (TH2D *)inf1->Get(Form("S1S2Mass_%d",a));
  S1S2Mass->SetMinimum(0);
  TH1D *S1S2hx = S1S2Mass->ProjectionX("S1S2hx");
  TH1D *S1S2hy = S1S2Mass->ProjectionY("S1S2hy");

  S1S2hx->Add(S1B2Mass->ProjectionX());
  S1S2hx->Add(S1SW2Mass->ProjectionX());
  S1S2hy->Add(B1S2Mass->ProjectionY());
  S1S2hy->Add(SW1S2Mass->ProjectionY());

  TH2D *SignalSwap12Mass = (TH2D *)inf1->Get(Form("SignalSwapMass_%d",a)); // Full Range of F1F2 Mass
  SignalSwap12Mass->SetMinimum(0.5);
  TH1D *SSWhx = SignalSwap12Mass->ProjectionX("SSWhx");
  TH1D *SSWhy = SignalSwap12Mass->ProjectionY("SSWhy");

  SSWhx->Add(SW1B2Mass->ProjectionX());
  SSWhx->Add(S1B2Mass->ProjectionX());
  SSWhy->Add(B1SW2Mass->ProjectionY());
  SSWhy->Add(B1S2Mass->ProjectionY());

  TH2D *SwapOnlyMass = (TH2D *)inf1->Get(Form("SW1SW2Mass_%d",a)); // Full Range of F1F2 Mass
  SwapOnlyMass->SetMinimum(0.5);
  SwapOnlyMass->RebinX(1);
  SwapOnlyMass->RebinY(1);
  SwapOnlyMass->SetTitle("SW1SW2 Mass");
  TH1D *SWhx = SwapOnlyMass->ProjectionX("SWhx");
  TH1D *SWhy = SwapOnlyMass->ProjectionY("SWhy");

  SWhx->Add(SW1B2Mass->ProjectionX());
  SWhx->Add(SW1S2Mass->ProjectionX());
  SWhy->Add(B1SW2Mass->ProjectionY());
  SWhy->Add(S1SW2Mass->ProjectionY());


  //step - fit the swap alone 


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
  // FitMassPlotsTF1(SwapSW1, Swap_M1);
  FitMassPlotsTF1(SwapSW1, SWhx);

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
  // FitMassPlotsTF1(SwapSW2, Swap_M2);  // first fix swap according to the MC
  FitMassPlotsTF1(SwapSW2, SWhy);  // first fix swap according to the MC

  // step - fit the kk pipi curves with the entire M1M2

  TF1 *kkpipi = new TF1("kk_pipi", "[0]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[1]), 1.97*(1+[2])) + 1.26*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);
  CB_fitting(kkpipi, kk_mass);

  TF1 *fit_signal = nullptr;
  TF1 *fit_swap = nullptr;
  TF1 *fit_background = nullptr;

  TF1 *F1_total = new TF1("F1_total", "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x + [17]*[0]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[6]), 1.97*(1+[18])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);
  F1_total->FixParameter(18, kkpipi->GetParameter(2));
  F1_total->FixParameter(19, kkpipi->GetParameter(3));
  TF1 *F2_total = new TF1("F2_total", "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x + [17]*[0]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[6]), 1.97*(1+[18])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);
  F2_total->FixParameter(18, kkpipi->GetParameter(2));
  F2_total->FixParameter(19, kkpipi->GetParameter(3));

  FitMassPlots_total(F1_total, S1S2hx, SSWhx, M1M2hx, SwapSW1, SWhx, M1_canvas, fit_signal, fit_swap, fit_background);

  FitMassPlots_total(F2_total, S1S2hy, SSWhy, M1M2hy, SwapSW2, SWhy, M2_canvas, fit_signal, fit_swap, fit_background);


  M1_canvas->SaveAs("M1_total.pdf");
  M2_canvas->SaveAs("M2_total.pdf");

  std::vector<TGraphErrors *> ratioGraphs;



  // step - begin analysis for each dphi bin 

  for (int dphi_bin = 1; dphi_bin < 6; dphi_bin++)
  // for (int dphi_bin = 1; dphi_bin < 2; dphi_bin++)
  {
    TH1D *Signal_M1 = (TH1D *)inf2->Get("signalMass_pt_6_1000_cent_0_90_y_0_1");
    TH1D *Signal_M2 = (TH1D *)inf3->Get("signalMass_pt_3_1000_cent_0_90_y_0_1");

    TH1D *Swap_M1 = (TH1D *)inf2->Get("swapMass_pt_6_1000_cent_0_90_y_0_1");
    TH1D *Swap_M2 = (TH1D *)inf3->Get("swapMass_pt_3_1000_cent_0_90_y_0_1");

    TH1D *SignalSwap_M1 = (TH1D *)inf2->Get("signalSwapMass_pt_6_1000_cent_0_90_y_0_1");
    TH1D *SignalSwap_M2 = (TH1D *)inf3->Get("signalSwapMass_pt_3_1000_cent_0_90_y_0_1");

    TH2D *B1SW2Mass = (TH2D *)inf1->Get(Form("B1SW2Mass_%d", dphi_bin));
    TH2D *B1S2Mass = (TH2D *)inf1->Get(Form("B1S2Mass_%d", dphi_bin));
    TH2D *SW1S2Mass = (TH2D *)inf1->Get(Form("SW1S2Mass_%d", dphi_bin)); // Full Range of F1F2 Mass
    SW1S2Mass->SetMinimum(0.5);
    SW1S2Mass->RebinX(1);
    SW1S2Mass->RebinY(1);
    TH2D *SW1B2Mass = (TH2D *)inf1->Get(Form("SW1B2Mass_%d", dphi_bin));
    TH2D *S1B2Mass = (TH2D *)inf1->Get(Form("S1B2Mass_%d", dphi_bin));
    TH2D *S1SW2Mass = (TH2D *)inf1->Get(Form("S1SW2Mass_%d", dphi_bin)); // Full Range of F1F2 Mass
    S1SW2Mass->SetMinimum(0.5);
    S1SW2Mass->RebinX(1);
    S1SW2Mass->RebinY(1);

    TH2D *M1M2Mass = (TH2D *)inf1->Get(Form("M1M2Mass_%d", dphi_bin));
    TH1D *M1M2hx = M1M2Mass->ProjectionX("M1M2hx");
    TH1D *M1M2hy = M1M2Mass->ProjectionY("M1M2hy");

    TH2D *S1S2Mass = (TH2D *)inf1->Get(Form("S1S2Mass_%d", dphi_bin));
    S1S2Mass->SetMinimum(0);
    TH1D *S1S2hx = S1S2Mass->ProjectionX("S1S2hx");
    TH1D *S1S2hy = S1S2Mass->ProjectionY("S1S2hy");

    S1S2hx->Add(S1B2Mass->ProjectionX());
    S1S2hx->Add(S1SW2Mass->ProjectionX());
    S1S2hy->Add(B1S2Mass->ProjectionY());
    S1S2hy->Add(SW1S2Mass->ProjectionY());

    TH2D *SignalSwap12Mass = (TH2D *)inf1->Get(Form("SignalSwapMass_%d", dphi_bin)); // Full Range of F1F2 Mass
    SignalSwap12Mass->SetMinimum(0.5);
    TH1D *SSWhx = SignalSwap12Mass->ProjectionX("SSWhx");
    TH1D *SSWhy = SignalSwap12Mass->ProjectionY("SSWhy");

    SSWhx->Add(SW1B2Mass->ProjectionX());
    SSWhx->Add(S1B2Mass->ProjectionX());
    SSWhy->Add(B1SW2Mass->ProjectionY());
    SSWhy->Add(B1S2Mass->ProjectionY());

    TH2D *SwapOnlyMass = (TH2D *)inf1->Get(Form("SW1SW2Mass_%d", dphi_bin)); // Full Range of F1F2 Mass
    SwapOnlyMass->SetMinimum(0.5);
    SwapOnlyMass->RebinX(1);
    SwapOnlyMass->RebinY(1);
    SwapOnlyMass->SetTitle("SW1SW2 Mass");
    TH1D *SWhx = SwapOnlyMass->ProjectionX("SWhx");
    TH1D *SWhy = SwapOnlyMass->ProjectionY("SWhy");

    SWhx->Add(SW1B2Mass->ProjectionX());
    SWhx->Add(SW1S2Mass->ProjectionX());
    SWhy->Add(B1SW2Mass->ProjectionY());
    SWhy->Add(S1SW2Mass->ProjectionY());

    // step - this part briefly gets the actual MC data points of each component S,SW,B for that bin so that we can plot them on top of the fit for each slice later....
    TH2D *signalMC, *swapMC, *bkgMC;
    GetComponentSlices(f_mc_comp, dphi_bin, signalMC, swapMC, bkgMC);
    signalMC->SetTitle("signalMC"); swapMC->SetTitle("swapMC"); bkgMC->SetTitle("bkgMC");

    TF1 *fitted_signal = nullptr;
    TF1 *fitted_swap = nullptr;
    TF1 *fitted_background = nullptr;
    TF1 *fitted_foreground_reference = nullptr;
    // step - this part fits the projection x and projection y of the M1M2 for this bin, gets an overall fit
    TString F1_bin_name = Form("F1_bin%d", dphi_bin);
    TF1 *F1_bin = new TF1(F1_bin_name, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x + [17]*[0]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[6]), 1.97*(1+[18])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);

    F1_bin->FixParameter(1, F1_total->GetParameter(1));
    F1_bin->FixParameter(2, F1_total->GetParameter(2));
    F1_bin->FixParameter(3, F1_total->GetParameter(3));
    F1_bin->FixParameter(4, F1_total->GetParameter(4));
    F1_bin->FixParameter(5, F1_total->GetParameter(5));
    F1_bin->FixParameter(6, F1_total->GetParameter(6));
    F1_bin->FixParameter(7, F1_total->GetParameter(7));

    F1_bin->FixParameter(12, F1_total->GetParameter(12));
    F1_bin->FixParameter(13, F1_total->GetParameter(13));
    F1_bin->FixParameter(14, F1_total->GetParameter(14));
    F1_bin->FixParameter(15, F1_total->GetParameter(15));
    F1_bin->FixParameter(16, F1_total->GetParameter(16));
    F1_bin->FixParameter(17, F1_total->GetParameter(17));
    F1_bin->FixParameter(18, kkpipi->GetParameter(2));
    F1_bin->FixParameter(19, kkpipi->GetParameter(3));
    FitMassPlots(F1_bin, Signal_M1, S1S2hx, M1M2hx, SwapSW1, SWhx, M1_canvas, fitted_signal, fitted_swap, fitted_background, signalMC->ProjectionX(), swapMC->ProjectionX(), bkgMC->ProjectionX());


    // results->cd();
    // F1_bin->Write();

    


    float scale1 = 0.0;
    float ratio = 0.0;

    TGraphErrors *gr_ratio = new TGraphErrors();
    gr_ratio->SetName(Form("gr_ratio_bin%d", dphi_bin));
    gr_ratio->SetTitle(Form("Signal/Swap ratio vs Slice - bin %d;Slice;Parameter 5", dphi_bin));
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->SetMarkerSize(0.9);

    int colors[] = {kRed, kBlue, kGreen + 2, kMagenta, kCyan + 2, kOrange};
    int color = colors[(dphi_bin - 1) % (sizeof(colors) / sizeof(colors[0]))];
    gr_ratio->SetMarkerColor(color);
    gr_ratio->SetLineColor(color);

    cout << "going to process slices for dphi bin " << dphi_bin << endl;
    ProcessSlicesForAxis(dphi_bin, inf1, results, F1_bin, fitted_signal, fitted_swap, fitted_background,
                         M1_canvas,
                         "M1_bin%d_slice%d", "M2_bin%d_slice%d",
                         false, // projectY = false -> use ProjectionX (M1)
                         "M1_fits_bin%d.pdf",
                         "F1", "S1", "SW1", "B1",
                         number_of_slices,
                        gr_ratio);
  ratioGraphs.push_back(gr_ratio);

    TString F2_bin_name = Form("F2_bin%d", dphi_bin);
    TF1 *F2_bin = new TF1(F2_bin_name, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x + [17]*[0]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[6]), 1.97*(1+[18])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);
    F2_bin->FixParameter(1, F2_total->GetParameter(1));
    F2_bin->FixParameter(2, F2_total->GetParameter(2));
    F2_bin->FixParameter(3, F2_total->GetParameter(3));
    F2_bin->FixParameter(4, F2_total->GetParameter(4));
    F2_bin->FixParameter(5, F2_total->GetParameter(5));
    F2_bin->FixParameter(6, F2_total->GetParameter(6));
    F2_bin->FixParameter(7, F2_total->GetParameter(7));

    F2_bin->FixParameter(12, F2_total->GetParameter(12));
    F2_bin->FixParameter(13, F2_total->GetParameter(13));
    F2_bin->FixParameter(14, F2_total->GetParameter(14));
    F2_bin->FixParameter(15, F2_total->GetParameter(15));
    F2_bin->FixParameter(16, F2_total->GetParameter(16));
    F2_bin->FixParameter(17, F2_total->GetParameter(17));
    F2_bin->FixParameter(18, kkpipi->GetParameter(2));
    F2_bin->FixParameter(19, kkpipi->GetParameter(3));
    FitMassPlots(F2_bin, Signal_M2, S1S2hy, M1M2hy, SwapSW2, SWhy, M2_canvas, fitted_signal, fitted_swap, fitted_background, signalMC->ProjectionY(), swapMC->ProjectionY(), bkgMC->ProjectionY());
    ProcessSlicesForAxis(dphi_bin, inf1, results, F2_bin, fitted_signal, fitted_swap, fitted_background,
                     M2_canvas,
                     "M1_bin%d_slice%d", "M2_bin%d_slice%d",
                     true,                        // projectY = true -> use ProjectionY (M2)
                     "M2_fits_bin%d.pdf",
                     "F2", "S2", "SW2", "B2",
                     number_of_slices);




  } // loop over each dphi bin
  TCanvas *c = new TCanvas("c", "Param5 vs Slice", 900, 700);

TH1F *frame = new TH1F("frame",
                       "Signal/Swap ratio vs Slice - pythiaMC;Slice;Parameter 5",
                       30, 1, 30);
frame->Draw();

TLegend *leg = new TLegend(0.15, 0.7, 0.35, 0.9);
int markers[5] = {20, 21, 22, 23, 33};

for (int b = 0; b < 5; b++)
{
    ratioGraphs[b]->SetMarkerStyle(markers[b]);
    ratioGraphs[b]->Draw("P SAME");
    leg->AddEntry(ratioGraphs[b], Form("Bin %d", b + 1), "p");
}

leg->Draw();

c->SaveAs("ratio_per_slice_per_bin.pdf");

  results->Close();
}
