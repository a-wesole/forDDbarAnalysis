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
  f->SetParNames("normalization", "mean mass", "sigma 1", "sigma 2", "gaus ratio", "signalfrac", "smearing", "sigma swap", "const", "x", "x*x");
  f->SetParName(11, "x*x*x");
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
  f->SetParLimits(7, 0.05, 0.2); // sigma1

  // f->SetParLimits(15, 0.5, 2.0);
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

  // signal = new TF1("signal_updated", "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[7]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[7])*TMath::Gaus(x,[1],[8]*(1.0 +[6]))/(sqrt(2*3.14159)*[8]*(1.0 +[6]))))", fit_range_low, fit_range_high);
  // signal->FixParameter(0, f->GetParameter(0));  // normailzation
  // signal->FixParameter(1, f->GetParameter(1));  // mean mass
  // signal->FixParameter(2, f->GetParameter(2));  // sigma 1
  // signal->FixParameter(3, f->GetParameter(3));  // sigma 2
  // signal->FixParameter(4, f->GetParameter(4));  // ratio btw gaus 1 & 2
  // signal->FixParameter(5, f->GetParameter(5));  // signal fraction
  // signal->FixParameter(6, f->GetParameter(6));  // smearing
  // signal->FixParameter(7, f->GetParameter(12)); // ratio btw gaus 2 & 3
  // signal->FixParameter(8, f->GetParameter(14)); // sigma 3

  /// the below signal is wronng!! i will come back and find out why exactly

  swap1 = new TF1("swap1", "[0]*(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6]))))", fit_range_low, fit_range_high);

  swap1->FixParameter(0, f->GetParameter(0));
  swap1->FixParameter(5, f->GetParameter(5));
  swap1->FixParameter(6, f->GetParameter(6));
  swap1->FixParameter(7, f->GetParameter(7));
  swap1->FixParameter(13, f->GetParameter(13));
  swap1->FixParameter(15, f->GetParameter(15));
  swap1->FixParameter(16, f->GetParameter(16));

  // swap1 = new TF1("swap", "[0]*(1-[5])*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 +[6])) + (1-[4])*(TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6]))))", fit_range_low, fit_range_high);
  // swap1->FixParameter(0, f->GetParameter(0));  // norm
  // swap1->FixParameter(1, f->GetParameter(16)); // mean mass
  // swap1->FixParameter(2, f->GetParameter(7));  // sigma1
  // swap1->FixParameter(3, f->GetParameter(15)); // sigma2
  // swap1->FixParameter(4, f->GetParameter(13)); // ratio
  // swap1->FixParameter(5, f->GetParameter(5));  // signal fraction
  // swap1->FixParameter(6, f->GetParameter(6));  // smearing
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

  // TPad *pad4 = (TPad*)canvas->GetPad(4);
  // TPad *pad4Clone = (TPad*)pad4->Clone("pad4Clone");

  // TCanvas *cup = new TCanvas("cup","", 800, 600);
  // pad4Clone->SetCanvas(cup);
  // pad4Clone->Draw();

  // c1->SaveAs("TF1_outputs/M2_check.pdf");
}

void CB_fitting(TF1 *f, TH1D *kk_mass)
{
  cout << "-------------------" << endl;
  cout << "fitting CB" << endl;
  cout << "-------------------" << endl;

  f->SetLineColor(kRed);
  f->SetLineWidth(2);

  // f->FixParameter(0,0);  //start with only kk
  // f->FixParameter(2,0);
  // f->FixParameter(1, 0.0170797); //fix to match f1 total
  f->FixParameter(1, 0);
  // f->ReleaseParameter(4);
  f->ReleaseParameter(3);
  f->ReleaseParameter(0);
  // f->ReleaseParameter(1);
  f->ReleaseParameter(2);

  // TF1 *kkpipi = new TF1("kk_pipi", "[0]*ROOT::Math::crystalball_function(x, 2.2, 17, 0.0267*(1+[1]), 1.96*(1+[2])) + 4*[4]*(ROOT::Math::crystalball_function(x, 0.34, 5, 0.0146*(1+[1]), 1.7734*(1+[3])) )", fit_range_low, fit_range_high);

  for (int i = 0; i < 20; i++)
  {
    kk_mass->Fit(f, "q R N 0", "", fit_range_low, fit_range_high);
  }
  //   // kk_mass->Fit(kkpipi, "L R", "", fit_range_low, fit_range_high);
  // /// fit pipi
  // f->FixParameter(4, f->GetParameter(4));
  // f->FixParameter(3, f->GetParameter(3));
  // // f->FixParameter(4, f->GetParameter(0));
  // f->ReleaseParameter(0);
  // f->ReleaseParameter(2);

  // for (int i =0; i< 20; i++){
  //   kk_mass->Fit(f, "q R", "", fit_range_low, fit_range_high);
  // }
  kk_mass->Fit(f, "L M R q N 0", "", fit_range_low, fit_range_high);

  // f->SetParameter(0, 100); //scaling

  // kk_mass->SetTitle("kk and pipi Fit");
  // kk_mass->SetMaximum(2 * kk_mass->GetMaximum());
  // kk_mass->SetLineColor(kBlack);
  // kk_mass->Draw();
  // // pipi_mass->SetLineColor(kBlack);
  // // pipi_mass->Draw("same");
  // f->Draw("same");
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
  f->SetParNames("normalization", "mean mass", "sigma 1", "sigma 2", "gaus ratio", "signalfrac", "smearing", "sigma swap", "const", "x", "x*x");
  f->SetParName(11, "x*x*x");
  f->SetParName(12, "raio guas 2/3");
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

void fitM1M2_slices_5bins_crystalball_pythiaMC()
{

  TCanvas *M1_canvas = new TCanvas("M1_canvas", "M1_canvas", 1000, 1000);
  TCanvas *M2_canvas = new TCanvas("M2_canvas", "M2_canvas", 1000, 1000);

  // TString outfile = TString("fited_M1M2_MC_with5bins_includingSharedDtrks.root");
  TString outfile = TString("fited_M1M2_MC_test.root");
  // TString outfile = TString("fited_M1M2_MC_with5bins_30slices.root");
  // TString outfile = TString("fited_M1M2_MC_with5bins_updatedRidge.root");
  TFile *results = new TFile(outfile, "recreate");

  TFile *f_mc_comp = TFile::Open("updatedM1M2_with30Slices_v1.root");

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

  // M1_canvas->SaveAs("TF1_outputs/30slices/5bins/MC_release_sig_swap_ratio/M1_total.pdf");
  // M2_canvas->SaveAs("TF1_outputs/30slices/5bins/MC_release_sig_swap_ratio/M2_total.pdf");
  // return;



  //begin analysis for each bin 
  for (int dphi_bin = 1; dphi_bin < 6; dphi_bin++)
  // for (int dphi_bin = 1; dphi_bin < 2; dphi_bin++)
  {
    if (dphi_bin ==2 || dphi_bin == 3 || dphi_bin ==5) continue;

    TH2D *signalMC, *swapMC, *bkgMC;
    GetComponentSlices(f_mc_comp, dphi_bin, signalMC, swapMC, bkgMC);
    signalMC->SetTitle("signalMC"); swapMC->SetTitle("swaMC"); bkgMC->SetTitle("blgMC");
    cout << "signalMC entries = " << signalMC->GetEntries() <<   "swapMC entries = " << swapMC->GetEntries() <<   "bkgMC entries = " << bkgMC->GetEntries() << endl;
    for (int i = 0; i < number_of_slices ; i++)
    // for (int i = 0; i < 1 ; i++)
    {

      TString F1_slice_name = Form("F1_slice%d", i + 0);
      TF1 *F1 = new TF1(F1_slice_name, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x + [17]*[0]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[6]), 1.97*(1+[18])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);

      F1->FixParameter(1, F1_total->GetParameter(1));
      F1->FixParameter(2, F1_total->GetParameter(2));
      F1->FixParameter(3, F1_total->GetParameter(3));
      F1->FixParameter(4, F1_total->GetParameter(4));
      if (i != 99) F1->FixParameter(5, F1_total->GetParameter(5));
      F1->FixParameter(6, F1_total->GetParameter(6));
      F1->FixParameter(7, F1_total->GetParameter(7));

      F1->FixParameter(12, F1_total->GetParameter(12));
      F1->FixParameter(13, F1_total->GetParameter(13));
      F1->FixParameter(14, F1_total->GetParameter(14));
      F1->FixParameter(15, F1_total->GetParameter(15));
      F1->FixParameter(16, F1_total->GetParameter(16));
      F1->FixParameter(17, F1_total->GetParameter(17));
      F1->FixParameter(18, kkpipi->GetParameter(2));
      F1->FixParameter(19, kkpipi->GetParameter(3));

      TString F2_slice_name = Form("F2_slice%d", i + 0);
      TF1 *F2 = new TF1(F2_slice_name, "[0]*([5]*([4]*TMath::Gaus(x,[1],[2]*(1.0 +[6]))/(sqrt(2*3.14159)*[2]*(1.0 + [6]))+(1-[4])*[12]*TMath::Gaus(x,[1],[3]*(1.0 +[6]))/(sqrt(2*3.14159)*[3]*(1.0 +[6])) + (1-[4])*(1-[12])*TMath::Gaus(x,[1],[14]*(1.0 +[6]))/(sqrt(2*3.14159)*[14]*(1.0 +[6])))+(1-[5])*(([13]*TMath::Gaus(x,[16],[7]*(1.0 +[6]))/(sqrt(2*3.14159)*[7]*(1.0 +[6]))) + (1-[13])*(TMath::Gaus(x,[16],[15]*(1.0 +[6]))/(sqrt(2*3.14159)*[15]*(1.0 +[6])))))+ [8] + [9]*x + [10]*x*x + [11]*x*x*x + [17]*[0]*ROOT::Math::crystalball_function(x, 2.2, 6, 0.025*(1+[6]), 1.97*(1+[18])) + 1.26*[17]*[0]*(ROOT::Math::crystalball_function(x, 0.34, 8, 0.023*(1+[6]), 1.7734*(1+[19])))", fit_range_low, fit_range_high);

      F2->FixParameter(1, F2_total->GetParameter(1));
      F2->FixParameter(2, F2_total->GetParameter(2));
      F2->FixParameter(3, F2_total->GetParameter(3));
      F2->FixParameter(4, F2_total->GetParameter(4));
      if (i != 99) F2->FixParameter(5, F2_total->GetParameter(5));
      F2->FixParameter(6, F2_total->GetParameter(6));
      F2->FixParameter(7, F2_total->GetParameter(7));

      F2->FixParameter(12, F2_total->GetParameter(12));
      F2->FixParameter(13, F2_total->GetParameter(13));
      F2->FixParameter(14, F2_total->GetParameter(14));
      F2->FixParameter(15, F2_total->GetParameter(15));
      F2->FixParameter(16, F2_total->GetParameter(16));
      F2->FixParameter(17, F2_total->GetParameter(17));
      F2->FixParameter(18, kkpipi->GetParameter(2));
      F2->FixParameter(19, kkpipi->GetParameter(3));

      cout << "M1 slice " << i << " in dphi bin " << dphi_bin << endl;
      TH2D *M1_slice = (TH2D *)inf1->Get(Form("M1_bin%d_slice%d", dphi_bin, i+1));
      TH2D *M2_slice = (TH2D *)inf1->Get(Form("M2_bin%d_slice%d", dphi_bin, i+1));

      TH1D *M1M2hx = M1_slice->ProjectionX("M1M2hx");
      TH1D *M1M2hy = M2_slice->ProjectionY("M1M2hy");
      M1M2hx->Rebin(10);
      M1M2hy->Rebin(10);

      M1M2hx->SetTitle(Form("M1 mass, bin %d, slice %d, fixed ratio", dphi_bin, i + 1));
      M1M2hy->SetTitle(Form("M2 mass, bin %d, slice %d, fixed ratio", dphi_bin, i + 1));

      // double gaus signal and double gaus swap
      TF1 *fitted_signal = nullptr;
      TF1 *fitted_swap = nullptr;
      TF1 *fitted_background = nullptr;
      TH1 *F_getHistogram = nullptr;
      TH1 *signal_getHistogram = nullptr;
      TH1 *swap_getHistogram = nullptr;
      TH1 *background_getHistogram = nullptr;

      FitMassPlots(F1, Signal_M1, S1S2hx, M1M2hx, SwapSW1, SWhx, M1_canvas, fitted_signal, fitted_swap, fitted_background, signalMC->ProjectionY("y", i, i), swapMC->ProjectionY("wy", i, i), bkgMC->ProjectionY("by", i, i));

      
      fitted_signal->SetName(Form("S1_bin%d_slice%d", dphi_bin, i + 1));
      fitted_swap->SetName(Form("SW1_bin%d_slice%d", dphi_bin, i + 1));
      fitted_background->SetName(Form("B1_bin%d_slice%d", dphi_bin, i + 1));

      F_getHistogram = (TH1 *)F1->GetHistogram()->Clone(Form("F1histo_bin%d_slice%d", dphi_bin, i + 1));
      signal_getHistogram = (TH1 *)fitted_signal->GetHistogram()->Clone(Form("S1histo_bin%d_slice%d", dphi_bin, i + 1));
      swap_getHistogram = (TH1 *)fitted_swap->GetHistogram()->Clone(Form("SW1histo_bin%d_slice%d", dphi_bin, i + 1));
      background_getHistogram = (TH1 *)fitted_background->GetHistogram()->Clone(Form("B1histo_bin%d_slice%d", dphi_bin, i + 1));

      signal_getHistogram->SetName(Form("S1histo_bin%d_slice%d", dphi_bin, i + 1));
      swap_getHistogram->SetName(Form("SW1histo_bin%d_slice%d", dphi_bin, i + 1));
      background_getHistogram->SetName(Form("B1histo_bin%d_slice%d", dphi_bin, i + 1));
      F_getHistogram->SetTitle(Form("F1histo_bin%d_slice%d", dphi_bin, i + 1));

      results->cd();
      fitted_signal->Write();
      fitted_swap->Write();
      fitted_background->Write();
      signal_getHistogram->SetDirectory(0);
      swap_getHistogram->SetDirectory(0);
      background_getHistogram->SetDirectory(0);
      F_getHistogram->SetDirectory(0);
      signal_getHistogram->Write();
      swap_getHistogram->Write();
      background_getHistogram->Write();
      F_getHistogram->Write();

      fitted_signal = nullptr;
      fitted_swap = nullptr;
      fitted_background = nullptr;
      // FitMassPlots(F2, Signal_M2, SignalSwap_M2, M1M2hy, SwapSW2, Swap_M2, M2_canvas, fitted_signal, fitted_swap, fitted_background);
      FitMassPlots(F2, Signal_M2, S1S2hy, M1M2hy, SwapSW1, SWhy, M2_canvas, fitted_signal, fitted_swap, fitted_background, signalMC->ProjectionY("y", i+1, i+1), swapMC->ProjectionY("wy", i+1, i+1), bkgMC->ProjectionY("by", i+1, i+1));
      // FitMassPlots(F2, S1S2hy, SSWhy, M1M2hy, SwapSW2, SWhy, M2_canvas, fitted_signal, fitted_swap, fitted_background);
      results->cd();

      fitted_signal->SetName(Form("S2_bin%d_slice%d", dphi_bin, i + 1));
      fitted_swap->SetName(Form("SW2_bin%d_slice%d", dphi_bin, i + 1));
      fitted_background->SetName(Form("B2_bin%d_slice%d", dphi_bin, i + 1));

      F_getHistogram = (TH1 *)F2->GetHistogram()->Clone(Form("F2histo_bin%d_slice%d", dphi_bin, i + 1));
      signal_getHistogram = (TH1 *)fitted_signal->GetHistogram()->Clone(Form("S2histo_bin%d_slice%d", dphi_bin, i + 1));
      swap_getHistogram = (TH1 *)fitted_swap->GetHistogram()->Clone(Form("SW2histo_bin%d_slice%d", dphi_bin, i + 1));
      background_getHistogram = (TH1 *)fitted_background->GetHistogram()->Clone(Form("B2histo_bin%d_slice%d", dphi_bin, i + 1));

      signal_getHistogram->SetName(Form("S2histo_bin%d_slice%d", dphi_bin, i + 1));
      swap_getHistogram->SetName(Form("SW2histo_bin%d_slice%d", dphi_bin, i + 1));
      background_getHistogram->SetName(Form("B2histo_bin%d_slice%d", dphi_bin, i + 1));
      F_getHistogram->SetTitle(Form("F1histo_bin%d_slice%d", dphi_bin, i + 1));

      fitted_signal->Write();
      fitted_swap->Write();
      fitted_background->Write();
      signal_getHistogram->SetDirectory(0);
      swap_getHistogram->SetDirectory(0);
      background_getHistogram->SetDirectory(0);
      F_getHistogram->SetDirectory(0);
      signal_getHistogram->Write();
      swap_getHistogram->Write();
      background_getHistogram->Write();
      F_getHistogram->Write();

      // TString m1_pdf_name = Form("TF1_outputs/30slices/5bins/MC_release_sig_swap_ratio/M1_dphiBin%d.pdf", dphi_bin);
      // TString m2_pdf_name = Form("TF1_outputs/30slices/5bins/MC_release_sig_swap_ratio/M2_dphiBin%d.pdf", dphi_bin);

      // if (i == 0)
      // {
      //   M1_canvas->SaveAs(m1_pdf_name + "[");
      //   M2_canvas->SaveAs(m2_pdf_name + "[");
      // }

      // M1_canvas->SaveAs(m1_pdf_name);
      // M2_canvas->SaveAs(m2_pdf_name);

      // if (i == number_of_slices -1)
      // {

      //   M1_canvas->SaveAs(m1_pdf_name + "]");
      //   M2_canvas->SaveAs(m2_pdf_name + "]");
      // }

    } // end of for loop over slices
  } // loop over each dphi bin
  results->Close();
}
