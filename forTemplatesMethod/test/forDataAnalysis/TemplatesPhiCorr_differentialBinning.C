// old version was called SignalSwapOnly.C
// updated to include all bins
#include <iostream>
#include "TRandom.h"
#include "TFile.h"
#include <iostream>
#include <random>
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include <TBranch.h>
#include <vector>
#include <TRandom3.h>
#include "TROOT.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"
#include <TPDF.h>
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TLegend.h"

double histo_entries = 3.0e+08;
TH2D *S1S2Template, *B1B2Template, *SW1SW2Template, *S1B2Template, *B1S2Template, *SW1B2Template, *B1SW2Template, *S1SW2Template, *SW1S2Template;
const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[5] = {in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val};
const bool check_templates = false;
const float fit_range_low = 1.75, fit_range_high = 1.86484 + 0.125;
const int nBins_X = 100, nBins_Y = 30;

using namespace std;

Double_t MyCustomFunction(Double_t *x, Double_t *par)
{
    Double_t xx = x[0]; // x-coordinate
    Double_t yy = x[1]; // y-coordinate

    // Find bin indices once using any of the histograms, as they have the same binning
    Int_t binX = B1B2Template->GetXaxis()->FindBin(xx);
    Int_t binY = B1B2Template->GetYaxis()->FindBin(yy);
    if (binX == nBins_X+1) binX = nBins_X;
    if (binY == nBins_Y+1) binY = nBins_Y; //so bin 101 doesnt overflow

    // Retrieve bin contents using the same bin indices for all templates
    Double_t content_S1S2 = S1S2Template->GetBinContent(binX, binY);
    Double_t content_B1B2 = B1B2Template->GetBinContent(binX, binY);
    Double_t content_SW1SW2 = SW1SW2Template->GetBinContent(binX, binY);
    Double_t content_S1B2 = S1B2Template->GetBinContent(binX, binY);
    Double_t content_B1S2 = B1S2Template->GetBinContent(binX, binY);
    Double_t content_SW1B2 = SW1B2Template->GetBinContent(binX, binY);
    Double_t content_B1SW2 = B1SW2Template->GetBinContent(binX, binY);
    Double_t content_S1SW2 = S1SW2Template->GetBinContent(binX, binY);
    Double_t content_SW1S2 = SW1S2Template->GetBinContent(binX, binY);

    // Linear combination of templates with parameters
    return par[0] * content_S1S2 + par[1] * content_B1B2 + par[2] * content_SW1SW2 + par[3] * content_S1B2 + par[4] * content_B1S2 + par[5] * content_SW1S2 + par[6] * content_S1SW2 + par[7] * content_B1SW2 + par[8] * content_SW1B2;
}



void TemplatesPhiCorr_differentialBinning(int bin)
{

    TString outfile = Form("bin%d_v1.root", bin);
    TFile *results = new TFile(outfile, "recreate");

    //~~First open the file nad retrieve all histograms & their x and y projections

    cout << "working to open file..." << endl;
    // TString inputfile1 = "TH2F_08Nov_yesskip.root";
    // TString inputfile1 = "/scratch/negishi/awesole/DDbar_analysis_17Oct/pT_6_pT_3/ROOT/with_20_100_bins/M1M2.root";
    // TString inputfile1 = "/scratch/negishi/awesole/DDbar_analysis_17Oct/pT_6_pT_3/ROOT/updated_separate_DD_DbarDbar/M1M2.root";
    TString inputfile1 = "/scratch/negishi/awesole/DDbar_Jan2026/ROOT/M1M2_includingSharedDtrk.root";
    TFile *inf1 = TFile::Open(inputfile1);
    if (!inf1)
        cout << "hey no inf1 !!!!!!!!!!!!!!!!" << endl;


    // TString templates_path = "test_templates_data_updated.root";
    // TString templates_path = "templates_out_for_6_3_updatedFit_20Jan.root";
    TString templates_path = "templates_out_for_6_3_30slices_includingSharedDtrk_range175.root";
    TFile *inf3 = TFile::Open(templates_path);

    TH2F *M1M2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("M1M2MassDDbar_bin%d", bin)));
    // TH2F *M1M2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("M1M2Mass_bin%d", bin)));
    if (!M1M2Mass) {
        cout << "hey no m1m2 mass .. !!!!!" << endl;
        return;
    }
    cout << "File Successfully Opened!" << endl;
    M1M2Mass->SetMinimum(0);

 

    S1S2Template = (TH2D *)inf3->Get(Form("S1S2Mass_Bin%d_Template", bin));
    if (!S1S2Template)
    {
        cout << "NO S1S2temp!!!!!!!!!!!!!!" << endl;
        return;
    }
    S1SW2Template = (TH2D *)inf3->Get(Form("S1SW2Mass_Bin%d_Template", bin));
    S1B2Template = (TH2D *)inf3->Get(Form("S1B2Mass_Bin%d_Template", bin));
    B1S2Template = (TH2D *)inf3->Get(Form("B1S2Mass_Bin%d_Template", bin));
    B1SW2Template = (TH2D *)inf3->Get(Form("B1SW2Mass_Bin%d_Template", bin));
    B1B2Template = (TH2D *)inf3->Get(Form("B1B2Mass_Bin%d_Template", bin));
    SW1S2Template = (TH2D *)inf3->Get(Form("SW1S2Mass_Bin%d_Template", bin));
    SW1SW2Template = (TH2D *)inf3->Get(Form("SW1SW2Mass_Bin%d_Template", bin));
    SW1B2Template = (TH2D *)inf3->Get(Form("SW1B2Mass_Bin%d_Template", bin));

    TCanvas *can = new TCanvas("can", "", 1000, 1000);
    can->Divide(2, 2);
    // can->cd(1);

    cout << "defining fit function" << endl;

    // TF2 *fitF1F2 = new TF2("fitF1F2", MyCustomFunction, fit_range_low, fit_range_high, fit_range_low, fit_range_high, 9);
    TF2 *fitF1F2 = new TF2("fitF1F2", MyCustomFunction, fit_range_low, fit_range_high, fit_range_low, fit_range_high, 9);
    for (int i = 0; i < 9; i++)
    {
        fitF1F2->SetParameter(i, 1.);
    }
    fitF1F2->SetParNames("s1s2", "b1b2", "sw1sw2", "s1b2", "b1s2", "sw1s2", "s1sw2", "b1sw2", "sw1b2");
    fitF1F2->SetNpx(80);
    fitF1F2->SetNpy(80);
    // fitF1F2->SetNpx(nBins_X);
    // fitF1F2->SetNpy(nBins_Y);

    for (int i = 0; i < 9; i++)
    {
        fitF1F2->SetParLimits(i, 000.0, M1M2Mass->GetEntries());
    }
    TH2D *M1M2MassClone = (TH2D *)M1M2Mass->Clone("M1M2MassClone");
    // TCanvas *cg = new TCanvas("cg", "cg", 800, 1200);
    // cg->Divide(2, 3);
    // cg->cd(1);

    for (int i = 0; i < 1; i++)
    {

        fitF1F2->SetParameter(0, 25000.); // s1s2
        fitF1F2->SetParameter(1, 386212); // 
        fitF1F2->SetParameter(2, 20e+03); // 
        fitF1F2->SetParameter(3, 300e+03); // 
        fitF1F2->SetParameter(4, 300e+03); // 
        fitF1F2->SetParameter(5, 100); // 
        fitF1F2->SetParameter(6, 100.); // 
        fitF1F2->SetParameter(7, 277e+03); // 
        fitF1F2->SetParameter(8, 277e+03); // 

        // fitF1F2->FixParameter(0, S1S2Mass->GetEntries());     // s1s2
        // fitF1F2->FixParameter(1, B1B2Mass->GetEntries());   // b1b2
        // fitF1F2->FixParameter(2, SW1SW2Mass->GetEntries()); // sw1sw2
        // fitF1F2->FixParameter(3, S1B2Mass->GetEntries());   // s1b2
        // fitF1F2->FixParameter(4, B1S2Mass->GetEntries());   // b1s2
        // fitF1F2->FixParameter(5, SW1S2Mass->GetEntries());  // sw1s2
        // fitF1F2->FixParameter(6, S1SW2Mass->GetEntries());  // s1sw2
        // fitF1F2->FixParameter(7, B1SW2Mass->GetEntries());  // b1sw2
        // fitF1F2->FixParameter(8, SW1B2Mass->GetEntries());  // sw1b2
        can->cd(1);
        M1M2Mass->Draw("lego2");

        can->cd(2);
        cout << "---------" << endl;
        cout << "iteration = " << endl;
        for (int i = 0; i < 10; i++)
        {
            cout << i << endl;
            M1M2MassClone->Fit(fitF1F2, "LQ");
        }
        fitF1F2->FixParameter(1, fitF1F2->GetParameter(1));
        fitF1F2->FixParameter(3, fitF1F2->GetParameter(3));
        fitF1F2->FixParameter(4, fitF1F2->GetParameter(4));
        fitF1F2->FixParameter(7, fitF1F2->GetParameter(7));
        fitF1F2->FixParameter(8, fitF1F2->GetParameter(8));

        // fitF1F2->SetParLimits(2, 0.7 * fitF1F2->GetParameter(0), 1.3 * fitF1F2->GetParameter(0));
        M1M2MassClone->Fit(fitF1F2, "LQ");
        cout << "___ M1M2 mass entries" << M1M2Mass->GetEntries() << endl;
        //  fitF1F2->SetParLimits(5, 0.008 * (fitF1F2->GetParameter(0) - fitF1F2->GetParameter(2)) , 0.012*(fitF1F2->GetParameter(0) - fitF1F2->GetParameter(2)));
        //  fitF1F2->SetParLimits(6, 0.008 * (fitF1F2->GetParameter(0) - fitF1F2->GetParameter(2)) , 0.012*(fitF1F2->GetParameter(0) - fitF1F2->GetParameter(2)));

        //  fitF1F2->SetParLimits(2, fitF1F2->GetParameter(5), M1M2Mass->GetEntries());

        for (int i = 0; i < 10; i++)
        {
            cout << i << endl;
            M1M2MassClone->Fit(fitF1F2, "LQ");
        }
        fitF1F2->ReleaseParameter(0);
        fitF1F2->ReleaseParameter(1);
        fitF1F2->ReleaseParameter(2);
        fitF1F2->ReleaseParameter(3);
        fitF1F2->ReleaseParameter(4);
        fitF1F2->ReleaseParameter(5);
        fitF1F2->ReleaseParameter(6);
        fitF1F2->ReleaseParameter(7);
        fitF1F2->ReleaseParameter(8);
        for (int i = 0; i < 9; i++)
        {
            fitF1F2->SetParLimits(i, 000.0, M1M2Mass->GetEntries());
        }
        M1M2MassClone->Fit(fitF1F2, "LM");
        can->cd(2);
        M1M2MassClone->Draw("lego1");

        can->cd(3);
        TH2D *FitCheck = new TH2D("FitCheck", "Data-Fit", nBins_X, fit_range_low, fit_range_high, nBins_Y, fit_range_low, fit_range_high);
        // FitCheck->Add(M1M2MassClone);
        cout << "FitCheck n bins x = " << FitCheck->GetNbinsX() << ", n bins y = " << FitCheck->GetNbinsY() << endl;
        cout << "M1M2Mass n bins x = " << M1M2Mass->GetNbinsX() << ", n bins y = " << M1M2Mass->GetNbinsY() << endl;
        FitCheck->Add(M1M2Mass);
        FitCheck->Add(fitF1F2, -1);
        FitCheck->SetMaximum(M1M2Mass->GetMaximum());
        FitCheck->Draw("lego1");

        can->cd(4);
        fitF1F2->Draw("lego1");

        can->Update();
        can->cd();

        cout << " <<<<<<<<<<<<" << endl;
        cout << " << BIN = " << bin << " << " << endl;
        cout << " <<<<<<<<<<<<" << endl;
        cout << "fit S1S2 yield = " << fitF1F2->GetParameter(0) << endl;
        cout << " <<<<<<<<<<<<" << endl;
    }

    // can->Draw();

    // fit_range_low -= 0.02;
    // fit_range_high +=0.02;

    // getchar();

    // /*
    results->cd();

    can->Write();
    M1M2MassClone->Write();
    M1M2Mass->Write();
    fitF1F2->Write();

    results->Close();
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    { // Ensure that exactly one argument is provided
        std::cerr << "Usage: " << argv[0] << " <bin>" << std::endl;
        return 1;
    }

    int bin = std::atoi(argv[1]);              // Convert the command-line argument to an integer
    TemplatesPhiCorr_differentialBinning(bin); // Call the function with the provided bin number

    return 0;
}
