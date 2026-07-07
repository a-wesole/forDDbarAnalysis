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
#include <iomanip>

double histo_entries = 3.0e+08;
TH2D *S1S2Template, *B1B2Template, *SW1SW2Template, *S1B2Template, *B1S2Template, *SW1B2Template, *B1SW2Template, *S1SW2Template, *SW1S2Template;
const int rebin_factor = 10;
const int number_of_slices = 30;
const bool check_templates = false;
// const float fit_range_low = 1.55,  fit_range_high = 2.2;
using namespace std;

Double_t MyCustomFunction(Double_t *x, Double_t *par)
{
    Double_t xx = x[0]; // x-coordinate
    Double_t yy = x[1]; // y-coordinate

    // Find bin indices once using any of the histograms, as they have the same binning
    Int_t binX = B1B2Template->GetXaxis()->FindBin(xx);
    Int_t binY = B1B2Template->GetYaxis()->FindBin(yy);

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

TH2D *CreateTemplate(TF1 *fit1, TF1 *fit2, TH2F *reference_histo)
{
    // your function here
    TString newTitle = TString(reference_histo->GetTitle()) + "Template";

    TH2D *Template = new TH2D(newTitle, newTitle, 1000 / rebin_factor, 1.55, 2.2, 1000 / rebin_factor, 1.55, 2.2);

    for (int i = 0; i < histo_entries; i++)
    {
        float x = fit1->GetRandom();
        float y = fit2->GetRandom();
        Template->Fill(x, y);
    }
    // Template->Scale(reference_histo->GetEntries() / histo_entries);
    Template->Scale(1.0 / Template->GetEntries());
    return Template;
}

TH2D *CreateClone(TH2D *inputHisto, TH2D *templato)
{
    TString newTitle = TString(inputHisto->GetTitle()) + " - Template";
    TH2D *histoClone = new TH2D("histoClone", newTitle,
                                inputHisto->GetNbinsX(), inputHisto->GetXaxis()->GetXmin(), inputHisto->GetXaxis()->GetXmax(),
                                inputHisto->GetNbinsY(), inputHisto->GetYaxis()->GetXmin(), inputHisto->GetYaxis()->GetXmax());
    inputHisto->Copy(*histoClone);
    histoClone->Scale(1.0 / histoClone->GetEntries());
    histoClone->SetTitle(newTitle);
    cout << "histo clone n bins x = " << histoClone->GetNbinsX() << ", n bins y = " << histoClone->GetNbinsY() << endl;
    cout << "templato n bins x = " << templato->GetNbinsX() << ", n bins y = " << templato->GetNbinsY() << endl;
    histoClone->Add(templato, -1);
    histoClone->SetLineColor(1);
    histoClone->SetMaximum(templato->GetMaximum());
    return histoClone;
}

void TemplatesPhiCorr_differentialBinning_MC(int bin)
{

   TString outfile = TString("test.root");
   TFile *results = new TFile(outfile, "recreate");

    //~~First open the file nad retrieve all histograms & their x and y projections

    cout << "working to open file..." << endl;
    // TString inputfile1 = "TH2F_08Nov_yesskip.root";
    // TString inputfile1 = "updatedM1M2_withSlices.root";
    TString inputfile1 = "updatedM1M2_with30Slices.root";
    // TString inputfile1 = "updatedM1M2_with30Slices_v1.root";
    TFile *inf1 = TFile::Open(inputfile1);
    if (!inf1)
        cout << "hey no inf1 !!!!!!!!!!!!!!!!" << endl;

    // TString fits_file_path = "TF1_outputs_12Nov_01_allbins.root";
    // TString fits_file_path = "TF1_outputs_12Nov_01_allbins.root";
    // TFile *inf2 = TFile::Open(fits_file_path);

    // TString templates_path = Form("templates_01_bin%d.root", bin); // actual one we need
    TString templates_path = "MC_templates_final_5July.root";
    // TString templates_path = "test_templates_MC_30slices.root";
    TFile *inf3 = TFile::Open(templates_path);

        TH2F *M1M2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("M1M2Mass_%d", bin)));
        if (!M1M2Mass)
            cout << "hey no m1m2 mass .. !!!!!" << endl;
        // M1M2Mass->RebinX(rebin_factor * 5);
        M1M2Mass->RebinY(rebin_factor);

        TH2F *S1S2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("S1S2Mass_%d", bin)));
        // S1S2Mass->RebinX(rebin_factor * 5);
        S1S2Mass->RebinY(rebin_factor);
        S1S2Mass->SetMinimum(0);

        TH2F *SignalSwap12Mass = dynamic_cast<TH2F *>(inf1->Get(Form("SignalSwapMass_%d", bin)));
        SignalSwap12Mass->SetMinimum(0);
        // SignalSwap12Mass->RebinX(rebin_factor * 5);
        SignalSwap12Mass->RebinY(rebin_factor);


        TH2F *SW1SW2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("SW1SW2Mass_%d", bin)));
        SW1SW2Mass->SetMinimum(0);
        // SW1SW2Mass->RebinX(rebin_factor * 5);
        SW1SW2Mass->RebinY(rebin_factor);

        TH2F *B1B2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("B1B2Mass_%d", bin)));
        B1B2Mass->SetMinimum(0);
        // B1B2Mass->RebinX(rebin_factor * 5);
        B1B2Mass->RebinY(rebin_factor);

        TH2F *S1SW2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("S1SW2Mass_%d", bin)));
        S1SW2Mass->SetMinimum(0);
        // S1SW2Mass->RebinX(rebin_factor * 5);
        S1SW2Mass->RebinY(rebin_factor);

        TH2F *SW1S2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("SW1S2Mass_%d", bin)));
        SW1S2Mass->SetMinimum(0);
        // SW1S2Mass->RebinX(rebin_factor * 5);
        SW1S2Mass->RebinY(rebin_factor);

        TH2F *S1B2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("S1B2Mass_%d", bin)));
        S1B2Mass->SetMinimum(0);
        // S1B2Mass->RebinX(rebin_factor * 5);
        S1B2Mass->RebinY(rebin_factor);

        TH2F *B1S2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("B1S2Mass_%d", bin)));
        B1S2Mass->SetMinimum(0);
        // B1S2Mass->RebinX(rebin_factor * 5);
        B1S2Mass->RebinY(rebin_factor);

        TH2F *SW1B2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("SW1B2Mass_%d", bin)));
        SW1B2Mass->SetMinimum(0);
        // SW1B2Mass->RebinX(rebin_factor * 5);
        SW1B2Mass->RebinY(rebin_factor);

        TH2F *B1SW2Mass = dynamic_cast<TH2F *>(inf1->Get(Form("B1SW2Mass_%d", bin)));
        B1SW2Mass->SetMinimum(0);
        // B1SW2Mass->RebinX(rebin_factor * 5);
        B1SW2Mass->RebinY(rebin_factor);

        cout << "File Successfully Opened!" << endl;
        /*
        M1M2Mass->Add(S1S2Mass);
        //M1M2Mass->Add(S1B2Mass);
        //M1M2Mass->Add(B1S2Mass);
        M1M2Mass->Add(BkgOnlyMass);
        M1M2Mass->Add(SwapOnlyMass);
        */
        // histo_entries = M1M2Mass->GetEntries();

        //~~Next define a signal swap and bkg function for each d0 and dbar candidiate
        //~~Signal is double gaus, swap is double gaus (see notes) bkg is linear pol

        double fit_range_low = 1.55, fit_range_high = 2.2;

        
        // S1S2Mass_Bin5_Template
        S1S2Template = (TH2D *)inf3->Get(Form("S1S2Mass_Bin%d_Template", bin));
        if (!S1S2Template) {
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

        TF2 *fitF1F2 = new TF2("fitF1F2", MyCustomFunction, 1.55, 2.2, 1.55, 2.2, 9);
        // for (int i = 0; i < 9; i++)
        // {
        //     fitF1F2->SetParameter(i, 1.);
        // }
        fitF1F2->SetParNames("s1s2", "b1b2", "sw1sw2", "s1b2", "b1s2", "sw1s2", "s1sw2", "b1sw2", "sw1b2");
        fitF1F2->SetNpx(80);
        fitF1F2->SetNpy(80);

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
            fitF1F2->SetParameter(0, 100.);                     // s1s2

            fitF1F2->SetParameter(0, 25000.);  // s1s2
            fitF1F2->SetParameter(1, 6.5e+06);   //
            fitF1F2->SetParameter(2, 25e+03);  //
            fitF1F2->SetParameter(3, 300e+03); //
            fitF1F2->SetParameter(4, 300e+03); //
            fitF1F2->SetParameter(5, 100);     //
            fitF1F2->SetParameter(6, 100.);    //
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
            M1M2Mass->Draw("lego1");

            // can->cd(2);
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

             fitF1F2->SetParLimits(2, 0.7 * fitF1F2->GetParameter(0), 1.3 * fitF1F2->GetParameter(0)); 
             //apply the same pT cuts as in data (6,3) that should give good estimate of sw1sw2 fraction to s1s2
             //draw M1(data) and M1(MC) on same plot ... same for M2
             M1M2MassClone->Fit(fitF1F2, "LQ");
             cout << "___ M1M2 mass entries" << M1M2Mass->GetEntries() << endl;


             for (int i = 0; i < 10; i++)
             {
                 cout << i << endl;
                 M1M2MassClone->Fit(fitF1F2, "LQM");
             }
             fitF1F2->ReleaseParameter(1);
             fitF1F2->ReleaseParameter(2);
             fitF1F2->ReleaseParameter(3);
             fitF1F2->ReleaseParameter(4);
             fitF1F2->ReleaseParameter(7);
             fitF1F2->ReleaseParameter(8);
             
             fitF1F2->SetParLimits(2, 0, fitF1F2->GetParameter(1));
            fitF1F2->SetParameter(2, 25e+03);  //
             M1M2MassClone->Fit(fitF1F2, "LM");
            can->cd(2);
            M1M2MassClone->Draw("lego1");

            can->cd(3);
            TH2D *FitCheck = new TH2D("FitCheck", "Data-Fit", number_of_slices, 1.55, 2.2, 1000 / rebin_factor, 1.55, 2.2);
            // FitCheck->Add(M1M2MassClone);
            cout << "FitCheck n bins x = " << FitCheck->GetNbinsX() << ", n bins y = " << FitCheck->GetNbinsY() << endl;
            cout << "M1M2Mass n bins x = " << M1M2Mass->GetNbinsX() << ", n bins y = " << M1M2Mass->GetNbinsY() << endl;
            FitCheck->Add(M1M2Mass);
            FitCheck->Add(fitF1F2, -1);
            // FitCheck->RebinX(2);
            // FitCheck->RebinY(2);
            FitCheck->SetMaximum(M1M2Mass->GetMaximum());
            FitCheck->Draw("lego1");

            can->cd(4);
            fitF1F2->Draw("lego1");

            // can->cd(5);
            // SW1B2Template->Draw("lego2");

            can->Update();
            can->cd();

            double fract = fitF1F2->GetParameter(0) / S1S2Mass->GetEntries();
            cout << " <<<<<<<<<<<<" << endl;
            cout << " << BIN = " << bin << " << " << endl;
            cout << " <<<<<<<<<<<<" << endl;
            cout << "fit S1S2 yield = " << fitF1F2->GetParameter(0) << endl;
            cout << "fit yield / data yield for S1S2 = " << fitF1F2->GetParameter(0) << " / " << S1S2Mass->GetEntries() << endl;
            cout << "= " << fract << endl;
            cout << " <<<<<<<<<<<<" << endl;
            cout << " " << endl << " " << endl;


            cout << "(0)S1S2    fit = " << setw(12) << fitF1F2->GetParameter(0)
                 << "   data = " << setw(12) << S1S2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(0) / S1S2Mass->GetEntries() << endl;

            cout << "(1)B1B2    fit = " << setw(12) << fitF1F2->GetParameter(1)
                 << "   data = " << setw(12) << B1B2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(1) / B1B2Mass->GetEntries() << endl;

            cout << "(2)SW1SW2  fit = " << setw(12) << fitF1F2->GetParameter(2)
                 << "   data = " << setw(12) << SW1SW2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(2) / SW1SW2Mass->GetEntries() << endl;

            cout << "(3)S1B2    fit = " << setw(12) << fitF1F2->GetParameter(3)
                 << "   data = " << setw(12) << S1B2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(3) / S1B2Mass->GetEntries() << endl;

            cout << "(4)B1S2    fit = " << setw(12) << fitF1F2->GetParameter(4)
                 << "   data = " << setw(12) << B1S2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(4) / B1S2Mass->GetEntries() << endl;

            cout << "(5)SW1S2   fit = " << setw(12) << fitF1F2->GetParameter(5)
                 << "   data = " << setw(12) << SW1S2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(5) / SW1S2Mass->GetEntries() << endl;

            cout << "(6)S1SW2   fit = " << setw(12) << fitF1F2->GetParameter(6)
                 << "   data = " << setw(12) << S1SW2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(6) / S1SW2Mass->GetEntries() << endl;

            cout << "(7)B1SW2   fit = " << setw(12) << fitF1F2->GetParameter(7)
                 << "   data = " << setw(12) << B1SW2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(7) / B1SW2Mass->GetEntries() << endl;

            cout << "(8)SW1B2   fit = " << setw(12) << fitF1F2->GetParameter(8)
                 << "   data = " << setw(12) << SW1B2Mass->GetEntries()
                 << "   closure = " << setw(12) << fitF1F2->GetParameter(8) / SW1B2Mass->GetEntries() << endl;
        }

        // can->Draw();

        // fit_range_low -= 0.02;
        // fit_range_high +=0.02;

        // getchar();

        // /*
        results->cd();

        S1S2Template->Write();
        S1SW2Template->Write();
        S1B2Template->Write();
        SW1S2Template->Write();
        SW1SW2Template->Write();
        SW1B2Template->Write();
        B1S2Template->Write();
        B1SW2Template->Write();
        B1B2Template->Write();
        can->Write();
    /*
    M1M2MassClone->Write();
    M1M2Mass->Write();
    fitF1F2->Write();
    */

    can->SaveAs("check_pdf.pdf");
    results->Close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) { // Ensure that exactly one argument is provided
        std::cerr << "Usage: " << argv[0] << " <bin>" << std::endl;
        return 1;
    }

    int bin = std::atoi(argv[1]); // Convert the command-line argument to an integer
    TemplatesPhiCorr_differentialBinning_MC(bin); // Call the function with the provided bin number

    return 0;
}

