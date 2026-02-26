// this creates all 9 templates
// how to run: root CreateTemplates.C

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
#include "TH1D.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"
#include <TPDF.h>
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"

const float histo_entries = 75.0e+07;
const float org_entries = 3.0e+06;
// const float nbins = 150.0;
const float nbins = 100.0;
TH2D *S1S2Template, *B1B2Template, *SW1SW2Template, *S1B2Template, *B1S2Template, *SW1B2Template, *B1SW2Template, *S1SW2Template, *SW1S2Template;
const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[5] = {in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val};
// const float fit_range_low = 1.86484-0.125, fit_range_high = 1.86484 + 0.125;
const float fit_range_low = 1.75, fit_range_high = 1.86484 + 0.125;
// const float fit_range_low = 1.55, fit_range_high = 2.2;

const int number_of_slices = 30;
const float slice_width = (fit_range_high - fit_range_low) / number_of_slices;

std::vector<std::vector<TF1 *>> F1_fits(6), S1_fits(6), SW1_fits(6), B1_fits(6);
std::vector<std::vector<TF1 *>> F2_fits(6), S2_fits(6), SW2_fits(6), B2_fits(6);
std::vector<std::vector<TH1 *>> S1_getHistograms(6), SW1_getHistograms(6), B1_getHistograms(6), F1_getHistograms(6);
std::vector<std::vector<TH1 *>> S2_getHistograms(6), SW2_getHistograms(6), B2_getHistograms(6), F2_getHistograms(6);

using namespace std;

void LoadFitsFromFile(const char *filename)
{
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file with slices fits!!!" << std::endl;
        return;
    }

    for (int dphi_bin = 1; dphi_bin < 6; dphi_bin++)
    {
        for (int i = 0; i < number_of_slices; ++i)
        {
            // M1 fits
            F1_fits[dphi_bin - 1].push_back((TF1 *)file->Get(Form("F1_bin%d_slice%d", dphi_bin, i + 1)));
            S1_fits[dphi_bin - 1].push_back((TF1 *)file->Get(Form("S1_bin%d_slice%d", dphi_bin, i + 1)));
            SW1_fits[dphi_bin - 1].push_back((TF1 *)file->Get(Form("SW1_bin%d_slice%d", dphi_bin, i + 1)));
            B1_fits[dphi_bin - 1].push_back((TF1 *)file->Get(Form("B1_bin%d_slice%d", dphi_bin, i + 1)));

            // M2 fits
            F2_fits[dphi_bin - 1].push_back((TF1 *)file->Get(Form("F2_bin%d_slice%d", dphi_bin, i + 1)));
            S2_fits[dphi_bin - 1].push_back((TF1 *)file->Get(Form("S2_bin%d_slice%d", dphi_bin, i + 1)));
            SW2_fits[dphi_bin - 1].push_back((TF1 *)file->Get(Form("SW2_bin%d_slice%d", dphi_bin, i + 1)));
            B2_fits[dphi_bin - 1].push_back((TF1 *)file->Get(Form("B2_bin%d_slice%d", dphi_bin, i + 1)));

            TH1 *s1_hist = dynamic_cast<TH1 *>(file->Get(Form("S1histo_bin%d_slice%d", dphi_bin, i + 1)));
            s1_hist->SetDirectory(0);
            S1_getHistograms[dphi_bin - 1].push_back(s1_hist);
            TH1 *sw1_hist = dynamic_cast<TH1 *>(file->Get(Form("SW1histo_bin%d_slice%d", dphi_bin, i + 1)));
            sw1_hist->SetDirectory(0);
            SW1_getHistograms[dphi_bin - 1].push_back(sw1_hist);
            TH1 *b1_hist = dynamic_cast<TH1 *>(file->Get(Form("B1histo_bin%d_slice%d", dphi_bin, i + 1)));
            b1_hist->SetDirectory(0);
            B1_getHistograms[dphi_bin - 1].push_back(b1_hist);

            TH1 *s2_hist = dynamic_cast<TH1 *>(file->Get(Form("S2histo_bin%d_slice%d", dphi_bin, i + 1)));
            if (s2_hist)
            {
                s2_hist->SetDirectory(0);
                S2_getHistograms[dphi_bin - 1].push_back(s2_hist);
            }
            else
            {
                std::cerr << "Missing or invalid histogram: S2histo_slice" << i + 1 << std::endl;
                return;
            }
            TH1 *sw2_hist = dynamic_cast<TH1 *>(file->Get(Form("SW2histo_bin%d_slice%d", dphi_bin, i + 1)));
            if (sw2_hist)
            {
                sw2_hist->SetDirectory(0);
                SW2_getHistograms[dphi_bin - 1].push_back(sw2_hist);
            }
            else
            {
                std::cerr << "Missing or invalid histogram: SW2histo_slice" << i + 1 << std::endl;
                return;
            }
            TH1 *b2_hist = dynamic_cast<TH1 *>(file->Get(Form("B2histo_bin%d_slice%d", dphi_bin, i + 1)));
            if (b2_hist)
            {
                b2_hist->SetDirectory(0);
                B2_getHistograms[dphi_bin - 1].push_back(b2_hist);
            }
            else
            {
                std::cerr << "Missing or invalid histogram: B2histo_slice" << i + 1 << std::endl;
                return;
            }

            TH1 *f1_hist = dynamic_cast<TH1 *>(file->Get(Form("F1histo_bin%d_slice%d", dphi_bin, i + 1)));
            f1_hist->SetDirectory(0);
            F1_getHistograms[dphi_bin - 1].push_back(f1_hist);
            TH1 *f2_hist = dynamic_cast<TH1 *>(file->Get(Form("F2histo_bin%d_slice%d", dphi_bin, i + 1)));
            f2_hist->SetDirectory(0);
            F2_getHistograms[dphi_bin - 1].push_back(f2_hist);

            // S2_getHistograms.push_back((TH1D*)file->Get(Form("S2histo_slice%d", i+1)));
            // SW2_getHistograms.push_back((TH1D*)file->Get(Form("SW2histo_slice%d", i+1)));
            // B2_getHistograms.push_back((TH1D*)file->Get(Form("B2histo_slice%d", i+1)));
        }
    }

    // file->Close();
}

TH2D *CreateTemplate(TF1 *fit1, TF1 *fit2, const TString &histoTitle)
{
    // your function here
    TString newTitle = histoTitle + "Template";

    TH2D *Template = new TH2D(newTitle, newTitle, nbins, fit_range_low, fit_range_high, number_of_slices, fit_range_low, fit_range_high);
    Template->SetOption("lego1");
    Template->SetXTitle("M1 Mass");
    Template->SetYTitle("M2 Mass");

    for (int i = 0; i < org_entries; i++)
    {
        float x = fit1->GetRandom();
        float y = fit2->GetRandom();
        Template->Fill(x, y);
    }
    Template->Scale(1.0 / Template->GetEntries());
    return Template;
}

TH2D *CreateTemplateNew(std::vector<TH1 *> &v1, std::vector<TH1 *> &v2, const TString &histoTitle)
{

    TH1 *leftM2_B2_histo = v2[0];
    TH1 *rightM2_B2_histo = v2[number_of_slices - 1];
    float count = 0;
    std::vector<float> scaling_factors;

    for (auto histogram : v1)
    {
        float point = fit_range_low + (0.5 + count) * slice_width; // choose the point in the middle of the slice
        float left_most_point = fit_range_low;
        float scaling_factor = leftM2_B2_histo->GetBinContent(leftM2_B2_histo->FindBin(point)) / histogram->GetBinContent(histogram->FindBin(left_most_point));
        // cout << "point = " << point << " && leftM2_B2_histo->GetBinContent(leftM2_B2_histo->FindBin(point)) = " << leftM2_B2_histo->GetBinContent(leftM2_B2_histo->FindBin(point)) << endl;
        // cout << " && histogram->GetBinContent(histogram->FindBin(left_most_point)) = " << histogram->GetBinContent(histogram->FindBin(left_most_point)) << endl;

        /// this is just to check the scaling//
        // result should be almost exactly 1//
        float point2 = point;
        float right_most_point = fit_range_high;
        float new_scaling_factor = rightM2_B2_histo->GetBinContent(rightM2_B2_histo->FindBin(point2)) / histogram->GetBinContent(histogram->FindBin(right_most_point));

        // cout << "-------" << endlo
        // cout << "histogram title = " << histogram->GetTitle() << endl;
        // cout << "for slice " << count + 1 << " scaling factor = " << scaling_factor << endl;
        // cout << "histogram title = " << histogram->GetTitle() << endl;
        // cout << "for slice " << count + 1 << " redone scaling factor = " << new_scaling_factor << endl;
        count += 1;
        histogram->Scale(scaling_factor);
        scaling_factors.push_back(scaling_factor);
    }

    int nbinsX = v2[0]->GetNbinsX();
    TString newTitle = histoTitle + "Template";
    TH2D *Template = new TH2D(newTitle, newTitle, nbins, fit_range_low, fit_range_high, number_of_slices, fit_range_low, fit_range_high);
    Template->SetOption("lego1");
    Template->SetXTitle("M1 Mass");
    Template->SetYTitle("M2 Mass");
    // cout << "nbinsX = " << nbinsX << " , number_of_slices = " << number_of_slices << endl;
    // int nbinsX = B2_getHistograms[0]->GetNbinsX();
    for (int M2_bin_number = 0; M2_bin_number < number_of_slices; M2_bin_number++) // M1 bins (0-30)
    // for (int M2_bin_number = 0; M2_bin_number < 1; M2_bin_number++) // M1 bins (0-20)
    {
        for (int M1_bin_number = 0; M1_bin_number < nbins; M1_bin_number++) // M2 bins (0-100)
        {
            // cout << "----" << endl;
            // cout << "M1_bin_number (0-30)= " << M1_bin_number << " , M2_bin_number(0-100) = " << M2_bin_number << endl;
            double content = v1[M2_bin_number]->GetBinContent(M1_bin_number + 1);
            double error = v1[M2_bin_number]->GetBinError(M1_bin_number + 1);
            // cout << "content = " << content << " , error = " << error << endl;
            Template->SetBinContent(M1_bin_number + 1, M2_bin_number + 1, content);
            Template->SetBinError(M1_bin_number + 1, M2_bin_number + 1, error);
        }
    }
    // cout << "template entries = " << Template->GetSumOfWeights() << endl;
    Template->Scale(1.0 / Template->GetSumOfWeights());
    Template->SetMinimum(0);
    // for (int i = 0; i < v2.size(); i++)
    // {
    //     TH1 *histo = v1[i];
    //     histo->Scale(1.0 / scaling_factors[i]);
    // }
    return Template;
}


// void CreateTemplates_wSlices(const TString &outfile)
void CreateTemplates_wSlices_5bins_M2Slices()
{

    // TString outfile = "test_templates_data.root";
    // TString outfile = "templates_out_for_6_3_30slices.root";
    TString outfile = "templates_out_for_6_3_30slices_includingSharedDtrk_range175.root";
    TFile *results = new TFile(outfile, "recreate");
    cout << "output_name = " << outfile << endl;

    // TString slices_fits_file_path = "/home/awesole/forDDbarAnalysis/forTemplatesMethod/test/TF1_outputs/20slices/pT_6_pT_3_20Slices_wGetHistogram_with175.root";
    // TString slices_fits_file_path = "fited_M1M2_with5bins_pT6_3_30slices.root";
    TString slices_fits_file_path = "fited_M1M2_with5bins_pT6_3_30slices_updated24Feb_forSharedDtrTrks_175range.root";

    LoadFitsFromFile(slices_fits_file_path);


    for (int dphi = 0; dphi < 5; dphi++)
    // for (int dphi = 0; dphi < 1; dphi++)
    {
        cout << "=================================" << endl;
        cout << "=== for phi bin = " << dphi + 1 << " ===" << endl;
        cout << "=================================" << endl;


        cout << "Generating template (1) - S1S2" << endl;
        TH2D *S1S2Template = CreateTemplate(S1_fits[dphi][0], S2_fits[dphi][0], Form("S1S2Mass_Bin%d_", dphi + 1));
        cout << "Generating template (2) - SW1SW2" << endl;
        TH2D *SW1SW2Template = CreateTemplate(SW1_fits[dphi][0], SW2_fits[dphi][0], Form("SW1SW2Mass_Bin%d_", dphi + 1));
        cout << "Generating template (3) - S1SW2" << endl;
        TH2D *S1SW2Template = CreateTemplate(S1_fits[dphi][0], SW2_fits[dphi][0], Form("S1SW2Mass_Bin%d_", dphi+1));
        cout << "Generating template (4) - SW1S2" << endl;
        TH2D *SW1S2Template = CreateTemplate(SW1_fits[dphi][0], S2_fits[dphi][0], Form("SW1S2Mass_Bin%d_", dphi+1));

        cout << "Generating template (5) - S1B2" << endl;
        TH2D *S1B2Template = CreateTemplateNew(S2_getHistograms[dphi], B2_getHistograms[dphi], Form("S1B2Mass_Bin%d_", dphi+1));
        cout << "Generating template (6) - B1S2" << endl;
        TH2D *B1S2Template = CreateTemplateNew(B2_getHistograms[dphi], S2_getHistograms[dphi], Form("B1S2Mass_Bin%d_", dphi+1));

        cout << "Generating template (7) - SW1B2" << endl;
        TH2D *SW1B2Template = CreateTemplateNew(SW1_getHistograms[dphi], B2_getHistograms[dphi], Form("SW1B2Mass_Bin%d_", dphi+1));
        cout << "Generating template (8) - B1SW2" << endl;
        TH2D *B1SW2Template = CreateTemplateNew(B1_getHistograms[dphi], SW2_getHistograms[dphi], Form("B1SW2Mass_Bin%d_", dphi+1));


    

        cout << "Generating template (9) - B1B2" << endl;
        TH2D *B1B2Template = CreateTemplateNew(B1_getHistograms[dphi], B2_getHistograms[dphi], Form("B1B2Mass_Bin%d_", dphi+1));
        B1B2Template->SetXTitle("M1 Mass");
        B1B2Template->SetYTitle("M2 Mass");
        B1B2Template->SetMinimum(0);
        // B1B2Template->Draw("lego1");

        TCanvas temp_canvas("temp_canvas", "temp_canvas", 1000, 1000);

        temp_canvas.Divide(3, 3);
        temp_canvas.cd(1);
        S1S2Template->Draw();
        temp_canvas.cd(2);
        SW1SW2Template->Draw();
        temp_canvas.cd(3);
        S1SW2Template->Draw();
        temp_canvas.cd(4);
        SW1S2Template->Draw();
        temp_canvas.cd(5);
        S1B2Template->Draw();
        temp_canvas.cd(6);
        B1S2Template->Draw();
        temp_canvas.cd(7);
        SW1B2Template->Draw();
        temp_canvas.cd(8);
        B1SW2Template->Draw();
        temp_canvas.cd(9);
        B1B2Template->Draw();


        results->cd();

        S1S2Template->Write();
        SW1SW2Template->Write();
        S1SW2Template->Write();
        SW1S2Template->Write();
        S1B2Template->Write();
        B1S2Template->Write();
        SW1B2Template->Write();
        B1SW2Template->Write();
        B1B2Template->Write();
        temp_canvas.SetName(Form("temp_canvas_bin%d", dphi+1));
        temp_canvas.Write();
    }
    results->Close();
}

