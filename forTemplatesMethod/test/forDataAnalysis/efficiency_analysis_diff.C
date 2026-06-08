#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TDirectory.h"
#include <iostream>

#include "BDTHandler.h"

int extract_pT_hat(const std::string &filename)
{ // this reads the ptHat value from the filename and reutrns the int to the main
    size_t pos = filename.find("PT-");
    if (pos == std::string::npos)
        return -1; // not found

    pos += 3; // move past PT-
    std::string num = "";

    while (pos < filename.size() && isdigit(filename[pos]))
    {
        num += filename[pos];
        pos++;
    }

    if (num.size() == 0)
        return -1; // no digits found
    return std::stoi(num);
}

bool pass_pT_hat_cuts(int pthat, float pt)
{
    if (pthat == 0)
        return (pt < 10);
    if (pthat == 1)
        return false; // skipping all pthat1 for now ....
    if (pthat == 8)
        return false; // skipping all pthat8 for now ....
    if (pthat == 10)
        return (pt >= 10 && pt < 20);
    if (pthat == 20)
        return (pt >= 20 && pt < 30);
    if (pthat == 30)
        return (pt >= 30 && pt < 40);
    if (pthat == 40)
        return (pt >= 40 && pt < 60);
    if (pthat == 60)
        return (pt >= 60 && pt < 80);
    if (pthat == 80)
        return (pt >= 80 && pt < 200);
    return false; // anything else return false (primarily pthat 8 )
}

void efficiency_analysis_diff(TString i_start = "", TString i_end = "", TString output_file = "test.root")
{

    //////// // ifstream file_stream("/home/awesole/forDDbarAnalysis/forTemplatesMethod/test/forDataAnalysis/lists/mcfiles_withGen_andBDT.list");
    // ifstream file_stream("ptHat_0_mc_files.list");
    ifstream file_stream("/home/awesole/forDDbarAnalysis/forTemplatesMethod/test/forDataAnalysis/lists/MC_v6_only2Dtrs.list");
    // TFile *results = new TFile("efficiency_D0_pthat0_kpi.root", "RECREATE");
    TFile *results = new TFile(output_file, "RECREATE");
    string filename;
    int ifile = 0;
    int istart = atoi(i_start);
    int iend = atoi(i_end);
    auto bdt = std::make_unique<BDTHandler>();
    // Define pt bin edges
    const int nBins = 12;
    const float pT_minimum = 0.5;

    // Histograms: denominator = gen, numerator = reco matched

    TH1F *hGen_raw = new TH1F("hGen_raw", "Gen D0; p_{T} [GeV]; Entries", 200, 0, 100);
    TH1F *hReco_raw = new TH1F("hReco_raw", "Reco D0 matched; p_{T} [GeV]; Entries", 200, 0, 100);
    hGen_raw->Sumw2();
    hReco_raw->Sumw2();

    // Final stitched histograms (after applying pThat weights)
    TH1F *hGen_stitched = new TH1F("hGen_stitched", "Weighted stitched GEN pT; p_{T} [GeV]; Entries", 200, 0, 100);
    hGen_stitched->Sumw2();

    TH1F *hReco_stitched = new TH1F("hReco_stitched", "Weighted stitched RECO pT; p_{T} [GeV]; Entries", 200, 0, 100);
    hReco_stitched->Sumw2();
    //------------------------------------------------------
    TH1F *hGen_stitched_L = new TH1F("hGen_stitched_L", "Weighted stitched GEN pT; p_{T} [GeV]; Entries", 200, 0, 100);
    hGen_stitched_L->Sumw2();

    TH1F *hReco_stitched_L = new TH1F("hReco_stitched_L", "Weighted stitched RECO pT; p_{T} [GeV]; Entries", 200, 0, 100);
    hReco_stitched_L->Sumw2();
    //------------------------------------------------------
    TH1F *hGen_stitched_S = new TH1F("hGen_stitched_S", "Weighted stitched GEN pT; p_{T} [GeV]; Entries", 200, 0, 100);
    hGen_stitched_S->Sumw2();

    TH1F *hReco_stitched_S = new TH1F("hReco_stitched_S", "Weighted stitched RECO pT; p_{T} [GeV]; Entries", 200, 0, 100);
    hReco_stitched_S->Sumw2();

    // pThat bins we will stitch
    std::vector<int> pthat_bins = {80, 60, 40, 30, 20, 10, 0};

    ///// for hGen pThat distributions iwth no cuts, so we can scale the tail
    std::map<int, TH1F *> hGen_ptHat;
    for (int p : pthat_bins)
    {
        hGen_ptHat[p] = new TH1F(Form("hGen_ptHat_%d", p),
                                 Form("Gen pT for pThat %d", p),
                                 200, 0, 100);
        hGen_ptHat[p]->Sumw2();
    }

    std::vector<float> gen_pts;
    std::vector<int> gen_pthats;

    std::vector<float> reco_pts;
    std::vector<int> reco_pthats;

    // --- Set up branches for GEN tree ---
    Int_t N_genD0s;
    Float_t genD0_pt[1000], genD0_y[1000]; // must match max_gen_size

    while (true)
    {
        file_stream >> filename;

        if (file_stream.eof())
            break;

        if (ifile < istart)
        {
            ifile++;
            continue;
        }

        if (ifile >= iend)
            break;
        if (ifile % 2 == 0)
        {
            ifile++;
            continue;
        } // for testing, only run on half the files (every other file)
        if (ifile % 3 == 0)
        {
            ifile++;
            continue;
        } // for testing, only run on half the files (every other file)
        if (ifile % 5 == 0)
        {
            ifile++;
            continue;
        } // for testing, only run on half the files (every other file)

        // cout << "file = " << filename.c_str() << endl;
        TFile *f = TFile::Open(filename.c_str());
        int pthat = extract_pT_hat(filename);
        if (pthat == 1 || pthat == 8)
        {
            ifile++;
            continue;
        }

        if (!f || f->IsZombie())
        {
            // std::cerr << "Cannot open your file" << std::endl;
            ifile++;
            continue;
        }

        // Open file
        // cout << "File opened successfully" << std::endl;

        // Get directory
        TDirectory *d = (TDirectory *)f->Get("d0Analyzer");
        if (!d)
        {
            // std::cerr << "Cannot find directory d0Analyzer" << std::endl;
            ifile++;
            continue;
        }

        // Get trees
        TTree *tReco = (TTree *)d->Get("VertexCompositeNtuple");
        TTree *tGen = (TTree *)d->Get("AllGensNtuple");
        if (!tReco || !tGen || (tGen->GetEntries() != tReco->GetEntries()))
        {
            // std::cerr << "Cannot find trees or they do not have the same number of entries" << std::endl;

            ifile++;
            continue;
        }
        // cout << "Trees retrieved successfully" << std::endl;
        cout << "ifile = " << ifile << endl;
        if (ifile % 100 == 0)
        {
            cout << "----- ------ ----- -----" << endl;
            cout << "progress = " << ifile << " / 5553 = " << (float(ifile) / 5553) * 100 << " %" << endl;
        }

        tGen->SetBranchStatus("*", 0);
        tGen->SetBranchStatus("N_genD0s", 1);
        tGen->SetBranchStatus("genD0_pt", 1);
        tGen->SetBranchStatus("genD0_y", 1);
        tGen->SetBranchAddress("N_genD0s", &N_genD0s);
        tGen->SetBranchAddress("genD0_pt", genD0_pt);
        tGen->SetBranchAddress("genD0_y", genD0_y);

        // --- Set up branches for RECO tree ---
        // Adjust names/types if your ntuple uses different ones
        Int_t D0_N, centrality;
        Float_t D0_pt[50000], y[50000], Dtrk1Pt[50000], Dtrk2Pt[50000], mva[50000], pT_gen[50000];
        bool matchGEN[50000], isSwap[50000];

        tReco->SetBranchStatus("*", 0);
        tReco->SetBranchStatus("candSize", 1);
        tReco->SetBranchStatus("pT", 1);
        tReco->SetBranchStatus("matchGEN", 1);
        tReco->SetBranchStatus("isSwap", 1);
        tReco->SetBranchStatus("pTD1", 1);
        tReco->SetBranchStatus("pTD2", 1);
        tReco->SetBranchStatus("y", 1);
        tReco->SetBranchStatus("pT_gen", 1);
        tReco->SetBranchStatus("mva", 1); // if needed

        tReco->SetBranchAddress("candSize", &D0_N);
        tReco->SetBranchAddress("centrality", &centrality);
        tReco->SetBranchAddress("pT", D0_pt);
        tReco->SetBranchAddress("matchGEN", matchGEN);
        tReco->SetBranchAddress("isSwap", isSwap);
        tReco->SetBranchAddress("pTD1", Dtrk1Pt);
        tReco->SetBranchAddress("pTD2", Dtrk2Pt);
        tReco->SetBranchAddress("mva", mva);
        tReco->SetBranchAddress("y", y);
        tReco->SetBranchAddress("pT_gen", pT_gen);

        float bdt_cut_value = -999.9;

        // Sanity check: same number of entries
        Long64_t nEntriesReco = tReco->GetEntries();
        Long64_t nEntriesGen = tGen->GetEntries();
        if (nEntriesReco != nEntriesGen)
        {
            std::cerr << "Reco and Gen trees have different number of entries!" << std::endl;
            std::cerr << "file: " << filename.c_str() << std::endl;
            return;
        }

        // Loop over events
        // cout << "Starting event loop...nEntries: " << tReco->GetEntries() << endl;
        for (Long64_t i = 0; i < tReco->GetEntries(); ++i)
        {
            tGen->GetEntry(i);
            tReco->GetEntry(i);

            // Fill denominator: all gen D0s
            // cout << "NgenD0s = " << N_genD0s << ", D0_N = " << D0_N << endl;
            for (int j = 0; j < N_genD0s; ++j)
            {

                if (genD0_pt[j] >= 100)
                    continue;
                if (fabs(genD0_y[j]) > 1)
                    continue;
                // cout << "filling gen " << endl;
                hGen_ptHat[pthat]->Fill(genD0_pt[j]);
                hGen_raw->Fill(genD0_pt[j]);

                // for raw distrobution (before scaling)
                //  cout << "filled gen " << endl;
                if (!pass_pT_hat_cuts(pthat, genD0_pt[j]))
                    continue;
                // cout << "Filling denominator..." << endl;

                gen_pts.push_back(genD0_pt[j]);
                gen_pthats.push_back(pthat);
            }

            // Fill numerator: reco D0 with matchGEN == 1
            // cout << "Filling numerator..." << endl;
            for (int k = 0; k < D0_N; ++k)
            {

                if (D0_pt[k] >= 100)
                    continue;
                if (!(matchGEN[k] == 1 && isSwap[k] == 0))
                    continue;
                if (fabs(y[k]) > 1)
                    continue;
                if (Dtrk1Pt[k] < pT_minimum || Dtrk2Pt[k] < pT_minimum)
                    continue;
                bdt_cut_value = bdt->getBDTCut(y[k], centrality, D0_pt[k]);
                if (mva[k] <= bdt_cut_value)
                    continue; // background cuts
                hReco_raw->Fill(D0_pt[k]);

                if (!pass_pT_hat_cuts(pthat, pT_gen[k]))
                    continue;
                double pt = D0_pt[k];

                reco_pts.push_back(pt);
                reco_pthats.push_back(pthat);
            }
        }
        f->Close();
        delete f;
        ifile++;
    } // while

    // begin calcualting the weights
    cout << "Calculating weights..." << std::endl;
    std::map<int, double> weights;
    weights[80] = 1.0; // keep pthat 80 fixed
    for (int i = 1; i < pthat_bins.size(); i++)
    {
        int higher = pthat_bins[i - 1];
        int lower = pthat_bins[i];

        double lower_bound = higher;
        double upper_bound = 200.0;

        int bin_low = hGen_ptHat[higher]->FindBin(lower_bound);
        int bin_high = hGen_ptHat[higher]->FindBin(upper_bound);

        double I_high = hGen_ptHat[higher]->Integral(bin_low, bin_high);
        double I_low = hGen_ptHat[lower]->Integral(bin_low, bin_high);

        weights[lower] = weights[higher] * (I_high / I_low);
        cout << "pThat " << lower << ": weight = " << weights[lower] << std::endl;
    }
    // Fill stitched GEN histogram
    for (size_t i = 0; i < gen_pts.size(); i++)
    {
        int pthat = gen_pthats[i];
        double w = weights[pthat];
        hGen_stitched->Fill(gen_pts[i], w);
        if (gen_pts[i] >= 6) hGen_stitched_L->Fill(gen_pts[i], w);
        if (gen_pts[i] >= 3) hGen_stitched_S->Fill(gen_pts[i], w);
    }

    // Fill stitched RECO histogram
    for (size_t i = 0; i < reco_pts.size(); i++)
    {
        int pthat = reco_pthats[i];
        double w = weights[pthat];
        hReco_stitched->Fill(reco_pts[i], w);
        if (reco_pts[i] >= 6) hReco_stitched_L->Fill(reco_pts[i], w);
        if (reco_pts[i] >= 3) hReco_stitched_S->Fill(reco_pts[i], w);
    }
    TH1F *hEff_stitched = (TH1F *)hReco_stitched->Clone("hEff_stitched");
    hEff_stitched->SetTitle("Stitched efficiency; p_{T} [GeV]; Efficiency");
    hEff_stitched->Divide(hGen_stitched);

    TH1F *hEff_stitched_L = (TH1F *)hReco_stitched_L->Clone("hEff_stitched_L");
    hEff_stitched_L->SetTitle("Stitched efficiency_L; p_{T} [GeV]; Efficiency");
    hEff_stitched_L->Divide(hGen_stitched_L);

    TH1F *hEff_stitched_S = (TH1F *)hReco_stitched_S->Clone("hEff_stitch_S");
    hEff_stitched_S->SetTitle("Stitched efficiency_S; p_{T} [GeV]; Efficiency");
    hEff_stitched_S->Divide(hGen_stitched_S);

    TH1F *hEff_raw = (TH1F *)hReco_raw->Clone("hEff_raw");
    hEff_raw->SetTitle("Raw efficiency; p_{T} [GeV]; Efficiency");
    hEff_raw->Divide(hGen_raw);

    // ------------------------------------------------------------
    // Draw summary canvas
    // ------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "Efficiency Summary", 1800, 1200);
    c1->Divide(3, 2);

    c1->cd(1);
    hGen_raw->SetTitle("hGen_raw");
    hGen_raw->Draw();
    c1->cd(2);
    hReco_raw->SetTitle("hReco_raw");
    hReco_raw->Draw();

    c1->cd(3);
    hEff_raw->SetTitle("hEff_raw");
    hEff_raw->SetMaximum(1.0);
    hEff_raw->Draw();

    c1->cd(4);
    hGen_stitched->SetTitle("hGen_weighted");
    hGen_stitched->GetXaxis()->SetRangeUser(0, 100);
    hGen_stitched->Draw();

    c1->cd(5);
    hReco_stitched->SetTitle("hReco_weighted");
    hReco_stitched->GetXaxis()->SetRangeUser(0, 100);
    hReco_stitched->Draw();

    c1->cd(6);
    hEff_stitched->SetTitle("hEff_weighted");
    hEff_stitched->SetMaximum(1.0);
    hEff_stitched->GetXaxis()->SetRangeUser(0, 100);
    hEff_stitched->Draw();

    // Draw

    results->cd();
    hReco_raw->Write();
    hGen_raw->Write();
    hEff_raw->Write();
    hGen_stitched_L->Write();
    hReco_stitched_L->Write();
    hEff_stitched_L->Write();
    hGen_stitched_S->Write();
    hReco_stitched_S->Write();
    hEff_stitched_S->Write();
    hGen_stitched->Write();
    hReco_stitched->Write();
    hEff_stitched->Write();
    c1->Write("eff_plots_summary");

    cout << "done!" << endl;
    cout << "pT_L efficiency = " <<  hReco_stitched_L->Integral() << " / " << hGen_stitched_L->Integral() << " = " <<  hReco_stitched_L->Integral() / hGen_stitched_L->Integral() << endl;
    cout << "pT_S efficiency = " <<  hReco_stitched_S->Integral() << " / " << hGen_stitched_S->Integral() << " = " <<  hReco_stitched_S->Integral() / hGen_stitched_S->Integral() << endl;
    results->Close();
    // c1->SaveAs("efficiency_D0_pthat0_kpi.pdf");
    // c1->SaveAs("efficiency_D0_pthat0_kpi.root");

    // std::cout << "Done. Efficiency plot saved as efficiency_D0_pt.png" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <istart> <iend> <output file>" << std::endl;
        return 1;
    }

    efficiency_analysis_diff(argv[1], argv[2], argv[3]);
    return 0;
}
