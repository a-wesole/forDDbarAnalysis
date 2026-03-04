// this code performs event mixing and calculates the phi and xT
// test file:     root://eos.cms.rcac.purdue.edu//store/user/junseok/D0_2023PbPb/HIPhysicsRawPrime31/HIPhysicsRawPrime31_D0_05Aug2025/250821_071550/0004/TTree_4569.root

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
#include <deque>
#include <fstream>

#include "phiTools.h"
#include "BDTHandler.h"

// ---config ---
const bool isTestRun = false;
const int number_of_slices = 30;

const int k_Num_PvtxZ_Bins = 40;
const int k_Num_Buffer_Size = 1; // number of previous events kept and calculate dphi with
const float PvtxZ_min = -15.0, PvtxZ_max = -1 * PvtxZ_min;

const float pT_minimum = 0.5;
const int max_size = 50000;

const float ptMin = 6.0;
const float ptMin_2 = 3.0;
const int centMin = 0, centMax = 90;
const float yMin = 0.0, yMax = 1.0;

const float nbinsX = 100;
const float nbinsY = 30;
const float fit_range_low = 1.75, fit_range_high = 1.86484 + 0.125;
const float slice_width = (fit_range_high - fit_range_low) / number_of_slices;
const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[6] = {in_val, in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val, in_val + 5 * sc_val};

//  -- lightweight storage of D0 cands to analyze--
struct DCand
{
    float phi;
    float pT;
    float eta;
    float mass;
    float y;
    float flavor; // 1 for d0 and -1 for dbar
    float Dtrk1Pt;
    float Dtrk2Pt;
};

struct EventPool
{
    std::vector<DCand> d;    // point is we keep only the d0 candidates for mixing in a vector and the Dcands are very light weight with only the variables above.
    std::vector<DCand> dbar; // point is we keep only the dbar candidates for mixing in a vector and the Dcands are very light weight with only the variables above.
};

// -- function to determine PvtxZ bin--
int getPvtxZBin(float vtx_z)
{
    if (vtx_z < PvtxZ_min || vtx_z >= PvtxZ_max)
        return -1;
    float width = (PvtxZ_max - PvtxZ_min) / k_Num_PvtxZ_Bins;
    int bin = static_cast<int>((vtx_z - PvtxZ_min) / width);
    if (bin < 0)
        bin = 0;
    if (bin >= k_Num_PvtxZ_Bins)
        bin = k_Num_PvtxZ_Bins - 1;
    return bin;
}

template <typename Cand>
constexpr std::pair<const Cand &, const Cand &> determineLeadingAndSubleading(const Cand &a, const Cand &b)
{
    if (a.pT >= b.pT)
    {
        return {a, b}; // a is leading, b is subleading
    }
    else
    {
        return {b, a}; // b is leading, a is subleading
    }
}

template <typename CandContainerA, typename CandContainerB>
void thatFunction(const CandContainerA &A,
                  const CandContainerB &B,
                  int bin,
                  double ptMin,
                  double ptMin_2,
                  TH1D *h_single_pT,
                  TH1D *h_single_phi_leading,
                  TH1D *h_single_phi_subleading,
                  TH1D *h_dphi_mixedEvent_bin,
                  TH1D *h_dphi_mixing_total,
                  TH2F *M1M2,
                  std::vector<TH2F *> &M1M2_bins,
                  std::vector<std::vector<TH2F *>> &M1_bins_slices,
                  std::vector<std::vector<TH2F *>> &M2_bins_slices,
                  TH1D *h_pT_leading_mixedEvt, TH1D *h_pT_subleading_mixedEvt, TH1D *h_pT_leading_sameEvt, TH1D *h_pT_subleading_sameEvt, TNtuple *tup)
{
    for (const auto &a : A)
    {
        for (const auto &b : B)
        {
            if (std::max(a.pT, b.pT) >= ptMin && std::min(a.pT, b.pT) >= ptMin_2)
            {
                auto [leading, subleading] = determineLeadingAndSubleading(a, b);

                if (abs(leading.pT - subleading.pT) < 1.0e-5 && abs(leading.eta - subleading.eta) < 1.0e-5 && abs(leading.phi - subleading.phi) < 1.0e-5)
                    continue;

                h_single_phi_leading->Fill(leading.phi);
                h_single_phi_subleading->Fill(subleading.phi);
                float raw_phi = leading.phi - subleading.phi;
                double dphi = transition_phi(raw_phi);
                h_dphi_mixedEvent_bin->Fill(dphi);
                h_dphi_mixing_total->Fill(dphi);

                //--get dphi bin [1-5]
                int dphi_bin = 0;
                if (dphi_bin != 0)
                    cout << "Error!!!!!!!  phi = " << dphi << " and dphi_bin =" << dphi_bin << endl;
                if (dphi >= phi_array[0] && dphi < phi_array[1])
                    dphi_bin = 1;
                if (dphi >= phi_array[1] && dphi < phi_array[2])
                    dphi_bin = 2;
                if (dphi >= phi_array[2] && dphi < phi_array[3])
                    dphi_bin = 3;
                if (dphi >= phi_array[3] && dphi < phi_array[4])
                    dphi_bin = 4;
                if (dphi >= phi_array[4] && dphi <= phi_array[5])
                    dphi_bin = 5;

                int sliceX = std::min(int(number_of_slices), int((leading.mass - fit_range_low) / slice_width));
                int sliceY = std::min(int(number_of_slices), int((subleading.mass - fit_range_low) / slice_width));

                M1M2->Fill(leading.mass, subleading.mass);
                M1_bins_slices[dphi_bin - 1][sliceY]->Fill(leading.mass, subleading.mass);
                M2_bins_slices[dphi_bin - 1][sliceX]->Fill(leading.mass, subleading.mass);
                M1M2_bins[dphi_bin - 1]->Fill(leading.mass, subleading.mass);
                h_pT_leading_mixedEvt->Fill(leading.pT); h_pT_subleading_mixedEvt->Fill(subleading.pT);
                // h_pT_leading_sameEvt->Fill(leading.pT); h_pT_subleading_sameEvt->Fill(subleading.pT);
                tup->Fill(leading.pT, subleading.pT);
                h_single_pT->Fill(leading.pT);
                h_single_pT->Fill(subleading.pT);
            }
        }
    }
}

void first_macro(TString i_start = "", TString i_end = "")
{
    int Dsize, centrality;
    float PvtxZ;
    float mass[max_size], phi[max_size], flavor[max_size], Dtrk1Pt[max_size], Dtrk2Pt[max_size], mva[max_size], eta[max_size];
    float y[max_size], pT[max_size];

    ifstream file_stream("/home/awesole/forDDbarAnalysis/forTemplatesMethod/dataFiles/Run_2023_19Sept.list");
    string filename;
    int ifile = 0;
    int istart = atoi(i_start);
    int iend = atoi(i_end);
    string output_name = "/scratch/negishi/awesole/DDbar_eventmixing/ROOT/M1M2_file_" + to_string(istart) + "_" + to_string(iend) + ".root";
    if (isTestRun)
        output_name = "oo.root";
    TFile *results = new TFile(output_name.c_str(), "recreate");
    // --histograms--
    // TH1D* h_dphi_sameEvent[k_Num_PvtxZ_Bins];
    TH1D *h_dphi_mixedEvent[k_Num_PvtxZ_Bins];
    for (int b = 0; b < k_Num_PvtxZ_Bins; b++)
    {
        // h_dphi_sameEvent[b] = new TH1D(Form("h_dphi_sameEvent_PvtxZbin_%d", b), Form("h_dphi_sameEvent_PvtxZbin_%d", b), 5, TMath::Pi()/2, 3*TMath::Pi()/2);
        h_dphi_mixedEvent[b] = new TH1D(Form("h_dphi_mixedEvent_PvtxZbin_%d", b), Form("h_dphi_mixedEvent_PvtxZbin_%d", b), 5, -1 * TMath::Pi() / 2, 3 * TMath::Pi() / 2);
        h_dphi_mixedEvent[b]->SetMinimum(0);
    }

    TH1D *h_dphi_mixing_total = new TH1D("h_dphi_mixing_total", "h_dphi_mixing_total", 5, -1 * TMath::Pi() / 2, 3 * TMath::Pi() / 2);
    h_dphi_mixing_total->SetMinimum(0);
    TH1D *h_single_pT = new TH1D("h_single_pT", "h_single_pT", 40, 0, 50);
    h_single_pT->SetMinimum(0);
    TH1D *h_single_phi_leading = new TH1D("h_single_phi_leading", "h_single_phi_leading", 20, -1 * TMath::Pi(), TMath::Pi());
    h_single_phi_leading->SetMinimum(0);
    TH1D *h_single_phi_subleading = new TH1D("h_single_phi_subleading", "h_single_phi_subleading", 20, -1 * TMath::Pi(), TMath::Pi());
    h_single_phi_subleading->SetMinimum(0);
    TH1D *h_pT_subleading_mixedEvt = new TH1D("h_pT_subleading_mixedEvt", "h_pT_subleading_mixedEvt", 40, 0, 50);
    h_pT_subleading_mixedEvt->SetMinimum(0);
    TH1D *h_pT_leading_mixedEvt = new TH1D("h_pT_leading_mixedEvt", "h_pT_leading_mixedEvt", 40, 0, 50);
    h_pT_leading_mixedEvt->SetMinimum(0);
    TH1D *h_pT_subleading_sameEvt = new TH1D("h_pT_subleading_sameEvt", "h_pT_subleading_sameEvt", 40, 0, 50);
    h_pT_subleading_sameEvt->SetMinimum(0);
    TH1D *h_pT_leading_sameEvt = new TH1D("h_pT_leading_sameEvt", "h_pT_leading_sameEvt", 40, 0, 50);
    h_pT_leading_sameEvt->SetMinimum(0);
    TNtuple *tup = new TNtuple("tup", "tup", "pT_lead:pT_sublead");

    auto createSlice = [&](const std::string &name) -> TH2F * // creates a TH2F with the same binning as M1M2, for each slice this will be filled later
    {
        TH2F *slice = new TH2F(name.c_str(), name.c_str(),
                               nbinsX,
                               fit_range_low,
                               fit_range_high,
                               nbinsY,
                               fit_range_low,
                               fit_range_high);
        slice->SetOption("surf1");
        return slice;
    };

    std::vector<TH2F *> M1M2_bins; // 6 bins
    std::vector<TH2F *> M1_slices;
    std::vector<TH2F *> M2_slices;
    for (int i = 0; i < number_of_slices; ++i)
    { // calls the abovve function and creates the histograms for the slices
        std::string m1_name = "M1_slice" + std::to_string(i);
        std::string m2_name = "M2_slice" + std::to_string(i);

        M1_slices.push_back(createSlice(m1_name));
        M2_slices.push_back(createSlice(m2_name));
    }

    std::vector<std::vector<TH2F *>> M1_bins_slices(6); // 6 bins
    std::vector<std::vector<TH2F *>> M2_bins_slices(6);

    TH2F *M1M2 = new TH2F("M1M2", "M1M2", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2->SetXTitle("M1 Mass");
    M1M2->SetYTitle("M2 Mass");
    M1M2->SetOption("SURF1");

    std::vector<TH2F *> M1M2Mass;
    for (int i = 1; i < 7; i++)
    {
        TString name1 = TString::Format("M1M2Mass_bin%d", i);
        TH2F *M1M2_bin = new TH2F(name1, name1, nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
        M1M2_bin->SetXTitle("M1 Mass");
        M1M2_bin->SetYTitle("M2 Mass");
        M1M2_bin->SetOption("lego2");
        M1M2_bins.push_back(M1M2_bin);

        for (int j = 0; j < number_of_slices; ++j)
        { // calls the abovve function and creates the histograms for the slices
            std::string m1_name = "M1_bin" + std::to_string(i) + "_slice" + std::to_string(j);
            std::string m2_name = "M2_bin" + std::to_string(i) + "_slice" + std::to_string(j);

            M1_bins_slices[i - 1].push_back(createSlice(m1_name));
            M2_bins_slices[i - 1].push_back(createSlice(m2_name));
        }
    }

    // --per-bin rolling buffers of prior events--
    std::deque<EventPool> buffers[k_Num_PvtxZ_Bins];

    auto bdt = std::make_unique<BDTHandler>();

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

        // cout << "file = " << filename.c_str() << endl;
        TFile *fin = TFile::Open(filename.c_str());
        if (fin->IsZombie() || fin->Get("d0Analyzer/VertexCompositeNtuple") == nullptr)
        {
            cout << "BIG PROBLEM" << endl;
            ifile++;
            continue;
        }

        TTree *t = (TTree *)fin->Get("d0Analyzer/VertexCompositeNtuple");

        t->SetBranchAddress("PvtxZ", &PvtxZ);
        t->SetBranchAddress("candSize", &Dsize);
        t->SetBranchAddress("centrality", &centrality);
        t->SetBranchAddress("mass", mass);
        t->SetBranchAddress("y", y);
        t->SetBranchAddress("eta", eta);
        t->SetBranchAddress("pT", pT);
        t->SetBranchAddress("phi", phi);
        t->SetBranchAddress("flavor", flavor);
        t->SetBranchAddress("mva", mva);
        t->SetBranchAddress("pTD1", Dtrk1Pt);
        t->SetBranchAddress("pTD2", Dtrk2Pt);

        if (ifile % 10 == 0)
            cout << "ifile = " << ifile << " ----- tree entries = " << t->GetEntries() << endl;

        // for each event
        for (int i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);
            // if (Dsize <= 1)
            if (Dsize < 1)
                continue;
            if (centrality < 2 * centMin || centrality >= 2 * centMax)
                continue;
            int bin = getPvtxZBin(PvtxZ);
            if (bin < 0)
                continue; // skip all Pvtx outside of the boundary

            //--set up for the buffer all the candidates in the events--
            std::vector<DCand> dCands_inEvent;
            std::vector<DCand> dbarCands_inEvent;
            dCands_inEvent.reserve(Dsize);
            dbarCands_inEvent.reserve(Dsize);

            // -- keep all D/Dbar candidates -- //

            for (int f = 0; f < Dsize; f++)
            {
                if (Dtrk1Pt[f] < pT_minimum || Dtrk2Pt[f] < pT_minimum)
                    continue; // apply pT cuts tp each daughter (iparticle and j particle)

                if (abs(y[f]) < yMin || abs(y[f]) >= yMax)
                    continue;

                double bdt_cut_value = bdt->getBDTCut(y[f], centrality, pT[f]);
                if (mva[f] <= bdt_cut_value)
                    continue; // background cuts
                if (pT[f] < ptMin_2)
                    continue;

                if (mass[f] <= fit_range_low || mass[f] >= fit_range_high)
                    continue;
                // determine D0 or Dbar candidate

                DCand c{phi[f], pT[f], eta[f], mass[f], y[f], static_cast<int>(flavor[f]), Dtrk1Pt[f], Dtrk2Pt[f]};
                // h_single_phi->Fill(phi[f]);
                if (flavor[f] == 1)
                    dCands_inEvent.push_back(c);
                else if (flavor[f] == -1)
                    dbarCands_inEvent.push_back(c);
            } // for first kpi pair
            if (dCands_inEvent.empty() && dbarCands_inEvent.empty())
                continue;

            //--mixxing ---
            // mix current d0 with all previous d0s in the same PvtxZ bin and all current d0bars with all previous d0bars in same pvtxZ bin
            for (const auto &pool : buffers[bin])
            {
                thatFunction(
                    dCands_inEvent,
                    pool.dbar,
                    // dbarCands_inEvent,
                    bin,
                    ptMin,
                    ptMin_2,
                    h_single_pT,
                    h_single_phi_leading,
                    h_single_phi_subleading,
                    h_dphi_mixedEvent[bin],
                    h_dphi_mixing_total,
                    M1M2,
                    M1M2_bins,
                    M1_bins_slices,
                    M2_bins_slices,
                    h_pT_leading_mixedEvt, h_pT_subleading_mixedEvt, h_pT_leading_sameEvt, h_pT_subleading_sameEvt, tup);
                thatFunction(
                    dbarCands_inEvent,
                    pool.d,
                    // dCands_inEvent,
                    bin,
                    ptMin,
                    ptMin_2,
                    h_single_pT,
                    h_single_phi_leading,
                    h_single_phi_subleading,
                    h_dphi_mixedEvent[bin],
                    h_dphi_mixing_total,
                    M1M2,
                    M1M2_bins,
                    M1_bins_slices,
                    M2_bins_slices,
                    h_pT_leading_mixedEvt, h_pT_subleading_mixedEvt, h_pT_leading_sameEvt, h_pT_subleading_sameEvt, tup);
            }
            EventPool currentPool;
            currentPool.d = std::move(dCands_inEvent);
            currentPool.dbar = std::move(dbarCands_inEvent);
            if (!currentPool.dbar.empty() || !currentPool.d.empty())
            {                                                   // only events with at least 1 d0 or dbar
                buffers[bin].push_back(std::move(currentPool)); // adds the current event to the end of the buffer for thiis PvtxZ Bin, we add after the mxing to make sure we do not mix events together
                if (static_cast<int>(buffers[bin].size()) > k_Num_Buffer_Size)
                {                             // checks if the buffer has grown larger than the allowed number of stored past events k_num_beffer_size
                    buffers[bin].pop_front(); // if the buffer is now too large it removes the oldest event bc newest evetns get added at the back
                }
            }

        } // tree entries
        ifile++;
    } // while
    results->cd();
    for (int b = 0; b < k_Num_PvtxZ_Bins; b++)
    {
        h_dphi_mixedEvent[b]->Write();
    }
    for (const auto &binVec : M1_bins_slices)
    {
        for (TH2F *hist : binVec)
        {
            hist->Write();
        }
    }

    for (const auto &binVec : M2_bins_slices)
    {
        for (TH2F *hist : binVec)
        {
            hist->Write();
        }
    }
    for (auto slice : M1_slices)
    {
        slice->Write();
    }
    for (auto slice : M2_slices)
    {
        slice->Write();
    }
    for (auto bin : M1M2_bins)
    {
        bin->Write();
    }
    h_dphi_mixing_total->Write();
    h_single_pT->Write();
    h_single_phi_leading->Write();
    h_single_phi_subleading->Write();
    h_pT_leading_mixedEvt->Write();
    h_pT_subleading_mixedEvt->Write();
    h_pT_leading_sameEvt->Write();
    h_pT_subleading_sameEvt->Write();
    tup->Write();

    M1M2->Write();
    results->Close();

} // first_macro
int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <istart> <iend>" << std::endl;
        return 1;
    }

    first_macro(argv[1], argv[2]);
    return 0;
}
