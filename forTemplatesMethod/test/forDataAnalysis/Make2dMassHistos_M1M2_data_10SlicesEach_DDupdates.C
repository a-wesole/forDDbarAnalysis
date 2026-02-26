#include <TFile.h>
#include <TH2F.h>
#include <vector>
#include <string>
#include "TFile.h"
#include "TNtuple.h"
#include <TTree.h>
#include "TTree.h"
#include <TBranch.h>

#include "BDTHandler.h"
#include "phiTools.h"

const bool isTestRun = false;
const int number_of_slices = 30;

const float nbinsX = 100;
const float nbinsY = 30;
const float fit_range_low = 1.75, fit_range_high = 1.86484 + 0.125;
// const float fit_range_low = 1.74, fit_range_high = 1.99;
const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / 5);
const float phi_array[6] = {in_val, in_val + sc_val, in_val + 2 * sc_val, in_val + 3 * sc_val, in_val + 4 * sc_val, in_val + 5 * sc_val};

const float ptMin = 6.0;
const float ptMin_2 = 3.0;
const int centMin = 0, centMax = 90;
const float yMin = 0.0, yMax = 1.0;
const float pT_minimum = 0.5;

bool fisD0candidate, fisDbarcandidate, gisD0candidate, gisDbarcandidate;

using namespace std;

template <typename T>
constexpr std::pair<int, int> determineLeadingAndSubleading(const T &pT_value, const int f, const int g)
{
    if (pT_value[0] >= pT_value[1])
    {
        return {f, g}; // f is leading, g is subleading
    }
    else
    {
        return {g, f}; // g is leading, f is subleading
    }
}

void Make2dMassHistos_M1M2_data_10SlicesEach_DDupdates(TString i_start = "", TString i_end = "")
{

    ifstream file_stream("/home/awesole/forDDbarAnalysis/forTemplatesMethod/dataFiles/Run_2023_19Sept.list");
    // ifstream file_stream("/home/awesole/forDDbarAnalysis/forTemplatesMethod/dataFiles/prime1.list");
    string filename;
    int ifile = 0;
    int istart = atoi(i_start);
    int iend = atoi(i_end);
    string output_name = "/scratch/negishi/awesole/DDbar_Jan2026/ROOT/M1M2_file_" + to_string(istart) + "_" + to_string(iend) + ".root";
    if(isTestRun) output_name = "oo.root";
    TFile *results = new TFile(output_name.c_str(), "recreate");
    cout << "output_name = " << output_name << endl;
    auto bdt = std::make_unique<BDTHandler>();

    // Prepare slices
    std::vector<TH2F *> M1_slices;
    std::vector<TH2F *> M2_slices;


    std::vector<TH2F *> M1M2_bins; // 6 bins
    std::vector<std::vector<TH2F *>> M1_bins_slices(6); // 6 bins
    std::vector<std::vector<TH2F *>> M2_bins_slices(6);
    std::vector<TH2F *> M1M2DDbar_bins; // 6 bins
    std::vector<std::vector<TH2F *>> M1DDbar_bins_slices(6); // 6 bins
    std::vector<std::vector<TH2F *>> M2DDbar_bins_slices(6);
    std::vector<TH2F *> M1M2DD_bins; // 6 bins
    std::vector<std::vector<TH2F *>> M1DD_bins_slices(6); // 6 bins
    std::vector<std::vector<TH2F *>> M2DD_bins_slices(6);
    std::vector<TH2F *> M1M2DbarDbar_bins; // 6 bins
    std::vector<std::vector<TH2F *>> M1DbarDbar_bins_slices(6); // 6 bins
    std::vector<std::vector<TH2F *>> M2DbarDbar_bins_slices(6);

    // prepare tuples
    TNtuple *M1M2Tuple = new TNtuple("M1M2Tuple", "M1M2Tuple", "m1:m2:phi1:phi2:pt1:pt2:eta1:eta2:DDbar:DD:DbarDbar");

    int bin = 0;
    const int max_size = 50000;
    float bdt_cut_value = -999.9;
    float delta_phi = -99;

    int Dsize, centrality;
    float mass[max_size], phi[max_size], flavor[max_size], Dtrk1Pt[max_size], Dtrk2Pt[max_size], mva[max_size], eta[max_size];
    float y[max_size], pT[max_size];
    std::vector<double> pT_vals = {};

    float slice_width = (fit_range_high - fit_range_low) / number_of_slices;

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

    TH2F *M1M2 = new TH2F("M1M2", "M1M2", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2->SetXTitle("M1 Mass");
    M1M2->SetYTitle("M2 Mass");
    M1M2->SetOption("SURF1");

    for (int i = 0; i < number_of_slices; ++i)
    { // calls the abovve function and creates the histograms for the slices
        std::string m1_name = "M1_slice" + std::to_string(i);
        std::string m2_name = "M2_slice" + std::to_string(i);

        M1_slices.push_back(createSlice(m1_name));
        M2_slices.push_back(createSlice(m2_name));
    }

    TH2F *M1M2_DDbar = new TH2F("M1M2_DDBar", "M1M2_DDBar", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2_DDbar->SetXTitle("M1 Mass");
    M1M2_DDbar->SetYTitle("M2 Mass");
    M1M2_DDbar->SetOption("SURF1");

    TH2F *M1M2_DD = new TH2F("M1M2_DD", "M1M2_DD", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2_DD->SetXTitle("M1 Mass");
    M1M2_DD->SetYTitle("M2 Mass");
    M1M2_DD->SetOption("SURF1");

    TH2F *M1M2_DbarDbar = new TH2F("M1M2_DbarDbar", "M1M2_DbarDbar", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2_DbarDbar->SetXTitle("M1 Mass");
    M1M2_DbarDbar->SetYTitle("M2 Mass");
    M1M2_DbarDbar->SetOption("SURF1");

    std::vector<TH2F *> M1M2Mass;
    for (int i = 1; i < 7; i++)
    {
        TString name1 = TString::Format("M1M2Mass_bin%d", i);
        TH2F *M1M2_bin = new TH2F(name1, name1, nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
        M1M2_bin->SetXTitle("M1 Mass");
        M1M2_bin->SetYTitle("M2 Mass");
        M1M2_bin->SetOption("lego2");
        M1M2_bins.push_back(M1M2_bin);

        TString name2 = TString::Format("M1M2MassDDbar_bin%d", i);
        TH2F *M1M2DDbar_bin = new TH2F(name2, name2, nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
        M1M2DDbar_bin->SetXTitle("M1 Mass");
        M1M2DDbar_bin->SetYTitle("M2 Mass");
        M1M2DDbar_bin->SetOption("lego2");
        M1M2DDbar_bins.push_back(M1M2DDbar_bin);

        TString name3 = TString::Format("M1M2MassDbarDbar_bin%d", i);
        TH2F *M1M2DbarDbar_bin = new TH2F(name3, name3, nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
        M1M2DbarDbar_bin->SetXTitle("M1 Mass");
        M1M2DbarDbar_bin->SetYTitle("M2 Mass");
        M1M2DbarDbar_bin->SetOption("lego2");
        M1M2DbarDbar_bins.push_back(M1M2DbarDbar_bin);

        TString name4 = TString::Format("M1M2MassDD_bin%d", i);
        TH2F *M1M2DD_bin = new TH2F(name4, name4, nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
        M1M2DD_bin->SetXTitle("M1 Mass");
        M1M2DD_bin->SetYTitle("M2 Mass");
        M1M2DD_bin->SetOption("lego2");
        M1M2DD_bins.push_back(M1M2DD_bin);


        for (int j = 0; j < number_of_slices; ++j)
        { // calls the abovve function and creates the histograms for the slices
            std::string m1_name = "M1_bin" + std::to_string(i) + "_slice" + std::to_string(j);
            std::string m2_name = "M2_bin" + std::to_string(i) + "_slice" + std::to_string(j);

            M1_bins_slices[i-1].push_back(createSlice(m1_name));
            M2_bins_slices[i-1].push_back(createSlice(m2_name));
        }
    }

    // fill the slices
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

        for (int i = 0; i < t->GetEntries(); i++)
        {
            // if (i != 42805) continue;
            t->GetEntry(i);
            if (Dsize <= 1)
                continue;
            if (centrality < 2 * centMin || centrality >= 2 * centMax)
                continue;

            for (int f = 0; f < Dsize - 1; f++)
            { // for first kpi pair

                fisD0candidate = fisDbarcandidate = false;

                if (Dtrk1Pt[f] < pT_minimum || Dtrk2Pt[f] < pT_minimum)
                    continue; // apply pT cuts tp each daughter (iparticle and j particle)

                if (abs(y[f]) < yMin || abs(y[f]) >= yMax)
                    continue;

                // determine D0 or Dbar candidate
                if (flavor[f] == 1)
                    fisD0candidate = true;
                if (flavor[f] == -1)
                    fisDbarcandidate = true;

                bdt_cut_value = bdt->getBDTCut(y[f], centrality, pT[f]);
                if (mva[f] <= bdt_cut_value)
                    continue; // background cuts

                // if ((mass[f] > 1.74 && mass[f] < 1.825) || (mass[f] > 1.91 && mass[f] < 1.99))
                if (mass[f] > fit_range_low && mass[f] < fit_range_high)
                {
                    for (int g = f + 1; g < Dsize; g++)
                    {
                        gisD0candidate = gisDbarcandidate = false;
                        bin = 0;

                        if (Dtrk1Pt[g] < pT_minimum || Dtrk2Pt[g] < pT_minimum)
                            continue; // apply pT cuts to each daughter (iparticle & jparticle)
                        if (abs(y[g]) < yMin || abs(y[g]) >= yMax)
                            continue;

                        bdt_cut_value = bdt->getBDTCut(y[g], centrality, pT[g]);
                        if (mva[g] <= bdt_cut_value)
                            continue; // background cuts

                        // float DDbar_phi = phi[f] - phi[g]; // store the dphi if hte kpi pair are from D0candidate then D0bar candidate
                        // float DbarD_phi = phi[g] - phi[f]; // store the dphi if hte kpi pair are from D0 bar candidate then D0candidate

                        // determine D0 and Dbar candidate:
                        if ((flavor[g] == 1))
                            gisD0candidate = true;
                        if ((flavor[g] == -1))
                            gisDbarcandidate = true;
                        // if ((fisD0candidate && gisD0candidate) || (fisDbarcandidate && gisDbarcandidate)) continue;

                        // if ((mass[g] > 1.74 && mass[g] < 1.825) || (mass[g] > 1.91 && mass[g] < 1.99))
                        if (mass[g] > fit_range_low && mass[g] < fit_range_high)
                        {
                            if (std::max(pT[f], pT[g]) >= ptMin && std::min(pT[f], pT[g]) >= ptMin_2)
                            {
                                pT_vals.clear();
                                pT_vals.push_back(pT[f]);
                                pT_vals.push_back(pT[g]);
                                auto [leading, subleading] = determineLeadingAndSubleading(pT_vals, f, g);
                                if (abs(pT[leading] - pT[subleading]) < 1.0e-5 && abs(eta[leading] - eta[subleading]) < 1.0e-5 && abs(phi[leading] - phi[subleading]) < 1.0e-5)
                                {
                                    continue; // skip for candidiate and its own swap
                                }

                                // if (!((abs(Dtrk1Pt[leading] - Dtrk1Pt[subleading]) < 1.0e-5) ||
                                //     (abs(Dtrk1Pt[leading] - Dtrk2Pt[subleading]) < 1.0e-5) ||
                                //     (abs(Dtrk2Pt[leading] - Dtrk1Pt[subleading]) < 1.0e-5) ||
                                //     (abs(Dtrk2Pt[leading] - Dtrk2Pt[subleading]) < 1.0e-5)))
                                // {
                                //     // cout << "skipping!!" << endl;
                                //     continue; // skip 2 candidiates that have a common daughter track.
                                // }

                                // if (abs(phi[leading] - phi[subleading] < 0.05)){
                                //     cout << "input file is = " << filename.c_str() << " and entry is " << i << endl;
                                //     cout << "mass 1 = " << mass[leading] << " mass 2 = " << mass[subleading] << endl;
                                //     cout << "fisD0candidate = " << fisD0candidate << " fisDbarcandidate = " << fisDbarcandidate << endl;
                                //     cout << "gisD0candidate = " << gisD0candidate << " gisDbarcandidate = " << gisDbarcandidate << endl;
                                // }
                                delta_phi = phi[leading] - phi[subleading];
                                transition_phi(delta_phi);

                                if (bin != 0)
                                    cout << "Error!!!!!!!  phi = " << delta_phi << " and bin =" << bin << endl;
                                if (delta_phi >= phi_array[0] && delta_phi < phi_array[1])
                                    bin = 1;
                                if (delta_phi >= phi_array[1] && delta_phi < phi_array[2])
                                    bin = 2;
                                if (delta_phi >= phi_array[2] && delta_phi < phi_array[3])
                                    bin = 3;
                                if (delta_phi >= phi_array[3] && delta_phi < phi_array[4])
                                    bin = 4;
                                if (delta_phi >= phi_array[4] && delta_phi <= phi_array[5])
                                    bin = 5;

                                // cout << "~~~~~~~ bin = " << bin << " ~~~~~~~" << endl;

                                int sliceX = std::min(int(number_of_slices), int((mass[leading] - fit_range_low) / slice_width));
                                int sliceY = std::min(int(number_of_slices), int((mass[subleading] - fit_range_low) / slice_width));
                                // cout << "slice X = " << sliceX << " sliceY = " << sliceY << endl;

                                ///first fill the overall M1M2 (w/o regard to dphi bin)
                                M1M2->Fill(mass[leading], mass[subleading]);
                                M1_slices[sliceY]->Fill(mass[leading], mass[subleading]); // Fill the corresponding slice with m1,m2
                                M2_slices[sliceX]->Fill(mass[leading], mass[subleading]); // Fill the corresponding slice with m1,m2

                                // cout << "filled overall " << endl;

                                //now fill the M1M2 wrt to apprpriate dphi bin
                                M1M2_bins[bin-1]->Fill(mass[leading], mass[subleading]);
                                if ((fisD0candidate && gisDbarcandidate) || (gisD0candidate && fisDbarcandidate))
                                {
                                    M1M2DDbar_bins[bin-1]->Fill(mass[leading], mass[subleading]);
                                }
                                else if (fisD0candidate && gisD0candidate)
                                {
                                    M1M2DD_bins[bin-1]->Fill(mass[leading], mass[subleading]);
                                }
                                else if (fisDbarcandidate && gisDbarcandidate)
                                {
                                    M1M2DbarDbar_bins[bin-1]->Fill(mass[leading], mass[subleading]);
                                }
                                // cout << "filled m1m2 bin" << endl;
                                M1_bins_slices[bin-1][sliceY]->Fill(mass[leading], mass[subleading]); // Fill the corresponding slice with m1,m2
                                M2_bins_slices[bin-1][sliceX]->Fill(mass[leading], mass[subleading]);

                                // cout << "filled for bins" << endl;


                                if ((fisD0candidate && gisDbarcandidate) || (gisD0candidate && fisDbarcandidate))
                                {
                                    M1M2_DDbar->Fill(mass[leading], mass[subleading]);
                                    M1M2Tuple->Fill(mass[leading], mass[subleading], phi[leading], phi[subleading], pT[leading], pT[subleading], eta[leading], eta[subleading], 1, 0, 0);
                                }
                                else if (fisD0candidate && gisD0candidate)
                                {
                                    M1M2_DD->Fill(mass[leading], mass[subleading]);
                                    M1M2Tuple->Fill(mass[leading], mass[subleading], phi[leading], phi[subleading], pT[leading], pT[subleading], eta[leading], eta[subleading], 0, 1, 0);
                                }
                                else if (fisDbarcandidate && gisDbarcandidate)
                                {
                                    M1M2_DbarDbar->Fill(mass[leading], mass[subleading]);
                                    M1M2Tuple->Fill(mass[leading], mass[subleading], phi[leading], phi[subleading], pT[leading], pT[subleading], eta[leading], eta[subleading], 0, 0, 1);
                                }
                            } // if one of the candidates has pT > max
                        } // if mass g
                    } // for g
                } // if mass f
            } // for f
        } // ttree entries
        ifile++;
    } // while

    results->cd();
    for (auto slice : M1_slices){slice->Write();}
    for (auto slice : M2_slices){slice->Write();}
    for (auto bin : M1M2_bins){ bin->Write();}
    for (auto bin : M1M2DDbar_bins){ bin->Write();}
    for (auto bin : M1M2DD_bins){ bin->Write();}
    for (auto bin : M1M2DbarDbar_bins){ bin->Write();}

    for (const auto &binVec : M1_bins_slices)
    {for (TH2F *hist : binVec){hist->Write();}}

    for (const auto &binVec : M2_bins_slices)
    {for (TH2F *hist : binVec){hist->Write();}}


    M1M2->Write();
    M1M2_DDbar->Write();
    M1M2_DD->Write();
    M1M2_DbarDbar->Write();
    M1M2Tuple->Write();

    results->Close();
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <istart> <iend>" << std::endl;
        return 1;
    }

    Make2dMassHistos_M1M2_data_10SlicesEach_DDupdates(argv[1], argv[2]);
    return 0;
}
