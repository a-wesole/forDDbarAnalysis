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
const int number_of_slices = 50;

const int number_of_phi_bins = 10;

const float nbinsX = 50;
const float nbinsY = 50;


const float fit_range_low = 1.86484 - 0.125, fit_range_high = 1.86484 + 0.125;
// const float fit_range_low = 1.74, fit_range_high = 1.99;
const float in_val = -0.5 * TMath::Pi(), sc_val = (2 * TMath::Pi() / number_of_phi_bins);
const float max_val = 2 * TMath::Pi() + in_val;

const float ptMin = 3.0;
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
    // string output_name = "/scratch/negishi/awesole/DDbar_July_2026_pT_2_1/ROOT/M1M2_file_" + to_string(istart) + "_" + to_string(iend) + ".root";
    string output_name = "/scratch/negishi/awesole/DDbar_July_2026_pT_" 
                   + to_string(static_cast<int>(ptMin)) + "_" 
                   + to_string(static_cast<int>(ptMin_2)) 
                   + "/ROOT/M1M2_skippingSharedTracks_" + to_string(istart) + "_" + to_string(iend) + ".root";
    if(isTestRun) output_name = "oo.root";
    TFile *results = new TFile(output_name.c_str(), "recreate");
    cout << "output_name = " << output_name << endl;
    auto bdt = std::make_unique<BDTHandler>();

    // Prepare slices
    std::vector<TH2F *> M1_slices;
    std::vector<TH2F *> M2_slices;


    std::vector<TH2F *> M1M2_bins; // 6 bins
    std::vector<TH2F *> M1M2DDbar_bins; // 6 bins
    std::vector<TH2F *> M1M2DD_bins; // 6 bins
    std::vector<TH2F *> M1M2DbarDbar_bins; // 6 bins

    // prepare tuples
    TNtuple *M1M2Tuple = new TNtuple("M1M2Tuple", "M1M2Tuple", "m1:m2:m1_slice:phi1:phi2:dphi:dphi_bin:pt1:pt2:eta1:eta2:DDbar:DD:DbarDbar");

    int bin = 0;
    const int max_size = 50000;
    float bdt_cut_value = -999.9;
    float delta_phi = -99;

    int Dsize, centrality;
    float mass[max_size], phi[max_size], flavor[max_size], Dtrk1Pt[max_size], Dtrk2Pt[max_size], mva[max_size], eta[max_size];
    float y[max_size], pT[max_size];
    std::vector<double> pT_vals = {};

    float slice_width = (fit_range_high - fit_range_low) / number_of_slices;


    TH2F *M1M2 = new TH2F("M1M2", "M1M2", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2->SetXTitle("M1 Mass");
    M1M2->SetYTitle("M2 Mass");
    M1M2->SetOption("lego2");
    M1M2->SetMinimum(0);


    TH2F *M1M2_DDbar = new TH2F("M1M2_DDBar", "M1M2_DDBar", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2_DDbar->SetXTitle("M1 Mass");
    M1M2_DDbar->SetYTitle("M2 Mass");
    M1M2_DDbar->SetOption("lego2");
    M1M2_DDbar->SetMinimum(0);

    TH2F *M1M2_DD = new TH2F("M1M2_DD", "M1M2_DD", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2_DD->SetXTitle("M1 Mass");
    M1M2_DD->SetYTitle("M2 Mass");
    M1M2_DD->SetOption("lego2");
    M1M2_DD->SetMinimum(0);

    TH2F *M1M2_DbarDbar = new TH2F("M1M2_DbarDbar", "M1M2_DbarDbar", nbinsX, fit_range_low, fit_range_high, nbinsY, fit_range_low, fit_range_high);
    M1M2_DbarDbar->SetXTitle("M1 Mass");
    M1M2_DbarDbar->SetYTitle("M2 Mass");
    M1M2_DbarDbar->SetOption("lego2");
    M1M2_DbarDbar->SetMinimum(0);

    std::vector<TH2F *> M1M2Mass;
    for (int i = 1; i < number_of_phi_bins + 2; i++)
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
                if (mass[f] >= fit_range_low && mass[f] <= fit_range_high)
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
                        if (mass[g] >= fit_range_low && mass[g] <= fit_range_high)
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

                                if (((abs(Dtrk1Pt[leading] - Dtrk1Pt[subleading]) < 1.0e-5) ||
                                    (abs(Dtrk1Pt[leading] - Dtrk2Pt[subleading]) < 1.0e-5) ||
                                    (abs(Dtrk2Pt[leading] - Dtrk1Pt[subleading]) < 1.0e-5) ||
                                    (abs(Dtrk2Pt[leading] - Dtrk2Pt[subleading]) < 1.0e-5)))
                                {
                                    // cout << "skipping!!" << endl;
                                    continue; // skip 2 candidiates that have a common daughter track.
                                    //July 2026 we need to remove this candidates as they make a large spike at dphi ~=0 
                                }

                                delta_phi = phi[leading] - phi[subleading];
                                transition_phi(delta_phi);

                                if (delta_phi == max_val) { bin = number_of_phi_bins;}
                                else {bin = static_cast<int>((delta_phi - in_val) / sc_val) + 1;}


                                int sliceX = std::min(int(number_of_slices), int((mass[leading] - fit_range_low) / slice_width));
                                int sliceY = std::min(int(number_of_slices), int((mass[subleading] - fit_range_low) / slice_width));

                                ///first fill the overall M1M2 (w/o regard to dphi bin)
                                M1M2->Fill(mass[leading], mass[subleading]);

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

                                // cout << "filled for bins" << endl;


                                if ((fisD0candidate && gisDbarcandidate) || (gisD0candidate && fisDbarcandidate))
                                {
                                    M1M2_DDbar->Fill(mass[leading], mass[subleading]);
                                    M1M2Tuple->Fill(mass[leading], mass[subleading], sliceY, phi[leading], phi[subleading], delta_phi, bin, pT[leading], pT[subleading], eta[leading], eta[subleading], 1, 0, 0);
                                }
                                else if (fisD0candidate && gisD0candidate)
                                {
                                    M1M2_DD->Fill(mass[leading], mass[subleading]);
                                    M1M2Tuple->Fill(mass[leading], mass[subleading], sliceY, phi[leading], phi[subleading], delta_phi, bin, pT[leading], pT[subleading], eta[leading], eta[subleading], 0, 1, 0);
                                }
                                else if (fisDbarcandidate && gisDbarcandidate)
                                {
                                    M1M2_DbarDbar->Fill(mass[leading], mass[subleading]);
                                    M1M2Tuple->Fill(mass[leading], mass[subleading],sliceY,  phi[leading], phi[subleading], delta_phi, bin, pT[leading], pT[subleading], eta[leading], eta[subleading], 0, 0, 1);
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
    // for (auto slice : M1_slices){slice->Write();}
    // for (auto slice : M2_slices){slice->Write();}
    for (auto bin : M1M2_bins){ bin->Write();}
    for (auto bin : M1M2DDbar_bins){ bin->Write();}
    for (auto bin : M1M2DD_bins){ bin->Write();}
    for (auto bin : M1M2DbarDbar_bins){ bin->Write();}


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
