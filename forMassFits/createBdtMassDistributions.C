//GENERATES HISTOGRAMS (MASS DISTRIBUTIONS) FOR ALL PT AND CENTRALITY RANGES
//all histograms will be written into one file
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <TMath.h>
#include <iostream>
# include <TF1.h>
# include <TH1F.h>
# include <TH1D.h>
# include <TCanvas.h>
#include <fstream>

using namespace std;
Bool_t istest = false;
int createBdtMassDistributions(TString i_start="", TString i_end="")
{


    ifstream file_stream("/home/awesole/forDDbarAnalysis/dataFiles/dataWithBdt.list");
    string filename;
    int ifile = 0;
    int istart = atoi(i_start);
    int iend = atoi(i_end);

    string output_name = "/scratch/negishi/awesole/BdtHistos/ROOT/outfile_histos_" + to_string(istart) + "_" + to_string(iend) + ".root";
    TFile *outf = new TFile(output_name.c_str(), "recreate");

    int candSize, centrality;
    float pT[10000], mass[10000], y[10000];
    double BDT_weight[10000];

    std::vector<double> pt_values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0};
    std::vector<int> cent_values = {0, 10, 30, 50, 90};
    std::vector<float> bdt_values = {-0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.00, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};
    std::vector<float> bdt_val = {-0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6};


    //aalist = centrlaity 0-10 pt 1-2
    std::vector<TH1D*> histograms_pt_1_2_cent_0_10;
    histograms_pt_1_2_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_1_2_cent_10_30;
    histograms_pt_1_2_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_1_2_cent_30_50;
    histograms_pt_1_2_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_1_2_cent_50_90;
    histograms_pt_1_2_cent_50_90.resize(bdt_values.size());

    std::vector<TH1D*> histograms_pt_2_3_cent_0_10;
    histograms_pt_2_3_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_2_3_cent_10_30;
    histograms_pt_2_3_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_2_3_cent_30_50;
    histograms_pt_2_3_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_2_3_cent_50_90;
    histograms_pt_2_3_cent_50_90.resize(bdt_values.size());

    std::vector<TH1D*> histograms_pt_3_4_cent_0_10;
    histograms_pt_3_4_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_3_4_cent_10_30;
    histograms_pt_3_4_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_3_4_cent_30_50;
    histograms_pt_3_4_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_3_4_cent_50_90;
    histograms_pt_3_4_cent_50_90.resize(bdt_values.size());

    std::vector<TH1D*> histograms_pt_4_5_cent_0_10;
    histograms_pt_4_5_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_4_5_cent_10_30;
    histograms_pt_4_5_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_4_5_cent_30_50;
    histograms_pt_4_5_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_4_5_cent_50_90;
    histograms_pt_4_5_cent_50_90.resize(bdt_values.size());

    std::vector<TH1D*> histograms_pt_5_6_cent_0_10;
    histograms_pt_5_6_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_5_6_cent_10_30;
    histograms_pt_5_6_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_5_6_cent_30_50;
    histograms_pt_5_6_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_5_6_cent_50_90;
    histograms_pt_5_6_cent_50_90.resize(bdt_values.size());

    std::vector<TH1D*> histograms_pt_6_8_cent_0_10;
    histograms_pt_6_8_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_6_8_cent_10_30;
    histograms_pt_6_8_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_6_8_cent_30_50;
    histograms_pt_6_8_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_6_8_cent_50_90;
    histograms_pt_6_8_cent_50_90.resize(bdt_values.size());

    std::vector<TH1D*> histograms_pt_8_10_cent_0_10;
    histograms_pt_8_10_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_8_10_cent_10_30;
    histograms_pt_8_10_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_8_10_cent_30_50;
    histograms_pt_8_10_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_8_10_cent_50_90;
    histograms_pt_8_10_cent_50_90.resize(bdt_values.size());

    std::vector<TH1D*> histograms_pt_10_15_cent_0_10;
    histograms_pt_10_15_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_10_15_cent_10_30;
    histograms_pt_10_15_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_10_15_cent_30_50;
    histograms_pt_10_15_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_10_15_cent_50_90;
    histograms_pt_10_15_cent_50_90.resize(bdt_values.size());

    std::vector<TH1D*> histograms_pt_15_20_cent_0_10;
    histograms_pt_15_20_cent_0_10.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_15_20_cent_10_30;
    histograms_pt_15_20_cent_10_30.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_15_20_cent_30_50;
    histograms_pt_15_20_cent_30_50.resize(bdt_values.size());
    std::vector<TH1D*> histograms_pt_15_20_cent_50_90;
    histograms_pt_15_20_cent_50_90.resize(bdt_values.size());


        for (int j = 0; j < bdt_values.size(); j++)
        {
            // you need a histo for each pt/cent bin and bdt value
            histograms_pt_1_2_cent_0_10[j] = new TH1D(Form("pt_1_2_cent_0_10_bdt_%d", j), Form("pt_1_2_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_1_2_cent_10_30[j] = new TH1D(Form("pt_1_2_cent_10_30_bdt_%d", j), Form("pt_1_2_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_1_2_cent_30_50[j] = new TH1D(Form("pt_1_2_cent_30_50_bdt_%d", j), Form("pt_1_2_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_1_2_cent_50_90[j] = new TH1D(Form("pt_1_2_cent_50_90_bdt_%d", j), Form("pt_1_2_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

            histograms_pt_2_3_cent_0_10[j] = new TH1D(Form("pt_2_3_cent_0_10_bdt_%d", j), Form("pt_2_3_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_2_3_cent_10_30[j] = new TH1D(Form("pt_2_3_cent_10_30_bdt_%d", j), Form("pt_2_3_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_2_3_cent_30_50[j] = new TH1D(Form("pt_2_3_cent_30_50_bdt_%d", j), Form("pt_2_3_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_2_3_cent_50_90[j] = new TH1D(Form("pt_2_3_cent_50_90_bdt_%d", j), Form("pt_2_3_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

            histograms_pt_3_4_cent_0_10[j] = new TH1D(Form("pt_3_4_cent_0_10_bdt_%d", j), Form("pt_3_4_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_3_4_cent_10_30[j] = new TH1D(Form("pt_3_4_cent_10_30_bdt_%d", j), Form("pt_3_4_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_3_4_cent_30_50[j] = new TH1D(Form("pt_3_4_cent_30_50_bdt_%d", j), Form("pt_3_4_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_3_4_cent_50_90[j] = new TH1D(Form("pt_3_4_cent_50_90_bdt_%d", j), Form("pt_3_4_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

            histograms_pt_4_5_cent_0_10[j] = new TH1D(Form("pt_4_5_cent_0_10_bdt_%d", j), Form("pt_4_5_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_4_5_cent_10_30[j] = new TH1D(Form("pt_4_5_cent_10_30_bdt_%d", j), Form("pt_4_5_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_4_5_cent_30_50[j] = new TH1D(Form("pt_4_5_cent_30_50_bdt_%d", j), Form("pt_4_5_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_4_5_cent_50_90[j] = new TH1D(Form("pt_4_5_cent_50_90_bdt_%d", j), Form("pt_4_5_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

            histograms_pt_5_6_cent_0_10[j] = new TH1D(Form("pt_5_6_cent_0_10_bdt_%d", j), Form("pt_5_6_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_5_6_cent_10_30[j] = new TH1D(Form("pt_5_6_cent_10_30_bdt_%d", j), Form("pt_5_6_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_5_6_cent_30_50[j] = new TH1D(Form("pt_5_6_cent_30_50_bdt_%d", j), Form("pt_5_6_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_5_6_cent_50_90[j] = new TH1D(Form("pt_5_6_cent_50_90_bdt_%d", j), Form("pt_5_6_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

            histograms_pt_6_8_cent_0_10[j] = new TH1D(Form("pt_6_8_cent_0_10_bdt_%d", j), Form("pt_6_8_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_6_8_cent_10_30[j] = new TH1D(Form("pt_6_8_cent_10_30_bdt_%d", j), Form("pt_6_8_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_6_8_cent_30_50[j] = new TH1D(Form("pt_6_8_cent_30_50_bdt_%d", j), Form("pt_6_8_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_6_8_cent_50_90[j] = new TH1D(Form("pt_6_8_cent_50_90_bdt_%d", j), Form("pt_6_8_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

            histograms_pt_8_10_cent_0_10[j] = new TH1D(Form("pt_8_10_cent_0_10_bdt_%d", j), Form("pt_8_10_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_8_10_cent_10_30[j] = new TH1D(Form("pt_8_10_cent_10_30_bdt_%d", j), Form("pt_8_10_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_8_10_cent_30_50[j] = new TH1D(Form("pt_8_10_cent_30_50_bdt_%d", j), Form("pt_8_10_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_8_10_cent_50_90[j] = new TH1D(Form("pt_8_10_cent_50_90_bdt_%d", j), Form("pt_8_10_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

            histograms_pt_10_15_cent_0_10[j] = new TH1D(Form("pt_10_15_cent_0_10_bdt_%d", j), Form("pt_10_15_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_10_15_cent_10_30[j] = new TH1D(Form("pt_10_15_cent_10_30_bdt_%d", j), Form("pt_10_15_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_10_15_cent_30_50[j] = new TH1D(Form("pt_10_15_cent_30_50_bdt_%d", j), Form("pt_10_15_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_10_15_cent_50_90[j] = new TH1D(Form("pt_10_15_cent_50_90_bdt_%d", j), Form("pt_10_15_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

            histograms_pt_15_20_cent_0_10[j] = new TH1D(Form("pt_15_20_cent_0_10_bdt_%d", j), Form("pt_15_20_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_15_20_cent_10_30[j] = new TH1D(Form("pt_15_20_cent_10_30_bdt_%d", j), Form("pt_15_20_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_15_20_cent_30_50[j] = new TH1D(Form("pt_15_20_cent_30_50_bdt_%d", j), Form("pt_15_20_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
            histograms_pt_15_20_cent_50_90[j] = new TH1D(Form("pt_15_20_cent_50_90_bdt_%d", j), Form("pt_15_20_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        }

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

        TFile *inf = TFile::Open(filename.c_str());
        cout << "file path = " << filename.c_str() << endl;
        if (inf->IsZombie() || inf->Get("mvaTree") == nullptr)
        {
            cout << "BIG PROBLEM" << endl;
            ifile++;
            continue;
        }

        TTree *t = (TTree *)inf->Get("mvaTree");
        t->SetMakeClass(1);
        for (int i = 0; i < t->GetEntries(); i++)
        // for(int i=0;i<1000;i++)
        {
            t->SetBranchAddress("candSize", &candSize);
            t->SetBranchAddress("pT", pT);
            t->SetBranchAddress("mass", mass);
            t->SetBranchAddress("centrality", &centrality);
            t->SetBranchAddress("y", y);
            t->SetBranchAddress("BDT_weight", BDT_weight);
            t->GetEntry(i);
            // cout << "candSize= " << candSize << endl;
            if (i % 100000 == 0)
                cout << i << " / " << t->GetEntries() << endl;

            //   for(int j=0; j<candSize; j++)
            for (int j = 0; j < candSize; j++)
            {
                if (y[j] > -1 && y[j] < 1)
                {
                    // cout << "y=" << y[j] << endl;
                    if (centrality >= 0 && centrality < 2 * 10)
                    {
                        if (pT[j] >= 1 && pT[j] < 2)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_1_2_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 2 && pT[j] < 3)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_2_3_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 3 && pT[j] < 4)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_3_4_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 4 && pT[j] < 5)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_4_5_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 5 && pT[j] < 6)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_5_6_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 6 && pT[j] < 8)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_6_8_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 8 && pT[j] < 10)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_8_10_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 10 && pT[j] < 15)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_10_15_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 15 && pT[j] < 20)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_15_20_cent_0_10[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                    } // if centrality
                    if (centrality >= 2 * 10 && centrality < 2 * 30)
                    {
                        if (pT[j] >= 1 && pT[j] < 2)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_1_2_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 2 && pT[j] < 3)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_2_3_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 3 && pT[j] < 4)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_3_4_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 4 && pT[j] < 5)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_4_5_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 5 && pT[j] < 6)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_5_6_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 6 && pT[j] < 8)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_6_8_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 8 && pT[j] < 10)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_8_10_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 10 && pT[j] < 15)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_10_15_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 15 && pT[j] < 20)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_15_20_cent_10_30[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                    } // if centrality
                    if (centrality >= 2 * 30 && centrality < 2 * 50)
                    {
                        if (pT[j] >= 1 && pT[j] < 2)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_1_2_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 2 && pT[j] < 3)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_2_3_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 3 && pT[j] < 4)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_3_4_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 4 && pT[j] < 5)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_4_5_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 5 && pT[j] < 6)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_5_6_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 6 && pT[j] < 8)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_6_8_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 8 && pT[j] < 10)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_8_10_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 10 && pT[j] < 15)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_10_15_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 15 && pT[j] < 20)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_15_20_cent_30_50[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                    } // if centrality
                    if (centrality >= 2 * 50 && centrality < 2 * 90)
                    {
                        if (pT[j] >= 1 && pT[j] < 2)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_1_2_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 2 && pT[j] < 3)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_2_3_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 3 && pT[j] < 4)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_3_4_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 4 && pT[j] < 5)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_4_5_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 5 && pT[j] < 6)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_5_6_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 6 && pT[j] < 8)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_6_8_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 8 && pT[j] < 10)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_8_10_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 10 && pT[j] < 15)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_10_15_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                        if (pT[j] >= 15 && pT[j] < 20)
                        {
                            for (int r = 0; r < bdt_values.size(); r++)
                            {
                                if (BDT_weight[j] > bdt_values[r])
                                {
                                    histograms_pt_15_20_cent_50_90[r]->Fill(mass[j]);
                                }
                            }
                        } // if pT
                    } // if centrality
                } // if y
            } // for candSize

        } // for entries
        ifile++;
        inf->Close();
    }//while 
        outf->cd();
        outf->Write();
        outf->Close();
        cout << "DONE --" << endl;
        return 0;
    } // function

int main(int argc, char *argv[])
{
    if(istest) argc=3;
    if(argc==3)
    {
        createBdtMassDistributions(argv[1], argv[2]);
    }
    else
    {
        std::cout << "Usage: <istart> <iend>" << std::endl;
        return 1;
    }
    return 0;
}
