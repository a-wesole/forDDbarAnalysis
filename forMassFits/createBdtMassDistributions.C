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

using namespace std;
Bool_t istest = true;
int createBdtMassDistributions(TString inputfile="", TString outfile1="")
{


    if(istest)
    {
        inputfile="/scratch/negishi/awesole/dataWithBdt/ROOT/outfile_ttree_2280_2310.root";
        //inputfile="/home/awesole/forDDbarAnalysis/BDT_training_output/BDT_signal.root";
        outfile1="test.root";
        cout << "TEST" << endl;
    }



    TFile* inf = TFile::Open(inputfile);
    TTree *t = (TTree*)inf->Get("mvaTree");
    t->SetMakeClass(1);
    TFile *outf = new TFile(outfile1,"recreate");

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

    for (int j=0; j<bdt_values.size(); j++)
    {
     //you need a histo for each pt/cent bin and bdt value 
        histograms_pt_1_2_cent_0_10[j] = new TH1D(Form("pt_1_2_cent_0_10_bdt_%d", j), Form("pt_1_2_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_1_2_cent_10_30[j] = new TH1D(Form("pt_1_2_cent_10_30_bdt_%d", j), Form("pt_1_2_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_1_2_cent_30_50[j] = new TH1D(Form("pt_1_2_cent_30_50_bdt_%d", j), Form("pt_1_2_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_1_2_cent_50_90[j] = new TH1D(Form("pt_1_2_cent_50_90_bdt_%d", j), Form("pt_1_2_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

        /*
        histograms_pt_2_3_cent_0_10[j] = new TH1D(Form("pt_2_3_cent_0_10_bdt_%f", bdt_values[j]), Form("pt_2_3_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_2_3_cent_10_30[j] = new TH1D(Form("pt_2_3_cent_10_30_bdt_%f", bdt_values[j]), Form("pt_2_3_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_2_3_cent_30_50[j] = new TH1D(Form("pt_2_3_cent_30_50_bdt_%f", bdt_values[j]), Form("pt_2_3_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_2_3_cent_50_90[j] = new TH1D(Form("pt_2_3_cent_50_90_bdt_%f", bdt_values[j]), Form("pt_2_3_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

        histograms_pt_3_4_cent_0_10[j] = new TH1D(Form("pt_3_4_cent_0_10_bdt_%f", bdt_values[j]), Form("pt_3_4_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_3_4_cent_10_30[j] = new TH1D(Form("pt_3_4_cent_10_30_bdt_%f", bdt_values[j]), Form("pt_3_4_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_3_4_cent_30_50[j] = new TH1D(Form("pt_3_4_cent_30_50_bdt_%f", bdt_values[j]), Form("pt_3_4_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_3_4_cent_50_90[j] = new TH1D(Form("pt_3_4_cent_50_90_bdt_%f", bdt_values[j]), Form("pt_3_4_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

        histograms_pt_4_5_cent_0_10[j] = new TH1D(Form("pt_4_5_cent_0_10_bdt_%f", bdt_values[j]), Form("pt_4_5_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_4_5_cent_10_30[j] = new TH1D(Form("pt_4_5_cent_10_30_bdt_%f", bdt_values[j]), Form("pt_4_5_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_4_5_cent_30_50[j] = new TH1D(Form("pt_4_5_cent_30_50_bdt_%f", bdt_values[j]), Form("pt_4_5_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_4_5_cent_50_90[j] = new TH1D(Form("pt_4_5_cent_50_90_bdt_%f", bdt_values[j]), Form("pt_4_5_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

        histograms_pt_5_6_cent_0_10[j] = new TH1D(Form("pt_5_6_cent_0_10_bdt_%f", bdt_values[j]), Form("pt_5_6_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_5_6_cent_10_30[j] = new TH1D(Form("pt_5_6_cent_10_30_bdt_%f", bdt_values[j]), Form("pt_5_6_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_5_6_cent_30_50[j] = new TH1D(Form("pt_5_6_cent_30_50_bdt_%f", bdt_values[j]), Form("pt_5_6_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_5_6_cent_50_90[j] = new TH1D(Form("pt_5_6_cent_50_90_bdt_%f", bdt_values[j]), Form("pt_5_6_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

        histograms_pt_6_8_cent_0_10[j] = new TH1D(Form("pt_6_8_cent_0_10_bdt_%f", bdt_values[j]), Form("pt_6_8_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_6_8_cent_10_30[j] = new TH1D(Form("pt_6_8_cent_10_30_bdt_%f", bdt_values[j]), Form("pt_6_8_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_6_8_cent_30_50[j] = new TH1D(Form("pt_6_8_cent_30_50_bdt_%f", bdt_values[j]), Form("pt_6_8_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_6_8_cent_50_90[j] = new TH1D(Form("pt_6_8_cent_50_90_bdt_%f", bdt_values[j]), Form("pt_6_8_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

        histograms_pt_8_10_cent_0_10[j] = new TH1D(Form("pt_8_10_cent_0_10_bdt_%f", bdt_values[j]), Form("pt_8_10_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_8_10_cent_10_30[j] = new TH1D(Form("pt_8_10_cent_10_30_bdt_%f", bdt_values[j]), Form("pt_8_10_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_8_10_cent_30_50[j] = new TH1D(Form("pt_8_10_cent_30_50_bdt_%f", bdt_values[j]), Form("pt_8_10_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_8_10_cent_50_90[j] = new TH1D(Form("pt_8_10_cent_50_90_bdt_%f", bdt_values[j]), Form("pt_8_10_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

        histograms_pt_10_15_cent_0_10[j] = new TH1D(Form("pt_10_15_cent_0_10_bdt_%f", bdt_values[j]), Form("pt_10_15_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_10_15_cent_10_30[j] = new TH1D(Form("pt_10_15_cent_10_30_bdt_%f", bdt_values[j]), Form("pt_10_15_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_10_15_cent_30_50[j] = new TH1D(Form("pt_10_15_cent_30_50_bdt_%f", bdt_values[j]), Form("pt_10_15_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_10_15_cent_50_90[j] = new TH1D(Form("pt_10_15_cent_50_90_bdt_%f", bdt_values[j]), Form("pt_10_15_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);

        histograms_pt_15_20_cent_0_10[j] = new TH1D(Form("pt_15_20_cent_0_10_bdt_%f", bdt_values[j]), Form("pt_15_20_cent_0_10_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_15_20_cent_10_30[j] = new TH1D(Form("pt_15_20_cent_10_30_bdt_%f", bdt_values[j]), Form("pt_15_20_cent_10_30_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_15_20_cent_30_50[j] = new TH1D(Form("pt_15_20_cent_30_50_bdt_%f", bdt_values[j]), Form("pt_15_20_cent_30_50_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        histograms_pt_15_20_cent_50_90[j] = new TH1D(Form("pt_15_20_cent_50_90_bdt_%f", bdt_values[j]), Form("pt_15_20_cent_50_90_bdt_%f", bdt_values[j]), 52, 1.74, 2.);
        */


    }

    for(int i=0;i<t->GetEntries();i++)
    //for(int i=0;i<1000;i++)
    {
        t->SetBranchAddress("candSize", &candSize); 
        t->SetBranchAddress("pT", pT);
        t->SetBranchAddress("mass", mass);
        t->SetBranchAddress("centrality", &centrality);
        t->SetBranchAddress("y", y);
        t->SetBranchAddress("BDT_weight", BDT_weight);
        t->GetEntry(i);
        // cout << "candSize= " << candSize << endl;
        if(i%100000==0) cout <<i<<" / "<<t->GetEntries()<<endl;

        //   for(int j=0; j<candSize; j++)
        for(int j=0; j<candSize; j++)
        { 
            if(y[j]>-1 && y[j]<1)
            {
                  //cout << "y=" << y[j] << endl;
                if (centrality>=0 && centrality<20 )
                {
                    if (pT[j]>=1 && pT[j]<2){
                        if(BDT_weight[j]>-0.6) histograms_pt_1_2_cent_0_10[0]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.55) histograms_pt_1_2_cent_0_10[1]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.5) histograms_pt_1_2_cent_0_10[2]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.45) histograms_pt_1_2_cent_0_10[3]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.4) histograms_pt_1_2_cent_0_10[4]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.35) histograms_pt_1_2_cent_0_10[5]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.3) histograms_pt_1_2_cent_0_10[6]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.25) histograms_pt_1_2_cent_0_10[7]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.2) histograms_pt_1_2_cent_0_10[8]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.15) histograms_pt_1_2_cent_0_10[9]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.1) histograms_pt_1_2_cent_0_10[10]->Fill(mass[j]);
                        if(BDT_weight[j]>-0.05) histograms_pt_1_2_cent_0_10[11]->Fill(mass[j]);
                        if(BDT_weight[j]>0.0) histograms_pt_1_2_cent_0_10[12]->Fill(mass[j]);
                        if(BDT_weight[j]>0.05) histograms_pt_1_2_cent_0_10[13]->Fill(mass[j]);
                        if(BDT_weight[j]>0.1) histograms_pt_1_2_cent_0_10[14]->Fill(mass[j]);
                        if(BDT_weight[j]>0.15) histograms_pt_1_2_cent_0_10[15]->Fill(mass[j]);
                        if(BDT_weight[j]>0.2) histograms_pt_1_2_cent_0_10[16]->Fill(mass[j]);
                        if(BDT_weight[j]>0.25) histograms_pt_1_2_cent_0_10[17]->Fill(mass[j]);
                        if(BDT_weight[j]>0.3) histograms_pt_1_2_cent_0_10[18]->Fill(mass[j]);
                        if(BDT_weight[j]>0.35) histograms_pt_1_2_cent_0_10[19]->Fill(mass[j]);
                        if(BDT_weight[j]>0.4) histograms_pt_1_2_cent_0_10[20]->Fill(mass[j]);
                        if(BDT_weight[j]>0.45) histograms_pt_1_2_cent_0_10[21]->Fill(mass[j]);
                        if(BDT_weight[j]>0.5) histograms_pt_1_2_cent_0_10[22]->Fill(mass[j]);
                        if(BDT_weight[j]>0.55) histograms_pt_1_2_cent_0_10[23]->Fill(mass[j]);
                        if(BDT_weight[j]>0.6) histograms_pt_1_2_cent_0_10[24]->Fill(mass[j]);
                   }//if pT
                }//if centrality
            }//if y
        }//for candSize

    }//for entries 
    outf->cd();
    outf->Write();
    outf->Close();
    cout << "DONE --" << endl;
    return 0;
}//function 


int main(int argc, char *argv[])
{
    if(istest) argc=3;
    if(argc==3)
    {
        createBdtMassDistributions(argv[1], argv[2]);
    }
    else
    {
        std::cout << "Usage: mergeForest <input_collection> <output_file>" << std::endl;
        return 1;
    }
    return 0;
}
