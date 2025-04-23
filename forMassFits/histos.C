// GENERATES HISTO GRAM FOR SIGNAL AND SIGNAL+SWAP FOR EACH PT RANGE.
// THE FIELS ARE NEEDED FOR FITTING
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <vector>
#include <TMath.h>
#include <iostream>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <fstream>

using namespace std;
void histos(float ptmin = 1, float ptmax = 2, int centmin = 10, int centmax = 20, float ymin = 0, float ymax = 1, int istart = 1, int iend = 2)
{

    ifstream file_stream("/home/awesole/forDDbarAnalysis/dataFiles/mcTTrees.list");
     TString outfile = TString::Format("/scratch/negishi/awesole/Apr22_25_BdtTrainings/SignalMassPlots/ROOT/massPlots_pt_%g_%g_cent_%d_%d_y_%g_%g.root", ptmin, ptmax, centmin, centmax, ymin, ymax);
    //TString outfile = TString::Format("test.root");
    TFile results(outfile, "recreate");

    string filename;

    int candSize, Dgen[50000], centrality, ifile = 0;
    bool matchGEN[50000], isSwap[50000];
    float pT[50000], mass[50000], y[50000], Dtrk1Chi2n[50000], Dtrk2Chi2n[50000];

    std::vector<TH1D *> vlist;
    // vlist.resize(10);
    vlist.resize(2);
    std::vector<float> massdata;
    TH1D *signalHisto = new TH1D(Form("signalMass_pt_%g_%g_cent_%d_%d_y_%g_%g", ptmin, ptmax, centmin, centmax, ymin, ymax), Form("signalMass_pt_%g_%g_cent_%d_%d_y_%g_%g", ptmin, ptmax, centmin, centmax, ymin, ymax), 52, 1.74, 2);
    TH1D *signalSwapHisto = new TH1D(Form("signalSwapMass_pt_%g_%g_cent_%d_%d_y_%g_%g", ptmin, ptmax, centmin, centmax, ymin, ymax), Form("signalSwapMass_pt_%g_%g_cent_%d_%d_y_%g_%g", ptmin, ptmax, centmin, centmax, ymin, ymax), 52, 1.74, 2);

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

        if (ifile % 200 ==0) {
        cout << "-- -- -- -- -- -- --" << endl;
        cout << "ifile = " << ifile << endl;
        cout << "file = " << filename.c_str() << endl;
        }

        TFile *fin = TFile::Open(filename.c_str());
        if (!fin)
            break;
        if (fin->IsZombie() || fin->Get("d0ana_newreduced/VertexCompositeNtuple") == nullptr)
        {
            ifile++;
            continue;
        }

        TTree *t = (TTree *)fin->Get("d0ana_newreduced/VertexCompositeNtuple");
        t->SetMakeClass(1);

        t->SetBranchAddress("candSize", &candSize);
        t->SetBranchAddress("centrality", &centrality);
        t->SetBranchAddress("pT", pT);
        t->SetBranchAddress("mass", mass);
        t->SetBranchAddress("matchGEN", matchGEN);
        t->SetBranchAddress("isSwap", isSwap);
        t->SetBranchAddress("y", y);
        for (int i = 0; i < t->GetEntries(); i++)
        {
            t->GetEntry(i);

            /*
            if (i % 1000 == 0)
                cout << i << " / " << t->GetEntries() << endl;
            //cout << "candSize = " << candSize << endl;
            */
            if (candSize == 0)
                continue;
            if (centrality >= 2 * centmax || centrality < 2 * centmin)
                continue;
            for (int j = 0; j < candSize; j++)
            {
                if (abs(y[j]) < ymin || abs(y[j]) >= ymax)
                    continue;
                if (matchGEN[j] != 1)
                    continue;
                if (pT[j] < ptmin || pT[j] >= ptmax)
                    continue;
                if (pT[j] >= ptmin && pT[j] < ptmax && centrality >= 2 * centmin && centrality < 2 * centmax && abs(y[j]) >= ymin && abs(y[j]) < ymax)
                {
                    if (matchGEN[j] == 1)
                    {
                        signalSwapHisto->Fill(mass[j]);
                        if (matchGEN[j] == 1 && isSwap[j] == 0)
                        {
                            signalHisto->Fill(mass[j]);
                        } // if dgen ==23333
                    } // if dgen==23333||23344
                }
            } // for nDsizpe
        } // for entries
        ifile++;
        fin->Close();
    } // while
    results.Write();
    results.Close();
} // function

int main(int argc, char *argv[])
{
    if (argc == 9) // Expecting 8 arguments + program name
    {
        float ptmin = atof(argv[1]);
        float ptmax = atof(argv[2]);
        int centmin = atoi(argv[3]);
        int centmax = atoi(argv[4]);
        float ymin = atoi(argv[5]);
        float ymax = atoi(argv[6]);
        int istart = std::stoi(argv[7]);
        int iend = std::stoi(argv[8]);

        histos(ptmin, ptmax, centmin, centmax, ymin, ymax, istart, iend);
    }
    else
    {
        std::cout << "Usage: mergeForest <ptmin> <ptmax> <centmin> <centmax> <ymin> <ymax> <istart> <iend>" << std::endl;
        return 1;
    }
    return 0;
}
